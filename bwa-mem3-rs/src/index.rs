//! Reference index handle.

use std::ffi::{CStr, CString};
use std::path::Path;

use crate::error::{shim_err, Error, Result};

/// On-disk files bwa-mem3 opens while loading an index from disk. Each is read
/// with an `exit()`-on-failure path inside the vendored C++, so a missing file
/// would abort the host process before [`BwaIndex::load`] could return (issue
/// #19). We check for them up front and surface a recoverable [`Error`]
/// instead. The shim builds its reference string from `.pac` in memory and
/// never reads `.0123`, so that suffix is deliberately absent here.
const REQUIRED_INDEX_SUFFIXES: &[&str] = &[".bwt.2bit.64", ".amb", ".ann", ".pac"];

/// The files a `--meth` dual index needs for its *original* reference. Seeding
/// runs against the converted seed index; the original is loaded by the shim's
/// `shim_meth_orig_ref_load` (`bwa_shim_align.cpp`), which reads only the
/// `bns_restore` files and `.pac`, never the original's FM-index
/// (`.bwt.2bit.64`). An index built with `bwa-mem3 index --meth --meth-only`,
/// which packs the original without an FM-index, therefore loads.
const ORIG_REF_SUFFIXES: &[&str] = &[".amb", ".ann", ".pac"];

/// [`check_files`] for a plain index's [`REQUIRED_INDEX_SUFFIXES`]. When the
/// FM-index is missing but `<prefix>.meth`'s exists, the prefix is most likely
/// the original reference of a `--meth --meth-only` dual index, so the error
/// points at [`BwaIndex::load_meth`].
fn check_index_files(prefix: &Path) -> Result<()> {
    check_files(prefix, REQUIRED_INDEX_SUFFIXES, "index").map_err(|err| {
        let Error::IndexLoad { path, msg } = err else {
            return err;
        };
        let suffixed = |suffix: &str| {
            let mut file = prefix.as_os_str().to_owned();
            file.push(suffix);
            std::path::PathBuf::from(file)
        };
        let msg = if !is_readable_file(&suffixed(".bwt.2bit.64"))
            && is_readable_file(&suffixed(".meth.bwt.2bit.64"))
        {
            format!(
                "{msg} (this looks like the original reference of a `bwa-mem3 index --meth \
                 --meth-only` index; load it with BwaIndex::load_meth(\"{}.meth\", \"{}\"))",
                prefix.display(),
                prefix.display()
            )
        } else {
            msg
        };
        Error::IndexLoad { path, msg }
    })
}

/// Whether `path` is a regular file that can be opened for reading. The
/// `is_file` guard short-circuits before [`File::open`](std::fs::File::open),
/// so a FIFO (whose `O_RDONLY` open would block) is rejected without ever
/// being opened.
fn is_readable_file(path: &Path) -> bool {
    path.is_file() && std::fs::File::open(path).is_ok()
}

/// Return [`Error::IndexLoad`] naming every `<prefix><suffix>` file in
/// `suffixes` that is not a readable regular file; `what` names the files in
/// the message (e.g. "index"). `Ok(())` means all are present and openable.
///
/// A bare existence check is not enough: a directory, FIFO, or unreadable
/// path satisfies [`Path::exists`] yet still trips the vendored loader's
/// `exit(1)` on the subsequent `fopen`/`fread`. We mirror what the loader
/// actually needs, a regular file it can open for reading (see
/// [`is_readable_file`]).
fn check_files(prefix: &Path, suffixes: &[&str], what: &str) -> Result<()> {
    let base = prefix.as_os_str();
    let unusable: Vec<String> = suffixes
        .iter()
        .map(|suffix| {
            let mut file = base.to_owned();
            file.push(suffix);
            file
        })
        .filter(|file| !is_readable_file(Path::new(file)))
        .map(|file| Path::new(&file).display().to_string())
        .collect();
    if unusable.is_empty() {
        Ok(())
    } else {
        Err(Error::IndexLoad {
            path: prefix.to_owned(),
            msg: format!(
                "missing or unreadable {what} file(s): {}",
                unusable.join(", ")
            ),
        })
    }
}

/// Validate `prefix` as UTF-8 and convert it to a `CString` for the FFI.
/// `what` names the argument for the error message (e.g. "prefix").
fn prefix_to_cstring(prefix: &Path, what: &str) -> Result<CString> {
    let s = prefix
        .to_str()
        .ok_or_else(|| Error::InvalidInput(format!("{what} must be valid UTF-8")))?;
    Ok(CString::new(s)?)
}

/// Handle to a loaded bwa-mem3 reference index.
///
/// The index is immutable after loading. Multiple threads may share a
/// `BwaIndex` via `&BwaIndex` or `Arc<BwaIndex>` for concurrent alignment.
pub struct BwaIndex {
    handle: *mut bwa_mem3_sys::BwaIndex,
    /// Unique per loaded index for the process lifetime (never reused, unlike
    /// the handle's address), so state built against one index can recognize
    /// a different one.
    id: u64,
}

/// Source of [`BwaIndex`] identities; starts at 1 so 0 can mean "none".
static NEXT_INDEX_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

fn next_index_id() -> u64 {
    NEXT_INDEX_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

impl BwaIndex {
    /// Load a prebuilt bwa-mem3 index.
    ///
    /// `prefix` is the path without extension; bwa-mem3 appends its own
    /// suffixes (`.bwt.2bit.64`, `.ann`, `.pac`, etc) internally.
    ///
    /// If a shared-memory segment for `prefix` is staged (see
    /// [`crate::shm::stage`]), bwa-mem3 transparently attaches to it
    /// instead of reading from disk. Callers that want fail-fast on a
    /// missing segment should probe first:
    ///
    /// ```ignore
    /// if !bwa_mem3_rs::shm::is_staged(&prefix)? {
    ///     anyhow::bail!("index not staged in shm");
    /// }
    /// let idx = bwa_mem3_rs::BwaIndex::load(&prefix)?;
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Error::IndexLoad`] if a required on-disk index file is
    /// missing or not a readable regular file (issue #19: the vendored loader
    /// would otherwise call `exit(1)` and abort the process). When a
    /// shared-memory segment is staged for `prefix`, the disk check is skipped
    /// since bwa-mem3 attaches to the segment instead of reading from disk.
    pub fn load(prefix: impl AsRef<Path>) -> Result<Self> {
        let path = prefix.as_ref();

        // Guard against the vendored loader's `exit()` paths: a staged shm
        // segment is attached in lieu of disk files, so only validate the
        // on-disk index when nothing is staged for this prefix. If the shm
        // probe itself fails, fall back to the disk check (the safe default).
        if !crate::shm::is_staged(path).unwrap_or(false) {
            check_index_files(path)?;
        }

        let c = prefix_to_cstring(path, "prefix")?;
        let handle = unsafe { bwa_mem3_sys::bwa_shim_idx_load(c.as_ptr()) };
        if handle.is_null() {
            return Err(Error::IndexLoad {
                path: path.to_owned(),
                msg: shim_err("idx load").to_string(),
            });
        }
        Ok(BwaIndex {
            handle,
            id: next_index_id(),
        })
    }

    /// As [`load`](Self::load), but loads the FM-index with `n_threads`
    /// threads (the CLI's `-t` behavior for index load). `0` is treated as `1`.
    ///
    /// # Errors
    ///
    /// Same as [`load`](Self::load).
    pub fn load_with_threads(prefix: impl AsRef<Path>, n_threads: usize) -> Result<Self> {
        let path = prefix.as_ref();

        if !crate::shm::is_staged(path).unwrap_or(false) {
            check_index_files(path)?;
        }

        let c = prefix_to_cstring(path, "prefix")?;
        let n = i32::try_from(n_threads.max(1)).unwrap_or(i32::MAX);
        // SAFETY: `c` is a valid, NUL-terminated C string for the duration of
        // this call.
        let handle = unsafe { bwa_mem3_sys::bwa_shim_idx_load_threads(c.as_ptr(), n) };
        if handle.is_null() {
            return Err(Error::IndexLoad {
                path: path.to_owned(),
                msg: shim_err("idx load").to_string(),
            });
        }
        Ok(BwaIndex {
            handle,
            id: next_index_id(),
        })
    }

    /// Load a bisulfite (BS-seq) **dual index** for `--meth` alignment.
    ///
    /// `seed_prefix` is the converted seed index (`<ref>.meth`, built by
    /// `bwa-mem3 index --meth`); `orig_prefix` is the un-converted original
    /// reference (`<ref>`). Both stay resident: seeding runs against the
    /// converted index, while chaining/extension/output run in original
    /// coordinates. Align with a [`MemOpts`](crate::MemOpts) whose
    /// [`set_meth`](crate::MemOpts::set_meth) is enabled; the emitted records
    /// carry Bismark `XG`/`XR`/`XM` tags and original-reference coordinates, so
    /// the contig accessors below report the original contigs.
    ///
    /// The original reference needs only its `.amb`, `.ann` and `.pac` files,
    /// not its FM-index, so an index built with `bwa-mem3 index --meth
    /// --meth-only` loads too.
    pub fn load_meth(seed_prefix: impl AsRef<Path>, orig_prefix: impl AsRef<Path>) -> Result<Self> {
        Self::load_meth_with_threads(seed_prefix, orig_prefix, 1)
    }

    /// As [`load_meth`](Self::load_meth), but loads the converted seed FM-index
    /// with `n_threads` threads, as `bwa-mem3 mem --meth -t` does. The original
    /// reference is read single-threaded either way. `0` is treated as `1`.
    ///
    /// # Errors
    ///
    /// Same as [`load_meth`](Self::load_meth).
    pub fn load_meth_with_threads(
        seed_prefix: impl AsRef<Path>,
        orig_prefix: impl AsRef<Path>,
        n_threads: usize,
    ) -> Result<Self> {
        let seed = seed_prefix.as_ref();
        let orig = orig_prefix.as_ref();

        if !crate::shm::is_staged(seed).unwrap_or(false) {
            check_index_files(seed)?;
        }
        check_files(orig, ORIG_REF_SUFFIXES, "original-reference")?;

        let seed_c = prefix_to_cstring(seed, "seed prefix")?;
        let orig_c = prefix_to_cstring(orig, "orig prefix")?;
        let n = i32::try_from(n_threads.max(1)).unwrap_or(i32::MAX);
        // SAFETY: `seed_c` and `orig_c` are valid, NUL-terminated C strings for
        // the duration of this call.
        let handle = unsafe {
            bwa_mem3_sys::bwa_shim_idx_load_meth_threads(seed_c.as_ptr(), orig_c.as_ptr(), n)
        };
        if handle.is_null() {
            return Err(Error::IndexLoad {
                path: seed.to_owned(),
                msg: shim_err("meth idx load").to_string(),
            });
        }
        Ok(BwaIndex {
            handle,
            id: next_index_id(),
        })
    }

    /// Whether this index was loaded as a bisulfite dual index via
    /// [`load_meth`](Self::load_meth). Aligning requires
    /// [`MemOpts::meth`](crate::MemOpts::meth) to agree with this.
    #[must_use]
    pub fn is_meth(&self) -> bool {
        unsafe { bwa_mem3_sys::bwa_shim_idx_is_meth(self.handle) != 0 }
    }

    #[must_use]
    pub fn n_contigs(&self) -> usize {
        unsafe { bwa_mem3_sys::bwa_shim_idx_n_contigs(self.handle) }
    }

    /// Contig `i`'s name, or `""` when there is no such contig or its name is
    /// not UTF-8.
    #[must_use]
    pub fn contig_name(&self, i: usize) -> &str {
        unsafe {
            let c = bwa_mem3_sys::bwa_shim_idx_contig_name(self.handle, i);
            if c.is_null() {
                ""
            } else {
                CStr::from_ptr(c).to_str().unwrap_or("")
            }
        }
    }

    /// Contig `i`'s name, or an error where [`contig_name`](Self::contig_name)
    /// would return `""`: no such contig, an empty name, or one that is not
    /// UTF-8.
    fn checked_contig_name(&self, i: usize) -> Result<&str> {
        // SAFETY: the handle is live; the shim returns NULL or a pointer to a
        // NUL-terminated name owned by the index, which outlives `&self`.
        let name = unsafe {
            let c = bwa_mem3_sys::bwa_shim_idx_contig_name(self.handle, i);
            (!c.is_null()).then(|| CStr::from_ptr(c))
        };
        match name.map(CStr::to_str) {
            Some(Ok(n)) if !n.is_empty() => Ok(n),
            Some(Ok(_)) => Err(Error::Shim(format!(
                "bam_header: contig {i} has an empty name"
            ))),
            Some(Err(e)) => Err(Error::Shim(format!(
                "bam_header: contig {i}'s name is not UTF-8: {e}"
            ))),
            None => Err(Error::Shim(format!(
                "bam_header: no contig {i} in the index"
            ))),
        }
    }

    #[must_use]
    pub fn contig_len(&self, i: usize) -> i64 {
        unsafe { bwa_mem3_sys::bwa_shim_idx_contig_len(self.handle, i) }
    }

    pub fn contigs(&self) -> impl Iterator<Item = (&str, i64)> + '_ {
        (0..self.n_contigs()).map(move |i| (self.contig_name(i), self.contig_len(i)))
    }

    /// The SAM header text `bwa-mem3 mem` writes for this index under `opts`,
    /// built by bwa-mem3's own `bwa_build_sam_hdr`: the compat target's `@HD`
    /// (see [`MemOpts::compat_hd_line`](crate::MemOpts::compat_hd_line)), then
    /// one `@SQ` per contig, carrying `AH:*` on an ALT contig. For a `--meth`
    /// dual index the contigs are the original reference's. Every line ends
    /// in `\n`.
    ///
    /// This is the generated header, with none of the CLI's overrides. Callers
    /// may append `@PG`, which names the program, and `@RG` (and other `-H`
    /// lines without `@HD`/`@SQ`), which the CLI writes from `-R`/`-H` and
    /// which records carry as `RG:Z` under
    /// [`MemOpts::set_read_group_id`](crate::MemOpts::set_read_group_id).
    /// `@HD` and `@SQ` records from elsewhere must **replace or merge** these,
    /// never be appended: `bwa-mem3 mem` takes a `-H` `@HD`/`@SQ` (or else the
    /// index sidecar's, `<prefix>.hdr` / `.dict`) in place of the generated
    /// ones, and under `--meth` merges the sidecar's identity tags
    /// (`M5`/`UR`/`AS`/`SP`) into the generated `@SQ` instead. With a sidecar
    /// present, then, these lines can differ from the CLI's.
    ///
    /// Packed records name contigs by their position in this `@SQ` block (the
    /// index's contig order), and [`bam_header`](Self::bam_header) writes the
    /// binary reference list in that same order. So a caller writing through
    /// `bam_header` may merge tags into these `@SQ` lines but must keep them
    /// in index order, leaving every record's `tid` and `mtid` as emitted. A
    /// caller that reorders the contigs must instead write its own binary
    /// reference list in the new order and remap every `tid` and `mtid` to
    /// match, as upstream's BAM writer does (`bam_writer.cpp`'s `tid_map`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Shim`] if the header could not be built (allocation
    /// failure) or is not UTF-8 (it is built from the index's contig names,
    /// not from anything the caller passed).
    pub fn sam_header(&self, opts: &crate::MemOpts) -> Result<String> {
        // SAFETY: both handles are live for the call; the shim returns NULL or
        // an owned NUL-terminated buffer, copied out here and then released
        // exactly once with its matching free.
        unsafe {
            let p = bwa_mem3_sys::bwa_shim_header_text(self.handle, opts.handle);
            if p.is_null() {
                return Err(shim_err("sam_header"));
            }
            let text = CStr::from_ptr(p).to_str().map(str::to_owned);
            bwa_mem3_sys::bwa_shim_header_text_free(p);
            text.map_err(|e| Error::Shim(format!("sam_header: SAM header is not UTF-8: {e}")))
        }
    }

    /// The binary BAM header for this index around SAM header `text` (usually
    /// [`sam_header`](Self::sam_header) plus the caller's `@PG`/`@RG` lines):
    /// the `BAM\1` magic, `l_text` and `text`, then the reference list,
    /// `n_ref` and each contig's `l_name`/name/NUL/`l_ref` in the index's
    /// order, which is the order packed records' `tid`/`mtid` index. Write it
    /// once, uncompressed bytes, at the start of the BGZF stream, before the
    /// first record.
    ///
    /// The reference list is always the index's, in index order, whatever
    /// `text` says, so any `@SQ` lines in `text` must list the same contigs in
    /// the same order (see [`sam_header`](Self::sam_header)); this is not
    /// checked.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if `text`, the contig count or a contig
    /// length is not below 2^31, as the BAM format requires of `l_text`,
    /// `n_ref` and `l_ref`, or a contig name does not fit `l_name`'s 32 bits,
    /// and [`Error::Shim`] if a contig name in the index is missing, empty or
    /// not UTF-8 (which [`contig_name`](Self::contig_name) reports as `""`).
    pub fn bam_header(&self, text: &str) -> Result<Vec<u8>> {
        let l_text = bam_header_field(text.len(), "the header text length", true)?;
        let n_ref = bam_header_field(self.n_contigs(), "the contig count", true)?;
        let mut out = Vec::with_capacity(12 + text.len() + 16 * self.n_contigs());
        out.extend_from_slice(b"BAM\x01");
        out.extend_from_slice(&l_text.to_le_bytes());
        out.extend_from_slice(text.as_bytes());
        out.extend_from_slice(&n_ref.to_le_bytes());
        for i in 0..self.n_contigs() {
            let (name, len) = (self.checked_contig_name(i)?, self.contig_len(i));
            let l_name = bam_header_field(name.len() + 1, "a contig name length", false)?;
            out.extend_from_slice(&l_name.to_le_bytes());
            out.extend_from_slice(name.as_bytes());
            out.push(0);
            let len = usize::try_from(len).unwrap_or(usize::MAX);
            out.extend_from_slice(&bam_header_field(len, "a contig length", true)?.to_le_bytes());
        }
        Ok(out)
    }

    /// This index's process-unique identity (see the field doc).
    pub(crate) fn id(&self) -> u64 {
        self.id
    }

    pub(crate) fn raw(&self) -> *mut bwa_mem3_sys::BwaIndex {
        self.handle
    }
}

impl Drop for BwaIndex {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { bwa_mem3_sys::bwa_shim_idx_free(self.handle) };
        }
    }
}

// SAFETY: the shim guarantees BwaIndex is read-only after construction; no
// mutable state is exposed via any accessor. Multiple threads may share a
// `&BwaIndex` or an `Arc<BwaIndex>` for concurrent alignment.
unsafe impl Send for BwaIndex {}
unsafe impl Sync for BwaIndex {}

/// One binary-BAM-header length field as its little-endian `u32` value.
/// `below_2_31` holds it to the format's `[0, 2^31)` range (`l_text`, `n_ref`,
/// `l_ref`); otherwise it need only fit 32 bits (`l_name`).
fn bam_header_field(n: usize, what: &str, below_2_31: bool) -> Result<u32> {
    let max = if below_2_31 {
        i32::MAX as u32
    } else {
        u32::MAX
    };
    u32::try_from(n).ok().filter(|&v| v <= max).ok_or_else(|| {
        Error::InvalidInput(format!(
            "bam_header: {what} ({n}) exceeds the BAM limit of {max}"
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::{bam_header_field, check_files, check_index_files, ORIG_REF_SUFFIXES};
    use crate::error::Error;

    /// A `--meth` original reference needs each of `.amb`/`.ann`/`.pac`, and
    /// not the FM-index `.bwt.2bit.64` that `bwa-mem3 index --meth-only` omits:
    /// with exactly one of the three absent, the error names exactly that file.
    #[rstest::rstest]
    #[case::amb(".amb")]
    #[case::ann(".ann")]
    #[case::pac(".pac")]
    fn orig_ref_check_requires_exactly_amb_ann_pac(#[case] absent: &str) {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("ref.fa");
        for suffix in ORIG_REF_SUFFIXES.iter().filter(|s| **s != absent) {
            std::fs::write(dir.path().join(format!("ref.fa{suffix}")), b"x").unwrap();
        }
        let Err(Error::IndexLoad { msg, .. }) =
            check_files(&prefix, ORIG_REF_SUFFIXES, "original-reference")
        else {
            panic!("expected an IndexLoad error");
        };
        assert_eq!(
            msg,
            format!(
                "missing or unreadable original-reference file(s): {}{absent}",
                prefix.display()
            )
        );

        std::fs::write(dir.path().join(format!("ref.fa{absent}")), b"x").unwrap();
        check_files(&prefix, ORIG_REF_SUFFIXES, "original-reference").unwrap();
    }

    /// A plain load of a `--meth-only` original prefix (no FM-index, but a
    /// `.meth` seed index beside it) points at `load_meth`; any other missing
    /// FM-index gets the plain message.
    #[rstest::rstest]
    #[case::meth_only_original(true)]
    #[case::no_seed_index(false)]
    fn plain_check_hints_load_meth_for_a_meth_only_original(#[case] seed_present: bool) {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("ref.fa");
        for suffix in ORIG_REF_SUFFIXES {
            std::fs::write(dir.path().join(format!("ref.fa{suffix}")), b"x").unwrap();
        }
        if seed_present {
            std::fs::write(dir.path().join("ref.fa.meth.bwt.2bit.64"), b"x").unwrap();
        }
        let Err(Error::IndexLoad { msg, .. }) = check_index_files(&prefix) else {
            panic!("expected an IndexLoad error");
        };
        assert!(msg.contains("ref.fa.bwt.2bit.64"), "{msg}");
        assert_eq!(msg.contains("BwaIndex::load_meth"), seed_present, "{msg}");
    }

    #[test]
    fn bam_header_fields_hold_the_format_limits() {
        let max = i32::MAX as usize;
        assert_eq!(
            bam_header_field(max, "l_text", true).unwrap(),
            i32::MAX as u32
        );
        let err = bam_header_field(max + 1, "l_text", true)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("bam_header: l_text (2147483648) exceeds the BAM limit of 2147483647"),
            "{err}"
        );
        // l_name need only fit 32 bits.
        assert_eq!(
            bam_header_field(max + 1, "l_name", false).unwrap(),
            1u32 << 31
        );
        #[cfg(target_pointer_width = "64")]
        assert!(bam_header_field(u32::MAX as usize + 1, "l_name", false).is_err());
    }
}
