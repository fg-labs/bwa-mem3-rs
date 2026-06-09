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

/// Return [`Error::IndexLoad`] naming every [`REQUIRED_INDEX_SUFFIXES`] file
/// that is not a readable regular file for `prefix`. `Ok(())` means all are
/// present and openable.
///
/// A bare existence check is not enough: a directory, FIFO, or unreadable
/// path satisfies [`Path::exists`] yet still trips the vendored loader's
/// `exit(1)` on the subsequent `fopen`/`fread`. We mirror what the loader
/// actually needs — a regular file it can open for reading. The `is_file`
/// guard short-circuits before [`File::open`], so a FIFO (whose `O_RDONLY`
/// open would block) is rejected without ever being opened.
fn check_index_files(prefix: &Path) -> Result<()> {
    let base = prefix.as_os_str();
    let unusable: Vec<String> = REQUIRED_INDEX_SUFFIXES
        .iter()
        .map(|suffix| {
            let mut file = base.to_owned();
            file.push(suffix);
            file
        })
        .filter(|file| {
            let path = Path::new(file);
            !path.is_file() || std::fs::File::open(path).is_err()
        })
        .map(|file| Path::new(&file).display().to_string())
        .collect();
    if unusable.is_empty() {
        Ok(())
    } else {
        Err(Error::IndexLoad {
            path: prefix.to_owned(),
            msg: format!(
                "missing or unreadable index file(s): {}",
                unusable.join(", ")
            ),
        })
    }
}

/// Handle to a loaded bwa-mem3 reference index.
///
/// The index is immutable after loading. Multiple threads may share a
/// `BwaIndex` via `&BwaIndex` or `Arc<BwaIndex>` for concurrent alignment.
pub struct BwaIndex {
    handle: *mut bwa_mem3_sys::BwaIndex,
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

        let s = path
            .to_str()
            .ok_or_else(|| Error::InvalidInput("prefix must be valid UTF-8".into()))?;
        let c = CString::new(s)?;
        let handle = unsafe { bwa_mem3_sys::bwa_shim_idx_load(c.as_ptr()) };
        if handle.is_null() {
            return Err(Error::IndexLoad {
                path: path.to_owned(),
                msg: shim_err("idx load").to_string(),
            });
        }
        Ok(BwaIndex { handle })
    }

    #[must_use]
    pub fn n_contigs(&self) -> usize {
        unsafe { bwa_mem3_sys::bwa_shim_idx_n_contigs(self.handle) }
    }

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

    #[must_use]
    pub fn contig_len(&self, i: usize) -> i64 {
        unsafe { bwa_mem3_sys::bwa_shim_idx_contig_len(self.handle, i) }
    }

    pub fn contigs(&self) -> impl Iterator<Item = (&str, i64)> + '_ {
        (0..self.n_contigs()).map(move |i| (self.contig_name(i), self.contig_len(i)))
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
