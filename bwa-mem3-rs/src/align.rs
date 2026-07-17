//! Read-pair input, seed/alignment batch handles, and the public align functions.

use crate::error::{shim_err, Error, Result};
use crate::index::BwaIndex;
use crate::opts::{MemOpts, MemPeStat};

/// A single read pair. Borrows its sequence and quality bytes; zero copies
/// across the FFI boundary.
#[derive(Debug, Clone, Copy)]
pub struct ReadPair<'a> {
    pub name_r1: &'a [u8],
    pub seq_r1: &'a [u8],
    pub qual_r1: Option<&'a [u8]>,
    pub name_r2: &'a [u8],
    pub seq_r2: &'a [u8],
    pub qual_r2: Option<&'a [u8]>,
}

impl<'a> ReadPair<'a> {
    pub fn validate(&self) -> Result<()> {
        if self.seq_r1.is_empty() || self.seq_r2.is_empty() {
            return Err(Error::InvalidInput("empty sequence".into()));
        }
        if let Some(q) = self.qual_r1 {
            if q.len() != self.seq_r1.len() {
                return Err(Error::InvalidInput(format!(
                    "R1 qual/seq length mismatch: {} vs {}",
                    q.len(),
                    self.seq_r1.len()
                )));
            }
        }
        if let Some(q) = self.qual_r2 {
            if q.len() != self.seq_r2.len() {
                return Err(Error::InvalidInput(format!(
                    "R2 qual/seq length mismatch: {} vs {}",
                    q.len(),
                    self.seq_r2.len()
                )));
            }
        }
        Ok(())
    }
}

fn validate_all(pairs: &[ReadPair<'_>]) -> Result<()> {
    for (i, p) in pairs.iter().enumerate() {
        p.validate().map_err(|e| match e {
            Error::InvalidInput(m) => Error::InvalidInput(format!("pair {i}: {m}")),
            other => other,
        })?;
    }
    Ok(())
}

/// Reject a bisulfite mode/index mismatch before aligning. `MemOpts::set_meth`
/// only sets a flag; it must be paired with a meth dual index
/// ([`BwaIndex::load_meth`]). Aligning meth-mode reads against a plain index
/// (or vice versa) would silently produce wrong output — projected reads
/// against an unconverted index, or spurious `XR`/`XG`/`XM` tags — so we fail
/// fast instead.
fn check_meth_consistency(idx: &BwaIndex, opts: &MemOpts) -> Result<()> {
    match (opts.meth(), idx.is_meth()) {
        (true, false) => Err(Error::InvalidInput(
            "opts have --meth enabled but the index is not a meth dual index \
             (load it with BwaIndex::load_meth)"
                .into(),
        )),
        (false, true) => Err(Error::InvalidInput(
            "the index is a meth dual index but opts do not have --meth enabled \
             (call MemOpts::set_meth(true))"
                .into(),
        )),
        _ => Ok(()),
    }
}

fn to_c_pairs(pairs: &[ReadPair<'_>]) -> Vec<bwa_mem3_sys::BwaReadPair> {
    pairs
        .iter()
        .map(|p| bwa_mem3_sys::BwaReadPair {
            r1_name: p.name_r1.as_ptr().cast::<std::ffi::c_char>(),
            r1_name_len: p.name_r1.len(),
            r1_seq: p.seq_r1.as_ptr(),
            r1_seq_len: p.seq_r1.len(),
            r1_qual: p.qual_r1.map_or(std::ptr::null(), <[u8]>::as_ptr),
            r2_name: p.name_r2.as_ptr().cast::<std::ffi::c_char>(),
            r2_name_len: p.name_r2.len(),
            r2_seq: p.seq_r2.as_ptr(),
            r2_seq_len: p.seq_r2.len(),
            r2_qual: p.qual_r2.map_or(std::ptr::null(), <[u8]>::as_ptr),
        })
        .collect()
}

/// Phase-1 opaque handle: the seeds (chains) for a batch of read pairs.
///
/// Send across threads; a single thread consumes it via [`extend_batch`].
pub struct Seeds {
    handle: *mut bwa_mem3_sys::BwaSeeds,
}

impl Drop for Seeds {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { bwa_mem3_sys::bwa_shim_seeds_free(self.handle) };
        }
    }
}

// SAFETY: Seeds owns all its memory with no aliasing back into the shared
// index or options. NOT Sync — consumed by extend_batch.
unsafe impl Send for Seeds {}

/// A single packed BAM record, per the BAM spec:
/// `[u32 le block_size][block_size bytes of record data]`.
#[derive(Debug, Clone, Copy)]
pub struct Record<'a> {
    /// Which input pair (index into the batch) produced this record.
    pub pair_idx: usize,
    pub bytes: &'a [u8],
}

/// Phase-2 output: packed BAM records from one batch of aligned pairs.
pub struct AlignmentBatch {
    handle: *mut bwa_mem3_sys::BwaBatch,
}

impl AlignmentBatch {
    pub fn len(&self) -> usize {
        unsafe { bwa_mem3_sys::bwa_shim_batch_n_records(self.handle) }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn get(&self, i: usize) -> Record<'_> {
        let pair_idx = unsafe { bwa_mem3_sys::bwa_shim_batch_pair_idx(self.handle, i) };
        let ptr = unsafe { bwa_mem3_sys::bwa_shim_batch_record_ptr(self.handle, i) };
        let len = unsafe { bwa_mem3_sys::bwa_shim_batch_record_len(self.handle, i) };
        let bytes = unsafe { std::slice::from_raw_parts(ptr, len) };
        Record { pair_idx, bytes }
    }

    pub fn iter(&self) -> impl Iterator<Item = Record<'_>> + '_ {
        (0..self.len()).map(move |i| self.get(i))
    }

    /// Group records by source pair index. Allocates one `Vec` per pair.
    #[must_use]
    pub fn by_pair(&self) -> Vec<Vec<Record<'_>>> {
        let n = self.iter().map(|r| r.pair_idx).max().map_or(0, |m| m + 1);
        let mut out: Vec<Vec<Record<'_>>> = (0..n).map(|_| Vec::new()).collect();
        for r in self.iter() {
            out[r.pair_idx].push(r);
        }
        out
    }
}

impl Drop for AlignmentBatch {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { bwa_mem3_sys::bwa_shim_batch_free(self.handle) };
        }
    }
}

unsafe impl Send for AlignmentBatch {}

/// Phase 1: seed the batch.
pub fn seed_batch(idx: &BwaIndex, opts: &MemOpts, pairs: &[ReadPair<'_>]) -> Result<Seeds> {
    validate_all(pairs)?;
    check_meth_consistency(idx, opts)?;
    let c_pairs = to_c_pairs(pairs);
    let handle = unsafe {
        bwa_mem3_sys::bwa_shim_seed_batch(idx.raw(), opts.as_ptr(), c_pairs.as_ptr(), c_pairs.len())
    };
    if handle.is_null() {
        return Err(shim_err("seed_batch"));
    }
    Ok(Seeds { handle })
}

/// Phase 2: extend seeds to full alignments. Consumes `seeds`.
pub fn extend_batch(
    idx: &BwaIndex,
    opts: &MemOpts,
    seeds: Seeds,
    pairs: &[ReadPair<'_>],
    pestat_in: Option<&MemPeStat>,
) -> Result<(AlignmentBatch, MemPeStat)> {
    validate_all(pairs)?;
    let c_pairs = to_c_pairs(pairs);
    let mut pestat_out = MemPeStat::zero()?;

    // Transfer seeds ownership to the shim; forget the wrapper so its Drop
    // doesn't double-free.
    let s_handle = seeds.handle;
    std::mem::forget(seeds);

    let batch = unsafe {
        bwa_mem3_sys::bwa_shim_extend_batch(
            idx.raw(),
            opts.as_ptr(),
            s_handle,
            c_pairs.as_ptr(),
            c_pairs.len(),
            pestat_in.map_or(std::ptr::null(), MemPeStat::as_ptr),
            pestat_out.as_mut_ptr(),
        )
    };
    if batch.is_null() {
        return Err(shim_err("extend_batch"));
    }
    Ok((AlignmentBatch { handle: batch }, pestat_out))
}

/// Convenience: seed + extend in one call.
pub fn align_batch(
    idx: &BwaIndex,
    opts: &MemOpts,
    pairs: &[ReadPair<'_>],
    pestat_in: Option<&MemPeStat>,
) -> Result<(AlignmentBatch, MemPeStat)> {
    validate_all(pairs)?;
    check_meth_consistency(idx, opts)?;
    let c_pairs = to_c_pairs(pairs);
    let mut pestat_out = MemPeStat::zero()?;

    let batch = unsafe {
        bwa_mem3_sys::bwa_shim_align_batch(
            idx.raw(),
            opts.as_ptr(),
            c_pairs.as_ptr(),
            c_pairs.len(),
            pestat_in.map_or(std::ptr::null(), MemPeStat::as_ptr),
            pestat_out.as_mut_ptr(),
        )
    };
    if batch.is_null() {
        return Err(shim_err("align_batch"));
    }
    Ok((AlignmentBatch { handle: batch }, pestat_out))
}

/// Estimate the insert-size model from a pilot batch. Discards alignments.
pub fn estimate_pestat(
    idx: &BwaIndex,
    opts: &MemOpts,
    pairs: &[ReadPair<'_>],
) -> Result<MemPeStat> {
    validate_all(pairs)?;
    check_meth_consistency(idx, opts)?;
    let c_pairs = to_c_pairs(pairs);
    let mut pestat_out = MemPeStat::zero()?;
    let rc = unsafe {
        bwa_mem3_sys::bwa_shim_estimate_pestat(
            idx.raw(),
            opts.as_ptr(),
            c_pairs.as_ptr(),
            c_pairs.len(),
            pestat_out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(shim_err("estimate_pestat"));
    }
    Ok(pestat_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_seq_rejected() {
        let p = ReadPair {
            name_r1: b"r1",
            seq_r1: b"",
            qual_r1: None,
            name_r2: b"r2",
            seq_r2: b"ACGT",
            qual_r2: None,
        };
        assert!(p.validate().is_err());
    }

    #[test]
    fn qual_length_mismatch_rejected() {
        let p = ReadPair {
            name_r1: b"r1",
            seq_r1: b"ACGT",
            qual_r1: Some(b"!!"),
            name_r2: b"r2",
            seq_r2: b"ACGT",
            qual_r2: None,
        };
        assert!(p.validate().is_err());
    }
}
