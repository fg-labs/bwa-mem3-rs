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

/// Longest read name a BAM record can carry: `l_read_name` is a `u8` that
/// counts the trailing NUL, so the name itself tops out at 254 bytes. The shim
/// aborts the process on a longer one when it builds the record, so the name is
/// rejected up front instead.
const MAX_READ_NAME_LEN: usize = 254;

/// Reject a read name too long for a BAM record (see [`MAX_READ_NAME_LEN`]).
fn validate_read_name(name: &[u8]) -> Result<()> {
    if name.len() > MAX_READ_NAME_LEN {
        return Err(Error::InvalidInput(format!(
            "read name is {} bytes; BAM allows at most {MAX_READ_NAME_LEN}",
            name.len()
        )));
    }
    Ok(())
}

impl<'a> ReadPair<'a> {
    pub fn validate(&self) -> Result<()> {
        if self.seq_r1.is_empty() || self.seq_r2.is_empty() {
            return Err(Error::InvalidInput("empty sequence".into()));
        }
        validate_read_name(self.name_r1)?;
        validate_read_name(self.name_r2)?;
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
pub(crate) fn check_meth_consistency(idx: &BwaIndex, opts: &MemOpts) -> Result<()> {
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

/// The C view of one pair, borrowing `p`'s bytes (no copy).
pub(crate) fn c_pair(p: &ReadPair<'_>) -> bwa_mem3_sys::BwaReadPair {
    bwa_mem3_sys::BwaReadPair {
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
    }
}

fn to_c_pairs(pairs: &[ReadPair<'_>]) -> Vec<bwa_mem3_sys::BwaReadPair> {
    pairs.iter().map(c_pair).collect()
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

impl MemPeStat {
    /// Insert-size model for one `-K` cohort: `mem_pestat` once over the
    /// paired reads of every batch in `regs` (singles are ignored), exactly as
    /// `bwa-mem3 mem` computes it per chunk (`bwamem.cpp:3030-3050`). Feed the
    /// result to every [`pair_emit`] of the same cohort.
    pub fn infer_cohort(idx: &BwaIndex, opts: &MemOpts, regs: &[AlnRegs]) -> Result<Self> {
        let mut out = Self::zero()?;
        let raw: Vec<*const bwa_mem3_sys::BwaRegs> = regs.iter().map(|r| r.raw()).collect();
        // SAFETY: `idx.raw()`/`opts.as_ptr()` are valid for the call; `raw`
        // holds `raw.len()` live `BwaRegs` pointers (each owned elsewhere by
        // its `AlnRegs`, still alive for the duration of this call since we
        // only borrowed them); `out.as_mut_ptr()` is a live `mem_pestat_t[4]`
        // owned by `out`.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_pestat_cohort(
                idx.raw(),
                opts.as_ptr(),
                raw.as_ptr(),
                raw.len(),
                out.as_mut_ptr(),
            )
        };
        if rc != 0 {
            return Err(shim_err("pestat_cohort"));
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// Three-phase API: seed_extend -> MemPeStat::infer_cohort -> pair_emit
// ---------------------------------------------------------------------------

/// A single-end read. Borrows its bytes; zero copies across the FFI boundary.
#[derive(Debug, Clone, Copy)]
pub struct SingleRead<'a> {
    pub name: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl SingleRead<'_> {
    pub fn validate(&self) -> Result<()> {
        if self.seq.is_empty() {
            return Err(Error::InvalidInput("empty sequence".into()));
        }
        validate_read_name(self.name)?;
        if let Some(q) = self.qual {
            if q.len() != self.seq.len() {
                return Err(Error::InvalidInput(format!(
                    "qual/seq length mismatch: {} vs {}",
                    q.len(),
                    self.seq.len()
                )));
            }
        }
        Ok(())
    }
}

/// One work item for [`seed_extend`]: pairs followed by singles. Either slice
/// may be empty. bwa-mem3 `-p` processes a cohort's single-end group and its
/// paired group separately (`fastmap.cpp:908-958`); a `ReadBatch` is one
/// slice of one of those groups, or of both.
#[derive(Debug, Clone, Copy)]
pub struct ReadBatch<'a> {
    pub pairs: &'a [ReadPair<'a>],
    pub singles: &'a [SingleRead<'a>],
}

impl ReadBatch<'_> {
    fn validate(&self) -> Result<()> {
        validate_all(self.pairs)?;
        for (i, s) in self.singles.iter().enumerate() {
            s.validate().map_err(|e| match e {
                Error::InvalidInput(m) => Error::InvalidInput(format!("single {i}: {m}")),
                other => other,
            })?;
        }
        Ok(())
    }
}

/// The C view of one single-end read, borrowing `s`'s bytes (no copy).
pub(crate) fn c_single(s: &SingleRead<'_>) -> bwa_mem3_sys::BwaSingleRead {
    bwa_mem3_sys::BwaSingleRead {
        name: s.name.as_ptr().cast::<std::ffi::c_char>(),
        name_len: s.name.len(),
        seq: s.seq.as_ptr(),
        seq_len: s.seq.len(),
        qual: s.qual.map_or(std::ptr::null(), <[u8]>::as_ptr),
    }
}

fn to_c_singles(singles: &[SingleRead<'_>]) -> Vec<bwa_mem3_sys::BwaSingleRead> {
    singles.iter().map(c_single).collect()
}

/// Per-thread reusable scratch (~24 MB of SIMD buffers after first use). Create
/// one per worker thread and pass it to every [`seed_extend`] / [`pair_emit`]
/// on that thread. `Send` so a pool can migrate it; not `Sync`.
pub struct AlignScratch {
    pub(crate) handle: *mut bwa_mem3_sys::BwaScratch,
}

/// Reads per bwa-mem3 kernel batch (`BATCH_SIZE` in the vendored `macro.h`:
/// 1024 on aarch64, 512 elsewhere). The shim runs each seed, extension and
/// mate-rescue kernel call on one such batch, as the CLI's workers do, so a
/// caller splitting a `-K` chunk into sub-batches should size them in
/// multiples of it. Not to be confused with a [`ReadBatch`] or `-K`.
#[must_use]
pub fn kernel_batch_size() -> usize {
    // SAFETY: a pure query with no preconditions.
    unsafe { bwa_mem3_sys::bwa_shim_kernel_batch_size() }
}

impl AlignScratch {
    /// The kernel thread slot this scratch runs bwa-mem3's kernels in. Test
    /// hook: up to 256 scratches alive together get distinct slots, spread so
    /// that concurrent workers do not write the same profiling-counter line;
    /// beyond 256 they share.
    #[doc(hidden)]
    #[must_use]
    pub fn tid_slot(&self) -> usize {
        // SAFETY: `handle` is a live scratch owned by `self`; the query reads
        // one field and has no other preconditions.
        let tid = unsafe { bwa_mem3_sys::bwa_shim_scratch_tid(self.handle) };
        usize::try_from(tid).expect("a live scratch has a slot")
    }

    pub fn new() -> Result<Self> {
        // SAFETY: `bwa_shim_scratch_new` is an opaque C allocator with no
        // preconditions; it returns null on failure, which we check below
        // before ever treating the result as a valid handle.
        let handle = unsafe { bwa_mem3_sys::bwa_shim_scratch_new() };
        if handle.is_null() {
            return Err(shim_err("scratch_new"));
        }
        Ok(Self { handle })
    }
}

impl Drop for AlignScratch {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            // SAFETY: `self.handle` is a valid non-null handle produced by
            // `bwa_shim_scratch_new`, owned exclusively by `self`, and this
            // is the only place it is freed -- `Drop` runs at most once.
            unsafe { bwa_mem3_sys::bwa_shim_scratch_free(self.handle) };
        }
    }
}

// SAFETY: the scratch is a bare opaque C heap allocation the shim only ever
// touches through the `&mut AlignScratch` FFI calls below (`seed_extend`,
// `pair_emit`); it holds no pointers into caller memory that must outlive a
// single call. Moving the whole allocation to another thread (Send) is sound;
// concurrent access from two threads (Sync) is not, since the C-side buffers
// are mutated in place with no internal locking -- callers must give each
// worker thread its own `AlignScratch`.
unsafe impl Send for AlignScratch {}

/// Phase-1 output for one [`ReadBatch`]: the copied reads plus every read's
/// candidate alignments. Consumed by [`pair_emit`]; contributes to
/// [`MemPeStat::infer_cohort`] by reference.
pub struct AlnRegs {
    handle: *mut bwa_mem3_sys::BwaRegs,
    n_pairs: usize,
    n_singles: usize,
}

impl AlnRegs {
    #[must_use]
    pub fn n_pairs(&self) -> usize {
        self.n_pairs
    }
    #[must_use]
    pub fn n_singles(&self) -> usize {
        self.n_singles
    }
    /// Bytes this handle holds on the C heap (read copies + alnreg arrays).
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        // SAFETY: `self.handle` is a valid non-null handle produced by
        // `bwa_shim_seed_extend`, owned by `self` for at least the duration
        // of this read-only call.
        unsafe { bwa_mem3_sys::bwa_shim_regs_heap_bytes(self.handle) }
    }

    pub(crate) fn raw(&self) -> *const bwa_mem3_sys::BwaRegs {
        self.handle
    }
}

impl Drop for AlnRegs {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            // SAFETY: `self.handle` is a valid non-null handle produced by
            // `bwa_shim_seed_extend`, owned exclusively by `self`, and this
            // is the only place it is freed -- `Drop` runs at most once.
            unsafe { bwa_mem3_sys::bwa_shim_regs_free(self.handle) };
        }
    }
}

// SAFETY: owns all of its memory (copied read bytes + alnreg arrays) with no
// aliasing back into the shared `BwaIndex` or caller buffers, so handing the
// whole allocation to another thread (Send) is sound. Not Sync: the C struct
// has no internal synchronization, so two threads must not share a `&AlnRegs`
// while one of them holds it long enough to matter (in practice access is
// always through an owning `&mut`/by-value handoff, e.g. into `pair_emit`).
unsafe impl Send for AlnRegs {}

/// Phase 1: seed + single-end-extend every read of `batch` (bwa-mem3's fused
/// `worker_bwt_aln`). Per-read independent, so a cohort may be split into any
/// number of batches on any number of threads without changing output.
pub fn seed_extend(
    idx: &BwaIndex,
    opts: &MemOpts,
    scratch: &mut AlignScratch,
    batch: &ReadBatch<'_>,
) -> Result<AlnRegs> {
    batch.validate()?;
    check_meth_consistency(idx, opts)?;
    let c_pairs = to_c_pairs(batch.pairs);
    let c_singles = to_c_singles(batch.singles);
    let c_batch = bwa_mem3_sys::BwaReadBatch {
        pairs: c_pairs.as_ptr(),
        n_pairs: c_pairs.len(),
        singles: c_singles.as_ptr(),
        n_singles: c_singles.len(),
    };
    // SAFETY: `idx.raw()`/`opts.as_ptr()` are valid for the call; `scratch.handle`
    // is a live `BwaScratch` owned by `scratch`; `c_batch` points at `c_pairs`/
    // `c_singles`, which outlive this call (they are dropped after it returns).
    let handle = unsafe {
        bwa_mem3_sys::bwa_shim_seed_extend(idx.raw(), opts.as_ptr(), scratch.handle, &c_batch)
    };
    if handle.is_null() {
        return Err(shim_err("seed_extend"));
    }
    Ok(AlnRegs {
        handle,
        n_pairs: batch.pairs.len(),
        n_singles: batch.singles.len(),
    })
}

/// Global read ordinals for one batch, reproducing bwa-mem3's `worker_sam`
/// ids (`bwamem.cpp:2795-2884`): pair `i` gets `first_pair_id + i` (the CLI's
/// `(n_processed >> 1) + pos`), single `i` gets `first_single_id + i` (the
/// CLI's `n_processed + i`). For a `-p` cohort whose first read is global read
/// `base` with `n_se` singles classified before the pairs
/// (`fastmap.cpp:924-944`): `first_single_id = base + singles_before_this_batch`,
/// `first_pair_id = ((base + n_se) >> 1) + pairs_before_this_batch`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdBases {
    pub first_single_id: u64,
    pub first_pair_id: u64,
}

/// Which input read a record came from (index into the batch's `pairs` or
/// `singles`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordOrigin {
    Pair(usize),
    Single(usize),
}

/// Receives packed BAM record **bodies** (no `u32 block_size` prefix) as they
/// are emitted, in input order: pairs (R1 side then R2 side, primary then
/// supplementary), then singles.
///
/// # Contract
///
/// `body` is a borrow from the shim's internal buffers, valid **only** for the
/// duration of this call -- copy it (as [`RecordVec`] does) if it must outlive
/// the call. `emit` is free to panic; [`pair_emit`] catches the unwind at the
/// FFI boundary and re-raises it once the underlying C call has returned, so a
/// panicking sink is sound (no unwind ever crosses the `extern "C"` frame) but
/// still surfaces to the `pair_emit` caller exactly as an ordinary panic would.
pub trait RecordSink {
    fn emit(&mut self, origin: RecordOrigin, body: &[u8]);
}

/// The simplest sink: owned copies in emission order.
#[derive(Default, Debug)]
pub struct RecordVec {
    pub records: Vec<(RecordOrigin, Vec<u8>)>,
}

impl RecordSink for RecordVec {
    fn emit(&mut self, origin: RecordOrigin, body: &[u8]) {
        self.records.push((origin, body.to_vec()));
    }
}

/// Trampoline state: the caller's sink plus, if `emit` unwound, the caught
/// panic payload to resume after the C call stack has unwound normally.
/// Shared by the packed ([`RecordSink`]) and structured
/// ([`AlignedFieldsSink`](crate::AlignedFieldsSink)) trampolines.
pub(crate) struct SinkCtx<'a, S: ?Sized> {
    pub(crate) sink: &'a mut S,
    pub(crate) panic: Option<Box<dyn std::any::Any + Send + 'static>>,
}

impl<S: ?Sized> SinkCtx<'_, S> {
    /// Run `f` against the sink unless an earlier call already panicked,
    /// catching (and stashing) any panic so it never unwinds across the
    /// calling `extern "C"` frame. The C call stack cannot unwind, so after a
    /// panic the trampoline keeps returning normally for the remaining records
    /// and the caller resumes the panic once the shim call itself returns.
    /// `AssertUnwindSafe` is fine because after a caught panic the sink is
    /// never touched again.
    pub(crate) fn emit_catching(&mut self, f: impl FnOnce(&mut S)) {
        if self.panic.is_some() {
            return;
        }
        let sink = &mut *self.sink;
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || f(sink)));
        if let Err(payload) = outcome {
            self.panic = Some(payload);
        }
    }
}

/// Map the shim's `origin_kind`/`origin_idx` pair to a [`RecordOrigin`].
pub(crate) fn record_origin(origin_kind: u32, origin_idx: usize) -> RecordOrigin {
    if origin_kind == bwa_mem3_sys::BWA_ORIGIN_SINGLE {
        RecordOrigin::Single(origin_idx)
    } else {
        RecordOrigin::Pair(origin_idx)
    }
}

/// C -> Rust trampoline for [`RecordSink`], passed to `bwa_shim_pair_emit` as
/// its `BwaRecordSinkFn`.
///
/// # Safety
///
/// The caller (`pair_emit`) must pass `ctx` as `&mut SinkCtx` cast to
/// `*mut c_void`, kept alive for the whole `bwa_shim_pair_emit` call, and must
/// not alias it elsewhere during that call. `body` must point at `body_len`
/// initialized, readable bytes that stay valid for the duration of this call
/// only (the shim's own contract for `BwaRecordSinkFn`).
pub(crate) unsafe extern "C" fn sink_trampoline(
    ctx: *mut std::ffi::c_void,
    origin_kind: u32,
    origin_idx: usize,
    body: *const u8,
    body_len: usize,
) {
    // SAFETY: see the function's `# Safety` section -- `ctx` is a live
    // `&mut SinkCtx` for the duration of this call.
    let ctx = unsafe { &mut *ctx.cast::<SinkCtx<'_, dyn RecordSink>>() };
    // SAFETY: see the function's `# Safety` section -- `body`/`body_len`
    // describe a valid slice for the duration of this call.
    let bytes = unsafe { std::slice::from_raw_parts(body, body_len) };
    let origin = record_origin(origin_kind, origin_idx);
    // Panics must not unwind across this `extern "C"` frame: with
    // `panic = "unwind"` that is UB, and with `panic = "abort"` it would take
    // down the process before the caller ever sees `pair_emit`'s `Result`.
    ctx.emit_catching(|sink| sink.emit(origin, bytes));
}

/// Phase 3: pairing + mate rescue + primary marking + emission (bwa-mem3's
/// `worker_sam`) for one batch. `pestat` is the cohort model from
/// [`MemPeStat::infer_cohort`]; it may be `None` only for a batch without
/// pairs. Consumes `regs`.
pub fn pair_emit(
    idx: &BwaIndex,
    opts: &MemOpts,
    scratch: &mut AlignScratch,
    regs: AlnRegs,
    pestat: Option<&MemPeStat>,
    ids: IdBases,
    sink: &mut dyn RecordSink,
) -> Result<()> {
    if regs.n_pairs > 0 && pestat.is_none() {
        return Err(Error::InvalidInput(
            "pair_emit: a batch with pairs requires the cohort pestat".into(),
        ));
    }
    let handle = regs.handle;
    // Ownership of the C allocation moves to the shim, which frees it on
    // every return path (see `bwa_shim_pair_emit`'s contract); forget the
    // Rust wrapper so its `Drop` doesn't double-free.
    std::mem::forget(regs);

    let mut sink_ctx = SinkCtx { sink, panic: None };
    let ctx_ptr = std::ptr::addr_of_mut!(sink_ctx).cast::<std::ffi::c_void>();

    // SAFETY: `idx`/`opts`/`scratch.handle` are valid for the call; `handle`
    // is the `BwaRegs` we just took ownership of and forgot above, consumed
    // by the shim on every path; `sink_trampoline`/`ctx_ptr` satisfy
    // `sink_trampoline`'s safety contract (`ctx_ptr` outlives this call and is
    // not aliased elsewhere during it).
    let rc = unsafe {
        bwa_mem3_sys::bwa_shim_pair_emit(
            idx.raw(),
            opts.as_ptr(),
            scratch.handle,
            handle,
            pestat.map_or(std::ptr::null(), MemPeStat::as_ptr),
            bwa_mem3_sys::BwaIdBases {
                first_single_id: ids.first_single_id,
                first_pair_id: ids.first_pair_id,
            },
            Some(sink_trampoline),
            ctx_ptr,
        )
    };

    if let Some(payload) = sink_ctx.panic {
        // The C call stack has already unwound normally (sink_trampoline
        // never let a panic cross it); re-raise it now so a panicking sink
        // still looks like an ordinary panic to `pair_emit`'s caller.
        std::panic::resume_unwind(payload);
    }
    if rc != 0 {
        return Err(shim_err("pair_emit"));
    }
    Ok(())
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

    #[test]
    fn single_read_validation() {
        let ok = SingleRead {
            name: b"s",
            seq: b"ACGT",
            qual: Some(b"IIII"),
        };
        assert!(ok.validate().is_ok());
        let bad = SingleRead {
            name: b"s",
            seq: b"ACGT",
            qual: Some(b"II"),
        };
        assert!(bad.validate().is_err());
        let empty = SingleRead {
            name: b"s",
            seq: b"",
            qual: None,
        };
        assert!(empty.validate().is_err());
    }

    #[test]
    fn record_vec_collects_in_order() {
        let mut v = RecordVec::default();
        v.emit(RecordOrigin::Pair(0), b"ab");
        v.emit(RecordOrigin::Single(3), b"c");
        assert_eq!(
            v.records,
            vec![
                (RecordOrigin::Pair(0), b"ab".to_vec()),
                (RecordOrigin::Single(3), b"c".to_vec())
            ]
        );
    }

    #[test]
    fn scratch_is_send_and_regs_is_send() {
        fn assert_send<T: Send>() {}
        assert_send::<AlignScratch>();
        assert_send::<AlnRegs>();
    }
}
