//! Resident-cohort alignment: one `-K` cohort's decoded reads and alignment
//! regions stay resident on the C heap from seed/extend through pair/emit, and
//! every sub-batch works on its own reserved range of them. See
//! [`ResidentCohort`].
//!
//! Also home to the structured-fields output ([`AlignedFields`],
//! [`AlignedFieldsSink`]) that [`ResidentCohort::pair_emit_fields`] reports
//! instead of packed BAM bodies.

use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, PoisonError, RwLock, RwLockReadGuard, RwLockWriteGuard};

use crate::align::{
    c_pair, c_single, check_meth_consistency, record_origin, sink_trampoline, AlignScratch,
    IdBases, ReadPair, RecordOrigin, RecordSink, SingleRead, SinkCtx,
};
use crate::error::{shim_err, shim_last_message, Error, Result};
use crate::index::BwaIndex;
use crate::opts::{MemOpts, MemPeStat};

// ---------------------------------------------------------------------------
// Structured-fields output
// ---------------------------------------------------------------------------

/// The structured fields of one emitted record: everything the packed-BAM
/// path serializes, as values.
/// Building a BAM record from these plus the input read reproduces the packed
/// body byte-for-byte; bwa-mem3 stays the sole authority on every value, the
/// consumer only decides where the bytes are assembled.
///
/// Every borrow is valid only for the duration of the
/// [`AlignedFieldsSink::emit`] call it is passed to, like
/// [`RecordSink::emit`]'s `body`, and points into memory that call owns
/// (never a process-wide or thread-local buffer), so nothing the sink does
/// can change a field while it reads it. Strings are the tag payloads without the
/// trailing NUL; `None` means the packed record omits the tag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignedFields<'a> {
    /// Reference id as serialized. A half-mapped pair's unmapped read carries
    /// its mate's placement (upstream's mate-coordinate copy).
    pub tid: i32,
    /// 0-based position as serialized (same placement rule as `tid`).
    pub pos: i32,
    /// RNEXT (`-1` without a mate).
    pub next_tid: i32,
    /// PNEXT (`-1` without a mate).
    pub next_pos: i32,
    pub tlen: i32,
    /// The final 16-bit FLAG. Its 0x10 bit is the emitted strand, see
    /// [`is_rev`](Self::is_rev).
    pub flag: u16,
    /// BAM `bin` as htslib's `bam_set1` computes it: over the reference span
    /// for a mapped record, over one base for an unmapped one placed at its
    /// mate, and 4680 for an unplaced one.
    pub bin: u16,
    pub mapq: u8,
    /// The emitted SEQ/QUAL is bases `[query_start, query_end)` of the input
    /// read in as-sequenced orientation, reverse-complemented (QUAL reversed)
    /// when [`is_rev`](Self::is_rev). The window is shorter than the read on a
    /// hard-clipped supplementary and empty on a secondary (`-a`) record, which
    /// carries no SEQ/QUAL. Base encoding follows the upstream writers: without
    /// `--meth` every base passes through bwa's 2-bit alphabet (case folded,
    /// anything but A/C/G/T becomes `N`); under `--meth` the original bases are
    /// emitted, forward ones upper-cased with IUPAC codes kept, reverse ones
    /// complemented with anything but A/C/G/T becoming `N`.
    pub query_start: usize,
    pub query_end: usize,
    /// CIGAR in BAM opcodes, clip rewrite (hard vs soft) already applied.
    pub cigar: &'a [u32],
    pub nm: Option<i32>,
    pub md: Option<&'a [u8]>,
    /// MC:Z, rendered with the emitting record's clip mode (upstream quirk).
    pub mc: Option<&'a [u8]>,
    pub mq: Option<i32>,
    /// AS:i.
    pub score: Option<i32>,
    /// XS:i.
    pub sub: Option<i32>,
    pub rg: Option<&'a [u8]>,
    pub sa: Option<&'a [u8]>,
    /// pa:f.
    pub pa: Option<f32>,
    pub xa: Option<&'a [u8]>,
    pub hn: Option<i32>,
    /// D3 (`--meth`) Bismark tags; `None` outside meth mode.
    pub meth: Option<MethFields<'a>>,
}

impl AlignedFields<'_> {
    /// The emitted strand (FLAG 0x10). Unlike bwa's raw `is_rev`, this is the
    /// effective strand: an unmapped read placed at its mate inherits the
    /// mate's strand, and its SEQ/QUAL are emitted reverse-complemented to it.
    pub fn is_rev(&self) -> bool {
        self.flag & 0x10 != 0
    }
}

/// D3 (`--meth`) tags of an emitted record. The packed record writes them
/// last, in the order XR, XG, XM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MethFields<'a> {
    /// XR:Z (read conversion), present on every meth-mode record.
    pub xr: &'a [u8],
    /// XG:Z (genome strand), mapped records only.
    pub xg: Option<&'a [u8]>,
    /// XM:Z (per-base methylation call), mapped records only.
    pub xm: Option<&'a [u8]>,
}

/// Which read of its template an emitted record belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mate {
    /// R1 of a pair, or a single-end read.
    R1,
    /// R2 of a pair.
    R2,
}

/// Receives each emitted record's [`AlignedFields`] instead of its packed BAM
/// body; see [`ResidentCohort::pair_emit_fields`]. Records arrive in the same
/// order a [`RecordSink`] would see them.
///
/// `mate` says which of a pair's two reads the record belongs to. `is_primary`
/// is true for the first record emitted for that read and false for each
/// supplementary (or, under `-a`, secondary) after it. The same borrow and
/// panic contract as [`RecordSink`] applies.
pub trait AlignedFieldsSink {
    fn emit(
        &mut self,
        origin: RecordOrigin,
        mate: Mate,
        is_primary: bool,
        fields: AlignedFields<'_>,
    );
}

/// Borrow a NUL-terminated shim string as bytes (`None` for NULL).
///
/// # Safety
/// `p` must be NULL or point at a NUL-terminated string that stays valid and
/// unmodified for `'a`.
unsafe fn opt_cstr<'a>(p: *const std::ffi::c_char) -> Option<&'a [u8]> {
    // SAFETY: forwarded from the caller.
    (!p.is_null()).then(|| unsafe { std::ffi::CStr::from_ptr(p) }.to_bytes())
}

impl<'a> AlignedFields<'a> {
    /// Borrow the shim's `BwaAlignedFields` as an [`AlignedFields`].
    ///
    /// # Safety
    /// Every pointer in `f` must be NULL or valid for `'a` as the shim's
    /// `BwaFieldSinkFn` contract states (strings NUL-terminated, `cigar`
    /// pointing at `n_cigar` ops when `n_cigar > 0`).
    unsafe fn from_raw(f: &'a bwa_mem3_sys::BwaAlignedFields) -> Self {
        let cigar = if f.n_cigar == 0 {
            &[][..]
        } else {
            // SAFETY: non-empty CIGAR -> `cigar` points at `n_cigar` ops (caller).
            unsafe { std::slice::from_raw_parts(f.cigar, f.n_cigar as usize) }
        };
        // SAFETY (all opt_cstr calls): string pointers are NULL or valid for 'a.
        let (md, mc, rg, sa, xa, xr, xg, xm) = unsafe {
            (
                opt_cstr(f.md),
                opt_cstr(f.mc),
                opt_cstr(f.rg),
                opt_cstr(f.sa),
                opt_cstr(f.xa),
                opt_cstr(f.xr),
                opt_cstr(f.xg),
                opt_cstr(f.xm),
            )
        };
        let flag_opt = |has: u8, v| (has != 0).then_some(v);
        Self {
            tid: f.tid,
            pos: f.pos,
            next_tid: f.next_tid,
            next_pos: f.next_pos,
            tlen: f.tlen,
            flag: f.flag,
            bin: f.bin,
            mapq: f.mapq,
            // Window bounds are never negative on real input; clamp defensively.
            query_start: usize::try_from(f.query_start).unwrap_or(0),
            query_end: usize::try_from(f.query_end).unwrap_or(0),
            cigar,
            nm: flag_opt(f.has_nm, f.nm),
            md,
            mc,
            mq: flag_opt(f.has_mq, f.mq),
            score: flag_opt(f.has_as, f.score),
            sub: flag_opt(f.has_xs, f.sub),
            rg,
            sa,
            pa: (f.has_pa != 0).then_some(f.pa),
            xa,
            hn: flag_opt(f.has_hn, f.hn),
            meth: xr.map(|xr| MethFields { xr, xg, xm }),
        }
    }
}

/// C -> Rust trampoline for [`AlignedFieldsSink`], passed to the shim's
/// `*_fields` pair-emit call as its `BwaFieldSinkFn`.
///
/// # Safety
///
/// `ctx` must be a `&mut SinkCtx<dyn AlignedFieldsSink>` cast to
/// `*mut c_void`, live and unaliased for the whole shim call. `f` must point
/// at a `BwaAlignedFields` whose pointers are valid for the duration of this
/// call only (the shim's `BwaFieldSinkFn` contract).
unsafe extern "C" fn field_sink_trampoline(
    ctx: *mut std::ffi::c_void,
    origin_kind: u32,
    origin_idx: usize,
    mate: u8,
    is_primary: std::ffi::c_int,
    f: *const bwa_mem3_sys::BwaAlignedFields,
) {
    // SAFETY: see `# Safety` -- `ctx` is a live `&mut SinkCtx` for this call.
    let ctx = unsafe { &mut *ctx.cast::<SinkCtx<'_, dyn AlignedFieldsSink>>() };
    // SAFETY: see `# Safety` -- `f` and everything it points at are valid for
    // this call, which outlives the borrow handed to `emit`.
    let fields = unsafe { AlignedFields::from_raw(&*f) };
    let origin = record_origin(origin_kind, origin_idx);
    let mate = if mate == 0 { Mate::R1 } else { Mate::R2 };
    ctx.emit_catching(|sink| sink.emit(origin, mate, is_primary != 0, fields));
}

// ---------------------------------------------------------------------------
// The resident cohort
// ---------------------------------------------------------------------------

/// Which region of a [`ResidentCohort`] a [`ResidentRange`] names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Region {
    Pairs,
    Singles,
}

/// Source of [`ResidentCohort`] identities. A range remembers its cohort's id,
/// not its address, so a range can never be accepted by a different cohort
/// that happens to be allocated where a dropped one used to be.
static NEXT_COHORT_ID: AtomicU64 = AtomicU64::new(0);

/// A reserved read-range of one [`ResidentCohort`]: a token for exactly one
/// sub-batch's reads.
///
/// Only [`ResidentCohort::reserve_pairs`]/[`reserve_singles`] create one, and
/// it is deliberately neither `Clone` nor `Copy`: every call that touches the
/// range's reads takes it by `&mut`, so the borrow checker proves no two calls
/// ever operate on the same reads at once. That, plus each range being its own
/// pointer-stable allocation inside the cohort, is what makes concurrent use
/// of one shared cohort sound. Offsets are in reads; a pair occupies two
/// consecutive reads.
///
/// [`reserve_singles`]: ResidentCohort::reserve_singles
#[derive(Debug)]
pub struct ResidentRange {
    cohort_id: u64,
    segment: *mut bwa_mem3_sys::BwaResidentSegment,
    region: Region,
    first: usize,
    n_reads: usize,
    /// Read bytes this range's writes copied into the cohort; handed back to
    /// the cohort's `heap_bytes` when the range is emitted.
    bytes: usize,
}

// SAFETY: the token is plain data plus a pointer that is only dereferenced
// (by the shim) through a call on the owning cohort, which checks `cohort_id`
// and serializes against whole-cohort calls. Moving it to another thread moves
// the exclusive right to use its reads, which is exactly what a pipeline hands
// from one step's worker to the next.
unsafe impl Send for ResidentRange {}
// SAFETY: every `&self` accessor reads immutable plain fields; the pointer is
// never dereferenced through a shared reference (range calls need `&mut`).
unsafe impl Sync for ResidentRange {}

impl ResidentRange {
    /// Read pairs in this range (0 for a single range).
    #[must_use]
    pub fn n_pairs(&self) -> usize {
        match self.region {
            Region::Pairs => self.n_reads / 2,
            Region::Singles => 0,
        }
    }
    /// Reads in this range (`2 * n_pairs` for a pair range, `n_singles` for a
    /// single range).
    #[must_use]
    pub fn n_reads(&self) -> usize {
        self.n_reads
    }
    /// Inclusive-start read offset of this range within its region, in
    /// reservation order. For a pair range `first() / 2` is the cohort-local
    /// index of its first pair, which the caller adds to the cohort's own global
    /// pair base for [`IdBases::first_pair_id`] (and may use as `origin_base` to
    /// keep record origins unique across ranges).
    #[must_use]
    pub fn first(&self) -> usize {
        self.first
    }
    /// Whether this is a pair range (else a single-end range).
    #[must_use]
    pub fn is_pairs(&self) -> bool {
        self.region == Region::Pairs
    }
}

/// A whole `-K` cohort's decoded reads + alnreg arrays, kept resident on the C
/// heap across `seed_extend -> infer_cohort -> pair_emit`. Every sub-batch
/// works on its own [`ResidentRange`] instead of owning a self-contained
/// [`AlnRegs`](crate::AlnRegs). Reads are copied in once, by the `write_*`
/// calls; each range's reads and alignment regions are released as soon as it
/// is emitted, and only the per-read headers stay until the cohort drops.
///
/// # Concurrency
///
/// `Send + Sync`, and every method is safe: share one cohort via `Arc` across a
/// worker pool. The guarantees, and what enforces each:
///
/// - **Range calls never alias.** [`write_pair`](Self::write_pair) (and the
///   other `write_*`), [`seed_extend`](Self::seed_extend),
///   [`pair_emit`](Self::pair_emit) and
///   [`pair_emit_fields`](Self::pair_emit_fields) take `&mut ResidentRange`, a
///   non-`Clone` token, so two of them can only run at once on different
///   ranges, which share no memory.
/// - **Reserving never moves a live range.** Each range is a separate C
///   allocation; [`reserve_pairs`](Self::reserve_pairs) only grows the table of
///   range pointers, which range calls never read. Reserves are serialized by a
///   mutex and may overlap in-flight range calls freely.
/// - **The barrier sees a quiescent cohort.** [`infer_cohort`](Self::infer_cohort)
///   reads every pair range, so it takes a read-write lock exclusively while
///   range calls hold it shared: it waits for in-flight range calls rather
///   than racing them. Because of that lock, a `pair_emit` sink must not call
///   back into its cohort except to reserve (see there), and no range call may
///   block waiting on another range call of the same cohort.
/// - **Lifecycle misuse is an error, not UB.** Writing a slot twice or after
///   `seed_extend`, seed-extending an unwritten or already-extended range,
///   emitting an unextended or already-emitted range, `infer_cohort` while a
///   pair range is unextended or after any was emitted, a range from another
///   cohort, a different [`BwaIndex`] than the one the cohort was first
///   extended against, and options whose `--meth` mode differs from the
///   cohort's all return `Err`.
///
/// A pair range reserved after `infer_cohort` is not part of that model; it is
/// the caller's to emit it with a model that fits.
pub struct ResidentCohort {
    handle: *mut bwa_mem3_sys::BwaResidentCohort,
    id: u64,
    /// Whether reads are written with the `--meth` projection; every call that
    /// aligns them must use options with the same mode.
    meth: bool,
    /// [`BwaIndex::id`] of the index the first `seed_extend` used (0 until
    /// then); every later call must use the same index, since the resident
    /// regions hold that index's coordinates.
    index_id: AtomicU64,
    /// C-heap bytes this cohort holds for its reads (see
    /// [`heap_bytes`](Self::heap_bytes)); raised by `reserve_*`/`write_*` and
    /// lowered when a range is emitted.
    heap_bytes: AtomicUsize,
    /// Serializes calls that read or grow the C segment table: `reserve_*`
    /// (writes it) and `infer_cohort` (walks it).
    table: Mutex<()>,
    /// Range calls hold this shared; `infer_cohort` holds it exclusively.
    ranges: RwLock<()>,
}

/// Map a resident shim status to an error: -3 is a lifecycle violation (caller
/// misuse, reported with the shim's complete message), anything else a shim
/// failure.
fn resident_status(rc: std::ffi::c_int, what: &str) -> Result<()> {
    match rc {
        0 => Ok(()),
        -3 => Err(Error::InvalidInput(shim_last_message())),
        _ => Err(shim_err(what)),
    }
}

impl ResidentCohort {
    /// Create an empty resident cohort. `meth` selects the per-read `--meth`
    /// projection applied at write time, so it must match the options every
    /// later call aligns with: a mismatch would seed projected reads against a
    /// plain index (or the reverse) and silently misalign, so those calls
    /// reject it with `Err`.
    pub fn new(meth: bool) -> Result<Self> {
        // SAFETY: an opaque C allocator with no preconditions; null on failure,
        // checked before use.
        let handle = unsafe { bwa_mem3_sys::bwa_shim_resident_cohort_new(meth.into()) };
        if handle.is_null() {
            return Err(shim_err("resident_cohort_new"));
        }
        Ok(Self {
            handle,
            id: NEXT_COHORT_ID.fetch_add(1, Ordering::Relaxed),
            meth,
            index_id: AtomicU64::new(0),
            heap_bytes: AtomicUsize::new(0),
            table: Mutex::new(()),
            ranges: RwLock::new(()),
        })
    }

    /// The locks guard no data (`()`), so a poisoned lock only means some other
    /// holder panicked; the C state it protects is still consistent (no C call
    /// panics), so recover the guard instead of propagating the poison.
    fn lock_table(&self) -> MutexGuard<'_, ()> {
        self.table.lock().unwrap_or_else(PoisonError::into_inner)
    }
    fn lock_ranges_shared(&self) -> RwLockReadGuard<'_, ()> {
        self.ranges.read().unwrap_or_else(PoisonError::into_inner)
    }
    fn lock_ranges_exclusive(&self) -> RwLockWriteGuard<'_, ()> {
        self.ranges.write().unwrap_or_else(PoisonError::into_inner)
    }

    /// Reject options whose `--meth` mode differs from the cohort's: the reads
    /// were written with (or without) the projection those options assume.
    fn check_cohort_meth(&self, opts: &MemOpts, what: &str) -> Result<()> {
        if self.meth == opts.meth() {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "{what}: cohort was created with meth={} but the options have meth={}",
                self.meth,
                opts.meth()
            )))
        }
    }

    /// Bind the cohort to `idx` on first use and reject any other index after:
    /// the resident regions hold coordinates into the index that produced
    /// them, and reading them against another would index its reference out
    /// of bounds.
    fn check_index(&self, idx: &BwaIndex, what: &str) -> Result<()> {
        let id = idx.id();
        match self
            .index_id
            .compare_exchange(0, id, Ordering::AcqRel, Ordering::Acquire)
        {
            Ok(_) => Ok(()),
            Err(bound) if bound == id => Ok(()),
            Err(_) => Err(Error::InvalidInput(format!(
                "{what}: this cohort was aligned against a different BwaIndex"
            ))),
        }
    }

    /// Reject a range reserved from a different cohort.
    fn check_owned(&self, range: &ResidentRange, what: &str) -> Result<()> {
        if range.cohort_id == self.id {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "{what}: range belongs to a different ResidentCohort"
            )))
        }
    }

    /// Reserve room for `n_pairs` read pairs and return the range naming them.
    /// A returned range is never invalidated by a later reserve.
    pub fn reserve_pairs(&self, n_pairs: usize) -> Result<ResidentRange> {
        let n_reads = n_pairs.checked_mul(2).ok_or_else(|| {
            Error::InvalidInput(format!(
                "reserve_pairs: {n_pairs} pairs overflow a read count"
            ))
        })?;
        self.reserve(Region::Pairs, n_reads)
    }

    /// Reserve room for `n_singles` single-end reads and return the range.
    pub fn reserve_singles(&self, n_singles: usize) -> Result<ResidentRange> {
        self.reserve(Region::Singles, n_singles)
    }

    fn reserve(&self, region: Region, n_reads: usize) -> Result<ResidentRange> {
        let mut first = 0usize;
        let segment = {
            let _table = self.lock_table();
            // SAFETY: `self.handle` is a valid cohort; the table mutex makes
            // this the only call reading or growing the segment table.
            unsafe {
                match region {
                    Region::Pairs => bwa_mem3_sys::bwa_shim_resident_reserve_pairs(
                        self.handle,
                        n_reads,
                        &mut first,
                    ),
                    Region::Singles => bwa_mem3_sys::bwa_shim_resident_reserve_singles(
                        self.handle,
                        n_reads,
                        &mut first,
                    ),
                }
            }
        };
        if segment.is_null() {
            return Err(shim_err("resident_reserve"));
        }
        // SAFETY: a pure query with no preconditions.
        let overhead = unsafe { bwa_mem3_sys::bwa_shim_resident_read_overhead() };
        self.heap_bytes
            .fetch_add(n_reads.saturating_mul(overhead), Ordering::Relaxed);
        Ok(ResidentRange {
            cohort_id: self.id,
            segment,
            region,
            first,
            n_reads,
            bytes: 0,
        })
    }

    /// Check that `range` is this cohort's and holds `region` reads.
    fn check_region(&self, range: &ResidentRange, region: Region, what: &str) -> Result<()> {
        self.check_owned(range, what)?;
        if range.region != region {
            let kind = match region {
                Region::Pairs => "pair",
                Region::Singles => "single",
            };
            return Err(Error::InvalidInput(format!("{what}: not a {kind} range")));
        }
        Ok(())
    }

    /// Validate a write to `range` of an item at local index `i` of `n_items`.
    fn check_write(
        &self,
        range: &ResidentRange,
        region: Region,
        i: usize,
        n_items: usize,
        what: &str,
    ) -> Result<()> {
        self.check_region(range, region, what)?;
        if i >= n_items {
            return Err(Error::InvalidInput(format!(
                "{what}: index {i} out of range ({n_items} in the range)"
            )));
        }
        Ok(())
    }

    /// Validate a whole-range batch write: owned range of `region` and exactly
    /// `expected` items.
    fn check_batch(
        &self,
        range: &ResidentRange,
        region: Region,
        n: usize,
        expected: usize,
        what: &str,
    ) -> Result<()> {
        self.check_region(range, region, what)?;
        if n != expected {
            return Err(Error::InvalidInput(format!(
                "{what}: a batch of {n} for a range of {expected}"
            )));
        }
        Ok(())
    }

    /// Record `added` bytes a write to `range` copied into the cohort. Bytes
    /// copied before a failure are held by the segment (a partly written range
    /// is never emitted), so they count either way.
    fn account_write(&self, range: &mut ResidentRange, added: usize) {
        self.heap_bytes.fetch_add(added, Ordering::Relaxed);
        range.bytes += added;
    }

    /// Decode pair `i` of a pair `range` from borrowed bytes into its resident
    /// slot. Each slot may be written once, before the range's `seed_extend`,
    /// and not at all in a range written with [`write_pairs`](Self::write_pairs).
    pub fn write_pair(
        &self,
        range: &mut ResidentRange,
        i: usize,
        pair: ReadPair<'_>,
    ) -> Result<()> {
        self.check_write(range, Region::Pairs, i, range.n_pairs(), "write_pair")?;
        pair.validate()?;
        let c = c_pair(&pair);
        let _ranges = self.lock_ranges_shared();
        let mut added = 0usize;
        // SAFETY: `range.segment` is a live segment of this cohort (checked
        // id; segments live until the cohort drops) that no other call is
        // touching (`&mut range`); `c` borrows `pair`'s bytes for the call
        // only; `added` is a live `usize`.
        let rc =
            unsafe { bwa_mem3_sys::bwa_shim_resident_write_pair(range.segment, i, &c, &mut added) };
        self.account_write(range, added);
        resident_status(rc, "resident_write_pair")
    }

    /// Decode every pair of pair `range` at once: `pairs.len()` must equal
    /// [`ResidentRange::n_pairs`], and no slot may have been written yet. All
    /// the reads' bytes go into one allocation (the CLI reader's per-chunk
    /// arena) instead of three per read, in one FFI call under one lock.
    pub fn write_pairs(&self, range: &mut ResidentRange, pairs: &[ReadPair<'_>]) -> Result<()> {
        self.check_batch(
            range,
            Region::Pairs,
            pairs.len(),
            range.n_pairs(),
            "write_pairs",
        )?;
        let cs = pairs
            .iter()
            .map(|pair| pair.validate().map(|()| c_pair(pair)))
            .collect::<Result<Vec<bwa_mem3_sys::BwaReadPair>>>()?;
        let _ranges = self.lock_ranges_shared();
        let mut added = 0usize;
        // SAFETY: `range.segment` is a live segment of this cohort (checked
        // id; segments live until the cohort drops) that no other call is
        // touching (`&mut range`); `cs` and the bytes it borrows outlive the
        // call, and the shim copies every byte before returning; `added` is a
        // live `usize`.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_write_pairs(
                range.segment,
                cs.as_ptr(),
                cs.len(),
                &mut added,
            )
        };
        self.account_write(range, added);
        resident_status(rc, "resident_write_pairs")
    }

    /// Decode single `i` of a single `range` from borrowed bytes. Each slot may
    /// be written once, before the range's `seed_extend`, and not at all in a
    /// range written with [`write_singles`](Self::write_singles).
    pub fn write_single(
        &self,
        range: &mut ResidentRange,
        i: usize,
        read: SingleRead<'_>,
    ) -> Result<()> {
        self.check_write(range, Region::Singles, i, range.n_reads, "write_single")?;
        read.validate()?;
        let c = c_single(&read);
        let _ranges = self.lock_ranges_shared();
        let mut added = 0usize;
        // SAFETY: as `write_pair`.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_write_single(range.segment, i, &c, &mut added)
        };
        self.account_write(range, added);
        resident_status(rc, "resident_write_single")
    }

    /// Decode every read of single `range` at once: `reads.len()` must equal
    /// [`ResidentRange::n_reads`], and no slot may have been written yet. See
    /// [`write_pairs`](Self::write_pairs).
    pub fn write_singles(&self, range: &mut ResidentRange, reads: &[SingleRead<'_>]) -> Result<()> {
        self.check_batch(
            range,
            Region::Singles,
            reads.len(),
            range.n_reads,
            "write_singles",
        )?;
        let cs = reads
            .iter()
            .map(|read| read.validate().map(|()| c_single(read)))
            .collect::<Result<Vec<bwa_mem3_sys::BwaSingleRead>>>()?;
        let _ranges = self.lock_ranges_shared();
        let mut added = 0usize;
        // SAFETY: as `write_pairs`.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_write_singles(
                range.segment,
                cs.as_ptr(),
                cs.len(),
                &mut added,
            )
        };
        self.account_write(range, added);
        resident_status(rc, "resident_write_singles")
    }

    /// Seed + SE-extend the reads of `range`. Every slot must be written, and a
    /// range is extended exactly once.
    pub fn seed_extend(
        &self,
        idx: &BwaIndex,
        opts: &MemOpts,
        scratch: &mut AlignScratch,
        range: &mut ResidentRange,
    ) -> Result<()> {
        self.check_owned(range, "seed_extend")?;
        self.check_cohort_meth(opts, "seed_extend")?;
        check_meth_consistency(idx, opts)?;
        self.check_index(idx, "seed_extend")?;
        let _ranges = self.lock_ranges_shared();
        // SAFETY: `range.segment` is a live segment of this cohort touched by no
        // other call (`&mut range`, checked id); `idx`/`opts`/`scratch.handle`
        // are valid for the call and `scratch` is exclusively borrowed.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_seed_extend(
                idx.raw(),
                opts.as_ptr(),
                scratch.handle,
                range.segment,
            )
        };
        resident_status(rc, "resident_seed_extend")
    }

    /// `mem_pestat` over the whole cohort's pair region. Every pair range must
    /// have been seed-extended and none emitted yet (else `Err`); waits for
    /// in-flight range calls. Feed the result to every
    /// [`pair_emit`](Self::pair_emit) of this cohort. The model is a histogram,
    /// so it does not depend on the order the ranges were reserved or extended.
    pub fn infer_cohort(&self, idx: &BwaIndex, opts: &MemOpts) -> Result<MemPeStat> {
        self.check_cohort_meth(opts, "infer_cohort")?;
        self.check_index(idx, "infer_cohort")?;
        let mut out = MemPeStat::zero()?;
        let _ranges = self.lock_ranges_exclusive();
        let _table = self.lock_table();
        // SAFETY: `idx`/`opts` valid for the call; `out.as_mut_ptr()` is a live
        // `mem_pestat_t[4]` owned by `out`; the exclusive lock plus the table
        // mutex make this the only call touching the cohort.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_pestat_cohort(
                idx.raw(),
                opts.as_ptr(),
                self.handle,
                out.as_mut_ptr(),
            )
        };
        resident_status(rc, "resident_pestat_cohort")?;
        Ok(out)
    }

    /// Bytes this cohort holds on the C heap for its reads: each reserved
    /// read's headers plus the name, bases and qualities copied into it (and,
    /// under `--meth`, its original bases). A range's copied bytes are held
    /// until it is emitted, when the shim releases its reads and regions and
    /// only the headers remain. It does not count the alignment
    /// regions, which `seed_extend` and mate rescue grow in place; that is the
    /// difference from [`AlnRegs::heap_bytes`](crate::AlnRegs::heap_bytes),
    /// which counts them. Lock-free: it never waits on in-flight range calls.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        self.heap_bytes.load(Ordering::Relaxed)
    }

    /// Whether `range`'s segment still holds any read string or alignment
    /// region on the C heap: `true` once written, `false` after it has been
    /// emitted. A test hook for the release-at-emit contract.
    #[doc(hidden)]
    #[must_use]
    pub fn range_holds_reads(&self, range: &ResidentRange) -> bool {
        if self.check_owned(range, "range_holds_reads").is_err() {
            return false;
        }
        let _ranges = self.lock_ranges_shared();
        // SAFETY: `range.segment` is a live segment of this cohort (checked id;
        // segments live until the cohort drops); the query only reads it, and
        // no call can mutate it concurrently because every mutating call needs
        // the range by `&mut`.
        unsafe { bwa_mem3_sys::bwa_shim_resident_segment_holds_reads(range.segment) != 0 }
    }

    /// Hand an emitted range's read bytes back: the shim released them at the
    /// end of the emit, even when the sink panicked partway.
    fn release_bytes(&self, range: &mut ResidentRange) {
        self.heap_bytes.fetch_sub(range.bytes, Ordering::Relaxed);
        range.bytes = 0;
    }

    /// Validate an emit call and take the shared range lock for it.
    fn begin_emit(
        &self,
        idx: &BwaIndex,
        opts: &MemOpts,
        range: &ResidentRange,
        pestat: Option<&MemPeStat>,
        what: &str,
    ) -> Result<RwLockReadGuard<'_, ()>> {
        self.check_owned(range, what)?;
        self.check_cohort_meth(opts, what)?;
        self.check_index(idx, what)?;
        if range.region == Region::Pairs && range.n_reads > 0 && pestat.is_none() {
            return Err(Error::InvalidInput(format!(
                "{what}: a range with pairs requires the cohort pestat"
            )));
        }
        Ok(self.lock_ranges_shared())
    }

    /// Pair / mate-rescue / emit the reads of a seed-extended `range`, then
    /// release its reads and alignment regions; each range is emitted once.
    /// `ids` gives the global read ordinal of the range's first pair/single;
    /// `origin_base` is added to the
    /// local index for each record's [`RecordOrigin`]. `pestat` is the cohort
    /// model from [`infer_cohort`](Self::infer_cohort); it may be `None` only
    /// for a single-end range.
    ///
    /// `sink` runs while this call holds the cohort's shared range lock, so it
    /// must not call any method of this cohort other than
    /// [`reserve_pairs`](Self::reserve_pairs)/[`reserve_singles`](Self::reserve_singles)
    /// and [`heap_bytes`](Self::heap_bytes). [`infer_cohort`](Self::infer_cohort)
    /// would wait on the lock the sink's own call holds, and a nested range
    /// call (`write_*`, `seed_extend`, `pair_emit*`) re-acquires it shared,
    /// which blocks forever once another thread is queued for it exclusively
    /// (the lock prefers writers).
    #[allow(clippy::too_many_arguments)]
    pub fn pair_emit(
        &self,
        idx: &BwaIndex,
        opts: &MemOpts,
        scratch: &mut AlignScratch,
        range: &mut ResidentRange,
        pestat: Option<&MemPeStat>,
        ids: IdBases,
        origin_base: usize,
        sink: &mut dyn RecordSink,
    ) -> Result<()> {
        let _ranges = self.begin_emit(idx, opts, range, pestat, "pair_emit")?;
        let mut sink_ctx = SinkCtx { sink, panic: None };
        let ctx_ptr = std::ptr::addr_of_mut!(sink_ctx).cast::<std::ffi::c_void>();
        // SAFETY: `range.segment` is a live segment of this cohort touched by no
        // other call (`&mut range`, checked id); `idx`/`opts`/`scratch.handle`
        // are valid; `sink_trampoline`/`ctx_ptr` satisfy the trampoline's
        // contract (`ctx_ptr` outlives the call and is not aliased during it).
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_pair_emit(
                idx.raw(),
                opts.as_ptr(),
                scratch.handle,
                range.segment,
                pestat.map_or(std::ptr::null(), MemPeStat::as_ptr),
                c_ids(ids),
                origin_base,
                Some(sink_trampoline),
                ctx_ptr,
            )
        };
        if rc == 0 {
            self.release_bytes(range);
        }
        if let Some(payload) = sink_ctx.panic {
            std::panic::resume_unwind(payload);
        }
        resident_status(rc, "resident_pair_emit")
    }

    /// [`pair_emit`](Self::pair_emit), reporting each record's
    /// [`AlignedFields`] to an [`AlignedFieldsSink`] instead of its packed BAM
    /// body. The same records arrive in the same order, and a record built
    /// from the reported fields plus the input read is byte-identical to the
    /// packed one. The same sink restriction as `pair_emit` applies.
    #[allow(clippy::too_many_arguments)]
    pub fn pair_emit_fields(
        &self,
        idx: &BwaIndex,
        opts: &MemOpts,
        scratch: &mut AlignScratch,
        range: &mut ResidentRange,
        pestat: Option<&MemPeStat>,
        ids: IdBases,
        origin_base: usize,
        sink: &mut dyn AlignedFieldsSink,
    ) -> Result<()> {
        let _ranges = self.begin_emit(idx, opts, range, pestat, "pair_emit_fields")?;
        let mut sink_ctx = SinkCtx { sink, panic: None };
        let ctx_ptr = std::ptr::addr_of_mut!(sink_ctx).cast::<std::ffi::c_void>();
        // SAFETY: as `pair_emit`, with `field_sink_trampoline`.
        let rc = unsafe {
            bwa_mem3_sys::bwa_shim_resident_pair_emit_fields(
                idx.raw(),
                opts.as_ptr(),
                scratch.handle,
                range.segment,
                pestat.map_or(std::ptr::null(), MemPeStat::as_ptr),
                c_ids(ids),
                origin_base,
                Some(field_sink_trampoline),
                ctx_ptr,
            )
        };
        if rc == 0 {
            self.release_bytes(range);
        }
        if let Some(payload) = sink_ctx.panic {
            std::panic::resume_unwind(payload);
        }
        resident_status(rc, "resident_pair_emit_fields")
    }
}

/// The C view of `ids`.
fn c_ids(ids: IdBases) -> bwa_mem3_sys::BwaIdBases {
    bwa_mem3_sys::BwaIdBases {
        first_single_id: ids.first_single_id,
        first_pair_id: ids.first_pair_id,
    }
}

impl Drop for ResidentCohort {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            // SAFETY: valid non-null handle owned exclusively by `self` (`&mut`
            // proves no call is in flight); this is the only place it is freed
            // and `Drop` runs at most once. Outstanding `ResidentRange`s become
            // unusable: every range call goes through a cohort, and no other
            // cohort accepts their id.
            unsafe { bwa_mem3_sys::bwa_shim_resident_cohort_free(self.handle) };
        }
    }
}

// SAFETY: the cohort is a C-owned allocation with no aliasing back into the
// shared index/opts and no thread affinity, so moving it between threads is
// sound.
unsafe impl Send for ResidentCohort {}
// SAFETY: shared (`&self`) access from several threads is sound because every
// path into the C state is synchronized (see the type's "Concurrency" doc):
// the segment table is only touched under `table`; a range's reads are only
// touched through an exclusive `&mut ResidentRange` whose segment is disjoint
// from every other range's and pointer-stable for the cohort's lifetime; and
// `infer_cohort`, which reads every pair range, holds `ranges` exclusively,
// excluding all range calls (which hold it shared). The only other shared
// state is atomics (`index_id`, `heap_bytes`). The C code keeps no other
// mutable shared state per cohort.
unsafe impl Sync for ResidentCohort {}
