//! Safe Rust API for [bwa-mem3] alignment.
//!
//! This crate exposes bwa-mem3 alignment through a blocking, reentrant
//! per-batch API with packed BAM output. Callers own all parallelism — every
//! function is synchronous on its calling thread and safe to call
//! concurrently from multiple threads sharing the same [`BwaIndex`].
//!
//! Phase-split API: [`seed_batch`] → [`extend_batch`] → packed BAM records.
//! [`align_batch`] is a thin wrapper = seed + extend. [`estimate_pestat`]
//! runs seed + SE extension + `mem_pestat` only, skipping pairing and
//! emission.
//!
//! The three-phase cohort API splits alignment further so a caller can share
//! one insert-size model across several sub-batches: [`seed_extend`] runs
//! seeding and SE extension for one sub-batch of reads, producing per-read
//! [`AlnRegs`]; [`MemPeStat::infer_cohort`] folds the `AlnRegs` from all of a
//! cohort's sub-batches into a single insert-size model; and [`pair_emit`]
//! takes that model plus one sub-batch's `AlnRegs` and performs pairing,
//! mate rescue, and emission, streaming packed BAM records to a caller-
//! supplied [`RecordSink`] as they're produced. Callers own reusable
//! [`AlignScratch`] buffers across calls (avoiding per-batch allocation) and
//! assign each read a global ordinal via [`IdBases`] so record IDs stay
//! stable regardless of how reads are split into sub-batches. Loading a
//! shared, thread-safe index for this API goes through
//! [`BwaIndex::load_with_threads`].
//!
//! [bwa-mem3]: https://github.com/fg-labs/bwa-mem3

pub mod align;
pub mod build_info;
pub mod error;
pub mod index;
pub mod opts;
pub mod resident;
pub mod shm;
pub mod xa;

pub use align::{
    align_batch, estimate_pestat, extend_batch, pair_emit, seed_batch, seed_extend, AlignScratch,
    AlignmentBatch, AlnRegs, IdBases, ReadBatch, ReadPair, Record, RecordOrigin, RecordSink,
    RecordVec, Seeds, SingleRead,
};
pub use build_info::{build_info, version, BuildInfo};
pub use error::{Error, Result};
pub use index::BwaIndex;
pub use opts::{MemOpts, MemPeStat, MethScoring, Mode, PeOrient, PeOrientation, SeedOrder};
pub use resident::{
    AlignedFields, AlignedFieldsSink, Mate, MethFields, ResidentCohort, ResidentRange,
};
pub use xa::{parse_xa, AuxHit};
