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
//! [bwa-mem3]: https://github.com/fg-labs/bwa-mem3

pub mod align;
pub mod build_info;
pub mod error;
pub mod index;
pub mod opts;
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
pub use xa::{parse_xa, AuxHit};
