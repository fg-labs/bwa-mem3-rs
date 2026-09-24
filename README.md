# bwa-mem3-rs

[![CI](https://github.com/fg-labs/bwa-mem3-rs/actions/workflows/check.yml/badge.svg)](https://github.com/fg-labs/bwa-mem3-rs/actions/workflows/check.yml)
[![crates.io: bwa-mem3-rs](https://img.shields.io/crates/v/bwa-mem3-rs.svg?label=bwa-mem3-rs)](https://crates.io/crates/bwa-mem3-rs)
[![crates.io: bwa-mem3-sys](https://img.shields.io/crates/v/bwa-mem3-sys.svg?label=bwa-mem3-sys)](https://crates.io/crates/bwa-mem3-sys)
[![crates.io: bwa-mem3-rs-cli](https://img.shields.io/crates/v/bwa-mem3-rs-cli.svg?label=bwa-mem3-rs-cli)](https://crates.io/crates/bwa-mem3-rs-cli)
[![docs.rs](https://img.shields.io/docsrs/bwa-mem3-rs?label=docs.rs)](https://docs.rs/bwa-mem3-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Rust FFI crate for [bwa-mem3](https://github.com/fg-labs/bwa-mem3) with packed-BAM output and caller-owned parallelism.

Hosted under [`fg-labs/bwa-mem3-rs`](https://github.com/fg-labs/bwa-mem3-rs); it vendors [`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3) at the `main` branch — our integration branch over upstream with Apple Silicon / NEON + Linux-arm64 / CI fixes applied.

## What it does

- Aligns read pairs via bwa-mem3, returning **packed BAM records** per the BAM spec (ready to stream to any BGZF writer).
- Blocking, reentrant, single-threaded per call — **caller drives all parallelism** via rayon, tokio `spawn_blocking`, or any work-stealing pool.
- Phase-split API: `seed_batch` → `extend_batch` lets work-stealing pools balance seeding and extension independently. `align_batch` is a thin convenience wrapper.
- Runs on **x86_64** (AVX2 by default) and **aarch64** (Apple Silicon + AWS Graviton).

## Quick start

```rust
use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

let idx = BwaIndex::load("hg38.fa")?;                 // bwa-mem3 index prefix
let mut opts = MemOpts::new()?;
opts.set_pe(true);

let pairs = vec![ReadPair {
    name_r1: b"read1",
    seq_r1:  b"ACGT...",
    qual_r1: Some(b"IIII..."),
    name_r2: b"read1",
    seq_r2:  b"TTTT...",
    qual_r2: Some(b"IIII..."),
}];

let (aln, _pestat) = align_batch(&idx, &opts, &pairs, None)?;
for record in aln.iter() {
    // record.bytes is a packed BAM record: [u32 le block_size][data].
    // Feed to noodles_bam, rust-htslib, or your own BGZF writer.
}
```

## Work-stealing parallelism

```rust
use bwa_mem3_rs::*;
use rayon::prelude::*;
use std::sync::Arc;

let idx = Arc::new(BwaIndex::load("hg38.fa")?);
let mut opts = MemOpts::new()?;
opts.set_pe(true);

let results: Vec<_> = batches
    .par_iter()
    .map(|batch| align_batch(&idx, &opts, batch, None))
    .collect::<Result<_, _>>()?;
```

The crate itself never spawns threads. Every public function blocks its calling thread and is safe to call concurrently from multiple threads sharing the same `BwaIndex`.

### Cohort-exact output from a work-stealing pool

`align_batch` estimates the insert-size model from whatever it is given, so output matches `bwa-mem3 mem` only when one call covers one `-K` cohort. The three-phase API keeps per-read work fine-grained and the model cohort-wide:

```rust
use bwa_mem3_rs::*;
let idx = BwaIndex::load_with_threads("hg38.fa", 8)?;
let opts = MemOpts::new()?;
let mut scratch = AlignScratch::new()?;               // one per worker thread
// cohort_sub_batches, cohort_first_pair_id, and my_sink are caller-supplied
let regs: Vec<AlnRegs> = cohort_sub_batches.iter()   // any size, any thread
    .map(|b| seed_extend(&idx, &opts, &mut scratch, &ReadBatch { pairs: b, singles: &[] }))
    .collect::<Result<_>>()?;
let pestat = MemPeStat::infer_cohort(&idx, &opts, &regs)?;
// Sub-batches may differ in length, so advance the pair id by each one's
// actual pair count -- a fixed stride would give later calls the wrong global
// ordinal and break the cohort-exact tie-break.
let mut next_pair_id = cohort_first_pair_id;
for (k, r) in regs.into_iter().enumerate() {
    let ids = IdBases { first_single_id: 0, first_pair_id: next_pair_id };
    pair_emit(&idx, &opts, &mut scratch, r, Some(&pestat), ids, &mut my_sink)?;
    next_pair_id += cohort_sub_batches[k].len() as u64;
}
```

`IdBases` carries the global read ordinal bwa-mem3 uses for tie-breaks; see the type's docs for the `-p` cohort formulas. `AlignScratch` is `Send`, `AlnRegs` is `Send`, `BwaIndex`/`MemOpts`/`MemPeStat` are `Send + Sync`.

### Resident cohorts: no per-sub-batch copies

Each `AlnRegs` above owns its own copy of its reads and alignment regions, allocated by `seed_extend` and freed by `pair_emit`. A `ResidentCohort` keeps a whole cohort's reads and regions resident instead: every sub-batch reserves a range of it, reads are copied in once, and each range's reads and regions are released as soon as it is emitted. Output is byte-identical to the three-phase API above.

```rust
use bwa_mem3_rs::*;
let cohort = ResidentCohort::new(opts.meth())?;       // one per -K cohort, shared via Arc
// Serially, in input order: reserve one range per sub-batch.
let mut ranges: Vec<ResidentRange> = cohort_sub_batches.iter()
    .map(|b| cohort.reserve_pairs(b.len()))
    .collect::<Result<_>>()?;
// On any thread: write each range's reads, then seed + extend them.
for (range, b) in ranges.iter_mut().zip(&cohort_sub_batches) {
    cohort.write_pairs(range, b)?;
    cohort.seed_extend(&idx, &opts, &mut scratch, range)?;
}
let pestat = cohort.infer_cohort(&idx, &opts)?;        // once every range is extended
for range in &mut ranges {
    let pair_base = range.first() / 2;                  // this range's first pair, cohort-local
    let ids = IdBases { first_single_id: 0, first_pair_id: cohort_first_pair_id + pair_base as u64 };
    // origin_base = pair_base keeps RecordOrigin::Pair(i) unique across ranges.
    cohort.pair_emit(&idx, &opts, &mut scratch, range, Some(&pestat), ids, pair_base, &mut my_sink)?;
}
```

`ResidentCohort` is `Send + Sync` and every method is safe: a `ResidentRange` is a move-only token (neither `Clone` nor `Copy`), and range calls take it by `&mut`, so two calls can only run at once on different ranges. Lifecycle misuse (writing a slot twice, emitting a range before extending it or twice, a range from another cohort, a different `BwaIndex` than the cohort was aligned against, options whose `--meth` mode differs from the cohort's) returns `Err`. `pair_emit_fields` reports each record as structured `AlignedFields` through an `AlignedFieldsSink` instead of a packed BAM body; a record rebuilt from those fields plus the input read is byte-identical to the packed one.

## CLI

The `bwa-mem3-rs-cli` crate ships a minimal `bwa-rs` binary that wraps the library:

```
bwa-rs mem <ref-prefix> <r1.fq[.gz]> <r2.fq[.gz]> [-o out.bam]
```

## Index

Use the upstream CLI to build an index:

```
bwa-mem3 index <prefix>.fa
```

Index building is intentionally out-of-scope for this crate.

## Architecture

- `bwa-mem3-sys` vendors upstream bwa-mem3 sources and compiles them alongside
  our custom C++ shim (`shim/bwa_shim*.cpp`). The shim bridges to upstream's
  public API via opaque pointers and emits packed BAM directly from bwa-mem3's
  `mem_aln_t` (no SAM round-trip).
- `bwa-mem3-rs` wraps the FFI crate in a safe Rust API with `Send`/`Sync`
  semantics that let callers drive parallelism externally.
- `bwa-mem3-rs-cli` is a thin `bwa-rs mem` binary.
- `bwa-mem3-py` is a [pyo3]-based Python wheel (`import bwa_mem3`) that
  exposes the same surface to Python — see `bwa-mem3-py/README.md`. It
  lives outside the cargo workspace and is built via `maturin` /
  `pixi run develop`.

See `CLAUDE.md` for the project-specific gotchas (CIGAR opcode remap, 2-bit
seq encoding, vendored `MATE_SORT=0` invariant, etc.).

## Testing

```bash
# Unit tests (no index required)
cargo test --workspace

# Integration tests (require a bwa-mem3 index prefix)
BWA_MEM3_RS_TEST_REF=/path/to/hg38.fa cargo test --workspace
```

The end-to-end regression test (PhiX reference + simulated reads + CLI round-trip) runs automatically in CI once `bwa-mem3` and `samtools` are available; locally set `BWA_MEM3_BIN=/path/to/bwa-mem3` if the CLI is not on `PATH`.

## License

MIT.

[pyo3]: https://pyo3.rs
