//! An empty batch is a valid batch: every `&[ReadPair]` entry point must
//! return an empty result rather than erroring or aborting.
//!
//! `bwa-rs` itself never reaches this path (`main.rs` breaks out of its loop on
//! `pending.is_empty()`), so only a library caller does — which is precisely
//! why it needs a test.
//!
//! What the shim does with `n_pairs == 0` is size *every* allocation on that
//! path from a zero read count: `copy_pairs_to_seqs`' `seqs` array and, because
//! this crate keeps the pre-0.9.0 nreads sizing (see the divergence note in
//! `worker_alloc`), all three of `regs` / `chain_scratch` / `seed_scratch`.
//! `calloc(0, n)` and `malloc(0)` may each legally return NULL, and the shim's
//! `xassert`s on those pointers do NOT compile out under `NDEBUG` — so an
//! unguarded zero-size return would abort a *release* build, and the `seqs` one
//! would surface as a spurious `Err`.
//!
//! Both are guarded now. Note this test does not *reproduce* the abort on every
//! platform: glibc and macOS both return a unique non-NULL pointer for a
//! zero-size request, so it passes with or without the guards there. It pins
//! the empty-batch contract at the API level; the guards are what make that
//! contract hold on an allocator that takes the other option.
//!
//! Gated on `bwa-mem3` on PATH (needed for `index`).

mod common;
mod phix_seq;

use bwa_mem3_rs::{align_batch, estimate_pestat, extend_batch, seed_batch, BwaIndex, MemOpts};

/// Build a PhiX index and return it alongside PE-enabled options.
fn phix_index_and_opts(dir: &std::path::Path, bwa: &str) -> (BwaIndex, MemOpts) {
    let ref_fa = common::setup_phix_index(dir, bwa, phix_seq::PHIX_SEQ);
    let idx = BwaIndex::load(&ref_fa).expect("load index");
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);
    (idx, opts)
}

#[test]
fn align_batch_on_empty_slice_returns_empty() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (idx, opts) = phix_index_and_opts(tmp.path(), &bwa);

    let (aln, _pes) = align_batch(&idx, &opts, &[], None).expect("align_batch on empty slice");
    assert_eq!(aln.len(), 0, "empty input must produce no records");
    assert_eq!(aln.iter().count(), 0, "iterator must yield nothing");
}

#[test]
fn seed_then_extend_on_empty_slice_returns_empty() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (idx, opts) = phix_index_and_opts(tmp.path(), &bwa);

    // The phase split is the case that matters most: it is the reason this
    // crate keeps the nreads-sized scratch upstream shrank to per-thread, so
    // it is the path with the most zero-size allocations to get wrong.
    let seeds = seed_batch(&idx, &opts, &[]).expect("seed_batch on empty slice");
    let (aln, _pes) = extend_batch(&idx, &opts, seeds, &[], None).expect("extend_batch");
    assert_eq!(aln.len(), 0, "empty input must produce no records");
}

#[test]
fn estimate_pestat_on_empty_slice_succeeds() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (idx, opts) = phix_index_and_opts(tmp.path(), &bwa);

    // No pairs means no insert-size evidence; mem_pestat leaves every
    // orientation failed. The contract under test is that it returns Ok at all.
    estimate_pestat(&idx, &opts, &[]).expect("estimate_pestat on empty slice");
}

#[test]
fn seeds_dropped_without_extending_on_empty_slice() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (idx, opts) = phix_index_and_opts(tmp.path(), &bwa);

    // Exercises worker_free's NULL-safety: Seeds' Drop runs shim_seeds_free on
    // a worker whose scratch pointers were never allocated.
    let seeds = seed_batch(&idx, &opts, &[]).expect("seed_batch on empty slice");
    drop(seeds);
}
