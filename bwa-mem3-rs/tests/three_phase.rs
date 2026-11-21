//! Wrapper-level three-phase tests. These build a small PhiX index at test
//! time (via `bwa-mem3 index`), so the safety-critical checks here -- byte
//! identity with `align_batch`, cross-thread `Send`, and panic-across-FFI
//! unwind -- run in CI wherever the `bwa-mem3` binary is available, rather than
//! being gated on a full hg38 index that CI never provides.
//!
//! A missing `bwa-mem3` makes each test skip (returning early with a message),
//! unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set, in which case it is a hard
//! failure -- the same skip-vs-panic convention the sys and cli crates use, so
//! a CI job that promises the tools cannot silently degrade into a no-op.
//!
//! PhiX (~5.5 kb) loads in milliseconds, so unlike the old hg38 path there is
//! no OOM risk from concurrent loads; the tests that only *read* the index
//! still share one `Arc<BwaIndex>` via [`shared_idx`] to avoid re-indexing per
//! test. [`load_with_threads_matches_load`] is the deliberate exception: it
//! exists to test loading itself, so it loads its own pair from the prefix.

use std::sync::Arc;

use bwa_mem3_rs::{
    align_batch, build_info, extend_batch, pair_emit, seed_batch, seed_extend, version,
    AlignScratch, AlignedFields, AlignedFieldsSink, AlnRegs, BwaIndex, Error, IdBases, Mate,
    MemOpts, MemPeStat, ReadBatch, ReadPair, RecordOrigin, RecordSink, RecordVec, ResidentCohort,
    ResidentRange, SingleRead,
};
use rstest::rstest;

mod support;
use support::*;

#[test]
fn version_and_build_info_are_populated() {
    assert_eq!(version(), bwa_mem3_sys::build_info::VERSION);
    let b = build_info();
    assert_eq!(b.bwa_mem3_version, version());
    assert!(!b.compiler.is_empty());
    if cfg!(target_arch = "x86_64") {
        assert_eq!(b.x86_tiers.len(), 5);
    } else {
        assert!(b.x86_tiers.is_empty());
    }
}

#[test]
fn three_phase_matches_align_batch_single_cohort() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let opts = MemOpts::new().unwrap();
    let (names, r1, r2) = pairs(64);
    let qual = vec![b'I'; 150];
    let pv: Vec<ReadPair<'_>> = (0..64)
        .map(|i| ReadPair {
            name_r1: names[i].as_bytes(),
            seq_r1: &r1[i],
            qual_r1: Some(&qual),
            name_r2: names[i].as_bytes(),
            seq_r2: &r2[i],
            qual_r2: Some(&qual),
        })
        .collect();

    let (legacy, _) = align_batch(&idx, &opts, &pv, None).unwrap();
    let legacy: Vec<(usize, Vec<u8>)> = legacy
        .iter()
        .map(|r| (r.pair_idx, r.bytes[4..].to_vec()))
        .collect();

    let mut scratch = AlignScratch::new().unwrap();
    let regs: Vec<AlnRegs> = pv
        .chunks(16)
        .map(|c| {
            seed_extend(
                &idx,
                &opts,
                &mut scratch,
                &ReadBatch {
                    pairs: c,
                    singles: &[],
                },
            )
            .unwrap()
        })
        .collect();
    let pestat = MemPeStat::infer_cohort(&idx, &opts, &regs).unwrap();
    let mut got = Vec::new();
    for (k, r) in regs.into_iter().enumerate() {
        let mut sink = RecordVec::default();
        let ids = IdBases {
            first_single_id: 0,
            first_pair_id: (k * 16) as u64,
        };
        pair_emit(&idx, &opts, &mut scratch, r, Some(&pestat), ids, &mut sink).unwrap();
        for (origin, body) in sink.records {
            let RecordOrigin::Pair(i) = origin else {
                panic!("unexpected single")
            };
            got.push((k * 16 + i, body));
        }
    }
    assert_eq!(got, legacy);
}

/// `ResidentCohort` driven in sub-chunks must reproduce `align_batch`
/// byte-for-byte at every sub-chunk size — the safe-wrapper equivalent of
/// `three_phase_matches_align_batch_single_cohort`, with the cohort's reads
/// staying resident (no per-sub-batch `AlnRegs`). The PhiX-derived pairs map
/// (FR/RF, indels, chimeras, half-mapped), so alignment, pairing, mate rescue
/// and the cohort model all run, not just the unmapped path.
#[test]
fn resident_cohort_matches_align_batch_at_every_sub_batch_size() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let opts = MemOpts::new().unwrap();
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();

    let (legacy, _) = align_batch(&idx, &opts, &pv, None).unwrap();
    let legacy: Vec<(usize, Vec<u8>)> = legacy
        .iter()
        .map(|r| (r.pair_idx, r.bytes[4..].to_vec()))
        .collect();

    let mut scratch = AlignScratch::new().unwrap();
    for sub in [1usize, 3, 16, 64, 128] {
        let cohort = ResidentCohort::new(false).unwrap();
        // Phase A: reserve + write + seed-extend each sub-chunk into a range.
        let mut chunks: Vec<(ResidentRange, usize)> = Vec::new();
        for (k, chunk) in pv.chunks(sub).enumerate() {
            let mut range = cohort.reserve_pairs(chunk.len()).unwrap();
            cohort.write_pairs(&mut range, chunk).unwrap();
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut range)
                .unwrap();
            chunks.push((range, k * sub));
        }
        // Phase B: one cohort pestat.
        let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
        // Phase C: pair-emit each sub-chunk in order.
        let mut got = Vec::new();
        for (range, p0) in &mut chunks {
            let ids = IdBases {
                first_single_id: 0,
                first_pair_id: *p0 as u64,
            };
            let mut sink = RecordVec::default();
            cohort
                .pair_emit(
                    &idx,
                    &opts,
                    &mut scratch,
                    range,
                    Some(&pestat),
                    ids,
                    *p0,
                    &mut sink,
                )
                .unwrap();
            for (origin, body) in sink.records {
                let RecordOrigin::Pair(i) = origin else {
                    panic!("unexpected single")
                };
                got.push((i, body));
            }
        }
        assert_eq!(
            got, legacy,
            "sub={sub}: resident cohort diverged from align_batch"
        );
    }
}

/// A resident cohort is shared across a worker pool (`Send + Sync`), and a
/// range token moves between workers (`Send`).
#[test]
fn resident_cohort_is_send_and_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<ResidentCohort>();
    assert_send_sync::<ResidentRange>();
    let cohort = ResidentCohort::new(false).unwrap();
    // Move it into another thread and drop it there.
    std::thread::spawn(move || drop(cohort)).join().unwrap();
}

#[test]
fn regs_and_scratch_cross_threads() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let opts = Arc::new(MemOpts::new().unwrap());
    let (names, r1, r2) = pairs(32);
    let qual = vec![b'I'; 150];
    let pv: Vec<ReadPair<'_>> = (0..32)
        .map(|i| ReadPair {
            name_r1: names[i].as_bytes(),
            seq_r1: &r1[i],
            qual_r1: Some(&qual),
            name_r2: names[i].as_bytes(),
            seq_r2: &r2[i],
            qual_r2: Some(&qual),
        })
        .collect();
    // Seed on one thread, pair on another: AlnRegs is Send.
    let regs = std::thread::scope(|s| {
        s.spawn(|| {
            let mut sc = AlignScratch::new().unwrap();
            seed_extend(
                &idx,
                &opts,
                &mut sc,
                &ReadBatch {
                    pairs: &pv,
                    singles: &[],
                },
            )
            .unwrap()
        })
        .join()
        .unwrap()
    });
    assert_eq!(regs.n_pairs(), 32);
    assert!(regs.heap_bytes() > 32 * 2 * 300);
    let pestat = MemPeStat::infer_cohort(&idx, &opts, std::slice::from_ref(&regs)).unwrap();
    let n = std::thread::scope(|s| {
        s.spawn(move || {
            let mut sc = AlignScratch::new().unwrap();
            let mut sink = RecordVec::default();
            pair_emit(
                &idx,
                &opts,
                &mut sc,
                regs,
                Some(&pestat),
                IdBases::default(),
                &mut sink,
            )
            .unwrap();
            sink.records.len()
        })
        .join()
        .unwrap()
    });
    assert_eq!(n, 64);
}

#[test]
fn pair_emit_without_pestat_is_an_error() {
    let Some(idx) = shared_idx() else {
        return;
    };
    let opts = MemOpts::new().unwrap();
    let (names, r1, r2) = pairs(2);
    let qual = vec![b'I'; 150];
    let pv: Vec<ReadPair<'_>> = (0..2)
        .map(|i| ReadPair {
            name_r1: names[i].as_bytes(),
            seq_r1: &r1[i],
            qual_r1: Some(&qual),
            name_r2: names[i].as_bytes(),
            seq_r2: &r2[i],
            qual_r2: Some(&qual),
        })
        .collect();
    let mut sc = AlignScratch::new().unwrap();
    let regs = seed_extend(
        &idx,
        &opts,
        &mut sc,
        &ReadBatch {
            pairs: &pv,
            singles: &[],
        },
    )
    .unwrap();
    let err = pair_emit(
        &idx,
        &opts,
        &mut sc,
        regs,
        None,
        IdBases::default(),
        &mut RecordVec::default(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("pestat"), "{err}");
}

#[test]
fn load_with_threads_matches_load() {
    let Some(prefix) = ref_prefix() else {
        return;
    };
    let a = BwaIndex::load(prefix).unwrap();
    let b = BwaIndex::load_with_threads(prefix, 4).unwrap();
    assert_eq!(
        a.contigs().collect::<Vec<_>>(),
        b.contigs().collect::<Vec<_>>()
    );
}

/// Loading the `--meth` seed index with threads must not change it: the dual
/// index from [`BwaIndex::load_meth_with_threads`] aligns a `--meth` batch
/// byte for byte like the single-threaded [`BwaIndex::load_meth`] one.
#[test]
fn load_meth_with_threads_matches_load_meth() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let seed = suffixed(&meth.prefix, ".meth");
    let threaded = BwaIndex::load_meth_with_threads(&seed, &meth.prefix, 4).unwrap();
    assert!(threaded.is_meth());
    assert_eq!(
        threaded.contigs().collect::<Vec<_>>(),
        meth.idx.contigs().collect::<Vec<_>>()
    );
    assert_eq!(meth_records(&threaded), meth_records(&meth.idx));
}

/// The field fixture aligned as a `--meth` batch (bwa-mem3's `--meth`
/// defaults) against `idx`, as `(pair_idx, record bytes)`. Panics if the
/// fixture yields no records, so a comparison of two results is never vacuous.
fn meth_records(idx: &BwaIndex) -> Vec<(usize, Vec<u8>)> {
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let mut opts = MemOpts::new().unwrap();
    opts.set_pe(true)
        .set_meth(true)
        .apply_meth_defaults()
        .unwrap();
    let (batch, _) = align_batch(idx, &opts, &pv, None).unwrap();
    let records: Vec<(usize, Vec<u8>)> = batch
        .iter()
        .map(|r| (r.pair_idx, r.bytes.to_vec()))
        .collect();
    assert!(!records.is_empty(), "the meth fixture produced no records");
    records
}

/// A `bwa-mem3 index --meth --meth-only` PhiX index, in its own temp dir (so a
/// test may remove files from it): `(dir, original prefix)`. `None` (skip)
/// without a `bwa-mem3` that supports `--meth-only`.
fn phix_meth_only() -> Option<(tempfile::TempDir, std::path::PathBuf)> {
    if !supports_meth_only() {
        return None;
    }
    let (dir, prefix) =
        build_reference(&[("phix", &phix_bases())], &[], &["--meth", "--meth-only"])?;
    assert!(
        !suffixed(&prefix, ".bwt.2bit.64").exists(),
        "--meth-only should not build the original's FM-index"
    );
    Some((dir, prefix))
}

/// A `bwa-mem3 index --meth --meth-only` index has no FM-index for the original
/// reference, which [`BwaIndex::load_meth`] never reads: it loads, and aligns a
/// `--meth` batch byte for byte like the full dual index.
#[test]
fn load_meth_accepts_a_meth_only_index() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let Some((_dir, prefix)) = phix_meth_only() else {
        return;
    };
    let meth_only = BwaIndex::load_meth(suffixed(&prefix, ".meth"), &prefix).unwrap();
    assert!(meth_only.is_meth());
    assert_eq!(meth_records(&meth_only), meth_records(&meth.idx));
}

/// On a `--meth-only` index, each original-reference file `load_meth` still
/// needs (`.amb`, `.ann`, `.pac`) is checked up front: without it, `load_meth`
/// returns a recoverable [`Error::IndexLoad`] naming that file, instead of the
/// vendored loader exiting the process.
#[rstest]
#[case::amb(".amb")]
#[case::ann(".ann")]
#[case::pac(".pac")]
fn load_meth_names_a_missing_original_reference_file(#[case] suffix: &str) {
    let Some((_dir, prefix)) = phix_meth_only() else {
        return;
    };
    let file = suffixed(&prefix, suffix);
    std::fs::remove_file(&file).unwrap();
    let Err(Error::IndexLoad { msg, .. }) =
        BwaIndex::load_meth(suffixed(&prefix, ".meth"), &prefix)
    else {
        panic!("expected an IndexLoad error");
    };
    assert_eq!(
        msg,
        format!(
            "missing or unreadable original-reference file(s): {}",
            file.display()
        )
    );
}

/// A `--meth-context cpg` combination `bwa-mem3 mem` refuses fails at
/// alignment too, not only in `apply_meth_defaults`: here em-seq is set after
/// the defaults validated a TAPS + CpG combination.
#[test]
fn align_rejects_a_cpg_context_set_out_of_order() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let fixture = field_fixture();
    let (r1, r2) = &fixture.pairs[0];
    let pv = [ReadPair {
        name_r1: &r1.name,
        seq_r1: &r1.seq,
        qual_r1: r1.qual.as_deref(),
        name_r2: &r2.name,
        seq_r2: &r2.seq,
        qual_r2: r2.qual.as_deref(),
    }];
    let mut opts = MemOpts::new().unwrap();
    opts.set_pe(true)
        .set_meth(true)
        .set_meth_chem(bwa_mem3_rs::MethChem::Taps)
        .set_meth_context(bwa_mem3_rs::MethContext::Cpg)
        .apply_meth_defaults()
        .unwrap();
    assert!(align_batch(&meth.idx, &opts, &pv, None).is_ok());
    opts.set_meth_chem(bwa_mem3_rs::MethChem::EmSeq);
    let Err(err) = align_batch(&meth.idx, &opts, &pv, None) else {
        panic!("em-seq + cpg aligned");
    };
    assert!(
        matches!(&err, bwa_mem3_rs::Error::InvalidOpts(m) if m.contains("requires --meth=taps")),
        "{err}"
    );
}

/// A mismatch penalty the 8-bit kernels cannot store (`bwa-mem3 mem`
/// refuses anything outside `0..=128`) fails at alignment rather than wrapping
/// to a free mismatch.
#[test]
fn align_rejects_a_mismatch_penalty_the_kernels_cannot_store() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(4).map(as_read_pair).collect();
    let mut opts = MemOpts::new().unwrap();
    opts.set_mismatch_penalty(128);
    assert!(align_batch(&idx, &opts, &pv, None).is_ok());
    opts.set_mismatch_penalty(256);
    let Err(err) = align_batch(&idx, &opts, &pv, None) else {
        panic!("a mismatch penalty of 256 aligned");
    };
    assert!(
        matches!(&err, bwa_mem3_rs::Error::InvalidOpts(m) if m.contains("outside 0..128")),
        "{err}"
    );
}

/// The seeding entry points have no earlier snapshot to compare against, so
/// their own option check is all that keeps a penalty the kernels cannot store
/// (`-B 256` wraps to a free mismatch in 8 bits) from reaching them. Each one
/// refuses it, and a refused resident seed binds nothing.
#[test]
fn seeding_entry_points_reject_options_the_kernels_cannot_take() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(4).map(as_read_pair).collect();
    let batch = ReadBatch {
        pairs: &pv,
        singles: &[],
    };
    let mut bad = MemOpts::new().unwrap();
    bad.set_mismatch_penalty(256);
    let good = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();
    let refused = |what: &str, err: Option<bwa_mem3_rs::Error>| match err {
        Some(bwa_mem3_rs::Error::InvalidOpts(m)) => {
            assert!(m.contains("outside 0..128"), "{what}: {m}");
        }
        other => panic!("{what}: expected InvalidOpts, got {other:?}"),
    };

    refused("seed_batch", seed_batch(&idx, &bad, &pv).err());
    refused(
        "seed_extend",
        seed_extend(&idx, &bad, &mut sc, &batch).err(),
    );
    refused(
        "estimate_pestat",
        bwa_mem3_rs::estimate_pestat(&idx, &bad, &pv).err(),
    );

    let cohort = ResidentCohort::new(false).unwrap();
    let mut range = cohort.reserve_pairs(pv.len()).unwrap();
    cohort.write_pairs(&mut range, &pv).unwrap();
    refused(
        "resident seed_extend",
        cohort.seed_extend(&idx, &bad, &mut sc, &mut range).err(),
    );
    let first_pair = (range.first() / 2) as u64;
    let marked = [None, Some(first_pair), None, None];
    refused(
        "seed_extend_with_reps",
        cohort
            .seed_extend_with_reps(&idx, &bad, &mut sc, &mut range, &marked)
            .err(),
    );
    // Neither refusal bound `bad`.
    cohort
        .seed_extend(&idx, &good, &mut sc, &mut range)
        .unwrap();
}

/// `Seeds` and `AlnRegs` hold coordinates into the index they were seeded
/// against; a later phase given another index is refused rather than reading
/// that index's reference out of bounds (it used to abort the process).
#[test]
fn split_phases_reject_a_different_index() {
    let (Some(idx), Some(other)) = (shared_idx(), alt_reference().map(|r| r.idx.clone())) else {
        eprintln!("skip: bwa-mem3 not available to build the PhiX and ALT indexes");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(4).map(as_read_pair).collect();
    let batch = ReadBatch {
        pairs: &pv,
        singles: &[],
    };
    let opts = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();
    let needle = "seeded against a different BwaIndex";

    let seeds = seed_batch(&idx, &opts, &pv).unwrap();
    assert_invalid(
        extend_batch(&other, &opts, seeds, &pv, None).map(|_| ()),
        needle,
        "extend_batch",
    );

    let regs = seed_extend(&idx, &opts, &mut sc, &batch).unwrap();
    assert_invalid(
        MemPeStat::infer_cohort(&other, &opts, std::slice::from_ref(&regs)).map(|_| ()),
        needle,
        "infer_cohort",
    );
    let pestat = MemPeStat::infer_cohort(&idx, &opts, std::slice::from_ref(&regs)).unwrap();
    assert_invalid(
        pair_emit(
            &other,
            &opts,
            &mut sc,
            regs,
            Some(&pestat),
            IdBases::default(),
            &mut RecordVec::default(),
        ),
        needle,
        "pair_emit",
    );
}

/// A `resolve_memo` refused for options that differ from the cohort's leaves
/// the memo open and the marks in place, so a later range can still add a
/// mark, a retry under the cohort's own options resolves every mark, and the
/// cohort emits normally.
#[test]
fn resolve_memo_refused_for_changed_options_can_be_retried() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(8).map(as_read_pair).collect();
    let opts = MemOpts::new().unwrap();
    let mut other = MemOpts::new().unwrap();
    other.set_mismatch_penalty(5);
    let mut sc = AlignScratch::new().unwrap();

    let cohort = ResidentCohort::new(false).unwrap();
    let mut first = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut first, &pv[..4]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut first)
        .unwrap();
    // The second range repeats the first range's first pair and marks it so.
    let repeat = [pv[0], pv[5], pv[6], pv[7]];
    let mut second = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut second, &repeat).unwrap();
    let rep_of = (first.first() / 2) as u64;
    cohort
        .seed_extend_with_reps(
            &idx,
            &opts,
            &mut sc,
            &mut second,
            &[Some(rep_of), None, None, None],
        )
        .unwrap();

    let Err(err) = cohort.resolve_memo(&idx, &other, &mut sc) else {
        panic!("resolve_memo accepted options that differ from the cohort's");
    };
    assert!(
        matches!(&err, bwa_mem3_rs::Error::InvalidOpts(m)
            if m.starts_with("resolve_memo: the options differ")),
        "{err}"
    );
    // The refusal left the memo open: a later range can still add a mark.
    let mut third = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut third, &[pv[1], pv[2]]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut sc, &mut third, &[Some(rep_of + 1), None])
        .unwrap();
    let stats = cohort.resolve_memo(&idx, &opts, &mut sc).unwrap();
    // Both the mark made before the refusal and the one made after it resolve.
    assert_eq!(stats.dup_pairs, 2, "{stats:?}");
    assert_eq!(stats.copied, 2, "{stats:?}");
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    for (p0, mut r) in [(0, first), (4, second), (8, third)] {
        let records = emit_range(&cohort, &idx, &opts, &mut sc, &mut r, &pestat, p0);
        assert!(!records.is_empty(), "range at pair {p0} emitted nothing");
    }
}

/// An empty pair range emitted early blocks neither the memo nor the cohort
/// pestat: it holds no pairs to copy into or to model, so a mark added after
/// it still resolves, and `infer_cohort` still runs.
#[test]
fn an_emitted_empty_pair_range_blocks_neither_the_memo_nor_the_pestat() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(6).map(as_read_pair).collect();
    let opts = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();

    let cohort = ResidentCohort::new(false).unwrap();
    let mut empty = cohort.reserve_pairs(0).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut empty)
        .unwrap();
    let mut records = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut sc,
            &mut empty,
            None,
            IdBases::default(),
            0,
            &mut records,
        )
        .unwrap();
    assert!(records.records.is_empty());

    let mut first = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut first, &pv[..4]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut first)
        .unwrap();
    let mut second = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut second, &[pv[0], pv[5]]).unwrap();
    let rep_of = (first.first() / 2) as u64;
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut sc, &mut second, &[Some(rep_of), None])
        .unwrap();

    let stats = cohort.resolve_memo(&idx, &opts, &mut sc).unwrap();
    assert_eq!(stats.dup_pairs, 1, "{stats:?}");
    assert_eq!(stats.copied, 1, "{stats:?}");
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    for (p0, mut r) in [(0, first), (4, second)] {
        let records = emit_range(&cohort, &idx, &opts, &mut sc, &mut r, &pestat, p0);
        assert!(!records.is_empty(), "range at pair {p0} emitted nothing");
    }
}

/// `infer_cohort` with no marks: an empty pair range emitted before it does
/// not count as an emitted range the model would have to include.
#[test]
fn infer_cohort_runs_after_an_empty_pair_range_is_emitted() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(4).map(as_read_pair).collect();
    let opts = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();

    let cohort = ResidentCohort::new(false).unwrap();
    let mut range = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut range, &pv).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut range)
        .unwrap();
    let mut empty = cohort.reserve_pairs(0).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut empty)
        .unwrap();
    let mut records = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut sc,
            &mut empty,
            None,
            IdBases::default(),
            0,
            &mut records,
        )
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let records = emit_range(&cohort, &idx, &opts, &mut sc, &mut range, &pestat, 0);
    assert!(!records.is_empty());
}

/// `BwaIndex::bam_header` lays out the magic, the caller's text, and the
/// index's contigs in index order (the order records' `tid` indexes).
#[test]
fn bam_header_encodes_the_text_and_the_index_contigs() {
    let Some(r) = alt_reference() else {
        eprintln!("skip: bwa-mem3 not available to build the ALT index");
        return;
    };
    let text = "@HD\tVN:1.5\n@PG\tID:t\n";
    let b = r.idx.bam_header(text).unwrap();
    let u32_at = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
    assert_eq!(&b[..4], b"BAM\x01");
    assert_eq!(u32_at(4), text.len());
    assert_eq!(&b[8..8 + text.len()], text.as_bytes());
    let mut o = 8 + text.len();
    assert_eq!(u32_at(o), r.contigs.len());
    o += 4;
    for (name, bases) in &r.contigs {
        let l_name = u32_at(o);
        assert_eq!(&b[o + 4..o + 4 + l_name], format!("{name}\0").as_bytes());
        o += 4 + l_name;
        assert_eq!(u32_at(o), bases.len(), "{name}");
        o += 4;
    }
    assert_eq!(o, b.len(), "trailing bytes after the reference list");
}

/// A contig name that is not UTF-8 reads back from `contig_name` as `""`, so
/// `bam_header` must refuse it rather than write an empty reference name.
#[test]
fn bam_header_refuses_a_contig_name_that_is_not_utf8() {
    let Some(bwa) = find_bwa_mem3() else {
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("ref.fa");
    let mut body = b">chr\xff1\n".to_vec();
    for chunk in phix_bases()[..2000].chunks(72) {
        body.extend_from_slice(chunk);
        body.push(b'\n');
    }
    std::fs::write(&fa, body).unwrap();
    let status = std::process::Command::new(bwa)
        .arg("index")
        .arg(&fa)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "bwa-mem3 index failed");
    let idx = BwaIndex::load(&fa).unwrap();
    assert_eq!(idx.contig_name(0), "", "fixture name decoded as UTF-8");
    let Err(err) = idx.bam_header("@HD\tVN:1.5\n") else {
        panic!("bam_header wrote a contig name that is not UTF-8");
    };
    assert!(
        err.to_string()
            .contains("bam_header: contig 0's name is not UTF-8"),
        "{err}"
    );
}

/// The split-phase entry points re-check the options too: `Seeds`, `AlnRegs`
/// and a resident range do not borrow `MemOpts`, so a caller can change them
/// between phases. Em-seq set after a valid TAPS + CpG seeding must fail every
/// later phase with `InvalidOpts` (and, for a resident range, leave the range
/// emittable once the options are put back).
#[test]
fn split_phases_reject_a_cpg_context_set_between_phases() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(8).map(as_read_pair).collect();
    let taps_cpg = || {
        let mut o = MemOpts::new().unwrap();
        o.set_pe(true)
            .set_meth(true)
            .set_meth_chem(bwa_mem3_rs::MethChem::Taps)
            .set_meth_context(bwa_mem3_rs::MethContext::Cpg)
            .apply_meth_defaults()
            .unwrap();
        o
    };
    let valid = taps_cpg();
    let mut invalid = taps_cpg();
    invalid.set_meth_chem(bwa_mem3_rs::MethChem::EmSeq);
    let is_cpg_err = |err: &bwa_mem3_rs::Error| matches!(err, bwa_mem3_rs::Error::InvalidOpts(m) if m.contains("requires --meth=taps"));
    let batch = ReadBatch {
        pairs: &pv,
        singles: &[],
    };

    // seed_batch -> extend_batch.
    let seeds = seed_batch(&meth.idx, &valid, &pv).unwrap();
    let Err(err) = extend_batch(&meth.idx, &invalid, seeds, &pv, None) else {
        panic!("extend_batch accepted em-seq + cpg");
    };
    assert!(is_cpg_err(&err), "{err}");
    let seeds = seed_batch(&meth.idx, &valid, &pv).unwrap();
    assert!(extend_batch(&meth.idx, &valid, seeds, &pv, None).is_ok());

    // seed_extend -> infer_cohort -> pair_emit.
    let mut sc = AlignScratch::new().unwrap();
    let regs = seed_extend(&meth.idx, &valid, &mut sc, &batch).unwrap();
    let Err(err) = MemPeStat::infer_cohort(&meth.idx, &invalid, std::slice::from_ref(&regs)) else {
        panic!("infer_cohort accepted em-seq + cpg");
    };
    assert!(is_cpg_err(&err), "{err}");
    let pestat = MemPeStat::infer_cohort(&meth.idx, &valid, std::slice::from_ref(&regs)).unwrap();
    let Err(err) = pair_emit(
        &meth.idx,
        &invalid,
        &mut sc,
        regs,
        Some(&pestat),
        IdBases::default(),
        &mut RecordVec::default(),
    ) else {
        panic!("pair_emit accepted em-seq + cpg");
    };
    assert!(is_cpg_err(&err), "{err}");

    // Resident cohort: infer_cohort and pair_emit, and the range survives.
    let cohort = ResidentCohort::new(true).unwrap();
    let mut range = cohort.reserve_pairs(pv.len()).unwrap();
    cohort.write_pairs(&mut range, &pv).unwrap();
    cohort
        .seed_extend(&meth.idx, &valid, &mut sc, &mut range)
        .unwrap();
    let Err(err) = cohort.infer_cohort(&meth.idx, &invalid) else {
        panic!("resident infer_cohort accepted em-seq + cpg");
    };
    assert!(is_cpg_err(&err), "{err}");
    let pestat = cohort.infer_cohort(&meth.idx, &valid).unwrap();
    let Err(err) = cohort.pair_emit(
        &meth.idx,
        &invalid,
        &mut sc,
        &mut range,
        Some(&pestat),
        IdBases::default(),
        0,
        &mut RecordVec::default(),
    ) else {
        panic!("resident pair_emit accepted em-seq + cpg");
    };
    assert!(is_cpg_err(&err), "{err}");
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(
            &meth.idx,
            &valid,
            &mut sc,
            &mut range,
            Some(&pestat),
            IdBases::default(),
            0,
            &mut sink,
        )
        .unwrap();
    assert!(
        !sink.records.is_empty(),
        "the meth fixture produced no records"
    );
}

/// A *valid* option change between phases is rejected too: regions seeded
/// under `--meth-context all` must not be paired, rescued or emitted under
/// `cpg` (CodeRabbit's All -> Cpg case), since the records would mix the two
/// scorings. The MD mode, which only shapes the tags written, may change.
#[test]
fn split_phases_reject_options_changed_between_phases() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(8).map(as_read_pair).collect();
    let taps = |ctx: bwa_mem3_rs::MethContext| {
        let mut o = MemOpts::new().unwrap();
        o.set_pe(true)
            .set_meth(true)
            .set_meth_chem(bwa_mem3_rs::MethChem::Taps)
            .set_meth_context(ctx)
            .apply_meth_defaults()
            .unwrap();
        o
    };
    let all = taps(bwa_mem3_rs::MethContext::All);
    let cpg = taps(bwa_mem3_rs::MethContext::Cpg);
    let mut all_md_none = taps(bwa_mem3_rs::MethContext::All);
    all_md_none.set_md_mode(bwa_mem3_rs::MdMode::None);
    let changed = |err: &bwa_mem3_rs::Error, what: &str| {
        matches!(err, bwa_mem3_rs::Error::InvalidOpts(m)
            if m.starts_with(&format!("{what}: the options differ")))
    };
    // The MD mode really does reach the records: every mapped record carries
    // MD:Z under FULL, and no record of any kind carries MD under NONE.
    let assert_md = |bodies: &[&[u8]], want: bool, what: &str| {
        let mapped: Vec<&&[u8]> = bodies.iter().filter(|b| has_cigar(b)).collect();
        assert!(!mapped.is_empty(), "{what}: no mapped record");
        if want {
            assert!(
                mapped.iter().all(|b| aux_type(b, b"MD") == Some(b'Z')),
                "{what}: expected MD:Z on every mapped record"
            );
        } else {
            assert!(
                bodies.iter().all(|b| aux_type(b, b"MD").is_none()),
                "{what}: expected no MD tag on any record"
            );
        }
    };
    let batch = ReadBatch {
        pairs: &pv,
        singles: &[],
    };

    // seed_batch -> extend_batch.
    let seeds = seed_batch(&meth.idx, &all, &pv).unwrap();
    let Err(err) = extend_batch(&meth.idx, &cpg, seeds, &pv, None) else {
        panic!("extend_batch accepted all -> cpg");
    };
    assert!(changed(&err, "extend_batch"), "{err}");
    for (opts, want) in [(&all, true), (&all_md_none, false)] {
        let seeds = seed_batch(&meth.idx, &all, &pv).unwrap();
        let (out, _) = extend_batch(&meth.idx, opts, seeds, &pv, None).unwrap();
        let bodies: Vec<&[u8]> = out.iter().map(|r| &r.bytes[4..]).collect();
        assert_md(&bodies, want, "extend_batch");
    }

    // seed_extend -> infer_cohort -> pair_emit.
    let mut sc = AlignScratch::new().unwrap();
    let regs = seed_extend(&meth.idx, &all, &mut sc, &batch).unwrap();
    let Err(err) = MemPeStat::infer_cohort(&meth.idx, &cpg, std::slice::from_ref(&regs)) else {
        panic!("infer_cohort accepted all -> cpg");
    };
    assert!(changed(&err, "infer_cohort"), "{err}");
    let pestat = MemPeStat::infer_cohort(&meth.idx, &all, std::slice::from_ref(&regs)).unwrap();
    let Err(err) = pair_emit(
        &meth.idx,
        &cpg,
        &mut sc,
        regs,
        Some(&pestat),
        IdBases::default(),
        &mut RecordVec::default(),
    ) else {
        panic!("pair_emit accepted all -> cpg");
    };
    assert!(changed(&err, "pair_emit"), "{err}");
    for (opts, want) in [(&all, true), (&all_md_none, false)] {
        let regs = seed_extend(&meth.idx, &all, &mut sc, &batch).unwrap();
        let mut sink = RecordVec::default();
        pair_emit(
            &meth.idx,
            opts,
            &mut sc,
            regs,
            Some(&pestat),
            IdBases::default(),
            &mut sink,
        )
        .unwrap();
        let bodies: Vec<&[u8]> = sink.records.iter().map(|(_, b)| b.as_slice()).collect();
        assert_md(&bodies, want, "pair_emit");
    }

    // Resident cohort: the first seed_extend binds the options.
    let cohort = ResidentCohort::new(true).unwrap();
    let mut ranges: Vec<ResidentRange> = pv
        .chunks(4)
        .map(|c| {
            let mut r = cohort.reserve_pairs(c.len()).unwrap();
            cohort.write_pairs(&mut r, c).unwrap();
            r
        })
        .collect();
    cohort
        .seed_extend(&meth.idx, &all, &mut sc, &mut ranges[0])
        .unwrap();
    let Err(err) = cohort.seed_extend(&meth.idx, &cpg, &mut sc, &mut ranges[1]) else {
        panic!("resident seed_extend accepted all -> cpg on a second range");
    };
    assert!(changed(&err, "seed_extend"), "{err}");
    cohort
        .seed_extend(&meth.idx, &all, &mut sc, &mut ranges[1])
        .unwrap();
    for (what, res) in [
        (
            "resolve_memo",
            cohort.resolve_memo(&meth.idx, &cpg, &mut sc).map(|_| ()),
        ),
        (
            "infer_cohort",
            cohort.infer_cohort(&meth.idx, &cpg).map(|_| ()),
        ),
    ] {
        let Err(err) = res else {
            panic!("resident {what} accepted all -> cpg");
        };
        assert!(changed(&err, what), "{err}");
    }
    let pestat = cohort.infer_cohort(&meth.idx, &all).unwrap();
    let Err(err) = cohort.pair_emit(
        &meth.idx,
        &cpg,
        &mut sc,
        &mut ranges[0],
        Some(&pestat),
        IdBases::default(),
        0,
        &mut RecordVec::default(),
    ) else {
        panic!("resident pair_emit accepted all -> cpg");
    };
    assert!(changed(&err, "pair_emit"), "{err}");
    // One range under FULL, the other under NONE: the MD mode may change.
    assert_eq!(ranges.len(), 2, "the zip below needs one range per MD mode");
    for (range, (opts, want)) in ranges.iter_mut().zip([(&all, true), (&all_md_none, false)]) {
        let mut sink = RecordVec::default();
        cohort
            .pair_emit(
                &meth.idx,
                opts,
                &mut sc,
                range,
                Some(&pestat),
                IdBases::default(),
                0,
                &mut sink,
            )
            .unwrap();
        let bodies: Vec<&[u8]> = sink.records.iter().map(|(_, b)| b.as_slice()).collect();
        assert_md(&bodies, want, "resident pair_emit");
    }
}

/// A seed-extend that fails binds no options: a native failure (an unwritten
/// range) after the tentative binding, through `seed_extend` and through
/// `seed_extend_with_reps` (with a valid mark, so every input check passes and
/// the call reaches the binding and the marks), keeps no later seed under
/// other options out, and the failed marked seed leaves no mark counted. And
/// `set_pe`, which the shim overrides per read group in every phase, may
/// differ between phases.
#[test]
fn failed_seeds_bind_no_options_and_pe_may_differ() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(8).map(as_read_pair).collect();
    let mut other = MemOpts::new().unwrap();
    other.set_mismatch_penalty(5);
    let opts = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();

    let cohort = ResidentCohort::new(false).unwrap();
    let mut unwritten = cohort.reserve_pairs(4).unwrap();
    let Err(err) = cohort.seed_extend(&idx, &other, &mut sc, &mut unwritten) else {
        panic!("seeding an unwritten range succeeded");
    };
    // The native call's refusal, reached after the tentative binding.
    assert!(
        err.to_string()
            .contains("resident_seed_extend: resident range lifecycle violated"),
        "{err}"
    );
    let mut unwritten_marked = cohort.reserve_pairs(4).unwrap();
    let first_pair = (unwritten_marked.first() / 2) as u64;
    let marked = [None, Some(first_pair), None, None];
    let Err(err) =
        cohort.seed_extend_with_reps(&idx, &other, &mut sc, &mut unwritten_marked, &marked)
    else {
        panic!("seeding an unwritten marked range succeeded");
    };
    assert!(
        err.to_string()
            .contains("resident_seed_extend_reps: resident range lifecycle violated"),
        "{err}"
    );
    // Neither failure bound `other`, so the default options seed fine.
    let mut range = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut range, &pv[..4]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut range)
        .unwrap();
    // The failed marked seed took its mark back: with none counted,
    // `resolve_memo` stays out of the native path, which would refuse the two
    // unextended ranges.
    let stats = cohort.resolve_memo(&idx, &opts, &mut sc).unwrap();
    assert_eq!(stats.dup_pairs, 0, "{stats:?}");

    let mut se = MemOpts::new().unwrap();
    se.set_pe(false);
    let mut pe = MemOpts::new().unwrap();
    pe.set_pe(true);
    let seeds = seed_batch(&idx, &se, &pv).unwrap();
    extend_batch(&idx, &pe, seeds, &pv, None).unwrap();
}

/// An empty range seeds and emits nothing, so it binds no options either way:
/// seeding it neither fixes the cohort's options for later ranges nor is held
/// to options an earlier range bound, and emitting it is not held to them
/// either.
#[test]
fn empty_range_seeds_bind_no_options() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().take(4).map(as_read_pair).collect();
    let mut other = MemOpts::new().unwrap();
    other.set_mismatch_penalty(5);
    let opts = MemOpts::new().unwrap();
    let mut sc = AlignScratch::new().unwrap();

    let cohort = ResidentCohort::new(false).unwrap();
    let mut empty = cohort.reserve_pairs(0).unwrap();
    cohort
        .seed_extend(&idx, &other, &mut sc, &mut empty)
        .unwrap();
    let mut range = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut range, &pv).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut sc, &mut range)
        .unwrap();
    // `opts` is now bound, and an empty range under `other` still seeds.
    let mut empty_after = cohort.reserve_singles(0).unwrap();
    cohort
        .seed_extend(&idx, &other, &mut sc, &mut empty_after)
        .unwrap();
    // Both empty ranges emit under `other`, through both sinks; the bound
    // range is still refused under it, and emits under `opts`.
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let mut records = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &other,
            &mut sc,
            &mut empty,
            Some(&pestat),
            IdBases::default(),
            0,
            &mut records,
        )
        .unwrap();
    let mut fields = CountingFieldSink::default();
    cohort
        .pair_emit_fields(
            &idx,
            &other,
            &mut sc,
            &mut empty_after,
            None,
            IdBases::default(),
            0,
            &mut fields,
        )
        .unwrap();
    assert!(records.records.is_empty());
    assert_eq!(fields.0, 0, "an empty range emitted structured records");
    let Err(err) = cohort.pair_emit(
        &idx,
        &other,
        &mut sc,
        &mut range,
        Some(&pestat),
        IdBases::default(),
        0,
        &mut records,
    ) else {
        panic!("a range emitted under options other than the cohort's");
    };
    assert!(
        err.to_string().contains("pair_emit: the options differ"),
        "{err}"
    );
    assert!(
        records.records.is_empty(),
        "a rejected pair_emit wrote records"
    );
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut sc,
            &mut range,
            Some(&pestat),
            IdBases::default(),
            0,
            &mut records,
        )
        .unwrap();
    assert!(!records.records.is_empty());

    // A non-empty range under `other` is still refused.
    let mut refused = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut refused, &pv).unwrap();
    let Err(err) = cohort.seed_extend(&idx, &other, &mut sc, &mut refused) else {
        panic!("a range seeded under options other than the cohort's");
    };
    assert!(
        err.to_string().contains("seed_extend: the options differ"),
        "{err}"
    );
}

/// A sink whose `emit` panics on its first call. Used to prove that a panic
/// inside `RecordSink::emit` propagates out of `pair_emit` as an ordinary
/// panic rather than aborting the process or corrupting state -- the
/// `sink_trampoline` catch/resume must be sound across the FFI boundary.
struct PanicOnFirstEmit {
    calls: usize,
}

impl RecordSink for PanicOnFirstEmit {
    fn emit(&mut self, _origin: RecordOrigin, _body: &[u8]) {
        self.calls += 1;
        panic!("boom from PanicOnFirstEmit");
    }
}

#[test]
fn pair_emit_propagates_a_panicking_sink_as_an_ordinary_panic() {
    let Some(idx) = shared_idx() else {
        return;
    };
    let opts = MemOpts::new().unwrap();
    let (names, r1, r2) = pairs(4);
    let qual = vec![b'I'; 150];
    let pv: Vec<ReadPair<'_>> = (0..4)
        .map(|i| ReadPair {
            name_r1: names[i].as_bytes(),
            seq_r1: &r1[i],
            qual_r1: Some(&qual),
            name_r2: names[i].as_bytes(),
            seq_r2: &r2[i],
            qual_r2: Some(&qual),
        })
        .collect();
    let mut sc = AlignScratch::new().unwrap();
    let regs = seed_extend(
        &idx,
        &opts,
        &mut sc,
        &ReadBatch {
            pairs: &pv,
            singles: &[],
        },
    )
    .unwrap();
    let pestat = MemPeStat::infer_cohort(&idx, &opts, std::slice::from_ref(&regs)).unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut sink = PanicOnFirstEmit { calls: 0 };
        pair_emit(
            &idx,
            &opts,
            &mut sc,
            regs,
            Some(&pestat),
            IdBases::default(),
            &mut sink,
        )
    }));
    let err = result.expect_err("panicking sink must unwind out of pair_emit");
    let msg = err
        .downcast_ref::<&str>()
        .map(|s| (*s).to_string())
        .or_else(|| err.downcast_ref::<String>().cloned())
        .unwrap_or_default();
    assert!(msg.contains("boom from PanicOnFirstEmit"), "{msg}");
}

/// Every record the structured-fields path reports, re-serialized, is
/// byte-identical to the packed body the packed path emits for the same input,
/// in the same order, with the same origin -- and `mate`/`is_primary` agree
/// with the record's own FLAG. Run under default, `-Y` and `-M` options at two
/// sub-chunk sizes, and each restructured record shape must actually occur.
#[test]
fn pair_emit_fields_matches_packed_pair_emit_byte_for_byte() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let mut coverage = Coverage::default();

    let default = MemOpts::new().unwrap();
    let mut soft_clip = MemOpts::new().unwrap();
    soft_clip.set_soft_clip_supplementary(true);
    let mut split_secondary = MemOpts::new().unwrap();
    split_secondary.set_mark_split_secondary(true);
    // RG:Z is covered in its own test binary (fields_rg.rs):
    // `set_read_group_id` writes bwa's process-global `bwa_rg_id`, which would
    // leak into every other test in this one.
    for (label, opts) in [
        ("default", default),
        ("-Y", soft_clip),
        ("-M", split_secondary),
    ] {
        for sub in [7usize, fixture.pairs.len()] {
            assert_fields_match_packed(&idx, &opts, &fixture, sub, label, &mut coverage);
        }
    }

    eprintln!("structured-sink coverage: {coverage:?}");
    for (what, n) in [
        ("reverse-strand", coverage.reverse),
        ("supplementary (0x800)", coverage.supplementary),
        ("-M split secondary (0x100)", coverage.split_secondary),
        ("hard-clipped SEQ window", coverage.hard_clipped_window),
        ("SA:Z", coverage.sa),
        ("MC:Z", coverage.mc),
        (
            "half-mapped unmapped read placed at its mate",
            coverage.placed_unmapped,
        ),
        ("unplaced unmapped", coverage.unmapped_unplaced),
        ("single-end", coverage.singles),
    ] {
        assert!(n > 0, "fixture never produced a {what} record");
    }
}

/// The same parity over a reference with an ALT contig, where reads hit both
/// the primary and the ALT copy: covers `XA:Z` and `pa:f`, which PhiX alone
/// never produces.
#[test]
fn pair_emit_fields_matches_packed_on_alt_hits() {
    let Some(alt) = alt_reference() else {
        eprintln!("skip: bwa-mem3 not available to build the ALT reference");
        return;
    };
    // Draw reads from the primary with the ALT copy's substitutions spliced into
    // the duplicated stretch: those reads match the ALT contig better, which is
    // what makes upstream record an ALT competitor (`alt_sc`) on the primary
    // hit and emit `pa:f` (bwamem.cpp mem_mark_primary_se).
    let mut source = alt.contigs[0].1.clone();
    source[500..1500].copy_from_slice(&alt.contigs[1].1);
    let fixture = field_fixture_from(&source);
    let opts = MemOpts::new().unwrap();
    let mut coverage = Coverage::default();
    for sub in [7usize, fixture.pairs.len()] {
        assert_fields_match_packed(&alt.idx, &opts, &fixture, sub, "alt", &mut coverage);
    }
    eprintln!("ALT coverage: {coverage:?}");
    assert!(coverage.xa > 0, "fixture never produced an XA:Z record");
    assert!(coverage.pa > 0, "fixture never produced a pa:f record");
}

/// The same parity under `--meth`: SEQ carries the original (unprojected)
/// bases, and every mapped record reports its XR/XG/XM tags.
#[test]
fn pair_emit_fields_matches_packed_under_meth() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let fixture = field_fixture();
    let mut opts = MemOpts::new().unwrap();
    opts.set_meth(true);
    let mut coverage = Coverage::default();
    for sub in [7usize, fixture.pairs.len()] {
        assert_fields_match_packed(&meth.idx, &opts, &fixture, sub, "meth", &mut coverage);
    }
    eprintln!("meth coverage: {coverage:?}");
    assert!(coverage.meth_mapped > 0, "no record carried XG:Z and XM:Z");
    assert!(coverage.reverse > 0, "no reverse-strand meth record");
}

/// A panicking structured sink surfaces as an ordinary panic from
/// `pair_emit_fields`, exactly like the packed sink (no unwind crosses FFI).
#[test]
fn pair_emit_fields_propagates_a_panicking_sink() {
    struct Boom;
    impl AlignedFieldsSink for Boom {
        fn emit(&mut self, _: RecordOrigin, _: Mate, _: bool, _: AlignedFields<'_>) {
            panic!("boom from structured sink");
        }
    }
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut range = cohort.reserve_singles(4).unwrap();
    for (i, r) in fixture.singles[..4].iter().enumerate() {
        let read = SingleRead {
            name: &r.name,
            seq: &r.seq,
            qual: r.qual.as_deref(),
        };
        cohort.write_single(&mut range, i, read).unwrap();
    }
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        cohort.pair_emit_fields(
            &idx,
            &opts,
            &mut scratch,
            &mut range,
            None,
            IdBases::default(),
            0,
            &mut Boom,
        )
    }));
    let payload = result.expect_err("a panicking sink must surface as a panic");
    assert_eq!(
        payload.downcast_ref::<&str>(),
        Some(&"boom from structured sink")
    );
}

/// The structured-sink types are borrowed POD: shareable and sendable.
#[test]
fn aligned_fields_are_send_and_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<AlignedFields<'static>>();
    assert_send_sync::<bwa_mem3_rs::MethFields<'static>>();
}

/// Many threads share one `Arc<ResidentCohort>`, each writing, seed-extending
/// and (after the pestat barrier) emitting only its own range, all released at
/// once by a barrier. The result must equal a serial run over the same reads:
/// concurrent disjoint-range access through the shared cohort is sound and
/// changes nothing.
#[test]
fn concurrent_disjoint_ranges_match_serial() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let (_, reference) = serial_reference(&idx, &opts, &fixture);

    let chunk = 10;
    let cohort = Arc::new(ResidentCohort::new(false).unwrap());
    // Reserve serially (the pipeline's Serial prepare step), then fan out.
    let ranges: Vec<(ResidentRange, usize)> = fixture
        .pairs
        .chunks(chunk)
        .enumerate()
        .map(|(k, _)| {
            let n = chunk.min(fixture.pairs.len() - k * chunk);
            (cohort.reserve_pairs(n).unwrap(), k * chunk)
        })
        .collect();

    // Phase A, concurrently: write + seed_extend each range on its own thread.
    let start = std::sync::Barrier::new(ranges.len());
    let ranges: Vec<(ResidentRange, usize)> = std::thread::scope(|s| {
        let handles: Vec<_> = ranges
            .into_iter()
            .map(|(mut range, p0)| {
                let (cohort, idx, opts, fixture, start) = (&cohort, &idx, &opts, &fixture, &start);
                s.spawn(move || {
                    let mut scratch = AlignScratch::new().unwrap();
                    start.wait();
                    for (i, p) in fixture.pairs[p0..p0 + range.n_pairs()].iter().enumerate() {
                        cohort.write_pair(&mut range, i, as_read_pair(p)).unwrap();
                    }
                    cohort
                        .seed_extend(idx, opts, &mut scratch, &mut range)
                        .unwrap();
                    (range, p0)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    // The barrier: one cohort pestat over every range.
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();

    // Phase C, concurrently: emit each range on its own thread.
    let start = std::sync::Barrier::new(ranges.len());
    let mut emitted: Vec<(usize, Records)> = std::thread::scope(|s| {
        let handles: Vec<_> = ranges
            .into_iter()
            .map(|(mut range, p0)| {
                let (cohort, idx, opts, pestat, start) = (&cohort, &idx, &opts, &pestat, &start);
                s.spawn(move || {
                    let mut scratch = AlignScratch::new().unwrap();
                    start.wait();
                    let records =
                        emit_range(cohort, idx, opts, &mut scratch, &mut range, pestat, p0);
                    (p0, records)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    emitted.sort_by_key(|(p0, _)| *p0);
    let got: Records = emitted.into_iter().flat_map(|(_, r)| r).collect();
    assert_eq!(
        got.len(),
        reference.len(),
        "record count differs from serial"
    );
    assert!(
        got == reference,
        "concurrent ranges diverged from the serial run"
    );
}

/// A range's reads are never moved by later reserves, even while the range is
/// mid-`seed_extend`/`pair_emit` on another thread and the cohort's segment
/// table is being regrown underneath it.
#[test]
fn reserves_during_in_flight_range_calls_never_move_a_range() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    // Take the pestat from the serial reference: this cohort will also hold
    // the extra ranges below, which must not feed its insert-size model.
    let (pestat, reference) = serial_reference(&idx, &opts, &fixture);

    // Whether the reserves actually overlap the worker's calls depends on
    // scheduling, so retry on a fresh cohort until one run overlaps enough to
    // regrow the table (starting at 8 slots and doubling, > 16 reserves means
    // at least two regrowths underneath the in-flight range). Every attempt's
    // output must match regardless.
    const ATTEMPTS: usize = 5;
    let mut overlapped = 0;
    for attempt in 0..ATTEMPTS {
        let cohort = ResidentCohort::new(false).unwrap();
        let mut range = written_pair_range(&cohort, &fixture.pairs);
        let done = std::sync::atomic::AtomicBool::new(false);
        let start = std::sync::Barrier::new(2);
        let (records, reserves) = std::thread::scope(|s| {
            let worker = s.spawn(|| {
                let mut scratch = AlignScratch::new().unwrap();
                start.wait();
                cohort
                    .seed_extend(&idx, &opts, &mut scratch, &mut range)
                    .unwrap();
                let records =
                    emit_range(&cohort, &idx, &opts, &mut scratch, &mut range, &pestat, 0);
                done.store(true, std::sync::atomic::Ordering::Release);
                records
            });
            // Grow the same pair region (and write into the new ranges) for as
            // long as the worker is busy, forcing table reallocations.
            start.wait();
            let mut extra = Vec::new();
            while !done.load(std::sync::atomic::Ordering::Acquire) {
                let mut r = cohort.reserve_pairs(1).unwrap();
                cohort
                    .write_pair(&mut r, 0, as_read_pair(&fixture.pairs[extra.len() % 8]))
                    .unwrap();
                extra.push(r);
            }
            (worker.join().unwrap(), extra.len())
        });
        assert!(
            records == reference,
            "attempt {attempt}: range output changed while the cohort grew \
             ({reserves} concurrent reserves)"
        );
        overlapped = reserves;
        if reserves > 16 {
            break;
        }
    }
    assert!(
        overlapped > 16,
        "no attempt overlapped the in-flight calls enough to regrow the table \
         (last: {overlapped} reserves)"
    );
}

/// A resident call that must be rejected, run against fresh state.
#[derive(Debug, Clone, Copy)]
enum Misuse {
    ForeignWritePair,
    ForeignWriteSingle,
    ForeignSeedExtend,
    ForeignPairEmit,
    ForeignPairEmitFields,
    SingleIntoPairRange,
    PairIntoSingleRange,
    IndexPastRange,
    NameTooLong,
    SeedExtendUnwritten,
    WriteSlotTwice,
    WriteAfterSeedExtend,
    SeedExtendTwice,
    SeedExtendOtherMeth,
    SeedExtendOtherIndex,
    PairEmitBeforeSeedExtend,
    PairEmitFieldsBeforeSeedExtend,
    PairEmitWithoutPestat,
    PairEmitOtherMeth,
    PairEmitFieldsOtherMeth,
    PairEmitTwice,
    InferBeforeSeedExtend,
    InferOtherMeth,
    InferAfterPairEmit,
    WritePairsPartial,
    WriteSinglesPartial,
    WritePairsAfterSlotWrite,
    PairsIntoSingleRange,
    SinglesIntoPairRange,
    WritePairAfterWritePairs,
    WriteSingleAfterWriteSingles,
    WritePairsTwice,
    PairEmitFieldsAfterPairEmit,
    SeedExtendAfterPairEmit,
}

/// Every lifecycle or contract misuse is reported as the `InvalidInput` naming
/// that misuse, rather than reaching the C kernels with a range in the wrong
/// state. Each case runs against its own fresh cohort.
#[rstest]
#[case::foreign_write_pair(Misuse::ForeignWritePair, "different ResidentCohort")]
#[case::foreign_write_single(Misuse::ForeignWriteSingle, "different ResidentCohort")]
#[case::foreign_seed_extend(Misuse::ForeignSeedExtend, "different ResidentCohort")]
#[case::foreign_pair_emit(Misuse::ForeignPairEmit, "different ResidentCohort")]
#[case::foreign_pair_emit_fields(Misuse::ForeignPairEmitFields, "different ResidentCohort")]
#[case::single_into_pair_range(Misuse::SingleIntoPairRange, "not a single range")]
#[case::pair_into_single_range(Misuse::PairIntoSingleRange, "not a pair range")]
#[case::index_past_range(Misuse::IndexPastRange, "out of range")]
#[case::name_too_long(Misuse::NameTooLong, "read name is 255 bytes")]
#[case::seed_extend_unwritten(Misuse::SeedExtendUnwritten, LIFECYCLE)]
#[case::write_slot_twice(Misuse::WriteSlotTwice, LIFECYCLE)]
#[case::write_after_seed_extend(Misuse::WriteAfterSeedExtend, LIFECYCLE)]
#[case::seed_extend_twice(Misuse::SeedExtendTwice, LIFECYCLE)]
#[case::seed_extend_other_meth(Misuse::SeedExtendOtherMeth, "meth=")]
#[case::seed_extend_other_index(Misuse::SeedExtendOtherIndex, "different BwaIndex")]
#[case::pair_emit_before_seed_extend(Misuse::PairEmitBeforeSeedExtend, LIFECYCLE)]
#[case::pair_emit_fields_before_seed_extend(Misuse::PairEmitFieldsBeforeSeedExtend, LIFECYCLE)]
#[case::pair_emit_without_pestat(Misuse::PairEmitWithoutPestat, "requires the cohort pestat")]
#[case::pair_emit_other_meth(Misuse::PairEmitOtherMeth, "meth=")]
#[case::pair_emit_fields_other_meth(Misuse::PairEmitFieldsOtherMeth, "meth=")]
#[case::pair_emit_twice(Misuse::PairEmitTwice, LIFECYCLE)]
#[case::infer_before_seed_extend(Misuse::InferBeforeSeedExtend, LIFECYCLE)]
#[case::infer_other_meth(Misuse::InferOtherMeth, "meth=")]
#[case::infer_after_pair_emit(Misuse::InferAfterPairEmit, LIFECYCLE)]
#[case::write_pairs_partial(Misuse::WritePairsPartial, "for a range of")]
#[case::write_singles_partial(Misuse::WriteSinglesPartial, "for a range of")]
#[case::write_pairs_after_slot_write(Misuse::WritePairsAfterSlotWrite, LIFECYCLE)]
#[case::pairs_into_single_range(Misuse::PairsIntoSingleRange, "not a pair range")]
#[case::singles_into_pair_range(Misuse::SinglesIntoPairRange, "not a single range")]
#[case::write_pair_after_write_pairs(Misuse::WritePairAfterWritePairs, LIFECYCLE)]
#[case::write_single_after_write_singles(Misuse::WriteSingleAfterWriteSingles, LIFECYCLE)]
#[case::write_pairs_twice(Misuse::WritePairsTwice, LIFECYCLE)]
#[case::pair_emit_fields_after_pair_emit(Misuse::PairEmitFieldsAfterPairEmit, LIFECYCLE)]
#[case::seed_extend_after_pair_emit(Misuse::SeedExtendAfterPairEmit, LIFECYCLE)]
fn resident_misuse_is_rejected(#[case] misuse: Misuse, #[case] needle: &str) {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let mut meth_opts = MemOpts::new().unwrap();
    meth_opts.set_meth(true);
    let mut scratch = AlignScratch::new().unwrap();
    let pair = as_read_pair(&fixture.pairs[0]);
    let single = SingleRead {
        name: &fixture.singles[0].name,
        seq: &fixture.singles[0].seq,
        qual: fixture.singles[0].qual.as_deref(),
    };
    let cohort = ResidentCohort::new(false).unwrap();
    let ids = IdBases::default();
    let mut sink = RecordVec::default();
    // A two-pair range, written in full.
    let written = || {
        let mut r = cohort.reserve_pairs(2).unwrap();
        cohort.write_pairs(&mut r, &[pair, pair]).unwrap();
        r
    };
    let foreign = || {
        ResidentCohort::new(false)
            .unwrap()
            .reserve_pairs(1)
            .unwrap()
    };

    let result: bwa_mem3_rs::Result<()> = match misuse {
        Misuse::WritePairsPartial => {
            cohort.write_pairs(&mut cohort.reserve_pairs(2).unwrap(), &[pair])
        }
        Misuse::WriteSinglesPartial => {
            cohort.write_singles(&mut cohort.reserve_singles(2).unwrap(), &[single])
        }
        Misuse::PairsIntoSingleRange => {
            cohort.write_pairs(&mut cohort.reserve_singles(1).unwrap(), &[pair])
        }
        Misuse::SinglesIntoPairRange => {
            cohort.write_singles(&mut cohort.reserve_pairs(1).unwrap(), &[single])
        }
        Misuse::WritePairAfterWritePairs => {
            let mut r = cohort.reserve_pairs(2).unwrap();
            cohort.write_pairs(&mut r, &[pair, pair]).unwrap();
            cohort.write_pair(&mut r, 0, pair)
        }
        Misuse::WriteSingleAfterWriteSingles => {
            let mut r = cohort.reserve_singles(1).unwrap();
            cohort.write_singles(&mut r, &[single]).unwrap();
            cohort.write_single(&mut r, 0, single)
        }
        Misuse::WritePairsTwice => {
            let mut r = written();
            cohort.write_pairs(&mut r, &[pair, pair])
        }
        Misuse::PairEmitFieldsAfterPairEmit | Misuse::SeedExtendAfterPairEmit => {
            let mut r = written();
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut r)
                .unwrap();
            let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
            cohort
                .pair_emit(
                    &idx,
                    &opts,
                    &mut scratch,
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut sink,
                )
                .unwrap();
            if matches!(misuse, Misuse::SeedExtendAfterPairEmit) {
                cohort.seed_extend(&idx, &opts, &mut scratch, &mut r)
            } else {
                cohort.pair_emit_fields(
                    &idx,
                    &opts,
                    &mut scratch,
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut NullFieldSink,
                )
            }
        }
        Misuse::WritePairsAfterSlotWrite => {
            let mut r = cohort.reserve_pairs(2).unwrap();
            cohort.write_pair(&mut r, 1, pair).unwrap();
            cohort.write_pairs(&mut r, &[pair, pair])
        }
        Misuse::ForeignWritePair => cohort.write_pair(&mut foreign(), 0, pair),
        Misuse::ForeignWriteSingle => {
            let mut r = ResidentCohort::new(false)
                .unwrap()
                .reserve_singles(1)
                .unwrap();
            cohort.write_single(&mut r, 0, single)
        }
        Misuse::ForeignSeedExtend => cohort.seed_extend(&idx, &opts, &mut scratch, &mut foreign()),
        Misuse::ForeignPairEmit => {
            let pestat = MemPeStat::zero().unwrap();
            let mut r = foreign();
            cohort.pair_emit(
                &idx,
                &opts,
                &mut scratch,
                &mut r,
                Some(&pestat),
                ids,
                0,
                &mut sink,
            )
        }
        Misuse::ForeignPairEmitFields => {
            let pestat = MemPeStat::zero().unwrap();
            let mut r = foreign();
            cohort.pair_emit_fields(
                &idx,
                &opts,
                &mut scratch,
                &mut r,
                Some(&pestat),
                ids,
                0,
                &mut NullFieldSink,
            )
        }
        Misuse::SingleIntoPairRange => {
            cohort.write_single(&mut cohort.reserve_pairs(1).unwrap(), 0, single)
        }
        Misuse::PairIntoSingleRange => {
            cohort.write_pair(&mut cohort.reserve_singles(1).unwrap(), 0, pair)
        }
        Misuse::IndexPastRange => cohort.write_pair(&mut cohort.reserve_pairs(2).unwrap(), 2, pair),
        Misuse::NameTooLong => {
            let name = vec![b'n'; 255];
            let long = ReadPair {
                name_r1: &name,
                ..pair
            };
            cohort.write_pair(&mut cohort.reserve_pairs(1).unwrap(), 0, long)
        }
        Misuse::SeedExtendUnwritten => {
            let mut r = cohort.reserve_pairs(2).unwrap();
            cohort.write_pair(&mut r, 0, pair).unwrap();
            cohort.seed_extend(&idx, &opts, &mut scratch, &mut r)
        }
        Misuse::WriteSlotTwice => {
            let mut r = cohort.reserve_pairs(2).unwrap();
            cohort.write_pair(&mut r, 0, pair).unwrap();
            cohort.write_pair(&mut r, 0, pair)
        }
        Misuse::WriteAfterSeedExtend | Misuse::SeedExtendTwice => {
            let mut r = cohort.reserve_pairs(1).unwrap();
            cohort.write_pair(&mut r, 0, pair).unwrap();
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut r)
                .unwrap();
            if matches!(misuse, Misuse::SeedExtendTwice) {
                cohort.seed_extend(&idx, &opts, &mut scratch, &mut r)
            } else {
                cohort.write_pair(&mut r, 0, pair)
            }
        }
        Misuse::SeedExtendOtherMeth => {
            cohort.seed_extend(&idx, &meth_opts, &mut scratch, &mut written())
        }
        Misuse::SeedExtendOtherIndex => {
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut written())
                .unwrap();
            let other = BwaIndex::load(ref_prefix().unwrap()).unwrap();
            cohort.seed_extend(&other, &opts, &mut scratch, &mut written())
        }
        Misuse::PairEmitBeforeSeedExtend => {
            let pestat = MemPeStat::zero().unwrap();
            let mut r = written();
            cohort.pair_emit(
                &idx,
                &opts,
                &mut scratch,
                &mut r,
                Some(&pestat),
                ids,
                0,
                &mut sink,
            )
        }
        Misuse::PairEmitFieldsBeforeSeedExtend => {
            let pestat = MemPeStat::zero().unwrap();
            let mut r = written();
            cohort.pair_emit_fields(
                &idx,
                &opts,
                &mut scratch,
                &mut r,
                Some(&pestat),
                ids,
                0,
                &mut NullFieldSink,
            )
        }
        Misuse::PairEmitWithoutPestat => {
            let mut r = written();
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut r)
                .unwrap();
            cohort.pair_emit(&idx, &opts, &mut scratch, &mut r, None, ids, 0, &mut sink)
        }
        Misuse::PairEmitOtherMeth
        | Misuse::PairEmitFieldsOtherMeth
        | Misuse::PairEmitTwice
        | Misuse::InferAfterPairEmit => {
            let mut r = written();
            cohort
                .seed_extend(&idx, &opts, &mut scratch, &mut r)
                .unwrap();
            let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
            let mut emit = |opts: &MemOpts, r: &mut ResidentRange, sink: &mut RecordVec| {
                cohort.pair_emit(&idx, opts, &mut scratch, r, Some(&pestat), ids, 0, sink)
            };
            match misuse {
                Misuse::PairEmitOtherMeth => emit(&meth_opts, &mut r, &mut sink),
                Misuse::PairEmitTwice => {
                    emit(&opts, &mut r, &mut sink).unwrap();
                    emit(&opts, &mut r, &mut sink)
                }
                Misuse::InferAfterPairEmit => {
                    emit(&opts, &mut r, &mut sink).unwrap();
                    cohort.infer_cohort(&idx, &opts).map(drop)
                }
                _ => cohort.pair_emit_fields(
                    &idx,
                    &meth_opts,
                    &mut AlignScratch::new().unwrap(),
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut NullFieldSink,
                ),
            }
        }
        Misuse::InferBeforeSeedExtend => {
            let _r = written();
            cohort.infer_cohort(&idx, &opts).map(drop)
        }
        Misuse::InferOtherMeth => cohort.infer_cohort(&idx, &meth_opts).map(drop),
    };
    assert_invalid(result, needle, &format!("{misuse:?}"));
    // A rejected pre-emit call must not have produced records.
    if !matches!(
        misuse,
        Misuse::PairEmitTwice
            | Misuse::InferAfterPairEmit
            | Misuse::PairEmitFieldsAfterPairEmit
            | Misuse::SeedExtendAfterPairEmit
    ) {
        assert!(
            sink.records.is_empty(),
            "{misuse:?}: a rejected call produced records"
        );
    }
}

/// A range reports its reservation: counts, region and reservation offset.
#[test]
fn resident_range_accessors_describe_the_reservation() {
    let cohort = ResidentCohort::new(false).unwrap();
    let a = cohort.reserve_pairs(3).unwrap();
    let b = cohort.reserve_pairs(2).unwrap();
    let singles = cohort.reserve_singles(4).unwrap();
    assert_eq!(
        (a.n_pairs(), a.n_reads(), a.first(), a.is_pairs()),
        (3, 6, 0, true)
    );
    assert_eq!(
        (b.n_pairs(), b.n_reads(), b.first(), b.is_pairs()),
        (2, 4, 6, true)
    );
    assert_eq!(
        (
            singles.n_pairs(),
            singles.n_reads(),
            singles.first(),
            singles.is_pairs()
        ),
        (0, 4, 0, false)
    );
    assert!(
        cohort.reserve_pairs(usize::MAX).is_err(),
        "a pair count whose read count overflows is rejected"
    );
}

/// `heap_bytes` is exact: each reserved read's headers plus every byte copied
/// in for its name, bases and qualities (the C side's own `+1`s included), and
/// alignment regions do not count, so seed_extend leaves it unchanged.
#[test]
fn heap_bytes_counts_reserved_headers_and_copied_reads() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let cohort = ResidentCohort::new(false).unwrap();
    // SAFETY: a pure query with no preconditions.
    let overhead = unsafe { bwa_mem3_sys::bwa_shim_resident_read_overhead() };
    let pairs = &fixture.pairs[..10];
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    assert_eq!(
        cohort.heap_bytes(),
        2 * pairs.len() * overhead,
        "headers only"
    );
    let pv: Vec<ReadPair<'_>> = pairs.iter().map(as_read_pair).collect();
    cohort.write_pairs(&mut range, &pv).unwrap();
    let copied: usize = pairs
        .iter()
        .flat_map(|(r1, r2)| [r1, r2])
        .map(|r| {
            (r.name.len() + 1) + (r.seq.len() + 1) + r.qual.as_ref().map_or(0, |q| q.len() + 1)
        })
        .sum();
    let expected = 2 * pairs.len() * overhead + copied;
    assert_eq!(cohort.heap_bytes(), expected, "headers + copied reads");
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    assert_eq!(
        cohort.heap_bytes(),
        expected,
        "alignment regions are not counted"
    );
}

/// Under `--meth` the resident path reproduces the legacy three-phase path
/// (`seed_extend` + `pair_emit`) byte-for-byte, singles included.
#[test]
fn resident_cohort_matches_legacy_under_meth() {
    let Some(meth) = phix_meth() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let idx = &meth.idx;
    let fixture = field_fixture();
    let mut opts = MemOpts::new().unwrap();
    opts.set_meth(true);
    let mut scratch = AlignScratch::new().unwrap();
    let pv: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let sv: Vec<SingleRead<'_>> = fixture
        .singles
        .iter()
        .map(|r| SingleRead {
            name: &r.name,
            seq: &r.seq,
            qual: r.qual.as_deref(),
        })
        .collect();
    let ids = IdBases {
        first_single_id: 2 * pv.len() as u64,
        first_pair_id: 0,
    };

    let regs = seed_extend(
        idx,
        &opts,
        &mut scratch,
        &ReadBatch {
            pairs: &pv,
            singles: &sv,
        },
    )
    .unwrap();
    let pestat = MemPeStat::infer_cohort(idx, &opts, std::slice::from_ref(&regs)).unwrap();
    let mut legacy = RecordVec::default();
    pair_emit(
        idx,
        &opts,
        &mut scratch,
        regs,
        Some(&pestat),
        ids,
        &mut legacy,
    )
    .unwrap();

    let cohort = ResidentCohort::new(true).unwrap();
    let mut pairs = cohort.reserve_pairs(pv.len()).unwrap();
    let mut singles = cohort.reserve_singles(sv.len()).unwrap();
    cohort.write_pairs(&mut pairs, &pv).unwrap();
    cohort.write_singles(&mut singles, &sv).unwrap();
    for range in [&mut pairs, &mut singles] {
        cohort.seed_extend(idx, &opts, &mut scratch, range).unwrap();
    }
    let resident_pestat = cohort.infer_cohort(idx, &opts).unwrap();
    let mut resident = RecordVec::default();
    for range in [&mut pairs, &mut singles] {
        cohort
            .pair_emit(
                idx,
                &opts,
                &mut scratch,
                range,
                Some(&resident_pestat),
                ids,
                0,
                &mut resident,
            )
            .unwrap();
    }
    assert!(!legacy.records.is_empty());
    assert_eq!(resident.records.len(), legacy.records.len(), "record count");
    assert!(
        resident.records == legacy.records,
        "resident --meth output diverged from the legacy path"
    );
}

/// A sink that panics mid-emission surfaces as a panic, and the cohort stays
/// fully usable afterwards: the unwind released the range lock (a whole-cohort
/// call does not deadlock) and another range still emits exactly the serial
/// reference's records.
#[test]
fn panicking_sink_leaves_the_cohort_usable() {
    struct Boom;
    impl RecordSink for Boom {
        fn emit(&mut self, _: RecordOrigin, _: &[u8]) {
            panic!("boom mid-emission");
        }
    }
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let (_, reference) = serial_reference(&idx, &opts, &fixture);

    let split = 80;
    let mut scratch = AlignScratch::new().unwrap();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut first = written_pair_range(&cohort, &fixture.pairs[..split]);
    let mut second = written_pair_range(&cohort, &fixture.pairs[split..]);
    for range in [&mut first, &mut second] {
        cohort
            .seed_extend(&idx, &opts, &mut scratch, range)
            .unwrap();
    }
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let ids = IdBases::default();
        cohort.pair_emit(
            &idx,
            &opts,
            &mut scratch,
            &mut first,
            Some(&pestat),
            ids,
            0,
            &mut Boom,
        )
    }));
    let payload = result.expect_err("the sink's panic must surface");
    assert_eq!(
        payload.downcast_ref::<&str>(),
        Some(&"boom mid-emission"),
        "the surfaced panic is the sink's own"
    );

    // `infer_cohort` takes the range lock exclusively, so reaching its answer
    // (the lifecycle rejection: a range has been emitted) proves the unwind
    // released the shared guard rather than leaking it.
    assert_invalid(
        cohort.infer_cohort(&idx, &opts).map(drop),
        LIFECYCLE,
        "infer_cohort after a panicked emit",
    );
    let got = emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut second,
        &pestat,
        split,
    );
    let expected = reference_slice(&reference, split, fixture.pairs.len() - split);
    assert!(!expected.is_empty());
    assert!(
        got == expected,
        "the other range's output changed after a panic"
    );
}

/// Which region a batch-write test writes.
#[derive(Clone, Copy, Debug)]
enum WriteRegion {
    Pairs,
    Singles,
}

/// How a test writes a range: one batch into the arena, or slot by slot.
#[derive(Clone, Copy, Debug)]
enum WriteMode {
    Batch,
    PerSlot,
}

/// The fixture's singles as borrowed reads.
fn fixture_singles(fixture: &Fixture) -> Vec<SingleRead<'_>> {
    fixture
        .singles
        .iter()
        .map(|r| SingleRead {
            name: &r.name,
            seq: &r.seq,
            qual: r.qual.as_deref(),
        })
        .collect()
}

/// Reserve and write a pair range holding `pairs`.
fn write_pair_range(
    cohort: &ResidentCohort,
    pairs: &[ReadPair<'_>],
    mode: WriteMode,
) -> ResidentRange {
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    match mode {
        WriteMode::Batch => cohort.write_pairs(&mut range, pairs).unwrap(),
        WriteMode::PerSlot => {
            for (i, p) in pairs.iter().enumerate() {
                cohort.write_pair(&mut range, i, *p).unwrap();
            }
        }
    }
    range
}

/// Reserve and write a single range holding `reads`.
fn write_single_range(
    cohort: &ResidentCohort,
    reads: &[SingleRead<'_>],
    mode: WriteMode,
) -> ResidentRange {
    let mut range = cohort.reserve_singles(reads.len()).unwrap();
    match mode {
        WriteMode::Batch => cohort.write_singles(&mut range, reads).unwrap(),
        WriteMode::PerSlot => {
            for (i, r) in reads.iter().enumerate() {
                cohort.write_single(&mut range, i, *r).unwrap();
            }
        }
    }
    range
}

/// Align `pairs` or `singles` in one range written with `mode`, returning the
/// packed records.
fn align_written_range(
    idx: &BwaIndex,
    opts: &MemOpts,
    meth: bool,
    region: WriteRegion,
    pairs: &[ReadPair<'_>],
    singles: &[SingleRead<'_>],
    mode: WriteMode,
) -> Records {
    let cohort = ResidentCohort::new(meth).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = match region {
        WriteRegion::Pairs => write_pair_range(&cohort, pairs, mode),
        WriteRegion::Singles => write_single_range(&cohort, singles, mode),
    };
    cohort
        .seed_extend(idx, opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = match region {
        WriteRegion::Pairs => Some(cohort.infer_cohort(idx, opts).unwrap()),
        WriteRegion::Singles => None,
    };
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(
            idx,
            opts,
            &mut scratch,
            &mut range,
            pestat.as_ref(),
            IdBases::default(),
            0,
            &mut sink,
        )
        .unwrap();
    sink.records
}

/// A range written in one batch (one arena) aligns byte-identically to the
/// same range written slot by slot, for pairs and singles, with and without
/// `--meth` (whose original bases and projection then live in the arena too).
#[rstest]
fn batch_write_matches_per_slot_write(
    #[values(WriteRegion::Pairs, WriteRegion::Singles)] region: WriteRegion,
    #[values(false, true)] meth: bool,
) {
    let Some((idx, opts)) = index_for(meth) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let singles = fixture_singles(&fixture);
    let run = |mode| align_written_range(&idx, &opts, meth, region, &pairs, &singles, mode);
    let (batched, per_slot) = (run(WriteMode::Batch), run(WriteMode::PerSlot));
    assert!(!batched.is_empty(), "the fixture must emit records");
    assert_eq!(batched, per_slot);
}

/// Under `--meth` the original bases are copied whole, not up to the first 0
/// byte: a read whose bases contain one aligns the same written either way.
#[rstest]
fn a_zero_base_byte_is_copied_whole_under_meth(
    #[values(WriteRegion::Pairs, WriteRegion::Singles)] region: WriteRegion,
) {
    let Some((idx, opts)) = index_for(true) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX meth index");
        return;
    };
    let fixture = field_fixture();
    let (r1, _) = &fixture.pairs[0];
    let mut seq_r1 = r1.seq.clone();
    seq_r1[1] = 0;
    let pair = ReadPair {
        seq_r1: &seq_r1,
        ..as_read_pair(&fixture.pairs[0])
    };
    let single = SingleRead {
        name: &r1.name,
        seq: &seq_r1,
        qual: r1.qual.as_deref(),
    };
    let run = |mode| align_written_range(&idx, &opts, true, region, &[pair], &[single], mode);
    let (batched, per_slot) = (run(WriteMode::Batch), run(WriteMode::PerSlot));
    assert!(!batched.is_empty());
    assert_eq!(batched, per_slot);
}

/// A batch write accounts exactly the bytes the per-slot writes do, with and
/// without `--meth` (which adds each read's original bases).
#[rstest]
fn batch_write_counts_the_same_heap_bytes(
    #[values(WriteRegion::Pairs, WriteRegion::Singles)] region: WriteRegion,
    #[values(false, true)] meth: bool,
) {
    let fixture = field_fixture();
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let singles = fixture_singles(&fixture);
    let bytes = |mode| {
        let cohort = ResidentCohort::new(meth).unwrap();
        let _range = match region {
            WriteRegion::Pairs => write_pair_range(&cohort, &pairs, mode),
            WriteRegion::Singles => write_single_range(&cohort, &singles, mode),
        };
        cohort.heap_bytes()
    };
    assert_eq!(bytes(WriteMode::Batch), bytes(WriteMode::PerSlot));
}

/// How a release test emits its range.
#[derive(Clone, Copy, Debug)]
enum EmitVia {
    Packed,
    Fields,
    PanickingSink,
}

/// Emitting a range releases its reads and regions in the shim and hands their
/// bytes back: the segment holds nothing afterwards and `heap_bytes` falls to
/// the headers-only figure, whichever emit call ran, however the range was
/// written, and even when the sink panics partway.
#[rstest]
fn emit_releases_the_ranges_reads(
    #[values(EmitVia::Packed, EmitVia::Fields, EmitVia::PanickingSink)] via: EmitVia,
    #[values(WriteMode::Batch, WriteMode::PerSlot)] mode: WriteMode,
) {
    struct Boom;
    impl RecordSink for Boom {
        fn emit(&mut self, _: RecordOrigin, _: &[u8]) {
            panic!("boom mid-emission");
        }
    }
    let Some((idx, opts)) = index_for(false) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    // SAFETY: a pure query with no preconditions.
    let overhead = unsafe { bwa_mem3_sys::bwa_shim_resident_read_overhead() };
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = write_pair_range(&cohort, &pairs, mode);
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let headers_only = 2 * pairs.len() * overhead;
    assert!(
        cohort.range_holds_reads(&range),
        "the reads are held before emit"
    );
    assert!(cohort.heap_bytes() > headers_only, "and counted");
    let ids = IdBases::default();
    match via {
        EmitVia::Packed => cohort
            .pair_emit(
                &idx,
                &opts,
                &mut scratch,
                &mut range,
                Some(&pestat),
                ids,
                0,
                &mut RecordVec::default(),
            )
            .unwrap(),
        EmitVia::Fields => cohort
            .pair_emit_fields(
                &idx,
                &opts,
                &mut scratch,
                &mut range,
                Some(&pestat),
                ids,
                0,
                &mut NullFieldSink,
            )
            .unwrap(),
        EmitVia::PanickingSink => {
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                cohort.pair_emit(
                    &idx,
                    &opts,
                    &mut scratch,
                    &mut range,
                    Some(&pestat),
                    ids,
                    0,
                    &mut Boom,
                )
            }));
            assert!(outcome.is_err(), "the sink's panic propagates");
        }
    }
    assert!(
        !cohort.range_holds_reads(&range),
        "the emit released the reads and regions"
    );
    assert_eq!(cohort.heap_bytes(), headers_only);
}

/// `kernel_batch_size` reports the kernel batch the vendored `macro.h` selects
/// for this build's architecture.
#[test]
fn kernel_batch_size_matches_the_build_target() {
    let expected = if cfg!(target_arch = "aarch64") {
        1024
    } else {
        512
    };
    assert_eq!(bwa_mem3_rs::kernel_batch_size(), expected);
}

/// An empty pair range writes, extends, infers and emits nothing, without error.
#[test]
fn an_empty_range_writes_and_emits_nothing() {
    let Some((idx, opts)) = index_for(false) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = cohort.reserve_pairs(0).unwrap();
    cohort.write_pairs(&mut range, &[]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut scratch,
            &mut range,
            Some(&pestat),
            IdBases::default(),
            0,
            &mut sink,
        )
        .unwrap();
    assert!(sink.records.is_empty());
}

/// The longest name the writers accept survives the arena intact in every
/// emitted record.
#[test]
fn the_longest_read_name_survives_the_arena() {
    let Some((idx, opts)) = index_for(false) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    // The longest name `ReadPair::validate` accepts (a BAM l_read_name of 255
    // with the NUL).
    let name = vec![b'n'; 254];
    let pair = ReadPair {
        name_r1: &name,
        name_r2: &name,
        ..as_read_pair(&fixture.pairs[0])
    };
    let records = align_written_range(
        &idx,
        &opts,
        false,
        WriteRegion::Pairs,
        &[pair],
        &[],
        WriteMode::Batch,
    );
    assert!(!records.is_empty());
    for (_, body) in &records {
        // BAM body: l_read_name at byte 8 (NUL included), name at 32.
        let l_name = usize::from(body[8]);
        assert_eq!(&body[32..32 + l_name - 1], name.as_slice());
    }
}
