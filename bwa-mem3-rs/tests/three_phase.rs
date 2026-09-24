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
    align_batch, build_info, pair_emit, seed_extend, version, AlignScratch, AlignedFields,
    AlignedFieldsSink, AlnRegs, BwaIndex, IdBases, Mate, MemOpts, MemPeStat, ReadBatch, ReadPair,
    RecordOrigin, RecordSink, RecordVec, ResidentCohort, ResidentRange, SingleRead,
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
}

/// A structured sink that discards every record.
struct NullFieldSink;
impl AlignedFieldsSink for NullFieldSink {
    fn emit(&mut self, _: RecordOrigin, _: Mate, _: bool, _: AlignedFields<'_>) {}
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
    if !matches!(misuse, Misuse::PairEmitTwice | Misuse::InferAfterPairEmit) {
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

/// A range written in one batch (one arena) aligns byte-identically to the same
/// range written slot by slot, with and without `--meth`.
#[rstest]
#[case::plain(false)]
#[case::meth(true)]
fn batch_write_matches_per_slot_write(#[case] meth: bool) {
    let idx = if meth {
        phix_meth().map(|r| r.idx.clone())
    } else {
        shared_idx()
    };
    let Some(idx) = idx else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let mut opts = MemOpts::new().unwrap();
    if meth {
        opts.set_meth(true);
    }
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let run = |batch: bool| -> Records {
        let cohort = ResidentCohort::new(meth).unwrap();
        let mut scratch = AlignScratch::new().unwrap();
        let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
        if batch {
            cohort.write_pairs(&mut range, &pairs).unwrap();
        } else {
            for (i, p) in pairs.iter().enumerate() {
                cohort.write_pair(&mut range, i, *p).unwrap();
            }
        }
        cohort
            .seed_extend(&idx, &opts, &mut scratch, &mut range)
            .unwrap();
        let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
        emit_range(&cohort, &idx, &opts, &mut scratch, &mut range, &pestat, 0)
    };
    let (batched, per_slot) = (run(true), run(false));
    assert!(!batched.is_empty(), "the fixture must emit records");
    assert_eq!(batched, per_slot);
}

/// A batch write accounts exactly the bytes the per-slot writes do.
#[rstest]
#[case::pairs(true)]
#[case::singles(false)]
fn batch_write_counts_the_same_heap_bytes(#[case] pairs_region: bool) {
    let fixture = field_fixture();
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let singles: Vec<SingleRead<'_>> = fixture
        .singles
        .iter()
        .map(|r| SingleRead {
            name: &r.name,
            seq: &r.seq,
            qual: r.qual.as_deref(),
        })
        .collect();
    let bytes = |batch: bool| {
        let cohort = ResidentCohort::new(false).unwrap();
        if pairs_region {
            let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
            if batch {
                cohort.write_pairs(&mut range, &pairs).unwrap();
            } else {
                for (i, p) in pairs.iter().enumerate() {
                    cohort.write_pair(&mut range, i, *p).unwrap();
                }
            }
        } else {
            let mut range = cohort.reserve_singles(singles.len()).unwrap();
            if batch {
                cohort.write_singles(&mut range, &singles).unwrap();
            } else {
                for (i, r) in singles.iter().enumerate() {
                    cohort.write_single(&mut range, i, *r).unwrap();
                }
            }
        }
        cohort.heap_bytes()
    };
    assert_eq!(bytes(true), bytes(false));
}

/// How a release test emits its range.
#[derive(Clone, Copy, Debug)]
enum EmitVia {
    Packed,
    Fields,
    PanickingSink,
}

/// Emitting a range hands its copied read bytes back: `heap_bytes` falls to the
/// headers-only figure, whichever emit call ran and even when the sink panics
/// (the shim has already released the reads by then).
#[rstest]
#[case::packed(EmitVia::Packed)]
#[case::fields(EmitVia::Fields)]
#[case::panicking_sink(EmitVia::PanickingSink)]
fn emit_releases_the_ranges_read_bytes(#[case] via: EmitVia) {
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
    // SAFETY: a pure query with no preconditions.
    let overhead = unsafe { bwa_mem3_sys::bwa_shim_resident_read_overhead() };
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    cohort.write_pairs(&mut range, &pairs).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let headers_only = 2 * pairs.len() * overhead;
    assert!(
        cohort.heap_bytes() > headers_only,
        "the reads are held before emit"
    );
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
    assert_eq!(cohort.heap_bytes(), headers_only);
}

/// After a range is emitted (and its reads released), every later call on it,
/// and the cohort pestat, is a lifecycle error rather than a read of freed
/// memory.
#[rstest]
#[case::emit_again("emit")]
#[case::emit_fields_again("emit_fields")]
#[case::extend_again("extend")]
#[case::infer_after_emit("infer")]
fn released_range_rejects_every_later_call(#[case] call: &str) {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let opts = MemOpts::new().unwrap();
    let pairs: Vec<ReadPair<'_>> = fixture.pairs.iter().map(as_read_pair).collect();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    cohort.write_pairs(&mut range, &pairs).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut range)
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let ids = IdBases::default();
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut scratch,
            &mut range,
            Some(&pestat),
            ids,
            0,
            &mut sink,
        )
        .unwrap();
    let result: bwa_mem3_rs::Result<()> = match call {
        "emit" => cohort.pair_emit(
            &idx,
            &opts,
            &mut scratch,
            &mut range,
            Some(&pestat),
            ids,
            0,
            &mut sink,
        ),
        "emit_fields" => cohort.pair_emit_fields(
            &idx,
            &opts,
            &mut scratch,
            &mut range,
            Some(&pestat),
            ids,
            0,
            &mut NullFieldSink,
        ),
        "extend" => cohort.seed_extend(&idx, &opts, &mut scratch, &mut range),
        "infer" => cohort.infer_cohort(&idx, &opts).map(drop),
        other => unreachable!("unknown case {other}"),
    };
    let err = result.unwrap_err();
    assert!(err.to_string().contains(LIFECYCLE), "{call}: {err}");
}
