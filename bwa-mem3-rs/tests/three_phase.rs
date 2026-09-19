//! Wrapper-level three-phase tests. Skipped unless `BWA_MEM3_RS_TEST_REF`
//! points at a bwa-mem3 index prefix (the sys crate covers the C behavior on
//! PhiX; this checks ownership, Send, and parity with `align_batch`).
//!
//! `BWA_MEM3_RS_TEST_REF` is typically a full hg38 index (~10 GB resident once
//! loaded). The harness runs `#[test]`s in one binary concurrently by default,
//! so tests that only need to *read* the index share one `Arc<BwaIndex>` via
//! [`shared_idx`] rather than each loading their own -- five independent
//! full-hg38 loads running at once is enough to trip the OOM killer on a
//! modest CI box. [`load_with_threads_matches_load`] is the deliberate
//! exception: it exists to test loading itself, so it loads its own pair.

use std::sync::{Arc, OnceLock};

use bwa_mem3_rs::{
    align_batch, build_info, pair_emit, seed_extend, version, AlignScratch, AlnRegs, BwaIndex,
    IdBases, MemOpts, MemPeStat, ReadBatch, ReadPair, RecordOrigin, RecordSink, RecordVec,
};

fn ref_prefix() -> Option<String> {
    let p = std::env::var("BWA_MEM3_RS_TEST_REF").ok()?;
    std::path::Path::new(&format!("{p}.bwt.2bit.64"))
        .exists()
        .then_some(p)
}

/// One `BwaIndex` shared by every test that only reads it, loaded at most
/// once per test binary run. See the module doc for why this matters.
fn shared_idx() -> Option<Arc<BwaIndex>> {
    static IDX: OnceLock<Option<Arc<BwaIndex>>> = OnceLock::new();
    IDX.get_or_init(|| {
        ref_prefix()
            .and_then(|p| BwaIndex::load(&p).ok())
            .map(Arc::new)
    })
    .clone()
}

fn pairs(n: usize) -> (Vec<String>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
    // Deterministic pseudo-reads: 150-mers over a fixed xorshift stream. Not
    // genomic, so most are unmapped -- enough to exercise ownership + parity.
    let mut x = 0x9E37_79B9u64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut names = Vec::new();
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for i in 0..n {
        names.push(format!("q{i}"));
        r1.push((0..150).map(|_| b"ACGT"[(next() % 4) as usize]).collect());
        r2.push((0..150).map(|_| b"ACGT"[(next() % 4) as usize]).collect());
    }
    (names, r1, r2)
}

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
        eprintln!("skip: set BWA_MEM3_RS_TEST_REF");
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

#[test]
fn regs_and_scratch_cross_threads() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: set BWA_MEM3_RS_TEST_REF");
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
    let a = BwaIndex::load(&prefix).unwrap();
    let b = BwaIndex::load_with_threads(&prefix, 4).unwrap();
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
