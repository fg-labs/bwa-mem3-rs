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

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, OnceLock};

use bwa_mem3_rs::{
    align_batch, build_info, pair_emit, seed_extend, version, AlignScratch, AlnRegs, BwaIndex,
    IdBases, MemOpts, MemPeStat, ReadBatch, ReadPair, RecordOrigin, RecordSink, RecordVec,
};

#[path = "../../bwa-mem3-rs-cli/tests/phix_seq.rs"]
mod phix_seq;

/// A PhiX index built once per test-binary run, kept alive for the process.
struct PhixRef {
    _dir: tempfile::TempDir, // holds the on-disk index files alive
    prefix: PathBuf,
    idx: Arc<BwaIndex>,
}

/// Locate `bwa-mem3` via `BWA_MEM3_BIN` or `PATH`. Returns `None` (with a skip
/// message) when it is absent, unless `BWA_MEM3_RS_REQUIRE_TOOLS` is set, which
/// turns the absence into a hard failure.
fn find_bwa_mem3() -> Option<String> {
    if let Ok(p) = std::env::var("BWA_MEM3_BIN") {
        if Path::new(&p).exists() {
            return Some(p);
        }
    }
    let out = Command::new("which").arg("bwa-mem3").output().ok();
    let found = out.and_then(|o| {
        o.status
            .success()
            .then(|| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .filter(|p| !p.is_empty())
    });
    if found.is_none() {
        assert!(
            std::env::var_os("BWA_MEM3_RS_REQUIRE_TOOLS").is_none(),
            "BWA_MEM3_RS_REQUIRE_TOOLS is set but bwa-mem3 was not found; \
             set BWA_MEM3_BIN or install it on PATH"
        );
        eprintln!("skip: bwa-mem3 not on PATH (set BWA_MEM3_BIN)");
    }
    found
}

/// Build+load a PhiX index once, shared across the whole test binary. `None`
/// (skip) when `bwa-mem3` is unavailable and not required.
fn phix() -> Option<&'static PhixRef> {
    static REF: OnceLock<Option<PhixRef>> = OnceLock::new();
    REF.get_or_init(|| {
        let bwa = find_bwa_mem3()?;
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("phix.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">phix").unwrap();
        for chunk in phix_seq::PHIX_SEQ.as_bytes().chunks(72) {
            f.write_all(chunk).unwrap();
            writeln!(f).unwrap();
        }
        drop(f);
        let status = Command::new(&bwa)
            .arg("index")
            .arg(&fa)
            .status()
            .expect("run bwa-mem3 index");
        assert!(status.success(), "bwa-mem3 index failed");
        let idx = Arc::new(BwaIndex::load(&fa).expect("load PhiX index"));
        Some(PhixRef {
            _dir: dir,
            prefix: fa,
            idx,
        })
    })
    .as_ref()
}

/// The PhiX index prefix (its FASTA path), or `None` (skip) without the tools.
fn ref_prefix() -> Option<&'static Path> {
    phix().map(|p| p.prefix.as_path())
}

/// One `BwaIndex` shared by every test that only reads it, loaded at most
/// once per test binary run. See the module doc for why this matters.
fn shared_idx() -> Option<Arc<BwaIndex>> {
    phix().map(|p| p.idx.clone())
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
