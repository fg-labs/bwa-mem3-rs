//! Wrapper-level tests for the read-pair memo: `PairMemo` +
//! `ResidentCohort::seed_extend_with_reps` + `ResidentCohort::resolve_memo`.
//! They build a small PhiX index at test time (see `support`), so they run in
//! CI wherever `bwa-mem3` is available and skip otherwise (a hard failure under
//! `BWA_MEM3_RS_REQUIRE_TOOLS`).

use std::sync::Arc;

use bwa_mem3_rs::{
    AlignScratch, AlignedFields, AlignedFieldsSink, BwaIndex, IdBases, Mate, MemOpts, MemPeStat,
    ReadPair, RecordOrigin, RecordVec, ResidentCohort,
};
use rstest::rstest;

mod support;
use support::*;

/// The fixture's pairs as borrowed reads.
fn fixture_pairs(fixture: &Fixture) -> Vec<ReadPair<'_>> {
    fixture.pairs.iter().map(as_read_pair).collect()
}

/// The PhiX index for `meth`, and options in that mode (as
/// `three_phase.rs::index_for`).
fn index_for(meth: bool) -> Option<(Arc<BwaIndex>, MemOpts)> {
    let idx = if meth {
        phix_meth().map(|r| r.idx.clone())
    } else {
        shared_idx()
    }?;
    let mut opts = MemOpts::new().unwrap();
    if meth {
        opts.set_meth(true);
    }
    Some((idx, opts))
}

/// A structured sink that discards every record.
struct NullFieldSink;
impl AlignedFieldsSink for NullFieldSink {
    fn emit(&mut self, _: RecordOrigin, _: Mate, _: bool, _: AlignedFields<'_>) {}
}

/// Align `pairs` as one cohort cut into ranges of `sub` pairs, each seeded by
/// `seed` (called with the range and its local pair slice), then pestat and
/// emit every range in order.
fn run_cohort_with(
    idx: &BwaIndex,
    opts: &MemOpts,
    meth: bool,
    pairs: &[ReadPair<'_>],
    sub: usize,
    mut seed: impl FnMut(
        &ResidentCohort,
        &mut AlignScratch,
        &mut bwa_mem3_rs::ResidentRange,
        &[ReadPair<'_>],
    ),
) -> Records {
    let cohort = ResidentCohort::new(meth).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut ranges = Vec::new();
    for (k, chunk) in pairs.chunks(sub).enumerate() {
        let mut range = cohort.reserve_pairs(chunk.len()).unwrap();
        cohort.write_pairs(&mut range, chunk).unwrap();
        seed(&cohort, &mut scratch, &mut range, chunk);
        ranges.push((range, k * sub));
    }
    let pestat = cohort.infer_cohort(idx, opts).unwrap();
    let mut records = Records::new();
    for (mut range, p0) in ranges {
        records.extend(emit_range(
            &cohort,
            idx,
            opts,
            &mut scratch,
            &mut range,
            &pestat,
            p0,
        ));
    }
    records
}

/// `seed_extend_with_reps` with no duplicate marks is `seed_extend`: same
/// records, with and without `--meth`, at several range sizes.
#[rstest]
fn all_none_reps_match_seed_extend(
    #[values(false, true)] meth: bool,
    #[values(1, 7, 64)] sub: usize,
) {
    let Some((idx, opts)) = index_for(meth) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pairs = fixture_pairs(&fixture);
    let plain = run_cohort_with(&idx, &opts, meth, &pairs, sub, |c, sc, r, _| {
        c.seed_extend(&idx, &opts, sc, r).unwrap();
    });
    let with_reps = run_cohort_with(&idx, &opts, meth, &pairs, sub, |c, sc, r, chunk| {
        c.seed_extend_with_reps(&idx, &opts, sc, r, &vec![None; chunk.len()])
            .unwrap();
    });
    assert!(
        with_reps == plain,
        "all-None reps diverged from seed_extend (meth={meth}, sub={sub})"
    );
}

/// A memo call that must be rejected, run against fresh state.
#[derive(Debug, Clone, Copy)]
enum MemoMisuse {
    RepsLenMismatch,
    RepsOnSingleRange,
    RepNotEarlier,
    RepsUnderMeth,
    EmitBeforeResolve,
    EmitFieldsBeforeResolve,
    InferBeforeResolve,
}

#[rstest]
#[case::reps_len_mismatch(MemoMisuse::RepsLenMismatch, "marks for a range of")]
#[case::reps_on_single_range(MemoMisuse::RepsOnSingleRange, "not a pair range")]
#[case::rep_not_earlier(MemoMisuse::RepNotEarlier, "not an earlier pair")]
#[case::reps_under_meth(MemoMisuse::RepsUnderMeth, "not supported under --meth")]
#[case::emit_before_resolve(MemoMisuse::EmitBeforeResolve, "call resolve_memo first")]
#[case::emit_fields_before_resolve(MemoMisuse::EmitFieldsBeforeResolve, "call resolve_memo first")]
#[case::infer_before_resolve(MemoMisuse::InferBeforeResolve, "call resolve_memo first")]
fn memo_misuse_is_rejected(#[case] misuse: MemoMisuse, #[case] needle: &str) {
    let meth = matches!(misuse, MemoMisuse::RepsUnderMeth);
    let Some((idx, opts)) = index_for(meth) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let pair = as_read_pair(&fixture.pairs[0]);
    let cohort = ResidentCohort::new(meth).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut sink = RecordVec::default();
    // A two-pair range, written, whose second pair duplicates the first.
    let written = || {
        let mut r = cohort.reserve_pairs(2).unwrap();
        cohort.write_pairs(&mut r, &[pair, pair]).unwrap();
        r
    };
    let result: bwa_mem3_rs::Result<()> = match misuse {
        MemoMisuse::RepsLenMismatch => {
            cohort.seed_extend_with_reps(&idx, &opts, &mut scratch, &mut written(), &[None])
        }
        MemoMisuse::RepsOnSingleRange => {
            let mut r = cohort.reserve_singles(0).unwrap();
            cohort.seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &[])
        }
        MemoMisuse::RepNotEarlier => cohort.seed_extend_with_reps(
            &idx,
            &opts,
            &mut scratch,
            &mut written(),
            &[None, Some(1)],
        ),
        MemoMisuse::RepsUnderMeth => cohort.seed_extend_with_reps(
            &idx,
            &opts,
            &mut scratch,
            &mut written(),
            &[None, Some(0)],
        ),
        MemoMisuse::EmitBeforeResolve
        | MemoMisuse::EmitFieldsBeforeResolve
        | MemoMisuse::InferBeforeResolve => {
            let mut r = written();
            cohort
                .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &[None, Some(0)])
                .unwrap();
            let pestat = MemPeStat::zero().unwrap();
            let ids = IdBases::default();
            match misuse {
                MemoMisuse::EmitBeforeResolve => cohort.pair_emit(
                    &idx,
                    &opts,
                    &mut scratch,
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut sink,
                ),
                MemoMisuse::EmitFieldsBeforeResolve => cohort.pair_emit_fields(
                    &idx,
                    &opts,
                    &mut scratch,
                    &mut r,
                    Some(&pestat),
                    ids,
                    0,
                    &mut NullFieldSink,
                ),
                _ => cohort.infer_cohort(&idx, &opts).map(drop),
            }
        }
    };
    assert_invalid(result, needle, &format!("{misuse:?}"));
    assert!(
        sink.records.is_empty(),
        "{misuse:?}: a rejected call produced records"
    );
}
