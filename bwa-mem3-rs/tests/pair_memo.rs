//! Wrapper-level tests for the read-pair memo: `PairMemo` +
//! `ResidentCohort::seed_extend_with_reps` + `ResidentCohort::resolve_memo`.
//! They build a small PhiX index at test time (see `support`), so they run in
//! CI wherever `bwa-mem3` is available and skip otherwise (a hard failure under
//! `BWA_MEM3_RS_REQUIRE_TOOLS`).

use std::sync::Arc;

use bwa_mem3_rs::{
    AlignScratch, BwaIndex, IdBases, MemOpts, MemPeStat, MemoStats, PairMemo, ReadPair, RecordVec,
    ResidentCohort, ResidentRange,
};
use rstest::rstest;

mod support;
use support::*;

/// The fixture's pairs as borrowed reads.
fn fixture_pairs(fixture: &Fixture) -> Vec<ReadPair<'_>> {
    fixture.pairs.iter().map(as_read_pair).collect()
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

/// Align `pairs` as one cohort cut into ranges of `sub` pairs. With `memo`, each
/// range is marked by a `PairMemo`, seeded with `seed_extend_with_reps`, and
/// the cohort is resolved before the pestat; without it, plain `seed_extend`.
/// Returns the records and, with `memo`, the resolve stats.
fn run_cohort_memo(
    idx: &BwaIndex,
    opts: &MemOpts,
    pairs: &[ReadPair<'_>],
    sub: usize,
    memo: bool,
) -> (Records, Option<MemoStats>) {
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut pair_memo = PairMemo::new();
    let mut ranges: Vec<(ResidentRange, usize)> = Vec::new();
    for (k, chunk) in pairs.chunks(sub.max(1)).enumerate() {
        let mut range = cohort.reserve_pairs(chunk.len()).unwrap();
        cohort.write_pairs(&mut range, chunk).unwrap();
        if memo {
            let reps = pair_memo.mark_range(&range, chunk).unwrap();
            cohort
                .seed_extend_with_reps(idx, opts, &mut scratch, &mut range, &reps)
                .unwrap();
        } else {
            cohort
                .seed_extend(idx, opts, &mut scratch, &mut range)
                .unwrap();
        }
        ranges.push((range, k * sub.max(1)));
    }
    let stats = memo.then(|| cohort.resolve_memo(idx, opts, &mut scratch).unwrap());
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
    (records, stats)
}

/// What [`copy_read`] does with the source read's QUAL.
#[derive(Clone, Copy)]
enum Qual {
    /// Keep it (a representative, exactly as the fixture built it).
    Keep,
    /// Reverse it, so a duplicate that kept its representative's QUAL shows.
    Reverse,
    /// Drop it when present, add one when absent.
    Flip,
}

/// An owned copy of `r` under `name`, its QUAL treated per `qual`.
fn copy_read(r: &FixtureRead, name: &str, qual: Qual) -> FixtureRead {
    let qual = match (&r.qual, qual) {
        (q, Qual::Keep) => q.clone(),
        (q, Qual::Reverse) => q.as_ref().map(|q| q.iter().rev().copied().collect()),
        (Some(_), Qual::Flip) => None,
        (None, Qual::Flip) => Some(vec![b'5'; r.seq.len()]),
    };
    FixtureRead {
        name: name.as_bytes().to_vec(),
        seq: r.seq.clone(),
        qual,
    }
}

/// The PhiX field fixture plus duplicates: every third pair gets a copy (new
/// name, reversed QUAL), every ninth a second copy (a triplicate) whose QUAL
/// presence is flipped, and near-misses that must NOT merge: swapped mates, a
/// pair with one mate replaced, and a lowercased copy. Copies are appended
/// after the originals, so they land in later ranges than their representatives
/// as well as in the same range.
fn dup_fixture(fixture: &Fixture) -> Vec<(FixtureRead, FixtureRead)> {
    let mut out: Vec<(FixtureRead, FixtureRead)> = fixture
        .pairs
        .iter()
        .map(|(a, b)| {
            (
                copy_read(a, &String::from_utf8_lossy(&a.name), Qual::Keep),
                copy_read(b, &String::from_utf8_lossy(&b.name), Qual::Keep),
            )
        })
        .collect();
    let n = out.len();
    for i in (0..n).step_by(3) {
        let name = format!("dup{i}");
        let (a, b) = (&fixture.pairs[i].0, &fixture.pairs[i].1);
        out.push((
            copy_read(a, &name, Qual::Reverse),
            copy_read(b, &name, Qual::Reverse),
        ));
        if i % 9 == 0 {
            let name = format!("tri{i}");
            out.push((
                copy_read(a, &name, Qual::Flip),
                copy_read(b, &name, Qual::Flip),
            ));
        }
    }
    let (a, b) = (&fixture.pairs[1].0, &fixture.pairs[1].1);
    out.push((
        copy_read(b, "swap", Qual::Keep),
        copy_read(a, "swap", Qual::Keep),
    ));
    out.push((
        copy_read(a, "onemate", Qual::Keep),
        copy_read(&fixture.pairs[2].1, "onemate", Qual::Keep),
    ));
    let mut lower = copy_read(a, "lower", Qual::Keep);
    lower.seq.make_ascii_lowercase();
    out.push((lower, copy_read(b, "lower", Qual::Keep)));
    out
}

/// Pairs of `pairs` whose bases (both mates) repeat an earlier pair's: the
/// exact duplicate count a `PairMemo` must find.
fn expected_dups(pairs: &[ReadPair<'_>]) -> u64 {
    let mut seen = std::collections::HashSet::new();
    pairs
        .iter()
        .filter(|p| !seen.insert((p.seq_r1.to_vec(), p.seq_r2.to_vec())))
        .count() as u64
}

/// All of `dups` duplicates copied, none fallen back.
fn all_copied(dups: u64) -> MemoStats {
    MemoStats {
        dup_pairs: dups,
        copied: dups,
        fallback_aligned: 0,
    }
}

/// Memo on is byte-identical to memo off at every range size, and the memo
/// actually fired: every marked duplicate was copied, none fell back.
#[rstest]
fn memo_on_matches_memo_off(#[values(1, 5, 32, 10_000)] sub: usize) {
    let Some((idx, opts)) = index_for(false) else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let owned = dup_fixture(&fixture);
    let pairs: Vec<ReadPair<'_>> = owned.iter().map(as_read_pair).collect();
    let (off, _) = run_cohort_memo(&idx, &opts, &pairs, sub, false);
    let (on, stats) = run_cohort_memo(&idx, &opts, &pairs, sub, true);
    let dups = expected_dups(&pairs);
    assert!(dups > 0, "the fixture has duplicates");
    assert_eq!(stats.unwrap(), all_copied(dups));
    assert_eq!(on.len(), off.len(), "record count (sub={sub})");
    assert!(on == off, "memo on diverged from memo off (sub={sub})");
}

/// Resolve a hand-built cohort: `pairs` in one range with marks `reps`; returns
/// the records and stats, plus the memo-off records of the same pairs.
fn resolve_hand_marked(
    idx: &BwaIndex,
    opts: &MemOpts,
    pairs: &[ReadPair<'_>],
    reps: &[Option<u64>],
) -> (Records, MemoStats, Records) {
    let (off, _) = run_cohort_memo(idx, opts, pairs, pairs.len().max(1), false);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut range = cohort.reserve_pairs(pairs.len()).unwrap();
    cohort.write_pairs(&mut range, pairs).unwrap();
    cohort
        .seed_extend_with_reps(idx, opts, &mut scratch, &mut range, reps)
        .unwrap();
    let stats = cohort.resolve_memo(idx, opts, &mut scratch).unwrap();
    let pestat = cohort.infer_cohort(idx, opts).unwrap();
    let got = emit_range(&cohort, idx, opts, &mut scratch, &mut range, &pestat, 0);
    (got, stats, off)
}

/// A pair marked against a pair with different bases is aligned normally at
/// resolve, on its own as a two-read batch: the output equals memo-off,
/// including its own name and QUAL, for every record shape the fixture has
/// (`field_fixture_from` builds shape `i % 8`).
#[rstest]
#[case::rf_pair(2)]
#[case::indels(3)]
#[case::chimeric(4)]
#[case::half_mapped(5)]
#[case::both_unmapped(6)]
#[case::iupac_and_lowercase(7)]
fn mismatched_mark_falls_back_to_normal_alignment(#[case] marked: usize) {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pairs: Vec<ReadPair<'_>> = fixture.pairs[..8].iter().map(as_read_pair).collect();
    let mut reps = vec![None; pairs.len()];
    reps[marked] = Some(0);
    let (got, stats, off) = resolve_hand_marked(&idx, &opts, &pairs, &reps);
    assert_eq!(
        stats,
        MemoStats {
            dup_pairs: 1,
            copied: 0,
            fallback_aligned: 1
        }
    );
    assert!(got == off, "fallback output diverged from memo off");
}

/// `N` and `R` both encode to 2-bit 4, so a hand-built mark
/// between them passes the 2-bit check and is copied, and the output still
/// equals memo off (SEQ is regenerated from the 2-bit bases).
#[test]
fn hand_marked_iupac_variant_is_copied_and_output_invariant() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let (a, b) = &fixture.pairs[0];
    let mut with_n = a.seq.clone();
    with_n[10] = b'N';
    let mut with_r = a.seq.clone();
    with_r[10] = b'R';
    let q = vec![b'I'; a.seq.len()];
    let qb = vec![b'I'; b.seq.len()];
    let p_n = ReadPair {
        name_r1: b"n",
        seq_r1: &with_n,
        qual_r1: Some(&q),
        name_r2: b"n",
        seq_r2: &b.seq,
        qual_r2: Some(&qb),
    };
    let p_r = ReadPair {
        name_r1: b"r",
        seq_r1: &with_r,
        name_r2: b"r",
        ..p_n
    };
    let (got, stats, off) = resolve_hand_marked(&idx, &opts, &[p_n, p_r], &[None, Some(0)]);
    assert_eq!(
        stats,
        MemoStats {
            dup_pairs: 1,
            copied: 1,
            fallback_aligned: 0
        }
    );
    assert!(got == off, "an IUPAC-variant copy changed the output");
}

/// A duplicate pair whose R2 is empty copies (empty) regs and
/// matches memo off.
#[test]
fn duplicate_with_an_empty_mate() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let a = &fixture.pairs[0].0;
    let mk = |name: &'static [u8]| ReadPair {
        name_r1: name,
        seq_r1: &a.seq,
        qual_r1: a.qual.as_deref(),
        name_r2: name,
        seq_r2: b"",
        qual_r2: None,
    };
    let pairs = [mk(b"e1"), mk(b"e2")];
    let cohort = ResidentCohort::new(false).unwrap();
    let reps = PairMemo::new()
        .mark_range(&cohort.reserve_pairs(2).unwrap(), &pairs)
        .unwrap();
    assert_eq!(reps, vec![None, Some(0)]);
    let (got, stats, off) = resolve_hand_marked(&idx, &opts, &pairs, &reps);
    assert_eq!(stats.copied, 1);
    assert!(got == off, "an empty-mate duplicate diverged");
}

/// A range spanning several kernel batches where whole
/// batches are duplicates (kernel1 skipped for them), with duplicates
/// straddling batch boundaries.
#[test]
fn multi_chunk_range_with_an_all_duplicate_chunk() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let per_batch = bwa_mem3_rs::kernel_batch_size() / 2; // pairs per kernel batch
    let distinct = fixture.pairs.len().min(per_batch);
    let mut owned: Vec<(FixtureRead, FixtureRead)> = Vec::new();
    for (i, (a, b)) in fixture.pairs[..distinct].iter().enumerate() {
        let name = format!("o{i}");
        owned.push((
            copy_read(a, &name, Qual::Keep),
            copy_read(b, &name, Qual::Keep),
        ));
    }
    for k in 0..(2 * per_batch + per_batch / 2) {
        let (a, b) = &fixture.pairs[k % distinct];
        let name = format!("c{k}");
        owned.push((
            copy_read(a, &name, Qual::Reverse),
            copy_read(b, &name, Qual::Reverse),
        ));
    }
    let pairs: Vec<ReadPair<'_>> = owned.iter().map(as_read_pair).collect();
    let (off, _) = run_cohort_memo(&idx, &opts, &pairs, pairs.len(), false);
    let (on, stats) = run_cohort_memo(&idx, &opts, &pairs, pairs.len(), true);
    assert_eq!(stats.unwrap(), all_copied(expected_dups(&pairs)));
    assert!(on == off, "multi-chunk memo diverged");
}

/// A duplicate whose representative's range was extended with plain
/// `seed_extend` (mixed API use) is copied like any other.
#[test]
fn representative_extended_with_plain_seed_extend() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let owned = dup_fixture(&fixture);
    let pairs: Vec<ReadPair<'_>> = owned.iter().map(as_read_pair).collect();
    let n0 = fixture.pairs.len(); // the originals: range 0; the copies: range 1
    let (off, _) = run_cohort_memo(&idx, &opts, &pairs, n0, false);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut memo = PairMemo::new();
    let mut r0 = cohort.reserve_pairs(n0).unwrap();
    cohort.write_pairs(&mut r0, &pairs[..n0]).unwrap();
    memo.mark_range(&r0, &pairs[..n0]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut r0)
        .unwrap();
    let mut r1 = cohort.reserve_pairs(pairs.len() - n0).unwrap();
    cohort.write_pairs(&mut r1, &pairs[n0..]).unwrap();
    let reps = memo.mark_range(&r1, &pairs[n0..]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r1, &reps)
        .unwrap();
    let stats = cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap();
    assert_eq!(stats, all_copied(expected_dups(&pairs)));
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let mut got = emit_range(&cohort, &idx, &opts, &mut scratch, &mut r0, &pestat, 0);
    got.extend(emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut r1,
        &pestat,
        n0,
    ));
    assert!(got == off, "mixed-API memo diverged");
}

/// Resolve is idempotent; a cohort with no marks resolves to zero stats; a
/// range reserved after resolve may not carry marks.
#[test]
fn resolve_idempotent_noop_and_late_range() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let owned = dup_fixture(&fixture);
    let pairs: Vec<ReadPair<'_>> = owned.iter().map(as_read_pair).collect();
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut memo = PairMemo::new();
    let mut r = cohort.reserve_pairs(pairs.len()).unwrap();
    cohort.write_pairs(&mut r, &pairs).unwrap();
    let reps = memo.mark_range(&r, &pairs).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &reps)
        .unwrap();
    let first = cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap();
    assert!(first.dup_pairs > 0);
    assert_eq!(
        cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap(),
        first
    );
    let mut late = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut late, &pairs[..2]).unwrap();
    assert_invalid(
        cohort.seed_extend_with_reps(&idx, &opts, &mut scratch, &mut late, &[None, Some(0)]),
        "memo is closed",
        "marks after resolve",
    );

    let plain = ResidentCohort::new(false).unwrap();
    let mut p = plain.reserve_pairs(pairs.len()).unwrap();
    plain.write_pairs(&mut p, &pairs).unwrap();
    plain
        .seed_extend(&idx, &opts, &mut scratch, &mut p)
        .unwrap();
    assert_eq!(
        plain.resolve_memo(&idx, &opts, &mut scratch).unwrap(),
        MemoStats::default()
    );
}

/// Resolve needs every pair range seed-extended, and says which range is not.
#[test]
fn resolve_with_an_unextended_range_names_it() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pair = as_read_pair(&fixture.pairs[0]);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut r = cohort.reserve_pairs(3).unwrap();
    cohort.write_pairs(&mut r, &[pair, pair, pair]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &[None, Some(0), Some(0)])
        .unwrap();
    let _unextended = cohort.reserve_pairs(1).unwrap();
    assert_invalid(
        cohort.resolve_memo(&idx, &opts, &mut scratch).map(drop),
        "is not seed-extended",
        "resolve with an unextended range",
    );
}

/// A `resolve_memo` that fails because a pair range is not seed-extended
/// leaves the memo open: that range can still be seeded with its marks, and a
/// retried `resolve_memo` then succeeds.
#[test]
fn a_failed_resolve_can_be_retried_once_the_range_is_extended() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pair = as_read_pair(&fixture.pairs[0]);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut first = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut first, &[pair, pair]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut first, &[None, Some(0)])
        .unwrap();
    let mut late = cohort.reserve_pairs(1).unwrap();
    assert!(cohort.resolve_memo(&idx, &opts, &mut scratch).is_err());
    cohort.write_pairs(&mut late, &[pair]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut late, &[Some(0)])
        .unwrap();
    assert_eq!(
        cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap(),
        all_copied(2)
    );
}

/// A mark naming a pair that is itself a duplicate (only a hand-built mark can)
/// resolves through it to the first copy, rather than failing: the output
/// equals memo off.
#[test]
fn chained_marks_resolve_to_the_first_copy() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pair = as_read_pair(&fixture.pairs[0]);
    let pairs = [pair, pair, pair];
    let (got, stats, off) = resolve_hand_marked(&idx, &opts, &pairs, &[None, Some(0), Some(1)]);
    assert_eq!(stats, all_copied(2));
    assert!(got == off, "chained marks diverged from memo off");
}

/// A fallback pair holding a `-` base aligns exactly as memo off does: the
/// fallback re-seeds bases that are already 2-bit, and `-` is the one code
/// (5) that a second conversion would change.
#[test]
fn fallback_keeps_a_dash_base() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let mut dashed = fixture.pairs[1].0.seq.clone();
    for i in [30, 31, 32, 60, 90] {
        dashed[i] = b'-';
    }
    let base = as_read_pair(&fixture.pairs[1]);
    let marked = ReadPair {
        name_r1: b"d",
        seq_r1: &dashed,
        name_r2: b"d",
        ..base
    };
    let first = as_read_pair(&fixture.pairs[0]);
    let (got, stats, off) = resolve_hand_marked(&idx, &opts, &[first, marked], &[None, Some(0)]);
    assert_eq!(
        stats,
        MemoStats {
            dup_pairs: 1,
            copied: 0,
            fallback_aligned: 1
        }
    );
    assert!(
        got == off,
        "a fallback pair with '-' diverged from memo off"
    );
}

/// Empty pair ranges around a representative's range do not confuse the
/// ordinal-to-range lookup at resolve.
#[test]
fn empty_ranges_around_the_representative() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pairs = fixture_pairs(&fixture);
    let dup = [pairs[3]];
    let chunks: [&[ReadPair<'_>]; 5] = [&[], &pairs[..4], &[], &[], &dup];
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut memo = PairMemo::new();
    let mut ranges = Vec::new();
    let mut p0 = 0;
    for chunk in chunks {
        let mut r = cohort.reserve_pairs(chunk.len()).unwrap();
        cohort.write_pairs(&mut r, chunk).unwrap();
        let reps = memo.mark_range(&r, chunk).unwrap();
        cohort
            .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &reps)
            .unwrap();
        ranges.push((r, p0));
        p0 += chunk.len();
    }
    assert_eq!(
        cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap(),
        all_copied(1)
    );
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    let mut got = Records::new();
    for (mut r, p0) in ranges {
        got.extend(emit_range(
            &cohort,
            &idx,
            &opts,
            &mut scratch,
            &mut r,
            &pestat,
            p0,
        ));
    }
    let all: Vec<ReadPair<'_>> = chunks.concat();
    let (off, _) = run_cohort_memo(&idx, &opts, &all, all.len(), false);
    assert!(got == off, "empty ranges changed the memo output");
}

/// A cohort with no marks resolves to zero stats whenever it is called, even
/// after its ranges were emitted.
#[test]
fn resolve_without_marks_is_a_no_op_even_after_emit() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pairs = fixture_pairs(&fixture);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut r = cohort.reserve_pairs(4).unwrap();
    cohort.write_pairs(&mut r, &pairs[..4]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut r)
        .unwrap();
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    emit_range(&cohort, &idx, &opts, &mut scratch, &mut r, &pestat, 0);
    assert_eq!(
        cohort.resolve_memo(&idx, &opts, &mut scratch).unwrap(),
        MemoStats::default()
    );
}

/// A single-end range does not depend on the memo, so it emits while pair
/// marks still await `resolve_memo`.
#[test]
fn singles_emit_while_pair_marks_are_unresolved() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pair = as_read_pair(&fixture.pairs[0]);
    let single = bwa_mem3_rs::SingleRead {
        name: &fixture.singles[0].name,
        seq: &fixture.singles[0].seq,
        qual: fixture.singles[0].qual.as_deref(),
    };
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut r = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut r, &[pair, pair]).unwrap();
    cohort
        .seed_extend_with_reps(&idx, &opts, &mut scratch, &mut r, &[None, Some(0)])
        .unwrap();
    let mut s = cohort.reserve_singles(1).unwrap();
    cohort.write_singles(&mut s, &[single]).unwrap();
    cohort
        .seed_extend(&idx, &opts, &mut scratch, &mut s)
        .unwrap();
    let mut sink = RecordVec::default();
    cohort
        .pair_emit(
            &idx,
            &opts,
            &mut scratch,
            &mut s,
            None,
            IdBases::default(),
            0,
            &mut sink,
        )
        .unwrap();
    assert!(!sink.records.is_empty());
}

/// An emit closes the memo just as `infer_cohort` does, even with a model from
/// elsewhere: a later mark is refused, and the cohort's other ranges still emit.
#[test]
fn marks_after_an_emit_are_refused_and_the_cohort_stays_usable() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pairs = fixture_pairs(&fixture);
    let (pestat, _) = serial_reference(&idx, &opts, &fixture);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut ranges = Vec::new();
    for chunk in pairs[..8].chunks(4) {
        let mut r = cohort.reserve_pairs(chunk.len()).unwrap();
        cohort.write_pairs(&mut r, chunk).unwrap();
        cohort
            .seed_extend(&idx, &opts, &mut scratch, &mut r)
            .unwrap();
        ranges.push(r);
    }
    emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut ranges[0],
        &pestat,
        0,
    );
    let mut late = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut late, &pairs[..2]).unwrap();
    assert_invalid(
        cohort.seed_extend_with_reps(&idx, &opts, &mut scratch, &mut late, &[None, Some(0)]),
        "memo is closed",
        "marks after an emit",
    );
    let records = emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut ranges[1],
        &pestat,
        4,
    );
    assert!(!records.is_empty());
}

/// Once `infer_cohort` has taken the cohort's model, no more duplicate marks
/// are accepted: otherwise a late range's marks would leave the cohort's
/// remaining ranges unemittable (emits wait for `resolve_memo`, which refuses
/// a cohort with an emitted range). The refusal leaves the cohort usable.
#[test]
fn marks_after_infer_cohort_are_refused_and_the_cohort_stays_usable() {
    let Some((idx, opts)) = index_for(false) else {
        return;
    };
    let fixture = field_fixture();
    let pairs = fixture_pairs(&fixture);
    let cohort = ResidentCohort::new(false).unwrap();
    let mut scratch = AlignScratch::new().unwrap();
    let mut ranges = Vec::new();
    for chunk in pairs[..8].chunks(4) {
        let mut r = cohort.reserve_pairs(chunk.len()).unwrap();
        cohort.write_pairs(&mut r, chunk).unwrap();
        cohort
            .seed_extend(&idx, &opts, &mut scratch, &mut r)
            .unwrap();
        ranges.push(r);
    }
    let pestat = cohort.infer_cohort(&idx, &opts).unwrap();
    emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut ranges[0],
        &pestat,
        0,
    );
    let mut late = cohort.reserve_pairs(2).unwrap();
    cohort.write_pairs(&mut late, &pairs[..2]).unwrap();
    assert_invalid(
        cohort.seed_extend_with_reps(&idx, &opts, &mut scratch, &mut late, &[None, Some(0)]),
        "memo is closed",
        "marks after infer_cohort",
    );
    // The remaining range still emits.
    let records = emit_range(
        &cohort,
        &idx,
        &opts,
        &mut scratch,
        &mut ranges[1],
        &pestat,
        4,
    );
    assert!(!records.is_empty());
}

