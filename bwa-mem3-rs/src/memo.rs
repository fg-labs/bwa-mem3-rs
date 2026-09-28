//! Exact-duplicate read-pair detection for one [`ResidentCohort`], feeding
//! [`ResidentCohort::seed_extend_with_reps`] and
//! [`ResidentCohort::resolve_memo`].
//!
//! bwa-mem3's CLI `--dedup-reads` aligns each distinct read pair once per `-K`
//! chunk and copies the result to its exact duplicates; the output is
//! byte-identical because every duplicate keeps its own name, qualities and
//! pair id through pairing and emission. [`PairMemo`] is the detection half
//! of that, run by the caller in cohort order: a pair is a duplicate of the
//! first earlier pair of the cohort with the same bases in both mates (names
//! and qualities are not compared). Keys are 64-bit hashes; the shim checks
//! the bases themselves at `resolve_memo`, so a collision costs one normal
//! alignment, never wrong output.
//!
//! [`ResidentCohort`]: crate::ResidentCohort
//! [`ResidentCohort::seed_extend_with_reps`]: crate::ResidentCohort::seed_extend_with_reps
//! [`ResidentCohort::resolve_memo`]: crate::ResidentCohort::resolve_memo

use std::collections::hash_map::{DefaultHasher, Entry, HashMap};
use std::hash::{BuildHasherDefault, Hash, Hasher};

use crate::align::ReadPair;
use crate::error::{Error, Result};
use crate::resident::ResidentRange;

/// Finds exact duplicate read pairs within one resident cohort. Create one per
/// cohort and mark every pair range in reservation order.
#[derive(Debug, Default)]
pub struct PairMemo {
    /// Cohort id of the first range marked; later ranges must match.
    cohort_id: Option<u64>,
    /// Hash of a pair's bases -> the cohort-local ordinal of its first occurrence.
    first: HashMap<u64, u64, BuildHasherDefault<PassThrough>>,
    /// The ordinal the next range must start at.
    next: u64,
    /// Pairs marked as duplicates so far.
    dups: u64,
}

impl PairMemo {
    /// An empty memo.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Marks for one pair range: element `i` is `Some(rep)` when pair `i` of
    /// `range` has the same bases in both mates as the earlier cohort pair
    /// `rep`, else `None`. Pair `i`'s ordinal is `range.first() / 2 + i`.
    ///
    /// # Errors
    ///
    /// `InvalidInput` when `range` is not a pair range, `pairs.len()` differs
    /// from `range.n_pairs()`, `range` belongs to a different cohort than the
    /// ranges marked before it, or ranges are not marked in reservation order
    /// with none skipped.
    pub fn mark_range(
        &mut self,
        range: &ResidentRange,
        pairs: &[ReadPair<'_>],
    ) -> Result<Vec<Option<u64>>> {
        if !range.is_pairs() {
            return Err(Error::InvalidInput("mark_range: not a pair range".into()));
        }
        if pairs.len() != range.n_pairs() {
            return Err(Error::InvalidInput(format!(
                "mark_range: {} pairs for a range of {}",
                pairs.len(),
                range.n_pairs()
            )));
        }
        if self.cohort_id.is_some_and(|id| id != range.cohort_id()) {
            return Err(Error::InvalidInput(
                "mark_range: range belongs to a different ResidentCohort".into(),
            ));
        }
        let start = (range.first() / 2) as u64;
        if start != self.next {
            return Err(Error::InvalidInput(format!(
                "mark_range: range starts at pair {start} but pair {} is next; ranges must be \
                 marked in reservation order, with none skipped",
                self.next
            )));
        }
        // Bind only once the call is accepted, so a rejected one changes nothing.
        self.cohort_id = Some(range.cohort_id());
        let reps = pairs
            .iter()
            .enumerate()
            .map(|(i, p)| match self.first.entry(pair_key(p)) {
                Entry::Occupied(e) => {
                    self.dups += 1;
                    Some(*e.get())
                }
                Entry::Vacant(v) => {
                    v.insert(start + i as u64);
                    None
                }
            })
            .collect();
        self.next = start + pairs.len() as u64;
        Ok(reps)
    }

    /// Pairs marked so far (duplicates or not).
    #[must_use]
    pub fn pairs_marked(&self) -> u64 {
        self.next
    }

    /// Pairs marked as duplicates so far.
    #[must_use]
    pub fn dup_pairs(&self) -> u64 {
        self.dups
    }
}

/// Hash of both mates' bases. `[u8]`'s `Hash` writes the length first, so the
/// boundary between the mates is part of the key. SipHash (std's default)
/// mixes long byte strings well; a collision is still safe (see the module doc).
fn pair_key(p: &ReadPair<'_>) -> u64 {
    let mut h = DefaultHasher::new();
    p.seq_r1.hash(&mut h);
    p.seq_r2.hash(&mut h);
    h.finish()
}

/// The map's hasher: its keys already are 64-bit SipHash values
/// ([`pair_key`]), so pass them through rather than hash them a second time.
#[derive(Debug, Default)]
struct PassThrough(u64);

impl Hasher for PassThrough {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write_u64(&mut self, key: u64) {
        self.0 = key;
    }
    /// Only `u64` keys reach this map, which `Hash` writes with `write_u64`;
    /// fold anything else in rather than drop it.
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(b);
        }
    }
}

/// What [`ResidentCohort::resolve_memo`](crate::ResidentCohort::resolve_memo)
/// did: every duplicate pair was either copied from its representative or,
/// when their bases differed, aligned normally. `dup_pairs == copied +
/// fallback_aligned`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MemoStats {
    /// Pairs marked as duplicates across the cohort.
    pub dup_pairs: u64,
    /// Duplicates whose alignment regions were copied from their representative.
    pub copied: u64,
    /// Duplicates whose bases differed from their representative's (a hash
    /// collision or a hand-built mark), aligned normally instead.
    pub fallback_aligned: u64,
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::resident::ResidentCohort;

    /// One pair from owned mates, with a fixed name and no QUAL: the memo
    /// never looks at either.
    fn pair<'a>(r1: &'a [u8], r2: &'a [u8]) -> ReadPair<'a> {
        ReadPair {
            name_r1: b"q",
            seq_r1: r1,
            qual_r1: None,
            name_r2: b"q",
            seq_r2: r2,
            qual_r2: None,
        }
    }

    /// Mark `mates` (each `"R1/R2"`) as one range of a fresh cohort and return
    /// the marks.
    fn mark_one_range(mates: &[&str]) -> Vec<Option<u64>> {
        let cohort = ResidentCohort::new(false).unwrap();
        let range = cohort.reserve_pairs(mates.len()).unwrap();
        let pairs: Vec<ReadPair<'_>> = mates
            .iter()
            .map(|m| {
                let (a, b) = m.split_once('/').unwrap();
                pair(a.as_bytes(), b.as_bytes())
            })
            .collect();
        PairMemo::new().mark_range(&range, &pairs).unwrap()
    }

    #[rstest]
    #[case::exact_duplicate(&["ACGT/TTGA", "ACGT/TTGA"], vec![None, Some(0)])]
    #[case::triplicate_points_at_first(
        &["ACGT/TTGA", "ACGT/TTGA", "ACGT/TTGA"],
        vec![None, Some(0), Some(0)]
    )]
    #[case::swapped_mates_differ(&["ACGT/TTGA", "TTGA/ACGT"], vec![None, None])]
    #[case::one_mate_identical_differs(&["ACGT/TTGA", "ACGT/TTGC"], vec![None, None])]
    #[case::length_differs(&["ACGT/TTGA", "ACGTA/TTGA"], vec![None, None])]
    #[case::boundary_shift_differs(&["ACG/TTTGA", "ACGT/TTGA"], vec![None, None])]
    #[case::case_differs(&["ACGT/TTGA", "acgt/TTGA"], vec![None, None])]
    #[case::n_differs(&["ACGT/TTGA", "ACNT/TTGA"], vec![None, None])]
    #[case::empty_mates_match(&["/TTGA", "/TTGA"], vec![None, Some(0)])]
    fn marks_within_one_range(#[case] mates: &[&str], #[case] expected: Vec<Option<u64>>) {
        assert_eq!(mark_one_range(mates), expected);
    }

    /// Ordinals are cohort-local pair indices from the range (`first() / 2 + i`),
    /// so a duplicate in the second range names its representative's ordinal
    /// in the first.
    #[test]
    fn ordinals_span_ranges() {
        let cohort = ResidentCohort::new(false).unwrap();
        let a = cohort.reserve_pairs(2).unwrap();
        let b = cohort.reserve_pairs(2).unwrap();
        let mut memo = PairMemo::new();
        let first = memo
            .mark_range(&a, &[pair(b"AAAA", b"CCCC"), pair(b"GGGG", b"TTTT")])
            .unwrap();
        let second = memo
            .mark_range(&b, &[pair(b"GGGG", b"TTTT"), pair(b"AAAA", b"CCCC")])
            .unwrap();
        assert_eq!(first, vec![None, None]);
        assert_eq!(second, vec![Some(1), Some(0)]);
        assert_eq!((memo.pairs_marked(), memo.dup_pairs()), (4, 2));
    }

    /// An empty range is legal and advances nothing.
    #[test]
    fn empty_range_is_a_no_op() {
        let cohort = ResidentCohort::new(false).unwrap();
        let empty = cohort.reserve_pairs(0).unwrap();
        let next = cohort.reserve_pairs(1).unwrap();
        let mut memo = PairMemo::new();
        assert_eq!(
            memo.mark_range(&empty, &[]).unwrap(),
            Vec::<Option<u64>>::new()
        );
        assert_eq!(
            memo.mark_range(&next, &[pair(b"A", b"C")]).unwrap(),
            vec![None]
        );
    }

    /// A call rejected for its order does not bind the memo to that range's
    /// cohort: the memo still takes another cohort's first range.
    #[test]
    fn a_rejected_mark_binds_nothing() {
        let a = ResidentCohort::new(false).unwrap();
        let b = ResidentCohort::new(false).unwrap();
        let _a0 = a.reserve_pairs(1).unwrap();
        let a1 = a.reserve_pairs(1).unwrap();
        let mut memo = PairMemo::new();
        let p = pair(b"ACGT", b"TTGA");
        assert!(memo.mark_range(&a1, &[p]).is_err());
        assert_eq!(
            memo.mark_range(&b.reserve_pairs(1).unwrap(), &[p]).unwrap(),
            vec![None]
        );
    }

    #[derive(Debug, Clone, Copy)]
    enum MarkMisuse {
        OutOfOrder,
        Gap,
        LengthMismatch,
        SingleRange,
        OtherCohort,
    }

    #[rstest]
    #[case::out_of_order(MarkMisuse::OutOfOrder, "reservation order")]
    #[case::gap(MarkMisuse::Gap, "reservation order")]
    #[case::length_mismatch(MarkMisuse::LengthMismatch, "for a range of")]
    #[case::single_range(MarkMisuse::SingleRange, "not a pair range")]
    #[case::other_cohort(MarkMisuse::OtherCohort, "different ResidentCohort")]
    fn misuse_is_rejected(#[case] misuse: MarkMisuse, #[case] needle: &str) {
        let cohort = ResidentCohort::new(false).unwrap();
        let other = ResidentCohort::new(false).unwrap();
        let p = pair(b"ACGT", b"TTGA");
        let mut memo = PairMemo::new();
        let result = match misuse {
            MarkMisuse::OutOfOrder => {
                let _a = cohort.reserve_pairs(1).unwrap();
                let b = cohort.reserve_pairs(1).unwrap();
                memo.mark_range(&b, &[p])
            }
            MarkMisuse::Gap => {
                let a = cohort.reserve_pairs(1).unwrap();
                let _b = cohort.reserve_pairs(1).unwrap();
                let c = cohort.reserve_pairs(1).unwrap();
                memo.mark_range(&a, &[p]).unwrap();
                memo.mark_range(&c, &[p])
            }
            MarkMisuse::LengthMismatch => memo.mark_range(&cohort.reserve_pairs(2).unwrap(), &[p]),
            MarkMisuse::SingleRange => memo.mark_range(&cohort.reserve_singles(1).unwrap(), &[]),
            MarkMisuse::OtherCohort => {
                memo.mark_range(&cohort.reserve_pairs(1).unwrap(), &[p])
                    .unwrap();
                // Start `other`'s range at pair 1, where this memo expects
                // the next range, so only the cohort check can reject it.
                let _pad = other.reserve_pairs(1).unwrap();
                memo.mark_range(&other.reserve_pairs(1).unwrap(), &[p])
            }
        };
        match result {
            Err(Error::InvalidInput(msg)) => {
                assert!(msg.contains(needle), "{misuse:?}: got {msg:?}");
            }
            other => panic!("{misuse:?}: expected InvalidInput naming {needle:?}, got {other:?}"),
        }
    }
}
