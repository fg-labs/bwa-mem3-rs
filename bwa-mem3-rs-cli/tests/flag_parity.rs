//! Flag-parity tests for `pair_and_emit`.
//!
//! Aligns simulated PhiX paired reads via `bwa_mem3_rs::align_batch`, parses
//! the packed BAM bytes, and asserts the SAM flag bits for each record match
//! what `bwa-mem3 mem` would have produced on the same input. Regression
//! cover for:
//!
//! - 0x2  (properly paired) — our shim called `mem_pair` but ignored its
//!   return value, leaving 0x2 unset on all records.
//! - 0x10 (read is reverse-strand) — our flag block set mate-reverse (0x20)
//!   but forgot the symmetric self-reverse bit.
//!
//! Gated on `bwa-mem3` on PATH (needed for `index`).

mod common;
mod phix_seq;

use bwa_mem3_rs::{align_batch, BwaIndex, MemOpts, ReadPair};

/// Read a 16-bit little-endian flag from a packed BAM record.
/// Layout: [u32 block_size][refID i32][pos i32][l_read_name u8][mapq u8]
///         [bin u16][n_cigar u16][flag u16]...
fn flag_of(bam_rec_bytes: &[u8]) -> u16 {
    assert!(bam_rec_bytes.len() >= 4 + 16);
    // Skip 4-byte block_size + 4 (refID) + 4 (pos) + 1 (l_read_name) +
    // 1 (mapq) + 2 (bin) + 2 (n_cigar_op) = 18 bytes, then read 2-byte flag.
    let flag_off = 4 + 4 + 4 + 1 + 1 + 2 + 2;
    u16::from_le_bytes([bam_rec_bytes[flag_off], bam_rec_bytes[flag_off + 1]])
}

/// Read is primary (neither secondary 0x100 nor supplementary 0x800).
fn is_primary(flag: u16) -> bool {
    flag & 0x900 == 0
}

#[test]
fn primary_records_have_proper_pair_and_self_reverse_flags() {
    let Some(bwa) = common::require_bwa_mem3() else {
        return;
    };

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let ref_fa = common::setup_phix_index(dir, &bwa, phix_seq::PHIX_SEQ);

    // Simulate 200 FR-orientation pairs: R1 forward slice, R2 revcomp.
    // Reference contains >>500 positions so pairs are unique per alignment,
    // mapq stays high, and pairs should be flagged 0x2 (properly paired).
    let n_pairs = 200;
    let read_len = 100;
    let insert = 250;
    let pairs =
        common::simulate_pairs(phix_seq::PHIX_SEQ.as_bytes(), n_pairs, read_len, insert, 42);

    let idx = BwaIndex::load(&ref_fa).expect("load index");
    let mut opts = MemOpts::new().expect("opts");
    opts.set_pe(true);

    let qual = vec![b'I'; read_len];
    let name_bufs: Vec<Vec<u8>> = pairs
        .iter()
        .map(|(n, _, _)| n.as_bytes().to_vec())
        .collect();
    let read_pairs: Vec<ReadPair<'_>> = pairs
        .iter()
        .zip(&name_bufs)
        .map(|((_, s1, s2), name)| ReadPair {
            name_r1: name,
            seq_r1: s1,
            qual_r1: Some(&qual),
            name_r2: name,
            seq_r2: s2,
            qual_r2: Some(&qual),
        })
        .collect();

    let (aln, _) = align_batch(&idx, &opts, &read_pairs, None).expect("align");
    assert!(
        aln.len() >= 2 * n_pairs,
        "expected at least {} records, got {}",
        2 * n_pairs,
        aln.len()
    );

    // Tally flag bits across primary records.
    let mut n_primary = 0u32;
    let mut n_proper_pair = 0u32; // 0x2
    let mut n_self_reverse = 0u32; // 0x10
    let mut n_r1_primary = 0u32; // 0x40 (first in pair)
    let mut n_r2_primary = 0u32; // 0x80 (last in pair)

    for rec in aln.iter() {
        let f = flag_of(rec.bytes);
        if !is_primary(f) {
            continue;
        }
        n_primary += 1;
        if f & 0x2 != 0 {
            n_proper_pair += 1;
        }
        if f & 0x40 != 0 {
            n_r1_primary += 1;
        }
        if f & 0x80 != 0 {
            n_r2_primary += 1;
            if f & 0x10 != 0 {
                n_self_reverse += 1;
            }
        }
    }

    eprintln!(
        "primaries={n_primary}, proper_pair={n_proper_pair}, self_reverse={n_self_reverse}, \
         r1_primaries={n_r1_primary}, r2_primaries={n_r2_primary}"
    );

    // Sanity: we should have exactly 2 primaries per pair (one R1, one R2).
    assert_eq!(n_primary, (2 * n_pairs) as u32);
    assert_eq!(n_r1_primary, n_pairs as u32);
    assert_eq!(n_r2_primary, n_pairs as u32);

    // Regression for the 0x2 bug: every FR-oriented pair with both mates
    // mapped and within the inferred insert-size band should be marked
    // properly-paired. With 200 perfect reads on PhiX, essentially all
    // primaries should be.
    assert!(
        n_proper_pair >= (2 * n_pairs as u32) * 9 / 10,
        "expected >=90% of primaries to have 0x2 (properly paired), got {n_proper_pair}/{n_primary}"
    );

    // Regression for the 0x10 bug: R2s are reverse-complement of the ref,
    // so they map to reverse strand and must carry 0x10. All (or nearly all)
    // R2 primaries should have self-reverse set.
    assert!(
        n_self_reverse >= n_pairs as u32 * 9 / 10,
        "expected >=90% of pairs to have an R2 with 0x10 (read reverse strand), \
         got {n_self_reverse}/{n_pairs}"
    );
}
