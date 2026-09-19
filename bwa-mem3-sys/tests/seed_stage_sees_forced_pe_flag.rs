//! Regression test for a fix-round finding on top of the const-correctness
//! change: `shim_seed_batch` copies the caller's `mem_opt_t` into
//! `s->opts_copy` and forces `n_threads = 1` / `flag |= MEM_F_PE` there (so
//! `bwa_shim_align_batch` never writes through the caller's pointer), and
//! `shim_extend_batch`/`shim_estimate_pestat` read that forced copy back via
//! `s->opts`. But `s->w.opt` — what the SEEDING kernel (`mem_kernel1_core`)
//! actually reads options through — was still set to the caller's raw,
//! un-forced `opts`. `mem_kernel1_core` reads `opt->flag & MEM_F_PE` to
//! decide whether the `--extend-mate-concordant` chain cap runs
//! (`bwamem.cpp:2501-2521`), so with the bug present the forced PE flag never
//! reached the seeding stage: a caller whose own `opts.flag` had `MEM_F_PE`
//! off got the non-mate-aware cap at seeding regardless of the shim's
//! always-paired semantics for `bwa_shim_align_batch`.
//!
//! This builds a small reference with one 150bp block duplicated at two
//! well-separated loci, so a 100bp read drawn from either copy matches both
//! equally well (an exact weight tie). With `max_extend_chains = 1` (so only
//! the top-1 chain survives ordinary capping) and `mate_concordant_window`
//! set, the *correct* disambiguator is which copy is concordant with the
//! mate — but with an exact tie, the non-mate-aware cap deterministically
//! keeps the lower-index (leftmost) chain regardless of concordance. The
//! mate is placed uniquely near the SECOND (rightmost) copy, so:
//!   - if the seeding kernel sees the forced PE flag, the mate-aware cap
//!     keeps (also) the concordant, rightmost chain, and read 1 lands next
//!     to its mate (small, correct TLEN);
//!   - if it does not (the bug), only the leftmost chain survives, and read 1
//!     lands far from its mate (large, wrong TLEN) even though the caller
//!     asked for paired alignment via `bwa_shim_align_batch`.
//!
//! The caller's `mem_opt_t.flag` has `MEM_F_PE` OFF throughout (mirroring
//! the bug scenario precisely): `bwa_shim_align_batch`'s contract is that it
//! always aligns as pairs internally regardless of the caller's flag, so this
//! also doubles as a check that the "always paired" semantic actually reaches
//! seeding, not just pairing/pestat.

mod common;

use bwa_mem3_sys as sys;
use common::Rng;

/// Body-offset layout of a packed BAM record (after the 4-byte block_size
/// prefix `collect_records` already strips into the `Vec<u8>` boundary —
/// note `three_phase_ffi.rs` indexes the SAME buffer with a `4 +` offset,
/// i.e. here `body` still includes that leading `[u32 block_size]`).
fn record_pos(body: &[u8]) -> i32 {
    i32::from_le_bytes(body[4 + 4..4 + 8].try_into().unwrap())
}
fn record_flag(body: &[u8]) -> u16 {
    u16::from_le_bytes([body[4 + 14], body[4 + 15]])
}
fn record_tlen(body: &[u8]) -> i32 {
    i32::from_le_bytes(body[4 + 28..4 + 32].try_into().unwrap())
}

const FIRST_IN_PAIR: u16 = 0x40;

/// One 150bp block duplicated at two well-separated loci, with unique random
/// flanking sequence everywhere else. Returns `(reference, repeat_copy1_start,
/// repeat_copy2_start)` (0-based, into `reference`).
fn build_repeat_reference() -> (Vec<u8>, usize, usize) {
    let mut rng = Rng(1);
    let gen = |rng: &mut Rng, n: usize| -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
    };
    let unique_a = gen(&mut rng, 300);
    let repeat = gen(&mut rng, 150);
    let unique_b = gen(&mut rng, 300);
    let unique_c = gen(&mut rng, 300);
    let copy1_start = unique_a.len();
    let copy2_start = copy1_start + repeat.len() + unique_b.len();
    let reference: Vec<u8> = [unique_a, repeat.clone(), unique_b, repeat, unique_c].concat();
    (reference, copy1_start, copy2_start)
}

#[test]
fn seeding_sees_the_forced_pe_flag_not_the_callers_raw_flag() {
    let (reference, copy1_start, copy2_start) = build_repeat_reference();
    let Some((_dir, prefix)) = common::index_custom_ref(&reference) else {
        return;
    };
    let idx = common::load_idx(&prefix);

    // Read 1: the 100bp window at offset 25 into the repeat block. Both
    // copies are byte-identical, so this read ties in weight between
    // `copy1_start + 25` and `copy2_start + 25`.
    let read_len = 100usize;
    let r1_start = copy2_start + 25;
    let r1 = reference[r1_start..r1_start + read_len].to_vec();
    // Read 2 (mate): uniquely anchored 100bp AFTER copy 2, so it is only ever
    // concordant with the rightmost repeat copy (small TLEN), never the
    // leftmost one (which would be ~450bp further away).
    let insert = 225usize;
    let r2_start = r1_start + insert - read_len;
    let r2 = common::revcomp(&reference[r2_start..r2_start + read_len]);

    let name = std::ffi::CString::new("r0").unwrap();
    let qual = vec![b'I'; read_len];
    let pair = sys::BwaReadPair {
        r1_name: name.as_ptr(),
        r1_name_len: name.as_bytes().len(),
        r1_seq: r1.as_ptr(),
        r1_seq_len: r1.len(),
        r1_qual: qual.as_ptr(),
        r2_name: name.as_ptr(),
        r2_name_len: name.as_bytes().len(),
        r2_seq: r2.as_ptr(),
        r2_seq_len: r2.len(),
        r2_qual: qual.as_ptr(),
    };

    let opts = common::new_opts();
    unsafe {
        // Mirror the bug scenario exactly: the caller's own flag has
        // MEM_F_PE off. `bwa_shim_align_batch` must still align as a pair
        // (that's its whole contract), including at the seeding stage.
        (*opts).flag &= !(sys::MEM_F_PE as i32);
        (*opts).max_extend_chains = 1;
        (*opts).mate_concordant_window = 300; // fixed window, bp
    }

    let recs = common::align_batch_records(idx, opts, &[pair]);
    let r1_rec = recs
        .iter()
        .map(|(_, b)| b)
        .find(|b| record_flag(b) & FIRST_IN_PAIR != 0)
        .expect("read 1's record must be present");

    assert_eq!(
        record_pos(r1_rec) as usize,
        r1_start,
        "read 1 must land at the mate-concordant (rightmost) repeat copy, \
         not the leftmost one the non-mate-aware cap would keep on an exact \
         weight tie -- if this is {copy1_start}, the seeding kernel is still \
         reading the caller's raw, un-forced MEM_F_PE instead of the shim's \
         forced copy"
    );
    assert_eq!(
        record_tlen(r1_rec).unsigned_abs() as usize,
        insert,
        "read 1's TLEN must reflect pairing with its mate at the concordant \
         locus, not the ~450bp-further discordant one"
    );

    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
