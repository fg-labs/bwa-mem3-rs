//! Differential test for the resident-cohort aligner: driving the public
//! resident C API (reserve a range per sub-chunk, write + seed-extend each, one
//! cohort `mem_pestat`, then pair-emit each range in order) must produce output
//! byte-identical to the legacy `bwa_shim_align_batch` over the same input, at
//! every sub-chunk granularity, and resident single-end emission must match the
//! legacy three-phase path. Each test builds a PhiX (or purpose-built) index and
//! skips when `bwa-mem3` is absent.
//!
//! ## Safety invariants shared by the `unsafe` blocks below
//! - `idx`/`opts`/`pes`/`sc`/`cohort` are valid, non-null handles from the
//!   matching shim constructor and stay alive for the whole call.
//! - `BwaReadPair`s/`BwaSingleRead`s borrow fixture-owned name/seq/qual buffers
//!   that outlive the call.
//! - A resident segment is only used through the cohort that reserved it, by one
//!   call at a time (these tests are single-threaded), and never after the
//!   cohort is freed.
//! - Each handle is freed exactly once, at the end of its test, and is not used
//!   afterward.

mod common;

use bwa_mem3_sys as sys;

/// Collects emitted records as `(origin_idx, [u32 block_size][body])`, the shape
/// `common::collect_records` returns for the legacy path.
#[derive(Default)]
struct Sink(Vec<(usize, Vec<u8>)>);

/// # Safety
/// `ctx` must be a valid `*mut Sink` that outlives the call, and `body`/`len`
/// must describe a readable buffer of `len` bytes -- both guaranteed by the
/// shim when it invokes this callback from a pair-emit call.
unsafe extern "C" fn sink_fn(
    ctx: *mut std::ffi::c_void,
    _kind: u32,
    idx: usize,
    body: *const u8,
    len: usize,
) {
    // SAFETY: `ctx` is the `&mut Sink` the test handed to pair-emit as context.
    let sink = unsafe { &mut *ctx.cast::<Sink>() };
    let mut rec = Vec::with_capacity(4 + len);
    rec.extend_from_slice(&u32::try_from(len).unwrap().to_le_bytes());
    // SAFETY: `body`/`len` describe a valid record buffer for this call.
    rec.extend_from_slice(unsafe { std::slice::from_raw_parts(body, len) });
    sink.0.push((idx, rec));
}

/// Run `pairs` through the public resident C API in sub-chunks of
/// `sub_batch_pairs` pairs (0 == one chunk) and collect the records, shaped like
/// `common::align_batch_records` so the two compare directly: each chunk's
/// origins and pair ids continue from the previous chunk's.
fn align_batch_records_resident(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
    sub_batch_pairs: usize,
) -> Vec<(usize, Vec<u8>)> {
    let chunk = if sub_batch_pairs == 0 {
        pairs.len().max(1)
    } else {
        sub_batch_pairs
    };
    // SAFETY: see the module's safety invariants for every call below.
    unsafe {
        let sc = sys::bwa_shim_scratch_new();
        let cohort = sys::bwa_shim_resident_cohort_new(0);
        assert!(
            !sc.is_null() && !cohort.is_null(),
            "{}",
            common::last_error()
        );
        let mut segments = Vec::new();
        for part in pairs.chunks(chunk) {
            let mut first = 0usize;
            let seg = sys::bwa_shim_resident_reserve_pairs(cohort, 2 * part.len(), &mut first);
            assert!(!seg.is_null(), "reserve failed: {}", common::last_error());
            let mut added = 0usize;
            for (i, pair) in part.iter().enumerate() {
                assert_eq!(
                    sys::bwa_shim_resident_write_pair(seg, i, pair, &mut added),
                    0,
                    "{}",
                    common::last_error()
                );
            }
            assert_eq!(
                sys::bwa_shim_resident_seed_extend(idx, opts, sc, seg),
                0,
                "{}",
                common::last_error()
            );
            segments.push((seg, first));
        }
        // The legacy path estimates the model over the whole batch; the
        // resident cohort's model covers the same pairs.
        let pes = common::new_pestat();
        assert_eq!(
            sys::bwa_shim_resident_pestat_cohort(idx, opts, cohort, pes),
            0,
            "{}",
            common::last_error()
        );
        let mut sink = Sink::default();
        for (seg, first) in segments {
            let pair_base = first / 2;
            let ids = sys::BwaIdBases {
                first_single_id: 0,
                first_pair_id: pair_base as u64,
            };
            let rc = sys::bwa_shim_resident_pair_emit(
                idx,
                opts,
                sc,
                seg,
                pes,
                ids,
                pair_base,
                Some(sink_fn),
                std::ptr::addr_of_mut!(sink).cast(),
            );
            assert_eq!(rc, 0, "{}", common::last_error());
        }
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_resident_cohort_free(cohort);
        sys::bwa_shim_scratch_free(sc);
        sink.0
    }
}

/// A range of sub-chunk sizes that brackets 1 read/chunk, the BATCH_SIZE/2
/// mate-rescue chunk boundary (256), a few multi-chunk sizes, the whole cohort,
/// and larger-than-cohort (one chunk) — plus 0 (one chunk of all pairs).
const SUB_BATCH_SIZES: [usize; 8] = [1, 3, 64, 256, 512, 1000, 4096, 0];

/// The core gate: unique-mapping FR pairs from PhiX. Every sub-chunk size must
/// reproduce the legacy single-call output byte-for-byte.
#[test]
fn resident_cohort_matches_legacy_at_every_sub_batch_size() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(1000, 150, 400, 41);
    let pairs = fx.pairs();

    let legacy = common::align_batch_records(idx, opts, &pairs);
    assert!(!legacy.is_empty(), "fixture produced no records");

    for sub in SUB_BATCH_SIZES {
        let got = align_batch_records_resident(idx, opts, &pairs, sub);
        assert_eq!(got.len(), legacy.len(), "sub={sub}: record count diverged");
        assert_eq!(
            got, legacy,
            "sub={sub}: resident-cohort bytes diverged from legacy align_batch"
        );
    }
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once, unused after.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Mate-rescue-heavy fixture: every 3rd R2 is mutated hard enough that seeding
/// finds nothing and only mate rescue from R1 can place it. This exercises the
/// resident `shim_resident_pair_emit_pairs` batched-rescue path, which reads
/// the resident 2-bit `seq` again at emit time — the reason the decode buffer
/// must stay resident from seed/extend to pair/emit rather than be per-worker
/// scratch.
/// The resident path must still match legacy at every sub-chunk size.
#[test]
fn resident_cohort_matches_legacy_on_rescue_heavy_fixture() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let mut fx = common::simulate(400, 100, 300, 61);
    for i in (0..fx.r2.len()).step_by(3) {
        for j in (0..fx.r2[i].len()).step_by(7) {
            fx.r2[i][j] = match fx.r2[i][j] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }
    let pairs = fx.pairs();
    let legacy = common::align_batch_records(idx, opts, &pairs);

    // Rescue must actually fire, or the comparison would pass vacuously.
    let mapped_r2 = legacy
        .iter()
        .filter(|(i, r)| {
            i % 3 == 0 && {
                let f = u16::from_le_bytes([r[4 + 14], r[4 + 15]]);
                f & 0x80 != 0 && f & 0x4 == 0
            }
        })
        .count();
    assert!(
        mapped_r2 > 50,
        "rescue did not fire ({mapped_r2} mutated R2 mapped)"
    );

    for sub in SUB_BATCH_SIZES {
        let got = align_batch_records_resident(idx, opts, &pairs, sub);
        assert_eq!(
            got, legacy,
            "sub={sub}: resident-cohort rescue path diverged from legacy"
        );
    }
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once, unused after.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Tandem-repeat reference: every pair maps to several EQUAL-score locations, so
/// which copy is primary (and the XA order of the rest) is decided by the
/// `id`-seeded tie-break hash (`id << 1 | i`). The resident driver derives each
/// sub-chunk's `first_pair_id` from its global pair offset; this fixture proves
/// those per-sub-chunk id bases reproduce the legacy `first_pair_id = 0..n`
/// sequence exactly — a bug in the id arithmetic would reorder primaries/XA and
/// this test would catch it where the unique-mapper fixtures cannot.
#[test]
fn resident_cohort_matches_legacy_on_tandem_repeat_ties() {
    let unit = &common::PHIX_SEQ.as_bytes()[500..800]; // 300 bp repeat unit
    let copies = 5usize;
    let mut refseq = Vec::with_capacity(unit.len() * copies);
    for _ in 0..copies {
        refseq.extend_from_slice(unit);
    }
    let Some((_dir, prefix)) = common::index_custom_ref(&refseq) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();

    let (read_len, insert, n) = (60usize, 120usize, 200usize);
    let mut names = Vec::new();
    let (mut r1, mut r2, mut qual) = (Vec::new(), Vec::new(), Vec::new());
    let mut rng = common::Rng(5);
    for i in 0..n {
        let start = (rng.next() as usize) % (unit.len() - insert + 1);
        r1.push(unit[start..start + read_len].to_vec());
        r2.push(common::revcomp(
            &unit[start + insert - read_len..start + insert],
        ));
        names.push(std::ffi::CString::new(format!("r{i}")).unwrap());
        qual.push(vec![b'I'; read_len]);
    }
    let pairs: Vec<sys::BwaReadPair> = (0..n)
        .map(|i| sys::BwaReadPair {
            r1_name: names[i].as_ptr(),
            r1_name_len: names[i].as_bytes().len(),
            r1_seq: r1[i].as_ptr(),
            r1_seq_len: r1[i].len(),
            r1_qual: qual[i].as_ptr(),
            r2_name: names[i].as_ptr(),
            r2_name_len: names[i].as_bytes().len(),
            r2_seq: r2[i].as_ptr(),
            r2_seq_len: r2[i].len(),
            r2_qual: qual[i].as_ptr(),
        })
        .collect();

    let legacy = common::align_batch_records(idx, opts, &pairs);
    for sub in SUB_BATCH_SIZES {
        let got = align_batch_records_resident(idx, opts, &pairs, sub);
        assert_eq!(
            got, legacy,
            "sub={sub}: resident-cohort tie-break output diverged from legacy \
             (per-sub-chunk first_pair_id arithmetic is wrong)"
        );
    }
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once, unused after.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Single-end reads: resident singles emission must reproduce the legacy
/// three-phase path (`bwa_shim_seed_extend` + `bwa_shim_pair_emit` over the same
/// singles) byte-for-byte, for mapped reads on both strands and unmapped ones.
#[test]
fn resident_singles_match_legacy_three_phase() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    // R1s and R2s of simulated FR pairs: forward and reverse mappers alike.
    let fx = common::simulate(300, 150, 400, 7);
    let pairs = fx.pairs();
    let singles: Vec<sys::BwaSingleRead> = pairs
        .iter()
        .flat_map(|p| {
            [
                sys::BwaSingleRead {
                    name: p.r1_name,
                    name_len: p.r1_name_len,
                    seq: p.r1_seq,
                    seq_len: p.r1_seq_len,
                    qual: p.r1_qual,
                },
                sys::BwaSingleRead {
                    name: p.r2_name,
                    name_len: p.r2_name_len,
                    seq: p.r2_seq,
                    seq_len: p.r2_seq_len,
                    qual: p.r2_qual,
                },
            ]
        })
        .collect();
    let ids = sys::BwaIdBases {
        first_single_id: 5,
        first_pair_id: 0,
    };

    // SAFETY: see the module's safety invariants for every call below.
    let (legacy, resident) = unsafe {
        let sc = sys::bwa_shim_scratch_new();
        let batch = sys::BwaReadBatch {
            pairs: std::ptr::null(),
            n_pairs: 0,
            singles: singles.as_ptr(),
            n_singles: singles.len(),
        };
        let regs = sys::bwa_shim_seed_extend(idx, opts, sc, &batch);
        assert!(!regs.is_null(), "{}", common::last_error());
        let mut legacy = Sink::default();
        let rc = sys::bwa_shim_pair_emit(
            idx,
            opts,
            sc,
            regs,
            std::ptr::null(),
            ids,
            Some(sink_fn),
            std::ptr::addr_of_mut!(legacy).cast(),
        );
        assert_eq!(rc, 0, "{}", common::last_error());

        let cohort = sys::bwa_shim_resident_cohort_new(0);
        let mut first = 0usize;
        let seg = sys::bwa_shim_resident_reserve_singles(cohort, singles.len(), &mut first);
        assert!(!seg.is_null(), "{}", common::last_error());
        let mut added = 0usize;
        for (i, read) in singles.iter().enumerate() {
            assert_eq!(
                sys::bwa_shim_resident_write_single(seg, i, read, &mut added),
                0
            );
        }
        assert_eq!(sys::bwa_shim_resident_seed_extend(idx, opts, sc, seg), 0);
        let mut resident = Sink::default();
        let rc = sys::bwa_shim_resident_pair_emit(
            idx,
            opts,
            sc,
            seg,
            std::ptr::null(),
            ids,
            0,
            Some(sink_fn),
            std::ptr::addr_of_mut!(resident).cast(),
        );
        assert_eq!(rc, 0, "{}", common::last_error());
        sys::bwa_shim_resident_cohort_free(cohort);
        sys::bwa_shim_scratch_free(sc);
        (legacy.0, resident.0)
    };
    assert!(!legacy.is_empty(), "fixture produced no records");
    assert_eq!(resident.len(), legacy.len(), "record count diverged");
    assert!(
        resident == legacy,
        "resident singles diverged from the legacy path"
    );
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once, unused after.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
