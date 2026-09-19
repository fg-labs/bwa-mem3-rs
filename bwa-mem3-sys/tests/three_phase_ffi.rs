//! C-level behavior tests for the shim, through the raw FFI. Each test builds
//! a PhiX index (skips when `bwa-mem3` is absent) and asserts on bytes.

mod common;

use bwa_mem3_sys as sys;

#[test]
fn version_and_build_info_are_exposed_from_c() {
    let v = unsafe { std::ffi::CStr::from_ptr(sys::bwa_shim_version()) };
    assert_eq!(v.to_str().unwrap(), sys::build_info::VERSION);
    let b = unsafe { std::ffi::CStr::from_ptr(sys::bwa_shim_build_info()) };
    assert!(b.to_str().unwrap().contains(sys::build_info::COMPILER));
}

/// `shim_seed_batch` wrote `n_threads = 1` and `flag |= MEM_F_PE` through the
/// shared `const mem_opt_t *` (bwa_shim_align.cpp:1105-1106). With one
/// `MemOpts` shared by many threads that is a data race. The alignment must
/// now run on a per-call copy and leave the caller's struct untouched — while
/// still producing paired output.
#[test]
fn align_batch_does_not_mutate_caller_opts() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    unsafe {
        (*opts).n_threads = 8;
        (*opts).flag &= !(sys::MEM_F_PE as i32);
    }
    let before = unsafe { std::ptr::read(opts) };

    let fx = common::simulate(64, 100, 300, 7);
    let recs = common::align_batch_records(idx, opts, &fx.pairs());

    let after = unsafe { std::ptr::read(opts) };
    let bytes = |o: &sys::mem_opt_t| unsafe {
        std::slice::from_raw_parts(
            (o as *const sys::mem_opt_t).cast::<u8>(),
            std::mem::size_of::<sys::mem_opt_t>(),
        )
        .to_vec()
    };
    assert_eq!(
        bytes(&before),
        bytes(&after),
        "mem_opt_t was mutated by align_batch"
    );
    assert_eq!(unsafe { (*opts).n_threads }, 8);

    // Still paired: every record carries FLAG 0x1 (offset 14 in the body after the u32 prefix).
    assert_eq!(recs.len(), 128, "one record per read for unique mappers");
    for (_, r) in &recs {
        let flag = u16::from_le_bytes([r[4 + 14], r[4 + 15]]);
        assert_ne!(flag & 0x1, 0, "record is not flagged paired");
    }
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

fn batch_of(pairs: &[sys::BwaReadPair]) -> sys::BwaReadBatch {
    sys::BwaReadBatch {
        pairs: pairs.as_ptr(),
        n_pairs: pairs.len(),
        singles: std::ptr::null(),
        n_singles: 0,
    }
}

/// The fused seed+extend must produce regs for every read, report the batch
/// shape, and account real heap bytes (the design's byte-bounded queues size
/// themselves from it).
#[test]
fn seed_extend_reports_shape_and_heap_bytes() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    assert!(!sc.is_null());

    let fx = common::simulate(50, 100, 300, 11);
    let pairs = fx.pairs();
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)) };
    assert!(
        !regs.is_null(),
        "seed_extend failed: {}",
        common::last_error()
    );
    unsafe {
        assert_eq!(sys::bwa_shim_regs_n_pairs(regs), 50);
        assert_eq!(sys::bwa_shim_regs_n_singles(regs), 0);
        // At least the copied names+seqs+quals (2 reads × (100 seq + 100 qual + ~3 name)).
        assert!(sys::bwa_shim_regs_heap_bytes(regs) >= 50 * 2 * 200);
        sys::bwa_shim_regs_free(regs);
        sys::bwa_shim_scratch_free(sc);
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Reusing one scratch across many calls must be byte-stable: run the legacy
/// path (which now goes through a scratch internally) and the new phase pair
/// twice each on one scratch, and require the four outputs identical.
#[test]
fn scratch_reuse_is_byte_stable_and_matches_legacy() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(300, 150, 400, 23);
    let pairs = fx.pairs();

    let legacy_a = common::align_batch_records(idx, opts, &pairs);
    let legacy_b = common::align_batch_records(idx, opts, &pairs);
    assert_eq!(legacy_a, legacy_b, "legacy path is not self-deterministic");

    for _ in 0..2 {
        let recs = three_phase(idx, opts, &pairs, 300, 0);
        assert_eq!(
            recs, legacy_a,
            "phase path diverged from legacy align_batch"
        );
    }
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Empty batch: valid, zero regs, no records, no error (matches the existing
/// `empty_batch.rs` contract of the legacy path).
#[test]
fn seed_extend_empty_batch_is_ok() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&[])) };
    assert!(!regs.is_null(), "{}", common::last_error());
    unsafe {
        assert_eq!(sys::bwa_shim_regs_n_pairs(regs), 0);
        sys::bwa_shim_regs_free(regs);
        sys::bwa_shim_scratch_free(sc);
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Sink that records `(origin_kind, origin_idx, [u32 block_size][body])` so
/// results compare directly with `common::collect_records`.
struct Sink(Vec<(u32, usize, Vec<u8>)>);

unsafe extern "C" fn sink_fn(
    ctx: *mut std::ffi::c_void,
    kind: u32,
    idx: usize,
    body: *const u8,
    len: usize,
) {
    let sink = &mut *ctx.cast::<Sink>();
    let mut rec = Vec::with_capacity(4 + len);
    rec.extend_from_slice(&(len as u32).to_le_bytes());
    rec.extend_from_slice(std::slice::from_raw_parts(body, len));
    sink.0.push((kind, idx, rec));
}

/// Run the three phases over `pairs` split into `sub` equal-size batches, one
/// cohort pestat, ids starting at `first_pair_id`. Returns legacy-shaped records.
fn three_phase(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
    sub: usize,
    first_pair_id: u64,
) -> Vec<(usize, Vec<u8>)> {
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let chunks: Vec<&[sys::BwaReadPair]> = pairs.chunks(sub.max(1)).collect();
    let regs: Vec<*mut sys::BwaRegs> = chunks
        .iter()
        .map(|c| {
            let r = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(c)) };
            assert!(!r.is_null(), "{}", common::last_error());
            r
        })
        .collect();
    let pes = common::new_pestat();
    let regs_const: Vec<*const sys::BwaRegs> = regs.iter().map(|r| r.cast_const()).collect();
    let rc = unsafe {
        sys::bwa_shim_pestat_cohort(idx, opts, regs_const.as_ptr(), regs_const.len(), pes)
    };
    assert_eq!(rc, 0, "{}", common::last_error());

    let mut out = Vec::new();
    let mut offset = 0usize;
    for (r, c) in regs.into_iter().zip(&chunks) {
        let mut sink = Sink(Vec::new());
        let ids = sys::BwaIdBases {
            first_single_id: 0,
            first_pair_id: first_pair_id + offset as u64,
        };
        let rc = unsafe {
            sys::bwa_shim_pair_emit(
                idx,
                opts,
                sc,
                r,
                pes,
                ids,
                Some(sink_fn),
                (&mut sink as *mut Sink).cast(),
            )
        };
        assert_eq!(rc, 0, "{}", common::last_error());
        for (kind, i, rec) in sink.0 {
            assert_eq!(kind, sys::BWA_ORIGIN_PAIR);
            out.push((offset + i, rec));
        }
        offset += c.len();
    }
    unsafe {
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_scratch_free(sc);
    }
    out
}

/// One cohort, one pestat, any sub-batch size: bytes identical to the legacy
/// single-call path (which is itself CLI-identical for one cohort).
#[test]
fn three_phase_is_invariant_to_sub_batch_size_and_matches_legacy() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(1000, 150, 400, 31);
    let pairs = fx.pairs();
    let legacy = common::align_batch_records(idx, opts, &pairs);
    for sub in [1usize, 7, 64, 512, 1000, 4096] {
        let got = three_phase(idx, opts, &pairs, sub, 0);
        assert_eq!(got.len(), legacy.len(), "sub={sub}: record count");
        assert_eq!(
            got, legacy,
            "sub={sub}: bytes diverged from legacy align_batch"
        );
    }
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// The pair id feeds hash_64 tie-breaks (mem_mark_primary_se / mem_pair,
/// seeded with `id<<1|i` at bwamem_pair.cpp:893-910). Offsetting `first_pair_id`
/// must therefore change *something* on a fixture with equal-score ties (proves
/// the id is wired), and legacy `align_batch` must equal `first_pair_id = 0`
/// (its own convention).
///
/// NB: `common::simulate` draws reads across the whole (repeat-free) PhiX
/// genome, where unique reads never tie, so the id would have no observable
/// effect there. This builds a tandem-repeat reference (a PhiX window repeated
/// several times) and draws reads from within one copy, so every pair maps to
/// several EQUAL-score locations; which copy is primary (and the XA order of
/// the rest) is exactly what the id-seeded hash decides.
#[test]
fn first_pair_id_reaches_the_tie_break_hash() {
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

    let (read_len, insert, n) = (60usize, 120usize, 50usize);
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

    let at_zero = three_phase(idx, opts, &pairs, n, 0);
    assert_eq!(at_zero, common::align_batch_records(idx, opts, &pairs));
    let shifted = three_phase(idx, opts, &pairs, n, 1 << 20);
    assert_eq!(shifted.len(), at_zero.len());
    assert_ne!(
        shifted, at_zero,
        "first_pair_id had no observable effect; the id is not plumbed"
    );
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Emission order and origin bookkeeping: for each pair, R1's records come
/// before R2's, primary before supplementary, and the origin index is the
/// batch-local pair index.
#[test]
fn pair_emit_orders_records_by_pair_then_side() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(20, 100, 300, 3);
    let pairs = fx.pairs();
    let recs = three_phase(idx, opts, &pairs, 20, 0);
    let mut last_pair = 0usize;
    for (pair_idx, _rec) in &recs {
        assert!(*pair_idx >= last_pair, "pair index went backwards");
        last_pair = *pair_idx;
    }
    // Each unique mapper yields exactly two records: flags 0x41 then 0x81 (+0x2/0x10/0x20).
    let mut per_pair = std::collections::BTreeMap::<usize, Vec<u16>>::new();
    for (pair_idx, rec) in &recs {
        per_pair
            .entry(*pair_idx)
            .or_default()
            .push(u16::from_le_bytes([rec[4 + 14], rec[4 + 15]]));
    }
    for (p, flags) in per_pair {
        assert_eq!(flags.len(), 2, "pair {p}: {flags:x?}");
        assert_ne!(
            flags[0] & 0x40,
            0,
            "pair {p}: first record is not R1: {flags:x?}"
        );
        assert_ne!(
            flags[1] & 0x80,
            0,
            "pair {p}: second record is not R2: {flags:x?}"
        );
    }
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// `pestat` may be omitted only for a pair-free batch; a paired batch without
/// a model is an error, not a silent per-batch re-estimate (gotcha #14).
#[test]
fn pair_emit_requires_pestat_for_pairs() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let fx = common::simulate(4, 100, 300, 9);
    let pairs = fx.pairs();
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)) };
    let mut sink = Sink(Vec::new());
    let ids = sys::BwaIdBases {
        first_single_id: 0,
        first_pair_id: 0,
    };
    let rc = unsafe {
        sys::bwa_shim_pair_emit(
            idx,
            opts,
            sc,
            regs,
            std::ptr::null(),
            ids,
            Some(sink_fn),
            (&mut sink as *mut Sink).cast(),
        )
    };
    assert_eq!(rc, -1);
    assert!(
        common::last_error().contains("pestat"),
        "{}",
        common::last_error()
    );
    assert!(sink.0.is_empty());
    unsafe {
        sys::bwa_shim_scratch_free(sc);
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
