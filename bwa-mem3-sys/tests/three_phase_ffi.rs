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

    let sc = unsafe { sys::bwa_shim_scratch_new() };
    for _ in 0..2 {
        let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)) };
        assert!(!regs.is_null(), "{}", common::last_error());
        // Task 4 replaces this with bwa_shim_pair_emit; until then the legacy
        // extend entry point consumes a BwaSeeds, so wrap the regs in one.
        let seeds = unsafe { sys::bwa_shim_seeds_from_regs(regs) };
        let pes = common::new_pestat();
        let b = unsafe {
            sys::bwa_shim_extend_batch(
                idx,
                opts,
                seeds,
                pairs.as_ptr(),
                pairs.len(),
                std::ptr::null(),
                pes,
            )
        };
        assert!(!b.is_null(), "{}", common::last_error());
        let recs = common::collect_records(b);
        unsafe {
            sys::bwa_shim_batch_free(b);
            sys::bwa_shim_pestat_free(pes);
        }
        assert_eq!(
            recs, legacy_a,
            "phase path diverged from legacy align_batch"
        );
    }
    unsafe {
        sys::bwa_shim_scratch_free(sc);
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
