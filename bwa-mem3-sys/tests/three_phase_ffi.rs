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
