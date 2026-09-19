//! C-level behavior tests for the shim, through the raw FFI. Each test builds
//! a PhiX index (skips when `bwa-mem3` is absent) and asserts on bytes.
//!
//! ## Safety invariants shared by the `unsafe` blocks below
//!
//! Every raw shim call here relies on the same handful of invariants, noted
//! tersely per block and spelled out once here (per CONTRIBUTING.md, each
//! `unsafe` block still carries its own comment):
//! - `idx`/`opts`/`sc`/`pes`/`regs` are valid, non-null handles returned by the
//!   matching shim constructor (`common::load_idx`/`common::new_opts`/
//!   `bwa_shim_scratch_new`/`common::new_pestat`/`bwa_shim_seed_extend`) and
//!   stay alive for the whole call.
//! - Read pointers -- `BwaReadPair`/`BwaReadBatch`/`BwaSingleRead` and the
//!   name/seq/qual buffers they reference -- borrow fixture-owned data that
//!   outlives the call, and every `n_*` count matches its backing slice.
//! - Each handle is freed exactly once, at the end of its test, and is not used
//!   afterward.
//! - The sink callback `sink_fn` receives a `ctx` that is a live `&mut Sink` for
//!   the call and a `body`/`len` describing a buffer valid for the callback.

mod common;

use bwa_mem3_sys as sys;

#[test]
fn version_and_build_info_are_exposed_from_c() {
    // SAFETY: both return pointers to static NUL-terminated C strings valid for
    // the whole program, so wrapping them in `CStr` is sound.
    let v = unsafe { std::ffi::CStr::from_ptr(sys::bwa_shim_version()) };
    assert_eq!(v.to_str().unwrap(), sys::build_info::VERSION);
    // SAFETY: static NUL-terminated C string valid for the program lifetime.
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
    // SAFETY: `opts` is a valid, uniquely-owned `mem_opt_t` from `new_opts`, so
    // writing its POD fields through the raw pointer is sound.
    unsafe {
        (*opts).n_threads = 8;
        (*opts).flag &= !(sys::MEM_F_PE as i32);
    }
    // SAFETY: `opts` points to an initialized `mem_opt_t`; `ptr::read` takes a
    // bitwise POD snapshot for the before/after byte comparison.
    let before = unsafe { std::ptr::read(opts) };

    let fx = common::simulate(64, 100, 300, 7);
    let recs = common::align_batch_records(idx, opts, &fx.pairs());

    // SAFETY: `opts` is still the valid initialized struct; bitwise POD snapshot.
    let after = unsafe { std::ptr::read(opts) };
    // SAFETY: `o` is a live `mem_opt_t`; viewing its `size_of::<mem_opt_t>()`
    // bytes as a `u8` slice is valid for reading a POD struct's bytes.
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
    // SAFETY: `opts` is still valid; this reads a single field.
    assert_eq!(unsafe { (*opts).n_threads }, 8);

    // Still paired: every record carries FLAG 0x1 (offset 14 in the body after the u32 prefix).
    assert_eq!(recs.len(), 128, "one record per read for unique mappers");
    for (_, r) in &recs {
        let flag = u16::from_le_bytes([r[4 + 14], r[4 + 15]]);
        assert_ne!(flag & 0x1, 0, "record is not flagged paired");
    }
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
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
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    assert!(!sc.is_null());

    let fx = common::simulate(50, 100, 300, 11);
    let pairs = fx.pairs();
    // SAFETY: valid `idx`/`opts`/`sc` handles; `batch_of(&pairs)` borrows the
    // live `pairs` slice (fixture-owned) for the duration of the call.
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)) };
    assert!(
        !regs.is_null(),
        "seed_extend failed: {}",
        common::last_error()
    );
    // SAFETY: `regs`/`sc`/`opts`/`idx` are the live owned handles; the accessors
    // only read `regs`, and each handle is then freed once and unused after.
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
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
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
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    // SAFETY: valid `idx`/`opts`/`sc`; `batch_of(&[])` is an empty (zero-count)
    // batch with null read pointers, which the shim accepts as a valid no-op.
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&[])) };
    assert!(!regs.is_null(), "{}", common::last_error());
    // SAFETY: live owned handles; the accessor reads `regs`, then each is freed
    // once and unused afterward.
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

/// # Safety
/// `ctx` must be a valid `*mut Sink` that outlives the call, and `body`/`len`
/// must describe a readable buffer of `len` bytes -- both guaranteed by the
/// shim when it invokes this callback from `bwa_shim_pair_emit`.
unsafe extern "C" fn sink_fn(
    ctx: *mut std::ffi::c_void,
    kind: u32,
    idx: usize,
    body: *const u8,
    len: usize,
) {
    // SAFETY: `ctx` is the `&mut Sink` the test handed to `pair_emit` as context.
    let sink = &mut *ctx.cast::<Sink>();
    let mut rec = Vec::with_capacity(4 + len);
    rec.extend_from_slice(&(len as u32).to_le_bytes());
    // SAFETY: `body`/`len` describe a valid record buffer for this call (see the
    // `# Safety` section); the bytes are copied out immediately.
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
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let chunks: Vec<&[sys::BwaReadPair]> = pairs.chunks(sub.max(1)).collect();
    let regs: Vec<*mut sys::BwaRegs> = chunks
        .iter()
        .map(|c| {
            // SAFETY: valid `idx`/`opts`/`sc`; `batch_of(c)` borrows the live
            // chunk slice `c` (fixture-owned) for the duration of the call.
            let r = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(c)) };
            assert!(!r.is_null(), "{}", common::last_error());
            r
        })
        .collect();
    let pes = common::new_pestat();
    let regs_const: Vec<*const sys::BwaRegs> = regs.iter().map(|r| r.cast_const()).collect();
    // SAFETY: valid `idx`/`opts`/`pes`; `regs_const` holds `regs_const.len()`
    // live `BwaRegs` pointers (owned by `regs`, still alive), and the length
    // matches the pointer array.
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
        // SAFETY: valid `idx`/`opts`/`sc`/`r`/`pes` handles; `sink_fn`'s context
        // is `&mut sink`, a live local for the whole call (see `sink_fn`'s
        // `# Safety`); `ids` is a plain POD value.
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
    // SAFETY: `pes`/`sc` are the live owned handles; each freed once, not used
    // afterward.
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
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
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
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
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
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Owned single-end reads borrowed as FFI structs.
struct Singles {
    names: Vec<std::ffi::CString>,
    seqs: Vec<Vec<u8>>,
    quals: Vec<Vec<u8>>,
}

fn simulate_singles(n: usize, read_len: usize, seed: u64) -> Singles {
    let reference = common::PHIX_SEQ.as_bytes();
    let mut rng = common::Rng(seed);
    let mut s = Singles {
        names: Vec::new(),
        seqs: Vec::new(),
        quals: Vec::new(),
    };
    for i in 0..n {
        let start = (rng.next() as usize) % (reference.len() - read_len);
        let mut seq = reference[start..start + read_len].to_vec();
        if rng.next() % 2 == 0 {
            seq = common::revcomp(&seq);
        }
        s.names
            .push(std::ffi::CString::new(format!("s{i}")).unwrap());
        s.seqs.push(seq);
        s.quals.push(vec![b'I'; read_len]);
    }
    s
}

impl Singles {
    fn ffi(&self) -> Vec<sys::BwaSingleRead> {
        (0..self.names.len())
            .map(|i| sys::BwaSingleRead {
                name: self.names[i].as_ptr(),
                name_len: self.names[i].as_bytes().len(),
                seq: self.seqs[i].as_ptr(),
                seq_len: self.seqs[i].len(),
                qual: self.quals[i].as_ptr(),
            })
            .collect()
    }
}

fn run_mixed(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
    singles: &[sys::BwaSingleRead],
    ids: sys::BwaIdBases,
) -> Vec<(u32, usize, Vec<u8>)> {
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let batch = sys::BwaReadBatch {
        pairs: pairs.as_ptr(),
        n_pairs: pairs.len(),
        singles: singles.as_ptr(),
        n_singles: singles.len(),
    };
    // SAFETY: valid `idx`/`opts`/`sc`; `batch` names live `pairs`/`singles`
    // slices with matching `n_pairs`/`n_singles` counts, alive for the call.
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch) };
    assert!(!regs.is_null(), "{}", common::last_error());
    // SAFETY: `regs` is the live handle just returned; the accessors only read it.
    unsafe {
        assert_eq!(sys::bwa_shim_regs_n_pairs(regs), pairs.len());
        assert_eq!(sys::bwa_shim_regs_n_singles(regs), singles.len());
    }
    let pes = common::new_pestat();
    let regs_c: [*const sys::BwaRegs; 1] = [regs.cast_const()];
    // SAFETY: valid `idx`/`opts`/`pes`; `regs_c` is a 1-element array of the live
    // `regs` pointer and the count `1` matches it.
    let rc = unsafe { sys::bwa_shim_pestat_cohort(idx, opts, regs_c.as_ptr(), 1, pes) };
    assert_eq!(rc, 0);
    let mut sink = Sink(Vec::new());
    let pes_arg = if pairs.is_empty() {
        std::ptr::null()
    } else {
        pes.cast_const()
    };
    // SAFETY: valid `idx`/`opts`/`sc`/`regs`; `pes_arg` is the live pestat for a
    // paired batch or null for a pair-free one (the shim's contract); `sink_fn`'s
    // context is `&mut sink`, live for the whole call (see `sink_fn`'s `# Safety`).
    let rc = unsafe {
        sys::bwa_shim_pair_emit(
            idx,
            opts,
            sc,
            regs,
            pes_arg,
            ids,
            Some(sink_fn),
            (&mut sink as *mut Sink).cast(),
        )
    };
    assert_eq!(rc, 0, "{}", common::last_error());
    // SAFETY: `pes`/`sc` are the live owned handles; each freed once, not used
    // afterward.
    unsafe {
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_scratch_free(sc);
    }
    sink.0
}

fn flag_of(rec: &[u8]) -> u16 {
    u16::from_le_bytes([rec[4 + 14], rec[4 + 15]])
}

/// Single-end only: one record per unique mapper, no pair bits, no mate fields.
#[test]
fn single_end_batch_emits_unpaired_records() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let s = simulate_singles(40, 100, 17);
    let ids = sys::BwaIdBases {
        first_single_id: 0,
        first_pair_id: 0,
    };
    let recs = run_mixed(idx, opts, &[], &s.ffi(), ids);
    assert_eq!(recs.len(), 40);
    for (i, (kind, origin, rec)) in recs.iter().enumerate() {
        assert_eq!(*kind, sys::BWA_ORIGIN_SINGLE);
        assert_eq!(*origin, i);
        let flag = flag_of(rec);
        assert_eq!(
            flag & (0x1 | 0x40 | 0x80 | 0x8 | 0x20),
            0,
            "read {i}: pair bits set: {flag:#x}"
        );
        let mtid = i32::from_le_bytes(rec[4 + 20..4 + 24].try_into().unwrap());
        let mpos = i32::from_le_bytes(rec[4 + 24..4 + 28].try_into().unwrap());
        let tlen = i32::from_le_bytes(rec[4 + 28..4 + 32].try_into().unwrap());
        assert_eq!(
            (mtid, mpos, tlen),
            (-1, -1, 0),
            "read {i}: mate fields must be unset"
        );
    }
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// Mixed batch: pairs then singles, each with its own id base; the single's
/// records are byte-identical to running the same singles alone with the same
/// `first_single_id`, and the pairs' records to running the pairs alone.
#[test]
fn mixed_batch_is_the_union_of_its_groups() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(30, 100, 300, 19);
    let pairs = fx.pairs();
    let s = simulate_singles(25, 100, 21);
    let singles = s.ffi();
    let ids = sys::BwaIdBases {
        first_single_id: 1000,
        first_pair_id: 5000,
    };

    let mixed = run_mixed(idx, opts, &pairs, &singles, ids);
    let only_pairs = run_mixed(idx, opts, &pairs, &[], ids);
    let only_singles = run_mixed(idx, opts, &[], &singles, ids);

    let (mixed_pairs, mixed_singles): (Vec<_>, Vec<_>) = mixed
        .iter()
        .cloned()
        .partition(|(k, _, _)| *k == sys::BWA_ORIGIN_PAIR);
    assert_eq!(mixed_pairs, only_pairs);
    assert_eq!(mixed_singles, only_singles);
    // Order: all pair records precede all single records.
    let first_single = mixed
        .iter()
        .position(|(k, _, _)| *k == sys::BWA_ORIGIN_SINGLE)
        .unwrap();
    assert!(mixed[first_single..]
        .iter()
        .all(|(k, _, _)| *k == sys::BWA_ORIGIN_SINGLE));
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// A read that maps nowhere still yields exactly one (0x4) record.
#[test]
fn unmapped_single_emits_one_unmapped_record() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let mut rng = common::Rng(0xdead_beef);
    let seq: Vec<u8> = (0..100)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    let name = std::ffi::CString::new("junk").unwrap();
    let qual = [b'I'; 100];
    let single = [sys::BwaSingleRead {
        name: name.as_ptr(),
        name_len: 4,
        seq: seq.as_ptr(),
        seq_len: seq.len(),
        qual: qual.as_ptr(),
    }];
    let ids = sys::BwaIdBases {
        first_single_id: 0,
        first_pair_id: 0,
    };
    let recs = run_mixed(idx, opts, &[], &single, ids);
    assert_eq!(recs.len(), 1);
    assert_ne!(flag_of(&recs[0].2) & 0x4, 0);
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// A parallel load must produce the same contig table and the same
/// alignments as the serial load.
#[test]
fn idx_load_threads_matches_serial_load() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let a = common::load_idx(&prefix);
    // SAFETY: `prefix` is a valid NUL-terminated C string alive for the call;
    // returns an owned index handle (null on failure, asserted below).
    let b = unsafe { sys::bwa_shim_idx_load_threads(prefix.as_ptr(), 4) };
    assert!(!b.is_null(), "{}", common::last_error());
    // SAFETY: `a`/`b` are live index handles; the accessors only read them and
    // `i` stays within `n_contigs`.
    unsafe {
        assert_eq!(
            sys::bwa_shim_idx_n_contigs(a),
            sys::bwa_shim_idx_n_contigs(b)
        );
        for i in 0..sys::bwa_shim_idx_n_contigs(a) {
            assert_eq!(
                sys::bwa_shim_idx_contig_len(a, i),
                sys::bwa_shim_idx_contig_len(b, i)
            );
        }
    }
    let opts = common::new_opts();
    let fx = common::simulate(100, 100, 300, 29);
    let pairs = fx.pairs();
    assert_eq!(
        common::align_batch_records(a, opts, &pairs),
        common::align_batch_records(b, opts, &pairs)
    );
    // SAFETY: `opts`/`a`/`b` are the live owned handles; each freed once, not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(a);
        sys::bwa_shim_idx_free(b);
    }
}

/// `n_threads < 1` must be clamped to 1, not passed through raw (which upstream's
/// FMI_search::load_index would treat as "no threads").
#[test]
fn idx_load_threads_clamps_non_positive_n_threads() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let a = common::load_idx(&prefix);
    for n in [0, -1, -100] {
        // SAFETY: `prefix` is a valid NUL-terminated C string alive for the call;
        // `n` is exercised precisely because the shim must clamp it to >= 1.
        let b = unsafe { sys::bwa_shim_idx_load_threads(prefix.as_ptr(), n) };
        assert!(!b.is_null(), "n_threads={n}: {}", common::last_error());
        // SAFETY: `a`/`b` are live index handles; the accessor reads them, then
        // `b` is freed once.
        unsafe {
            assert_eq!(
                sys::bwa_shim_idx_n_contigs(a),
                sys::bwa_shim_idx_n_contigs(b)
            );
            sys::bwa_shim_idx_free(b);
        }
    }
    // SAFETY: `a` is the live owned handle; freed once, not used afterward.
    unsafe {
        sys::bwa_shim_idx_free(a);
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
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let fx = common::simulate(4, 100, 300, 9);
    let pairs = fx.pairs();
    // SAFETY: valid `idx`/`opts`/`sc` handles; `batch_of(&pairs)` borrows the
    // live `pairs` slice (fixture-owned) for the duration of the call.
    let regs = unsafe { sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)) };
    let mut sink = Sink(Vec::new());
    let ids = sys::BwaIdBases {
        first_single_id: 0,
        first_pair_id: 0,
    };
    // SAFETY: valid `idx`/`opts`/`sc`/`regs`; the null `pestat` argument is the
    // point of the test (a paired batch without a model must be rejected, not
    // UB); `sink_fn`'s context is `&mut sink`, live for the call.
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
    // SAFETY: `sc`/`opts`/`idx` are the live owned handles; each freed once, not
    // used afterward.
    unsafe {
        sys::bwa_shim_scratch_free(sc);
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// On AVX2/NEON builds the CLI resolves pairs through the batched kswv mate
/// rescue; the shim used the scalar path. The two must agree byte-for-byte on
/// a fixture where rescue actually fires (one mate mutated hard enough that
/// seeding finds nothing and only mate rescue places it).
#[test]
fn rescue_heavy_fixture_matches_between_batched_and_legacy_paths() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let mut fx = common::simulate(400, 100, 300, 61);
    // Mutate every 3rd R2 at every 7th base: ~14 substitutions in 100 bp
    // leaves no 19-mer seed, so R2 is placeable only by rescue from R1.
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
    let phased = three_phase(idx, opts, &pairs, 64, 0);
    assert_eq!(phased, legacy);
    // Rescue must have fired: some R2 records are mapped (no 0x4) despite the
    // mutation load. Counted on the legacy output so a silent "everything
    // unmapped" regression cannot pass vacuously.
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
    // SAFETY: `opts`/`idx` are the live owned handles; each freed once and not
    // used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
