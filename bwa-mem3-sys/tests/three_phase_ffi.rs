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

/// A tandem-repeat reference (a 300 bp PhiX window repeated five times) and
/// reads drawn from within one copy, so every read maps to several
/// EQUAL-score locations and the tie-breaks decide which copy is primary (and
/// the XA order of the rest). `common::simulate` draws across the repeat-free
/// PhiX genome, where unique reads never tie and an id has no observable
/// effect.
struct RepeatReads {
    r1: Vec<Vec<u8>>,
    r2: Vec<Vec<u8>>,
    singles: Vec<Vec<u8>>,
    qual: Vec<u8>,
}

const REPEAT_READ_LEN: usize = 60;
const REPEAT_INSERT: usize = 120;
const REPEAT_N: usize = 50;

fn repeat_unit() -> &'static [u8] {
    &common::PHIX_SEQ.as_bytes()[500..800]
}

fn repeat_reference() -> Vec<u8> {
    repeat_unit().repeat(5)
}

impl RepeatReads {
    fn new() -> Self {
        let unit = repeat_unit();
        let (len, insert) = (REPEAT_READ_LEN, REPEAT_INSERT);
        let (mut r1, mut r2, mut singles) = (Vec::new(), Vec::new(), Vec::new());
        // Pairs and singles draw from separate generators, so the pairs are
        // the same whether or not a test uses the singles.
        let (mut rng, mut single_rng) = (common::Rng(5), common::Rng(6));
        for _ in 0..REPEAT_N {
            let start = (rng.next() as usize) % (unit.len() - insert + 1);
            r1.push(unit[start..start + len].to_vec());
            r2.push(common::revcomp(&unit[start + insert - len..start + insert]));
            let at = (single_rng.next() as usize) % (unit.len() - len + 1);
            singles.push(unit[at..at + len].to_vec());
        }
        RepeatReads {
            r1,
            r2,
            singles,
            qual: vec![b'I'; len],
        }
    }

    /// The pairs, R1 named from `r1_names` and R2 from `r2_names`.
    fn pairs(
        &self,
        r1_names: &[std::ffi::CString],
        r2_names: &[std::ffi::CString],
    ) -> Vec<sys::BwaReadPair> {
        (0..REPEAT_N)
            .map(|i| sys::BwaReadPair {
                r1_name: r1_names[i].as_ptr(),
                r1_name_len: r1_names[i].as_bytes().len(),
                r1_seq: self.r1[i].as_ptr(),
                r1_seq_len: self.r1[i].len(),
                r1_qual: self.qual.as_ptr(),
                r2_name: r2_names[i].as_ptr(),
                r2_name_len: r2_names[i].as_bytes().len(),
                r2_seq: self.r2[i].as_ptr(),
                r2_seq_len: self.r2[i].len(),
                r2_qual: self.qual.as_ptr(),
            })
            .collect()
    }

    fn singles(&self, names: &[std::ffi::CString]) -> Vec<sys::BwaSingleRead> {
        (0..REPEAT_N)
            .map(|i| sys::BwaSingleRead {
                name: names[i].as_ptr(),
                name_len: names[i].as_bytes().len(),
                seq: self.singles[i].as_ptr(),
                seq_len: self.singles[i].len(),
                qual: self.qual.as_ptr(),
            })
            .collect()
    }
}

/// `{prefix}0`, `{prefix}1`, ... as read names.
fn read_names(prefix: &str) -> Vec<std::ffi::CString> {
    (0..REPEAT_N)
        .map(|i| std::ffi::CString::new(format!("{prefix}{i}")).unwrap())
        .collect()
}

/// A packed BAM record (`block_size` included) without its QNAME, so records
/// of reads that differ only in name compare equal.
fn without_qname(rec: &[u8]) -> Vec<u8> {
    let l_read_name = rec[4 + 8] as usize;
    let mut out = rec[4..4 + 32].to_vec();
    out[8] = 0;
    out.extend_from_slice(&rec[4 + 32 + l_read_name..]);
    out
}

/// The pair id feeds hash_64 tie-breaks (mem_mark_primary_se / mem_pair,
/// seeded with `id<<1|i` at bwamem_pair.cpp:893-910). Offsetting `first_pair_id`
/// must therefore change *something* on a fixture with equal-score ties (proves
/// the id is wired), and legacy `align_batch` must equal `first_pair_id = 0`
/// (its own convention).
#[test]
fn first_pair_id_reaches_the_tie_break_hash() {
    let Some((_dir, prefix)) = common::index_custom_ref(&repeat_reference()) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let reads = RepeatReads::new();
    let names = read_names("r");
    let pairs = reads.pairs(&names, &names);

    let at_zero = three_phase(idx, opts, &pairs, REPEAT_N, 0);
    assert_eq!(at_zero, common::align_batch_records(idx, opts, &pairs));
    let shifted = three_phase(idx, opts, &pairs, REPEAT_N, 1 << 20);
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

/// `--tie-break name` (v0.18.0, `mem_opt_t::tie_break_mode`): upstream's
/// worker passes every tie-break seed through `mem_tie_break_id`, which under
/// `MEM_TIE_BREAK_NAME` (1) hashes the read name (mate 1's, for a pair)
/// instead of using the input ordinal. The shim derives its ids from
/// `BwaIdBases` itself, so it must apply the same helper on the pair and the
/// single path, or a raw `mem_opt_t` caller selecting name mode silently keeps
/// index tie-breaks.
///
/// The emit runs under a fixed FR insert-size model: the repeat's own estimate
/// fails every orientation (no read is unique), which would leave `mem_pair`,
/// and its id-seeded choice among equal-score pairings, unexercised.
#[test]
fn tie_break_name_hashes_mate_one_name_not_the_id() {
    let Some((_dir, prefix)) = common::index_custom_ref(&repeat_reference()) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let reads = RepeatReads::new();
    let (r, x, q) = (read_names("r"), read_names("x"), read_names("q"));
    let (s, t) = (read_names("s"), read_names("t"));

    let failed = sys::mem_pestat_t {
        low: 0,
        high: 0,
        failed: 1,
        avg: 0.0,
        std: 0.0,
    };
    let fr = sys::mem_pestat_t {
        low: (REPEAT_INSERT - 40) as i32,
        high: (REPEAT_INSERT + 40) as i32,
        failed: 0,
        avg: REPEAT_INSERT as f64,
        std: 10.0,
    };
    let pestat = [failed, fr, failed, failed]; // FF, FR, RF, RR
    type Records = Vec<(u32, usize, Vec<u8>)>;
    // The closure reads the options through its own copy of the pointer, so
    // the tie-break mode can still be switched through `opts` between runs.
    let opts_ro = opts.cast_const();
    let run = |r1: &[std::ffi::CString],
               r2: &[std::ffi::CString],
               sn: &[std::ffi::CString],
               base: u64|
     -> (Records, Records) {
        let ids = sys::BwaIdBases {
            first_single_id: base,
            first_pair_id: base,
        };
        let recs = run_mixed_with(
            idx,
            opts_ro,
            &reads.pairs(r1, r2),
            &reads.singles(sn),
            ids,
            Some(&pestat),
        );
        recs.into_iter()
            .partition(|(k, _, _)| *k == sys::BWA_ORIGIN_PAIR)
    };
    let sans = |recs: &Records| -> Vec<Vec<u8>> {
        recs.iter().map(|(_, _, rec)| without_qname(rec)).collect()
    };

    // Index mode: the fixture ties on both paths, and the pairs pair.
    let (pairs0, singles0) = run(&r, &r, &s, 0);
    let (pairs1, singles1) = run(&r, &r, &s, 1 << 20);
    assert_ne!(pairs0, pairs1, "the pairs never tie");
    assert_ne!(singles0, singles1, "the singles never tie");
    assert!(
        pairs0.iter().any(|(_, _, rec)| flag_of(rec) & 0x2 != 0),
        "no pair was paired, so mem_pair's tie-break never ran"
    );

    // SAFETY: `opts` is a valid, uniquely-owned `mem_opt_t`; POD field write.
    unsafe { (*opts).tie_break_mode = 1 }; // MEM_TIE_BREAK_NAME
    let (name_pairs, name_singles) = run(&r, &r, &s, 0);
    let (shifted_pairs, shifted_singles) = run(&r, &r, &s, 1 << 20);
    assert_eq!(
        name_pairs, shifted_pairs,
        "pair ids still reach the tie-break hash under --tie-break name"
    );
    assert_eq!(
        name_singles, shifted_singles,
        "single ids still reach the tie-break hash under --tie-break name"
    );
    // The name, and for a pair mate 1's name only, is what is hashed.
    let (r1_renamed, _) = run(&q, &r, &s, 0);
    assert_ne!(
        sans(&name_pairs),
        sans(&r1_renamed),
        "renaming R1 changed no tie-break: the pair seed ignores mate 1's name"
    );
    let (r2_renamed, _) = run(&r, &x, &s, 0);
    assert_eq!(
        sans(&name_pairs),
        sans(&r2_renamed),
        "renaming only R2 changed a tie-break: the pair seed reads mate 2's name"
    );
    let (_, singles_renamed) = run(&r, &r, &t, 0);
    assert_ne!(
        sans(&name_singles),
        sans(&singles_renamed),
        "renaming the singles changed no tie-break: the single seed ignores the name"
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
    run_mixed_with(idx, opts, pairs, singles, ids, None)
}

/// `run_mixed`, emitting under `fixed` (FF, FR, RF, RR) in place of the
/// batch's own insert-size estimate when given, as a CLI `-I` run would.
fn run_mixed_with(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
    singles: &[sys::BwaSingleRead],
    ids: sys::BwaIdBases,
    fixed: Option<&[sys::mem_pestat_t; 4]>,
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
    if let Some(fixed) = fixed {
        // SAFETY: `pes` is the live `mem_pestat_t[4]` handle; `fixed` holds 4.
        unsafe { std::ptr::copy_nonoverlapping(fixed.as_ptr(), pes, 4) };
    }
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

/// Set by `no_fields` if the shim ever calls it.
static FIELDS_EMITTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// A field sink the refusal test passes but the shim must never call. It
/// records the call rather than panicking: unwinding out of an `extern "C"`
/// function is undefined behaviour on the MSRV.
///
/// # Safety
/// Never dereferences its arguments.
unsafe extern "C" fn no_fields(
    _ctx: *mut std::os::raw::c_void,
    _origin_kind: u32,
    _origin_idx: usize,
    _mate: u8,
    _is_primary: std::os::raw::c_int,
    _f: *const sys::BwaAlignedFields,
) {
    FIELDS_EMITTED.store(true, std::sync::atomic::Ordering::SeqCst);
}

/// `--alt-project` (v0.18.0) is not supported: the shim's record writer emits
/// neither the lift-source `XL:Z` tag nor the demoted read's `OA:Z`, and under
/// `alt_project_primary` `mem_kernel2_core` asserts a side array the shim does
/// not allocate. `MemOpts` cannot set any of these knobs, but a raw
/// `mem_opt_t` caller can, so every C entry point refuses each one rather than
/// emit records that silently differ from the CLI's (or abort).
#[test]
fn c_entry_points_refuse_alt_project() {
    // A non-NULL table pointer is all the gate reads; the shim refuses before
    // anything could dereference it.
    static TABLE: u8 = 1;
    type Set = fn(*mut sys::mem_opt_t);
    let knobs: [(&str, Set); 4] = [
        ("alt_lift", |o| {
            // SAFETY: a valid, uniquely-owned `mem_opt_t`; one POD field write.
            unsafe { (*o).alt_lift = std::ptr::addr_of!(TABLE).cast() }
        }),
        ("alt_project_min_gap", |o| {
            // SAFETY: as above.
            unsafe { (*o).alt_project_min_gap = 1 }
        }),
        ("alt_project_pair_gap", |o| {
            // SAFETY: as above.
            unsafe { (*o).alt_project_pair_gap = 1 }
        }),
        ("alt_project_primary", |o| {
            // SAFETY: as above.
            unsafe { (*o).alt_project_primary = 1 }
        }),
    ];
    for (field, set) in knobs {
        assert_every_entry_point_refuses(set, |what| {
            format!("{what}: --alt-project is not supported (opts->{field} is set;")
        });
    }
}

/// `--tie-break` (v0.18.0) takes `index` (0) or `name` (1); the CLI's parser
/// accepts nothing else. The shim reads `tie_break_mode` at its three id sites
/// (through `mem_tie_break_id`, which treats any value but 1 as index), so a
/// raw `mem_opt_t` caller's out-of-range value is refused rather than silently
/// meaning index.
#[test]
fn c_entry_points_refuse_an_unknown_tie_break_mode() {
    assert_every_entry_point_refuses(
        |o| {
            // SAFETY: a valid, uniquely-owned `mem_opt_t`; one POD field write.
            unsafe { (*o).tie_break_mode = 2 }
        },
        |what| format!("{what}: opts->tie_break_mode is 2, not 0 (index) or 1 (name)"),
    );
}

/// A raw `mem_opt_t` caller skips `MemOpts::validate_for_alignment`, so every
/// C alignment entry point refuses an unstorable mismatch penalty itself:
/// `bwa_shim_opts_fill_scmat` leaves the matrices at the last storable `b`,
/// and aligning would score with them while MAPQ and primary marking read the
/// new `b`. Handles each call consumes are released on the refusal too.
#[test]
fn c_entry_points_refuse_an_unstorable_mismatch_penalty() {
    assert_every_entry_point_refuses(
        |o| {
            // A COLLAPSED-style `b = 2 * a` for `a = 100`: past what the
            // matrices hold.
            // SAFETY: `o` is a valid, uniquely-owned `mem_opt_t`; POD field
            // write, then a refill the shim skips for an unstorable `b`.
            unsafe {
                (*o).b = 200;
                sys::bwa_shim_opts_fill_scmat(o);
            }
        },
        |what| format!("{what}: the effective mismatch penalty is 200, outside 0..128"),
    );
}

/// Builds every input under clean options, applies `mutate` to them, then
/// calls each of the 13 C alignment entry points and asserts that each refuses
/// with an error starting `expect(<entry point>)`, and that nothing reached a
/// record or field sink. Every entry point gates on `shim_opts_alignable`, so
/// one refusal cause failing here at one call names the call that lost it.
fn assert_every_entry_point_refuses(
    mutate: impl FnOnce(*mut sys::mem_opt_t),
    expect: impl Fn(&str) -> String,
) {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let pes = common::new_pestat();
    // SAFETY: constructor with no arguments; returns an owned scratch handle.
    let sc = unsafe { sys::bwa_shim_scratch_new() };
    let fx = common::simulate(4, 100, 300, 11);
    let pairs = fx.pairs();
    // Inputs built under storable options, so only the penalty is wrong below.
    // SAFETY: valid `idx`/`opts`/`sc` handles; `pairs` (fixture-owned) outlives
    // both calls and `pairs.len()` matches it.
    let (seeds, regs) = unsafe {
        (
            sys::bwa_shim_seed_batch(idx, opts, pairs.as_ptr(), pairs.len()),
            sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)),
        )
    };
    assert!(
        !seeds.is_null() && !regs.is_null(),
        "{}",
        common::last_error()
    );
    // SAFETY: constructor taking a plain int; returns an owned cohort handle.
    let cohort = unsafe { sys::bwa_shim_resident_cohort_new(0) };
    let mut first = 0usize;
    // SAFETY: valid cohort handle; `first` is a live out-parameter.
    let seg = unsafe { sys::bwa_shim_resident_reserve_pairs(cohort, 2, &mut first) };
    assert!(
        !cohort.is_null() && !seg.is_null(),
        "{}",
        common::last_error()
    );

    mutate(opts);
    let refused = |what: &str| {
        let err = common::last_error();
        assert!(err.starts_with(&expect(what)), "{what}: {err}");
    };
    let ids = sys::BwaIdBases {
        first_single_id: 0,
        first_pair_id: 0,
    };
    let mut sink = Sink(Vec::new());
    let mut memo = sys::BwaMemoStats {
        dup_pairs: 0,
        copied: 0,
        fallback_aligned: 0,
    };
    // SAFETY: every handle is live and valid; `pairs` outlives the calls; the
    // sink context is `&mut sink`, live for each call; `extend_batch` and
    // `pair_emit` consume `seeds` / `regs` on every path, so neither is used or
    // freed afterward.
    unsafe {
        let n = pairs.len();
        assert!(sys::bwa_shim_seed_batch(idx, opts, pairs.as_ptr(), n).is_null());
        refused("seed_batch");
        let b = sys::bwa_shim_extend_batch(idx, opts, seeds, pairs.as_ptr(), n, pes, pes);
        assert!(b.is_null());
        refused("extend_batch");
        assert!(sys::bwa_shim_align_batch(idx, opts, pairs.as_ptr(), n, pes, pes).is_null());
        refused("align_batch");
        assert_eq!(
            sys::bwa_shim_estimate_pestat(idx, opts, pairs.as_ptr(), n, pes),
            -1
        );
        refused("estimate_pestat");
        assert!(sys::bwa_shim_seed_extend(idx, opts, sc, &batch_of(&pairs)).is_null());
        refused("seed_extend");
        let rs: [*const sys::BwaRegs; 1] = [regs];
        assert_eq!(
            sys::bwa_shim_pestat_cohort(idx, opts, rs.as_ptr(), 1, pes),
            -1
        );
        refused("pestat_cohort");
        let ctx = (&mut sink as *mut Sink).cast();
        let rc = sys::bwa_shim_pair_emit(idx, opts, sc, regs, pes, ids, Some(sink_fn), ctx);
        assert_eq!(rc, -1);
        refused("pair_emit");
        assert_eq!(sys::bwa_shim_resident_seed_extend(idx, opts, sc, seg), -1);
        refused("resident_seed_extend");
        let reps = [-1i64];
        let rc = sys::bwa_shim_resident_seed_extend_reps(idx, opts, sc, seg, reps.as_ptr(), 1);
        assert_eq!(rc, -1);
        refused("resident_seed_extend_reps");
        assert_eq!(
            sys::bwa_shim_resident_resolve_memo(idx, opts, sc, cohort, &mut memo),
            -1
        );
        refused("resident_resolve_memo");
        assert_eq!(
            sys::bwa_shim_resident_pestat_cohort(idx, opts, cohort, pes),
            -1
        );
        refused("resident_pestat_cohort");
        let rc =
            sys::bwa_shim_resident_pair_emit(idx, opts, sc, seg, pes, ids, 0, Some(sink_fn), ctx);
        assert_eq!(rc, -1);
        refused("resident_pair_emit");
        let rc = sys::bwa_shim_resident_pair_emit_fields(
            idx,
            opts,
            sc,
            seg,
            pes,
            ids,
            0,
            Some(no_fields),
            std::ptr::null_mut(),
        );
        assert_eq!(rc, -1);
        refused("resident_pair_emit_fields");
    }
    assert!(sink.0.is_empty());
    assert!(!FIELDS_EMITTED.load(std::sync::atomic::Ordering::SeqCst));
    // SAFETY: the remaining live owned handles, each freed once and unused
    // afterward (`seg` is owned by `cohort`).
    unsafe {
        sys::bwa_shim_resident_cohort_free(cohort);
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_scratch_free(sc);
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
