//! Aux-tag output switches on `mem_opt_t`, exercised through the raw FFI:
//! bwa-mem3 v0.17.0's `md_mode` (`--md`; the safe API exposes it too, as
//! `MemOpts::set_md_mode`) and the `MEM_F_XB` flag (`-u`), whose hit list
//! upstream names `XB:Z` as of v0.17.0 (fg-labs/bwa-mem3#690) and which only a
//! raw `mem_opt_t` caller can set. Each test builds an index (skips when
//! `bwa-mem3` is absent) and asserts on the aux-tag keys of the packed BAM
//! records.
//!
//! ## Safety invariants shared by the `unsafe` blocks below
//!
//! `opts` is a valid, uniquely-owned `mem_opt_t` from `common::new_opts`, alive
//! until it is freed exactly once at the end of its test; writing its POD
//! fields through the raw pointer is sound.

mod common;

use bwa_mem3_sys as sys;
use bwa_mem3_test_util::AuxTag;
use std::ffi::CString;

/// The aux tags of one `[u32 block_size][body]` record, in order (the shared
/// decoder, which panics on a truncated tag rather than reporting it present).
fn aux_tags(rec: &[u8]) -> Vec<AuxTag<'_>> {
    bwa_mem3_test_util::aux_tags(&rec[4..])
}

fn tag<'a>(rec: &'a [u8], key: &[u8; 2]) -> Option<AuxTag<'a>> {
    aux_tags(rec).into_iter().find(|t| &t.key == key)
}

/// A `Z` tag's string, or `None` if the record has no such tag of type `Z`.
fn z_tag(rec: &[u8], key: &[u8; 2]) -> Option<String> {
    tag(rec, key)
        .filter(|t| t.ty == b'Z')
        .map(|t| String::from_utf8(t.value.to_vec()).unwrap())
}

/// An integer tag's value, whatever width/signedness it was written with.
fn int_tag(rec: &[u8], key: &[u8; 2]) -> Option<i64> {
    let t = tag(rec, key)?;
    let v = t.value;
    Some(match t.ty {
        b'c' => i64::from(v[0] as i8),
        b'C' => i64::from(v[0]),
        b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
        b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
        b'i' => i64::from(i32::from_le_bytes(v.try_into().unwrap())),
        b'I' => i64::from(u32::from_le_bytes(v.try_into().unwrap())),
        _ => return None,
    })
}

fn n_cigar(rec: &[u8]) -> usize {
    u16::from_le_bytes([rec[4 + 12], rec[4 + 13]]) as usize
}

fn has_key(rec: &[u8], key: &[u8; 2]) -> bool {
    tag(rec, key).is_some()
}

/// A record's identity across two runs that differ only in aux tags: pair
/// index, `flag` (which mate, and primary/secondary/supplementary), `refID`,
/// `pos` and `n_cigar` (so a record that turned unmapped no longer matches).
type RecordId = (usize, u16, i32, i32, usize);

fn record_id(pair_idx: usize, rec: &[u8]) -> RecordId {
    let body = &rec[4..];
    let i32_at = |o: usize| i32::from_le_bytes(body[o..o + 4].try_into().unwrap());
    let flag = u16::from_le_bytes([body[14], body[15]]);
    (pair_idx, flag, i32_at(0), i32_at(4), n_cigar(rec))
}

/// The sorted identities of the records in `recs` that satisfy `keep`, as a
/// multiset: a dropped or duplicated record changes it.
fn ids_where(recs: &[(usize, Vec<u8>)], keep: impl Fn(&[u8]) -> bool) -> Vec<RecordId> {
    let mut ids: Vec<RecordId> = recs
        .iter()
        .filter(|(_, r)| keep(r))
        .map(|(p, r)| record_id(*p, r))
        .collect();
    ids.sort_unstable();
    ids
}

#[test]
fn md_mode_none_drops_md_and_full_keeps_it() {
    let Some((_dir, prefix)) = common::phix_index() else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();
    let fx = common::simulate(64, 100, 300, 11);

    // SAFETY: see the module docs; reads one POD field.
    assert_eq!(
        unsafe { (*opts).md_mode },
        1,
        "mem_opt_init defaults md_mode to MEM_MD_FULL"
    );
    let full = common::align_batch_records(idx, opts, &fx.pairs());
    let mapped: Vec<_> = full.iter().filter(|(_, r)| n_cigar(r) > 0).collect();
    assert!(!mapped.is_empty(), "fixture produced no mapped records");
    assert!(
        mapped.iter().all(|(_, r)| z_tag(r, b"MD").is_some()),
        "MEM_MD_FULL: every mapped record carries MD:Z"
    );

    // SAFETY: see the module docs; writes one POD field (MEM_MD_NONE).
    unsafe { (*opts).md_mode = 0 };
    let none = common::align_batch_records(idx, opts, &fx.pairs());
    assert_eq!(
        ids_where(&none, |_| true),
        ids_where(&full, |_| true),
        "MEM_MD_NONE changes no record's identity or placement"
    );
    assert!(
        none.iter().all(|(_, r)| !has_key(r, b"MD")),
        "MEM_MD_NONE: no record carries MD:Z"
    );
    // Each mapped record's NM, keyed by identity: dropping MD must leave the
    // edit distance itself unchanged, not merely the tag present.
    let mapped_nm = |recs: &[(usize, Vec<u8>)]| {
        let mut v: Vec<(RecordId, Option<i64>)> = recs
            .iter()
            .filter(|(_, r)| n_cigar(r) > 0)
            .map(|(p, r)| (record_id(*p, r), int_tag(r, b"NM")))
            .collect();
        v.sort_unstable();
        v
    };
    let full_nm = mapped_nm(&full);
    assert!(
        full_nm.iter().all(|(_, nm)| nm.is_some()),
        "MEM_MD_FULL: every mapped record carries NM"
    );
    assert_eq!(
        mapped_nm(&none),
        full_nm,
        "MEM_MD_NONE drops MD only: every record mapped under FULL is still mapped, with the same NM"
    );

    // SAFETY: each handle freed exactly once, never used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// A reference with one 1 kb segment present twice, and FR pairs drawn from
/// inside it, so both mates map to two places and carry an alternative-hit list.
fn duplicated_segment_fixture() -> (Vec<u8>, common::Fixture) {
    let mut rng = common::Rng(0x0dd5_eed5);
    let mut random =
        |n: usize| -> Vec<u8> { (0..n).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect() };
    let segment = random(1000);
    let mut reference = random(3000);
    reference.extend_from_slice(&segment);
    reference.extend(random(3000));
    reference.extend_from_slice(&segment);
    reference.extend(random(3000));

    let (read_len, insert) = (100, 300);
    let mut fx = common::Fixture {
        names: Vec::new(),
        r1: Vec::new(),
        r2: Vec::new(),
        qual: Vec::new(),
    };
    for i in 0..16 {
        let start = i * 40;
        fx.names.push(CString::new(format!("d{i}")).unwrap());
        fx.r1.push(segment[start..start + read_len].to_vec());
        fx.r2.push(common::revcomp(
            &segment[start + insert - read_len..start + insert],
        ));
        fx.qual.push(vec![b'I'; read_len]);
    }
    (reference, fx)
}

#[test]
fn mem_f_xb_names_the_hit_list_xb() {
    let (reference, fx) = duplicated_segment_fixture();
    let Some((_dir, prefix)) = common::index_custom_ref(&reference) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();

    let xa = common::align_batch_records(idx, opts, &fx.pairs());
    let with_xa = ids_where(&xa, |r| has_key(r, b"XA"));
    assert!(
        !with_xa.is_empty(),
        "fixture produced no XA:Z; it no longer tests -u"
    );
    assert!(xa.iter().all(|(_, r)| !has_key(r, b"XB")));

    // SAFETY: see the module docs; sets one bit of a POD field.
    unsafe { (*opts).flag |= sys::MEM_F_XB as i32 };
    let xb = common::align_batch_records(idx, opts, &fx.pairs());
    assert!(
        xb.iter().all(|(_, r)| !has_key(r, b"XA")),
        "MEM_F_XB: the hit list is XB:Z, never XA:Z"
    );
    assert_eq!(
        ids_where(&xb, |r| has_key(r, b"XB")),
        with_xa,
        "MEM_F_XB renames the tag on exactly the records that carried XA:Z"
    );

    // Each XB:Z entry is the matching XA:Z entry (`chr,pos,CIGAR,NM`) plus the
    // hit's alignment score and MAPQ (bwamem_extra.cpp, mem_gen_alt).
    let hit_lists = |recs: &[(usize, Vec<u8>)], key: &[u8; 2]| {
        let mut v: Vec<(RecordId, String)> = recs
            .iter()
            .filter_map(|(p, r)| Some((record_id(*p, r), z_tag(r, key)?)))
            .collect();
        v.sort_unstable();
        v
    };
    let (xa_lists, xb_lists) = (hit_lists(&xa, b"XA"), hit_lists(&xb, b"XB"));
    // Every tagged record's list must be type Z (an `XA:i` / `XB:i` would
    // drop out of `hit_lists` and leave the entry checks below with nothing).
    let ids = |lists: &[(RecordId, String)]| lists.iter().map(|(id, _)| *id).collect::<Vec<_>>();
    assert_eq!(ids(&xa_lists), with_xa, "every XA tag is XA:Z");
    assert_eq!(ids(&xb_lists), with_xa, "every XB tag is XB:Z");
    for ((id, xa_s), (_, xb_s)) in xa_lists.iter().zip(&xb_lists) {
        let xa_hits: Vec<&str> = xa_s.split_terminator(';').collect();
        let xb_hits: Vec<&str> = xb_s.split_terminator(';').collect();
        assert!(!xa_hits.is_empty(), "{id:?}: empty XA:Z list");
        assert_eq!(xa_hits.len(), xb_hits.len(), "{id:?}: {xa_s} vs {xb_s}");
        for (a, b) in xa_hits.iter().zip(&xb_hits) {
            let (a, b): (Vec<&str>, Vec<&str>) = (a.split(',').collect(), b.split(',').collect());
            assert_eq!(a.len(), 4, "{id:?}: XA entry {a:?}");
            assert_eq!(b.len(), 6, "{id:?}: XB entry {b:?} lacks AS and MAPQ");
            assert_eq!(
                b[..4],
                a[..],
                "{id:?}: XB entry does not extend its XA entry"
            );
            for field in &b[4..] {
                field
                    .parse::<i32>()
                    .unwrap_or_else(|_| panic!("{id:?}: XB field {field:?} is not an integer"));
            }
        }
    }

    // SAFETY: each handle freed exactly once, never used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// One record emitted by the resident structured-fields sink: who it is
/// (origin index, mate, primary-ness, FLAG, POS) and its alternative-hit list,
/// if any, as `(xa_is_xb, number of hits)`.
type HitList = ((usize, u8, bool, u16, i32), Option<(bool, usize)>);

#[derive(Default)]
struct HitListTags(Vec<HitList>);

/// # Safety
/// `ctx` must be a valid `*mut HitListTags` that outlives the call and `f` a
/// valid `BwaAlignedFields` for its duration -- both guaranteed by the shim
/// when it invokes this callback from `bwa_shim_resident_pair_emit_fields`.
unsafe extern "C" fn hit_list_sink(
    ctx: *mut std::ffi::c_void,
    _kind: u32,
    idx: usize,
    mate: u8,
    is_primary: std::ffi::c_int,
    f: *const sys::BwaAlignedFields,
) {
    // SAFETY: see the function's safety contract.
    let (sink, f) = unsafe { (&mut *ctx.cast::<HitListTags>(), &*f) };
    let list = (!f.xa.is_null()).then(|| {
        // SAFETY: a non-NULL `xa` is a NUL-terminated string valid for the call.
        let l = unsafe { std::ffi::CStr::from_ptr(f.xa) }.to_bytes();
        (
            f.xa_is_xb != 0,
            l.split(|&b| b == b';').filter(|e| !e.is_empty()).count(),
        )
    });
    sink.0
        .push(((idx, mate, is_primary != 0, f.flag, f.pos), list));
}

/// Run `pairs` as one resident cohort through the structured-fields sink and
/// return each record's [`HitListTags`] entry, in emission order.
fn resident_hit_list_tags(
    idx: *const sys::BwaIndex,
    opts: *const sys::mem_opt_t,
    pairs: &[sys::BwaReadPair],
) -> Vec<HitList> {
    // SAFETY: `idx`/`opts` are live handles; `sc`/`cohort`/`pes` are created,
    // used by this single thread through the one segment they own, and freed
    // exactly once below; each `BwaReadPair` borrows fixture bytes that outlive
    // the calls.
    unsafe {
        let sc = sys::bwa_shim_scratch_new();
        let cohort = sys::bwa_shim_resident_cohort_new(0);
        assert!(
            !sc.is_null() && !cohort.is_null(),
            "{}",
            common::last_error()
        );
        let mut first = 0usize;
        let seg = sys::bwa_shim_resident_reserve_pairs(cohort, 2 * pairs.len(), &mut first);
        assert!(!seg.is_null(), "reserve failed: {}", common::last_error());
        let mut added = 0usize;
        for (i, pair) in pairs.iter().enumerate() {
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
        let pes = common::new_pestat();
        assert_eq!(
            sys::bwa_shim_resident_pestat_cohort(idx, opts, cohort, pes),
            0,
            "{}",
            common::last_error()
        );
        let mut sink = HitListTags::default();
        let ids = sys::BwaIdBases {
            first_single_id: 0,
            first_pair_id: 0,
        };
        let rc = sys::bwa_shim_resident_pair_emit_fields(
            idx,
            opts,
            sc,
            seg,
            pes,
            ids,
            0,
            Some(hit_list_sink),
            std::ptr::addr_of_mut!(sink).cast(),
        );
        assert_eq!(rc, 0, "{}", common::last_error());
        sys::bwa_shim_pestat_free(pes);
        sys::bwa_shim_resident_cohort_free(cohort);
        sys::bwa_shim_scratch_free(sc);
        sink.0
    }
}

/// The structured-fields sink must say which tag its hit list is, or a caller
/// rebuilding the record writes the six-field XB list under `XA:Z`.
#[test]
fn mem_f_xb_is_reported_to_the_fields_sink() {
    let (reference, fx) = duplicated_segment_fixture();
    let Some((_dir, prefix)) = common::index_custom_ref(&reference) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();

    let xa = resident_hit_list_tags(idx, opts, &fx.pairs());
    // SAFETY: see the module docs; sets one bit of a POD field.
    unsafe { (*opts).flag |= sys::MEM_F_XB as i32 };
    let xb = resident_hit_list_tags(idx, opts, &fx.pairs());

    // -u changes only the tag: the same records, in the same order, and the
    // same ones carry a hit list, with the same number of hits.
    let ids = |r: &[HitList]| r.iter().map(|(id, _)| *id).collect::<Vec<_>>();
    assert_eq!(ids(&xb), ids(&xa), "-u changed which records are emitted");
    let with_list = |r: &[HitList]| {
        r.iter()
            .filter_map(|(id, l)| l.map(|(_, n)| (*id, n)))
            .collect::<Vec<_>>()
    };
    let (xa_lists, xb_lists) = (with_list(&xa), with_list(&xb));
    assert!(!xa_lists.is_empty(), "fixture produced no hit list");
    assert!(
        xa_lists.iter().all(|(_, n)| *n > 0),
        "an empty XA:Z list: {xa_lists:?}"
    );
    assert_eq!(
        xb_lists, xa_lists,
        "-u moved, dropped or resized a hit list"
    );
    // `map_or(true, ..)`, not `is_none_or`: the latter is stable since 1.82
    // and this workspace's MSRV is 1.75.
    #[allow(clippy::unnecessary_map_or)]
    let all_lists = |lists: &[HitList], want_xb: bool| {
        lists
            .iter()
            .all(|(_, l)| l.map_or(true, |(xb, _)| xb == want_xb))
    };
    assert!(all_lists(&xa, false), "XA:Z reported as XB without -u");
    assert!(all_lists(&xb, true), "MEM_F_XB: every hit list is XB");

    // SAFETY: each handle freed exactly once, never used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}

/// A NULL `compat` (a raw caller that memsets `mem_opt_t`) reads as the `off`
/// row everywhere, as upstream's header builder takes it: the records keep the
/// `MQ:i` / `HN:i` the default target emits, rather than the header saying one
/// target and the records another.
#[test]
fn null_compat_reads_as_the_off_row() {
    let (reference, fx) = duplicated_segment_fixture();
    let Some((_dir, prefix)) = common::index_custom_ref(&reference) else {
        return;
    };
    let idx = common::load_idx(&prefix);
    let opts = common::new_opts();

    let off = common::align_batch_records(idx, opts, &fx.pairs());
    assert!(
        off.iter().any(|(_, r)| has_key(r, b"MQ")) && off.iter().any(|(_, r)| has_key(r, b"HN")),
        "fixture emits no MQ:i or no HN:i; it no longer tests the compat switches"
    );
    // SAFETY: see the module docs; writes one POD pointer field.
    unsafe { (*opts).compat = std::ptr::null() };
    let null = common::align_batch_records(idx, opts, &fx.pairs());
    assert!(null == off, "a NULL compat changed the records");

    // SAFETY: each handle freed exactly once, never used afterward.
    unsafe {
        sys::bwa_shim_opts_free(opts);
        sys::bwa_shim_idx_free(idx);
    }
}
