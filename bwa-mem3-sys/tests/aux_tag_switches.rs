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
use std::ffi::CString;

/// One aux tag: its key, its type byte, and its raw value bytes (for `Z`/`H`,
/// without the trailing NUL).
struct AuxTag<'a> {
    key: [u8; 2],
    ty: u8,
    value: &'a [u8],
}

/// The aux tags of one `[u32 block_size][body]` record, in order. Panics on a
/// tag whose value runs past the end of the record, rather than reporting a
/// truncated tag as present.
fn aux_tags(rec: &[u8]) -> Vec<AuxTag<'_>> {
    let body = &rec[4..];
    assert!(
        body.len() >= 32,
        "record body shorter than its fixed fields"
    );
    let u16_at = |o: usize| u16::from_le_bytes([body[o], body[o + 1]]) as usize;
    let u32_at = |o: usize| u32::from_le_bytes(body[o..o + 4].try_into().unwrap()) as usize;
    let l_read_name = body[8] as usize;
    let n_cigar = u16_at(12);
    let l_seq = u32_at(16);
    // Where the aux tags start. A body shorter than its declared fields is
    // truncated, not a record without tags (`o == len` is a complete record
    // with none).
    let mut o = [l_read_name, 4 * n_cigar, l_seq.div_ceil(2), l_seq]
        .into_iter()
        .try_fold(32usize, usize::checked_add)
        .expect("record field lengths overflow usize");
    assert!(
        o <= body.len(),
        "record truncated: its fields need {o} bytes but the body has {}",
        body.len()
    );
    let mut tags = Vec::new();
    while o < body.len() {
        assert!(o + 3 <= body.len(), "aux tag header truncated at byte {o}");
        let (key, ty) = ([body[o], body[o + 1]], body[o + 2]);
        o += 3;
        let rest = &body[o..];
        let (value_len, skip) = match ty {
            b'A' | b'c' | b'C' => (1, 1),
            b's' | b'S' => (2, 2),
            b'i' | b'I' | b'f' => (4, 4),
            b'Z' | b'H' => {
                let nul = rest
                    .iter()
                    .position(|&b| b == 0)
                    .unwrap_or_else(|| panic!("{key:?}: string value has no NUL"));
                (nul, nul + 1)
            }
            b'B' => {
                assert!(rest.len() >= 5, "{key:?}: B-array header truncated");
                let width = match rest[0] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                let n = 5 + width * u32_at(o + 1);
                (n, n)
            }
            other => panic!("unknown aux type {:?}", other as char),
        };
        assert!(
            skip <= rest.len(),
            "{key:?}: value needs {skip} bytes but {} remain",
            rest.len()
        );
        tags.push(AuxTag {
            key,
            ty,
            value: &rest[..value_len],
        });
        o += skip;
    }
    tags
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
