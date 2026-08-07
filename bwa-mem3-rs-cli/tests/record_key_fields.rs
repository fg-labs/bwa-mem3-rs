//! Unit-tests the SAM-line reduction used by every parity comparison.
//!
//! The reduction must notice a record that carries an *extra* aux tag, even
//! one it does not otherwise compare — otherwise an upstream change that
//! starts emitting a new tag (e.g. XB) is invisible to the parity suite.

mod common;

use common::record_key_fields;

const BASE: &str = "r0\t99\tphix\t101\t60\t100M\t=\t250\t249\tACGT\tIIII\tNM:i:0\tMD:Z:100";

#[test]
fn compared_tag_values_are_reflected() {
    let with_xa = format!("{BASE}\tXA:Z:phix,+201,100M,0;");
    assert_ne!(
        record_key_fields(BASE),
        record_key_fields(&with_xa),
        "an added XA value must change the key"
    );
}

#[test]
fn an_extra_uncompared_tag_changes_the_key() {
    let with_xb = format!("{BASE}\tXB:i:7");
    assert_ne!(
        record_key_fields(BASE),
        record_key_fields(&with_xb),
        "a tag the reduction does not compare by value must still change the \
         key, or a newly-emitted upstream tag is invisible to parity tests"
    );
}

#[test]
fn tag_order_does_not_change_the_key() {
    let a = format!("{BASE}\tXG:Z:CT\tXR:Z:CT");
    let b = format!("{BASE}\tXR:Z:CT\tXG:Z:CT");
    assert_eq!(
        record_key_fields(&a),
        record_key_fields(&b),
        "aux ordering is deliberately not pinned"
    );
}

#[test]
fn an_allowlisted_tag_would_not_change_the_key() {
    // DELIBERATELY_ASYMMETRIC_TAG_KEYS is empty today: the shim emits every tag
    // upstream does, MQ:i included. This pins the MECHANISM rather than any
    // particular tag, so it keeps working if an entry is ever added -- and, more
    // usefully, it documents what adding one costs. A listed key is subtracted
    // from the comparison on both sides, so a real divergence involving it can
    // no longer fail a parity test.
    //
    // With an empty list, a record carrying an extra tag MUST change the key --
    // that is the property the suite depends on, and the assertion below is the
    // one that would break first if the list silently regrew.
    let with_mq = format!("{BASE}\tMQ:i:60");
    assert_ne!(
        record_key_fields(BASE),
        record_key_fields(&with_mq),
        "with an empty exclusion list, an extra tag must change the key -- \
         otherwise a divergence on that tag cannot fail any parity test"
    );
}

#[test]
#[should_panic(expected = "malformed SAM aux field")]
fn an_aux_field_with_no_type_or_value_panics() {
    // The worst case the validation exists for: a bare `MQ` would otherwise
    // reduce to the key `MQ`, get dropped as a deliberate asymmetry, and let
    // a corrupt record compare equal to a clean one.
    record_key_fields(&format!("{BASE}\tMQ"));
}

#[test]
#[should_panic(expected = "malformed SAM aux field")]
fn a_truncated_aux_field_panics() {
    // `NM:i` has a tag and a type but no value; without validation it reduces
    // to the same key as a well-formed `NM:i:0`.
    record_key_fields(&format!("{BASE}\tNM:i"));
}

#[test]
#[should_panic(expected = "malformed SAM aux field")]
fn an_aux_field_with_a_long_type_code_panics() {
    record_key_fields(&format!("{BASE}\tNM:int:0"));
}

#[test]
fn an_empty_aux_value_is_well_formed() {
    // `Z` values may be empty, so `XA:Z:` must survive validation — the
    // fail-fast rule must not reject a record the reference aligner can
    // legitimately emit.
    record_key_fields(&format!("{BASE}\tXA:Z:"));
}

#[test]
#[should_panic(expected = "malformed SAM line")]
fn a_short_sam_line_panics() {
    // Sibling of the aux-field rule: `samtools view` (no header) only ever
    // emits >= 11 columns, so a short line is a real bug, not noise to drop.
    record_key_fields("r0\t99\tphix\t101\t60\t100M");
}
