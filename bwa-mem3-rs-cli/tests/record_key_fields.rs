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
fn deliberately_asymmetric_tags_do_not_change_the_key() {
    // Upstream's SAM writer emits MQ:i whenever a mate pointer is present
    // (i.e. on essentially every record of a properly-paired alignment);
    // the shim never emits it. That asymmetry is by design, not a
    // regression, so a record with MQ:i and an otherwise-identical record
    // without it must still reduce to the same key -- otherwise every
    // paired parity test would read as a mismatch on a machine that has the
    // real bwa-mem3 CLI installed.
    let with_mq = format!("{BASE}\tMQ:i:60");
    assert_eq!(
        record_key_fields(BASE),
        record_key_fields(&with_mq),
        "MQ:i differs between the shim and upstream by design and must not change the key"
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
