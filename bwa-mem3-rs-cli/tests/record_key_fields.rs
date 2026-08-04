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
