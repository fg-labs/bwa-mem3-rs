//! The committed bindings must declare every `bwa_shim_*` prototype in the
//! header (a cheap textual guard; the authoritative check is the
//! `regenerate-bindings` build in CI).

#[test]
fn every_header_prototype_has_a_binding() {
    let header =
        std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/shim/bwa_shim.h")).unwrap();
    let bindings =
        std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/bindings.rs")).unwrap();
    let mut missing = Vec::new();
    for line in header.lines() {
        // Skip `#include "bwa_shim_types.h"` etc.: the filename tokenizes to
        // `bwa_shim_types`, which matches the `bwa_shim_` prefix but is not a
        // prototype.
        if line.trim_start().starts_with("#include") {
            continue;
        }
        for tok in line.split(|c: char| !(c.is_alphanumeric() || c == '_')) {
            if tok.starts_with("bwa_shim_") && !bindings.contains(&format!("pub fn {tok}(")) {
                missing.push(tok.to_string());
            }
        }
    }
    missing.sort();
    missing.dedup();
    assert!(
        missing.is_empty(),
        "prototypes without a committed binding: {missing:?}"
    );
}
