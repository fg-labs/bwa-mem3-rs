//! `RG:Z` through the structured-fields sink. A separate test binary on
//! purpose: `MemOpts::set_read_group_id` writes bwa's process-global read
//! group, so setting it in a shared binary would leak into every other test.

use bwa_mem3_rs::MemOpts;

mod support;
use support::*;

/// With a read group set, every record carries `RG:Z` in both emission paths,
/// and the structured fields still rebuild the packed record byte-for-byte.
#[test]
fn pair_emit_fields_matches_packed_with_a_read_group() {
    let Some(idx) = shared_idx() else {
        eprintln!("skip: bwa-mem3 not available to build a PhiX index");
        return;
    };
    let fixture = field_fixture();
    let mut opts = MemOpts::new().unwrap();
    opts.set_read_group_id(Some("rg1")).unwrap();
    let mut coverage = Coverage::default();
    assert_fields_match_packed(&idx, &opts, &fixture, 7, "RG", &mut coverage);
    let total = coverage.singles + fixture.pairs.len();
    assert!(coverage.rg > 0, "no record carried RG:Z");
    assert!(
        coverage.rg >= total,
        "only {} of at least {total} records carried RG:Z",
        coverage.rg
    );
}
