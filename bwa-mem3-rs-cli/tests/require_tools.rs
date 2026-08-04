//! Unit-tests the require-tools guard used by every CLI-parity test.
//!
//! `BWA_MEM3_RS_REQUIRE_TOOLS` is what turns a "tool missing → skip" into a
//! "tool missing → fail". CI sets it in the one job that installs `bwa-mem3`
//! and `samtools`, so a regression that removes the tools turns the job red
//! instead of silently green-with-skips. Never gate this on `CI`: `cargo
//! ci-test` runs these targets in a job that deliberately has no `bwa-mem3`.

mod common;

/// Both directions are asserted in a single test because the env var is
/// process-global and `cargo test` runs tests in parallel threads.
#[test]
fn require_tools_var_controls_skip_vs_panic() {
    // Point the lookup at a path that cannot exist, so the tool is "missing"
    // regardless of what is installed on the machine running the suite.
    let missing = concat!(env!("CARGO_MANIFEST_DIR"), "/does-not-exist-bwa-mem3");

    // Unset: missing tool yields None (caller skips).
    std::env::set_var("BWA_MEM3_BIN", missing);
    std::env::remove_var("BWA_MEM3_RS_REQUIRE_TOOLS");
    assert!(
        common::require_bwa_mem3().is_none(),
        "without BWA_MEM3_RS_REQUIRE_TOOLS a missing binary must skip, not panic"
    );

    // Set: missing tool panics.
    std::env::set_var("BWA_MEM3_RS_REQUIRE_TOOLS", "1");
    let panicked = std::panic::catch_unwind(common::require_bwa_mem3).is_err();
    std::env::remove_var("BWA_MEM3_RS_REQUIRE_TOOLS");
    std::env::remove_var("BWA_MEM3_BIN");
    assert!(
        panicked,
        "with BWA_MEM3_RS_REQUIRE_TOOLS set, a missing binary must panic"
    );
}
