//! Unit-tests the require-tools guard used by every CLI-parity test.
//!
//! `BWA_MEM3_RS_REQUIRE_TOOLS` is what turns a "tool missing → skip" into a
//! "tool missing → fail". CI sets it in the one job that installs `bwa-mem3`
//! and `samtools`, so a regression that removes the tools turns the job red
//! instead of silently green-with-skips. Never gate this on `CI`: `cargo
//! ci-test` runs these targets in a job that deliberately has no `bwa-mem3`.

mod common;

/// Restores a set of environment variables to their prior values on drop.
///
/// `require_tools_var_controls_skip_vs_panic` mutates process-global env vars
/// (including `PATH`), and an `assert!` failure partway through unwinds the
/// stack rather than running any explicit cleanup code that follows it. Tying
/// the restore to a `Drop` guard means the destructor still runs during that
/// unwind, so a failing assertion can't leave the mutated `PATH` /
/// `BWA_MEM3_BIN` / `BWA_MEM3_RS_REQUIRE_TOOLS` values in place for whichever
/// other test thread runs next in this process.
struct EnvRestore(Vec<(&'static str, Option<std::ffi::OsString>)>);

impl EnvRestore {
    fn capture(vars: &[&'static str]) -> Self {
        Self(vars.iter().map(|&k| (k, std::env::var_os(k))).collect())
    }
}

impl Drop for EnvRestore {
    fn drop(&mut self) {
        for (k, v) in &self.0 {
            match v {
                Some(v) => std::env::set_var(k, v),
                None => std::env::remove_var(k),
            }
        }
    }
}

/// Both directions are asserted in a single test because the env var is
/// process-global and `cargo test` runs tests in parallel threads.
#[test]
fn require_tools_var_controls_skip_vs_panic() {
    // Point BWA_MEM3_BIN at a path that cannot exist AND neutralize PATH via
    // an empty scratch directory. Both are required: find_bwa_mem3() falls
    // back to `which bwa-mem3` whenever the BWA_MEM3_BIN path doesn't exist,
    // so overriding BWA_MEM3_BIN alone leaves the "missing" precondition
    // false whenever a real `bwa-mem3` happens to be on PATH already. With
    // both overridden, the tool is "missing" regardless of what is installed
    // on the machine running the suite.
    let missing = concat!(env!("CARGO_MANIFEST_DIR"), "/does-not-exist-bwa-mem3");
    let empty_path_dir = tempfile::tempdir().expect("create empty PATH scratch dir");

    // Captured before any mutation and restored on drop — including on the
    // unwind from a failed assert! below — so this test can't leak a
    // clobbered PATH/BWA_MEM3_BIN/BWA_MEM3_RS_REQUIRE_TOOLS into whatever
    // other test runs next in this process.
    let _restore = EnvRestore::capture(&["PATH", "BWA_MEM3_BIN", "BWA_MEM3_RS_REQUIRE_TOOLS"]);
    std::env::set_var("PATH", empty_path_dir.path());
    std::env::set_var("BWA_MEM3_BIN", missing);
    std::env::remove_var("BWA_MEM3_RS_REQUIRE_TOOLS");

    // Unset: missing tool yields None (caller skips).
    assert!(
        common::require_bwa_mem3().is_none(),
        "without BWA_MEM3_RS_REQUIRE_TOOLS a missing binary must skip, not panic"
    );

    // Set: missing tool panics.
    std::env::set_var("BWA_MEM3_RS_REQUIRE_TOOLS", "1");
    let panicked = std::panic::catch_unwind(common::require_bwa_mem3).is_err();
    assert!(
        panicked,
        "with BWA_MEM3_RS_REQUIRE_TOOLS set, a missing binary must panic"
    );
}
