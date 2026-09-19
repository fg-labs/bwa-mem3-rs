//! The build script must export what it compiled: the vendored bwa-mem3
//! version (for `@PG VN:`), the C++ compiler line (so a benchmark can record
//! whether it ran on a floor-compliant toolchain), and the x86 kernel tiers.

#[test]
fn version_matches_vendor_file() {
    let vendored = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/bwa-mem3/version.txt"
    ))
    .expect("read vendor/bwa-mem3/version.txt");
    assert_eq!(bwa_mem3_sys::build_info::VERSION, vendored.trim());
}

#[test]
fn compiler_line_is_non_empty() {
    assert!(!bwa_mem3_sys::build_info::COMPILER.trim().is_empty());
}

#[test]
fn x86_tiers_are_reported_per_arch() {
    if cfg!(target_arch = "x86_64") {
        assert_eq!(
            bwa_mem3_sys::build_info::X86_TIERS,
            "sse41,sse42,avx,avx2,avx512bw"
        );
    } else {
        assert_eq!(bwa_mem3_sys::build_info::X86_TIERS, "");
    }
}
