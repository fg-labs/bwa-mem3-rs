//! What the build script compiled. Read by `bwa_mem3_rs::build_info()`.

/// Vendored bwa-mem3 version (`vendor/bwa-mem3/version.txt`), e.g. `"0.9.0"`.
pub const VERSION: &str = env!("BWA_MEM3_SYS_VERSION");
/// First line of `<cxx> --version` for the compiler that built the C++.
pub const COMPILER: &str = env!("BWA_MEM3_SYS_COMPILER");
/// Comma-separated x86_64 kernel tiers compiled in; empty on other arches.
pub const X86_TIERS: &str = env!("BWA_MEM3_SYS_X86_TIERS");
