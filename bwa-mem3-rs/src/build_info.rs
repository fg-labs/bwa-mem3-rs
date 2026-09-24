//! Provenance of the linked bwa-mem3: for `@PG VN:` and benchmark records.

/// Vendored bwa-mem3 version, e.g. `"0.9.0"`.
#[must_use]
pub fn version() -> &'static str {
    bwa_mem3_sys::build_info::VERSION
}

/// Snapshot of what `bwa-mem3-sys` compiled against. See [`build_info`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BuildInfo {
    /// Vendored bwa-mem3 version, e.g. `"0.9.0"`.
    pub bwa_mem3_version: &'static str,
    /// First line of `<cxx> --version` for the compiler that built the C++.
    pub compiler: &'static str,
    /// x86_64 kernel tiers compiled in (runtime-dispatched); empty elsewhere.
    pub x86_tiers: &'static [&'static str],
}

/// Provenance of the linked bwa-mem3, for `@PG VN:` lines and benchmark
/// records that need to know exactly what was compiled.
#[must_use]
pub fn build_info() -> BuildInfo {
    static TIERS: std::sync::OnceLock<Vec<&'static str>> = std::sync::OnceLock::new();
    let tiers = TIERS.get_or_init(|| {
        bwa_mem3_sys::build_info::X86_TIERS
            .split(',')
            .filter(|t| !t.is_empty())
            .collect()
    });
    BuildInfo {
        bwa_mem3_version: version(),
        compiler: bwa_mem3_sys::build_info::COMPILER,
        x86_tiers: tiers,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_matches_sys_crate() {
        assert_eq!(version(), bwa_mem3_sys::build_info::VERSION);
    }

    #[test]
    fn build_info_is_populated() {
        let b = build_info();
        assert_eq!(b.bwa_mem3_version, version());
        assert!(!b.compiler.is_empty());
        if cfg!(target_arch = "x86_64") {
            assert_eq!(b.x86_tiers.len(), 5);
        } else {
            assert!(b.x86_tiers.is_empty());
        }
    }
}
