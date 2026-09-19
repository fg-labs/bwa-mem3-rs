// Compiler-floor detection shared by `build.rs` (via `include!`) and the
// crate's unit tests (via `#[path]`). Pure functions only: no `cargo:` output.
//
// NB: this is a plain `//` comment, not a `//!` inner doc comment, because
// `build.rs` splices this file's tokens in via `include!` after its own
// `use` statements -- an inner doc comment there is a syntax error
// (E0753: expected outer doc comment), even though the same file compiles
// fine as a `//!`-documented module when pulled in via `#[path]` in `lib.rs`.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Toolchain {
    Clang,
    Gcc,
    Other,
}

/// Floors from the vendored bwa-mem3 `Makefile` (`unsupported compiler:
/// clang 15 < 19; set ALLOW_UNSUPPORTED_COMPILER=1 to override`). The crate
/// compiles below them — `cc` enforces nothing — but the result is slower
/// than the CLI and must not be benchmarked as representative.
pub const CLANG_FLOOR: u32 = 19;
pub const GCC_FLOOR: u32 = 15;

/// Classify the first line of `<cxx> --version` and pull the major version.
pub fn parse_version_line(line: &str) -> (Toolchain, Option<u32>) {
    let lower = line.to_ascii_lowercase();
    let toolchain = if lower.contains("clang") {
        Toolchain::Clang
    } else if lower.contains("gcc") || lower.starts_with("g++") || lower.starts_with("c++") {
        Toolchain::Gcc
    } else {
        return (Toolchain::Other, None);
    };
    // For clang the version follows the literal "version"; for gcc it is the
    // last `N.N.N` token on the line (the parenthesised vendor string may hold
    // another one, e.g. "(Ubuntu 15.1.0-1ubuntu1) 15.1.0").
    let major = match toolchain {
        Toolchain::Clang => lower
            .split_whitespace()
            .skip_while(|t| *t != "version")
            .nth(1)
            .and_then(|v| v.split('.').next())
            .and_then(|m| m.parse().ok()),
        Toolchain::Gcc => lower
            .split_whitespace()
            .rev()
            .find(|t| {
                t.split('.').count() >= 2 && t.chars().next().is_some_and(|c| c.is_ascii_digit())
            })
            .and_then(|v| v.split('.').next())
            .and_then(|m| m.parse().ok()),
        Toolchain::Other => None,
    };
    (toolchain, major)
}

/// `Some(warning)` when the compiler is below the floor or could not be
/// identified; `None` when it meets the floor.
pub fn floor_warning(toolchain: Toolchain, major: Option<u32>) -> Option<String> {
    let (name, floor) = match toolchain {
        Toolchain::Clang => ("clang", CLANG_FLOOR),
        Toolchain::Gcc => ("gcc", GCC_FLOOR),
        Toolchain::Other => {
            return Some(
                "bwa-mem3-sys: could not identify the C++ compiler; bwa-mem3 requires \
                 clang >= 19 or gcc >= 15 for benchmark-representative code"
                    .to_string(),
            )
        }
    };
    match major {
        Some(m) if m >= floor => None,
        Some(m) => Some(format!(
            "bwa-mem3-sys: unsupported compiler: {name} {m} < {floor}; the build will \
             succeed but is not benchmark-representative (old compilers emit redundant \
             SIMD that a modern clang elides). Install clang >= 19 or gcc >= 15."
        )),
        None => Some(format!(
            "bwa-mem3-sys: could not parse the {name} version; bwa-mem3 requires \
             {name} >= {floor} for benchmark-representative code"
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_clang_major() {
        assert_eq!(
            parse_version_line("Apple clang version 17.0.0 (clang-1700.0.13.3)"),
            (Toolchain::Clang, Some(17))
        );
        assert_eq!(
            parse_version_line("clang version 19.1.7"),
            (Toolchain::Clang, Some(19))
        );
    }

    #[test]
    fn parses_gcc_major() {
        assert_eq!(
            parse_version_line("g++ (GCC) 11.5.0 20240719 (Red Hat 11.5.0-5)"),
            (Toolchain::Gcc, Some(11))
        );
        assert_eq!(
            parse_version_line("g++ (Ubuntu 15.1.0-1ubuntu1) 15.1.0"),
            (Toolchain::Gcc, Some(15))
        );
        // GCC invoked through the `c++` alternative (Ubuntu's default C++ driver
        // symlink) banners as "c++ (...)"; a supported GCC must still classify
        // as Gcc, not Other. (Clang-as-c++ banners contain "clang" and are
        // caught by the clang branch first, so this cannot capture clang.)
        assert_eq!(
            parse_version_line("c++ (Ubuntu 15.1.0-1ubuntu1) 15.1.0"),
            (Toolchain::Gcc, Some(15))
        );
    }

    #[test]
    fn unknown_compiler_has_no_major() {
        assert_eq!(parse_version_line("icpx 2025.0"), (Toolchain::Other, None));
    }

    #[test]
    fn floor_warning_fires_below_floor_only() {
        assert!(floor_warning(Toolchain::Clang, Some(18)).is_some());
        assert!(floor_warning(Toolchain::Clang, Some(19)).is_none());
        assert!(floor_warning(Toolchain::Gcc, Some(14)).is_some());
        assert!(floor_warning(Toolchain::Gcc, Some(15)).is_none());
        // Unknown toolchain or unparseable version: warn that we could not check.
        assert!(floor_warning(Toolchain::Other, None).is_some());
    }

    #[test]
    fn floor_warning_names_the_floor() {
        let w = floor_warning(Toolchain::Clang, Some(15)).unwrap();
        assert!(w.contains("clang 15 < 19"), "{w}");
        assert!(w.contains("not benchmark-representative"), "{w}");
    }
}
