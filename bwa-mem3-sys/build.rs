use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Kernel TUs: bwa-mem3 v0.2.0 compiles these once per SIMD tier on x86_64
/// (`sse41`, `sse42`, `avx`, `avx2`, `avx512bw`) with `-DKERNEL_VARIANT=_<tier>`
/// so each per-tier compile emits mangled symbols (`make_kswv_kernel_avx2`,
/// `ksw_extend2_avx2`, …). The dispatcher wrappers in `simd_dispatch.cpp` (no
/// `KERNEL_VARIANT`) provide unmangled entry points that pick a tier at runtime.
///
/// On arm64 there is only the NEON tier; `simd_dispatch.cpp` calls the unmangled
/// kernels directly, and the kernels are compiled with `KERNEL_VARIANT` unset.
const KERNEL_SRCS: &[&str] = &["bandedSWA.cpp", "kswv.cpp", "ksw.cpp", "sam_encode.cpp"];

/// Per-tier ISA flags for the x86_64 kernel multi-build. Mirrors upstream's
/// `KERNEL_FLAGS_<tier>` groups in `vendor/bwa-mem3/Makefile`, with one
/// addition: `avx512bw` also carries `-mprefer-vector-width=256`. Upstream
/// applies that flag only on its baseline auto-vectorized path, not the kernel
/// TUs (`KERNEL_FLAGS_avx512bw` is just `-mavx -mavx2 -mavx512f -mavx512bw
/// -mpopcnt`). It is harmless here — the kernels use intrinsics, not
/// auto-vec — but capping vector width matches upstream's intent on Skylake-X
/// where 512-bit ops can downclock.
const KERNEL_TIERS_X86: &[(&str, &[&str])] = &[
    ("sse41", &["-msse4.1"]),
    ("sse42", &["-msse4.2", "-mpopcnt"]),
    ("avx", &["-mavx", "-mpopcnt"]),
    ("avx2", &["-mavx", "-mavx2", "-mpopcnt"]),
    (
        "avx512bw",
        &[
            "-mavx",
            "-mavx2",
            "-mavx512f",
            "-mavx512bw",
            "-mpopcnt",
            "-mprefer-vector-width=256",
        ],
    ),
];

#[allow(clippy::too_many_lines)]
fn main() {
    // Match cc's deployment target to cargo's to avoid SIGBUS on macOS arm64.
    // Without this, cc builds with the host SDK's deployment version (e.g.,
    // macOS 26) while cargo links at a lower version (e.g., macOS 11),
    // causing the dynamic linker to fault on binary start.
    if env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos")
        && env::var("MACOSX_DEPLOYMENT_TARGET").is_err()
    {
        env::set_var("MACOSX_DEPLOYMENT_TARGET", "11.0");
    }

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=shim/bwa_shim.cpp");
    println!("cargo:rerun-if-changed=shim/bwa_shim_align.cpp");
    println!("cargo:rerun-if-changed=shim/bwa_shim_layout_assert.cpp");
    println!("cargo:rerun-if-changed=shim/bwa_shim.h");
    println!("cargo:rerun-if-changed=shim/bwa_shim_types.h");
    println!("cargo:rerun-if-changed=vendor/COMMIT");
    println!("cargo:rerun-if-changed=patches");

    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let build_dir = out.join("build");

    // 1. Verify MATE_SORT=0 default in vendored Makefile (shim semantics depend on it),
    // and warn about any other -D define upstream's Makefile grows that build.rs
    // neither mirrors nor explicitly waives.
    let vendor_root = manifest.join("vendor/bwa-mem3");
    let makefile_path = vendor_root.join("Makefile");
    let makefile = fs::read_to_string(&makefile_path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", makefile_path.display(), e));
    assert!(
        makefile
            .lines()
            .any(|l| l.contains("CPPFLAGS") && l.contains("-DMATE_SORT=0")),
        "vendored bwa-mem3 Makefile must retain `-DMATE_SORT=0` default; shim pairing logic depends on it",
    );
    report_unmirrored_makefile_defines(&vendor_root);

    // 2. Copy vendor -> OUT_DIR/build (idempotent; vendor tree is never mutated in place).
    if build_dir.exists() {
        fs::remove_dir_all(&build_dir).unwrap();
    }
    fs::create_dir_all(&build_dir).unwrap();
    copy_dir(
        &manifest.join("vendor/bwa-mem3"),
        &build_dir.join("bwa-mem3"),
    );

    // 3. Apply any patches in patches/ lexicographic order. Expected empty in v1.
    let patches_dir = manifest.join("patches");
    if patches_dir.is_dir() {
        let mut patches: Vec<_> = fs::read_dir(&patches_dir)
            .unwrap()
            .filter_map(Result::ok)
            .filter(|e| {
                e.file_name()
                    .to_string_lossy()
                    .to_lowercase()
                    .ends_with(".patch")
            })
            .collect();
        patches.sort_by_key(std::fs::DirEntry::file_name);
        for p in patches {
            // Vendored source may have CRLF line endings; patch files we ship
            // are LF. `--ignore-whitespace` tolerates that mismatch.
            let status = Command::new("patch")
                .current_dir(build_dir.join("bwa-mem3"))
                .args(["-p1", "--ignore-whitespace", "-i"])
                .arg(p.path())
                .status()
                .expect("patch command not found");
            assert!(
                status.success(),
                "patch {} failed to apply",
                p.path().display()
            );
        }
    }

    let vendor_src = build_dir.join("bwa-mem3/src");
    // bwa-mem3 v0.2.2 (#105) dropped the intel/safestringlib dependency, sweeping
    // every `*_s` call back to plain libc. There is no longer a safestringlib to
    // build or include.
    let s2n = build_dir.join("bwa-mem3/ext/sse2neon");

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let is_x86 = target_arch == "x86_64";

    // 4a. Per-tier kernel TUs (x86_64 only). bwa-mem3 v0.2.0 picks the
    // matching tier at runtime in `simd_dispatch.cpp`; we must supply all
    // five mangled tier copies of bandedSWA/kswv/ksw/sam_encode so the
    // dispatcher's `make_bsw_kernel_<tier>` / `ksw_extend2_<tier>` /
    // `sam_encode_seq_fwd_<tier>` references resolve at link time.
    if is_x86 {
        for (tier, isa_flags) in KERNEL_TIERS_X86 {
            let mut k_build = cc::Build::new();
            k_build.cpp(true);
            for src in KERNEL_SRCS {
                k_build.file(vendor_src.join(src));
            }
            if s2n.is_dir() {
                k_build.include(&s2n);
            }
            k_build.include(manifest.join("shim"));
            k_build.include(&vendor_src);
            for flag in *isa_flags {
                k_build.flag_if_supported(flag);
            }
            k_build.std("c++17");
            k_build.define("ENABLE_PREFETCH", None);
            k_build.define("V17", Some("1"));
            k_build.define("MATE_SORT", Some("0"));
            // kernel_dispatch.h mangles every exported kernel symbol to
            // `<name><KERNEL_VARIANT>` (e.g. `_avx2`). The dispatcher in
            // simd_dispatch.cpp expects this exact suffix per tier.
            k_build.define("KERNEL_VARIANT", Some(format!("_{}", tier).as_str()));
            apply_common_warning_silencing(&mut k_build);
            k_build.compile(&format!("bwa-mem3-kernel-{}", tier));
        }
    }

    // 4b. bwa-mem3 + shim: C++.
    let mut build = cc::Build::new();
    build.cpp(true);

    let skip_common: &[&str] = &[
        "main.cpp",            // CLI entry point
        "bwtindex.cpp",        // index builder; out of scope (users run bwa-mem3 index)
        "bam_writer.cpp",      // bwa-mem3 CLI's htslib-based BAM writer; shim has its own
        "meth_bam.cpp",        // bisulfite BAM writer; htslib-dependent, out of scope
        "fm_index_writer.cpp", // index builder
        "index_prelude.cpp",   // index builder helper
        "libsais_build.cpp",   // index builder (libsais)
        "fastmap.cpp",         // CLI batch driver; see note below
    ];
    // fastmap.cpp is excluded because as of bwa-mem3 0.6.0 it pulls in the new
    // fast_reader FASTQ path (libdeflate + zlib-ng), numa, and htslib that this
    // crate does not vendor. The worker_alloc/worker_free it used to provide are
    // now carried in shim/bwa_shim_align.cpp — see the rationale comment there.
    //
    // The new C reader TUs (fast_reader.c, fast_reader_bseq.c, fr_fastq.c) are
    // not `.cpp`, so the *.cpp glob below never picks them up; they sit unused
    // in the vendor tree.
    //
    // On x86_64 the kernel TUs are compiled separately per tier in 4a
    // (`KERNEL_TIERS_X86`); skip them here so the baseline build doesn't
    // also emit unmangled copies (and so this build inherits no `-mavx512bw`
    // surprises from the baseline ISA). On arm64 the kernels compile once
    // here with KERNEL_VARIANT unset — `simd_dispatch.cpp`'s arm64 branch
    // calls the unmangled symbols directly.
    for entry in fs::read_dir(&vendor_src).unwrap() {
        let e = entry.unwrap();
        let name = e.file_name().to_string_lossy().into_owned();
        if e.path().extension().is_some_and(|x| x == "cpp")
            && !skip_common.iter().any(|s| *s == name)
            && !(is_x86 && KERNEL_SRCS.iter().any(|s| *s == name))
        {
            build.file(e.path());
        }
    }

    if s2n.is_dir() {
        build.include(&s2n);
    }

    build.file(manifest.join("shim/bwa_shim.cpp"));
    build.file(manifest.join("shim/bwa_shim_align.cpp"));
    // Symbol-free TU: static_asserts the POD mem_opt_t/mem_pestat_t layout
    // against upstream's real bwamem.h so a future vendor refresh can't
    // silently misalign the offsets Rust reads/writes (maintenance docs gotcha #2).
    build.file(manifest.join("shim/bwa_shim_layout_assert.cpp"));
    build.include(manifest.join("shim"));
    build.include(&vendor_src);

    apply_simd_flags(&mut build);

    build.std("c++17");
    build.define("ENABLE_PREFETCH", None);
    build.define("V17", Some("1"));
    build.define("MATE_SORT", Some("0"));
    apply_common_warning_silencing(&mut build);

    build.compile("bwa-mem3");

    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
    let target = env::var("TARGET").unwrap();
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }

    // 5. Generate Rust bindings for the shim header.
    generate_bindings(&manifest, &vendor_src, &out);
}

fn apply_common_warning_silencing(build: &mut cc::Build) {
    build.flag_if_supported("-Wno-unused-parameter");
    build.flag_if_supported("-Wno-sign-compare");
    build.flag_if_supported("-Wno-unused-variable");
    build.flag_if_supported("-Wno-unused-function");
    build.flag_if_supported("-Wno-unused-but-set-variable");
    build.flag_if_supported("-Wno-deprecated-declarations");
    build.flag_if_supported("-Wno-format");
    build.flag_if_supported("-Wno-format-truncation");
}

fn apply_simd_flags(build: &mut cc::Build) {
    let arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    match arch.as_str() {
        "x86_64" => {
            if cfg!(feature = "avx512") {
                build.flag_if_supported("-mavx512bw");
                build.flag_if_supported("-mavx512vl");
            } else if cfg!(feature = "sse42") {
                build.flag_if_supported("-msse4.2");
            } else if cfg!(feature = "native") {
                build.flag_if_supported("-march=native");
            } else {
                // default: AVX2 baseline (matches upstream's BASELINE_ARCH).
                build.flag_if_supported("-mavx2");
                build.flag_if_supported("-msse4.1");
            }
        }
        "aarch64" => {
            build.flag_if_supported("-march=armv8-a+simd");
            // Mirror upstream's Makefile SSE2NEON_FLAGS so sse2neon-translated
            // code paths (incl. SSSE3-gated banded-SW kernels) compile in.
            // TODO: if a future bwa-mem3 TU starts using _mm_crc32_*
            // (gated on __SSE4_2__), aarch64 will need `-march=armv8-a+crc`.
            // Not required by the current vendored snapshot.
            build.define("__SSE__", Some("1"));
            build.define("__SSE2__", Some("1"));
            build.define("__SSE3__", Some("1"));
            build.define("__SSSE3__", Some("1"));
            build.define("__SSE4_1__", Some("1"));
            build.define("__SSE4_2__", Some("1"));
        }
        other => panic!("unsupported target arch: {other}"),
    }
}

/// Every `-DNAME[=VALUE]` token in a Makefile's text, conditional blocks
/// included. Deliberately line-oriented and blind to `ifeq` nesting: we want
/// the union of every define upstream can pass, because `build.rs` has to
/// decide about each one whether it mirrors it, deliberately omits it, or has
/// simply never noticed it. Trims a stray leading/trailing `"`: upstream
/// forwards flags through a recursive `$(MAKE)` call as a single quoted shell
/// argument (e.g. `EXTRA_CXXFLAGS="$(EXTRA_CXXFLAGS) -DDISABLE_OUTPUT"`), and a
/// whitespace split leaves the closing quote glued to the last flag. No
/// legitimate `-D` value in this Makefile is itself a quoted string, so the
/// trim never eats a real character.
fn extract_makefile_defines(makefile: &str) -> Vec<String> {
    let mut out: Vec<String> = makefile
        .lines()
        .flat_map(str::split_whitespace)
        .filter_map(|tok| tok.strip_prefix("-D"))
        .map(|def| def.trim_matches('"').to_owned())
        .collect();
    out.sort();
    out.dedup();
    out
}

/// Defines upstream's Makefile passes that this crate deliberately does not.
/// Each entry needs a reason, because an unexplained omission is
/// indistinguishable from an oversight.
const DEFINES_DELIBERATELY_OMITTED: &[(&str, &str)] = &[
    (
        "LIBSAIS_OPENMP",
        "libsais is pruned from the vendor tree (refresh-bwa-mem3.sh DROP_SUBTREES); \
         the index builder it belongs to is out of scope",
    ),
    (
        "CACHE_LINE_BYTES",
        "no mirror in build.rs -- see the tracking issue; changing it is a perf change \
         needing its own benchmark",
    ),
    (
        "DISABLE_BATCHED_MATESW",
        "upstream sets this only for the proto-neon-kswv CI's on/off A/B test of the \
         batched mate-rescue SW port, never for a normal build",
    ),
    (
        "KERNEL_VARIANT",
        "set per tier by build.rs itself (KERNEL_TIERS_X86), not copied from the Makefile",
    ),
    (
        "APPLE_SILICON",
        "every use site ORs it with __aarch64__ / __ARM_NEON (compiler-builtin macros, \
         always true on this target), and simd_compat.h self-defines it under the same \
         condition; mirroring it in build.rs would be a no-op",
    ),
    (
        "USE_MIMALLOC",
        "only read by main.cpp, which build.rs excludes entirely (CLI entry point, out \
         of scope)",
    ),
    (
        "BWAMEM3_TESTING",
        "an opt-in hook in simd_dispatch.cpp for upstream's own regression tooling \
         (TESTING_BUILD=1); off by default upstream, which build.rs matches by never \
         setting it",
    ),
    (
        "STAGE_PROF",
        "opt-in profiling instrumentation (`make STAGE_PROF=1`); the default upstream \
         build, which build.rs mirrors, leaves it unset so the hooks compile to no-ops",
    ),
    (
        "DISABLE_OUTPUT",
        "only read by fastmap.cpp, which build.rs excludes entirely (CLI batch driver, \
         out of scope); used by the PGO profile-build target's compute-only binary",
    ),
    (
        "MI_BUILD_SHARED",
        "a CMake cache variable for configuring ext/mimalloc's own cmake project, not a \
         compiler define for bwa-mem3 sources; mimalloc is pruned from the vendor tree",
    ),
    (
        "MI_BUILD_STATIC",
        "ditto (ext/mimalloc CMake configure flag; mimalloc is pruned from the vendor tree)",
    ),
    (
        "MI_BUILD_OBJECT",
        "ditto (ext/mimalloc CMake configure flag; mimalloc is pruned from the vendor tree)",
    ),
    (
        "MI_BUILD_TESTS",
        "ditto (ext/mimalloc CMake configure flag; mimalloc is pruned from the vendor tree)",
    ),
    (
        "MI_OVERRIDE",
        "ditto (ext/mimalloc CMake configure flag; mimalloc is pruned from the vendor tree)",
    ),
    (
        "CMAKE_BUILD_TYPE",
        "a CMake cache variable shared by the ext/mimalloc and ext/zlib-ng cmake \
         sub-builds, not a compiler define for bwa-mem3 sources; both dependencies are \
         pruned from the vendor tree",
    ),
    (
        "BUILD_SHARED_LIBS",
        "a CMake cache variable for configuring ext/zlib-ng's own cmake project, not a \
         compiler define for bwa-mem3 sources; zlib-ng is pruned from the vendor tree",
    ),
    (
        "ZLIB_COMPAT",
        "ditto (ext/zlib-ng CMake configure flag; zlib-ng is pruned from the vendor tree)",
    ),
    (
        "ZLIB_ENABLE_TESTS",
        "ditto (ext/zlib-ng CMake configure flag; zlib-ng is pruned from the vendor tree)",
    ),
    (
        "ZLIBNG_ENABLE_TESTS",
        "ditto (ext/zlib-ng CMake configure flag; zlib-ng is pruned from the vendor tree)",
    ),
    (
        "WITH_GTEST",
        "ditto (ext/zlib-ng CMake configure flag; zlib-ng is pruned from the vendor tree)",
    ),
    (
        "CMAKE_POSITION_INDEPENDENT_CODE",
        "ditto (ext/zlib-ng CMake configure flag; zlib-ng is pruned from the vendor tree)",
    ),
];

/// Warn (do not fail) when upstream's Makefile grows a `-D` this crate neither
/// mirrors nor explicitly waives, or changes the *value* of one this crate
/// does mirror. A hard panic would be wrong here: the honest answer for a new
/// define is usually "a human must decide", and failing the build of an
/// otherwise-working refresh helps nobody. `cargo build` surfaces
/// `cargo:warning=` lines prominently, and the vendor-bump drift report
/// re-reports them.
fn report_unmirrored_makefile_defines(vendor_root: &Path) {
    let makefile = vendor_root.join("Makefile");
    let Ok(text) = fs::read_to_string(&makefile) else {
        println!(
            "cargo:warning=could not read {} to check define drift",
            makefile.display()
        );
        return;
    };
    // Names build.rs mirrors via a `.define()` call somewhere in this file,
    // paired with the value it compiles with: ENABLE_PREFETCH/V17/MATE_SORT
    // are the constant trio baked into every compile (main's steps 4a/4b);
    // the six SSE feature macros are reproduced one-for-one by
    // apply_simd_flags's aarch64 branch, which mirrors upstream's
    // SSE2NEON_FLAGS bridge.
    //
    // The value is carried here, not just the name, because a name-only
    // comparison passes an upstream flip of `-DV17=1` to `-DV17=0` as
    // "mirrored" while build.rs goes on compiling with the old value — the
    // exact silent divergence this check exists to surface. `ENABLE_PREFETCH`
    // is recorded as "1" because a bare `-DNAME` is what both the
    // preprocessor and `cc`'s `.define(name, None)` treat as 1.
    let mirrored = [
        ("ENABLE_PREFETCH", "1"),
        ("V17", "1"),
        ("MATE_SORT", "0"),
        ("__SSE__", "1"),
        ("__SSE2__", "1"),
        ("__SSE3__", "1"),
        ("__SSSE3__", "1"),
        ("__SSE4_1__", "1"),
        ("__SSE4_2__", "1"),
    ];
    for def in extract_makefile_defines(&text) {
        let (name, value) = canonical_define(&def);
        if let Some((_, mirrored_value)) = mirrored.iter().find(|(n, _)| *n == name) {
            if value != *mirrored_value {
                println!(
                    "cargo:warning=vendored Makefile now passes -D{def}, but build.rs \
                     compiles with {name}={mirrored_value}; update the .define() call \
                     and this list together, or waive the difference in \
                     DEFINES_DELIBERATELY_OMITTED"
                );
            }
            continue;
        }
        if DEFINES_DELIBERATELY_OMITTED.iter().any(|(n, _)| *n == name) {
            continue;
        }
        println!(
            "cargo:warning=vendored Makefile passes -D{def} which build.rs \
             neither mirrors nor waives; decide and add it to build.rs's \
             .define() calls or to DEFINES_DELIBERATELY_OMITTED"
        );
    }
}

/// Split one extracted `-D` body into `(name, value)`, giving a valueless
/// define the `1` the preprocessor gives it — so `-DV17` and `-DV17=1` compare
/// equal, and a mirrored define's value can be checked at all.
fn canonical_define(def: &str) -> (&str, &str) {
    def.split_once('=').unwrap_or((def, "1"))
}

fn copy_dir(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap();
    for entry in fs::read_dir(src).unwrap() {
        let e = entry.unwrap();
        let ft = e.file_type().unwrap();
        let from = e.path();
        let to = dst.join(e.file_name());
        if ft.is_dir() {
            copy_dir(&from, &to);
        } else {
            fs::copy(&from, &to).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{canonical_define, extract_makefile_defines};

    /// A Makefile excerpt in upstream's shape: unconditional CPPFLAGS plus the
    /// arch-conditional blocks that build.rs mirrors by hand.
    const MAKEFILE: &str = "\
CPPFLAGS+= -DENABLE_PREFETCH -DV17=1 -DMATE_SORT=0 -DLIBSAIS_OPENMP
ifeq ($(arch),arm64)
\tSSE2NEON_FLAGS= -D__SSE2NEON_H -DSSE2NEON_SUPPRESS_WARNINGS
\tCPPFLAGS+= -DCACHE_LINE_BYTES=128
endif
";

    #[test]
    fn extracts_unconditional_defines() {
        let d = extract_makefile_defines(MAKEFILE);
        assert!(d.contains(&"ENABLE_PREFETCH".to_string()));
        assert!(d.contains(&"V17=1".to_string()));
        assert!(d.contains(&"MATE_SORT=0".to_string()));
    }

    #[test]
    fn extracts_conditional_defines_too() {
        let d = extract_makefile_defines(MAKEFILE);
        assert!(
            d.contains(&"CACHE_LINE_BYTES=128".to_string()),
            "arch-conditional defines must be extracted; they are the ones \
             build.rs mirrors by hand"
        );
        assert!(d.contains(&"__SSE2NEON_H".to_string()));
    }

    #[test]
    fn ignores_non_define_flags() {
        let d = extract_makefile_defines("CPPFLAGS+= -O3 -Wall -DFOO=1\n");
        assert_eq!(d, vec!["FOO=1".to_string()]);
    }

    /// Regression test for a real false positive found in the vendored
    /// Makefile (v0.8.0-ish, commit a887e36): line 943 forwards
    /// `EXTRA_CXXFLAGS="$(EXTRA_CXXFLAGS) -DDISABLE_OUTPUT"` to a recursive
    /// `$(MAKE)` invocation. A naive whitespace split glues the closing quote
    /// to the flag, producing the bogus name `DISABLE_OUTPUT"` that no waiver
    /// entry could ever match by name.
    #[test]
    fn strips_quote_artifacts_from_forwarded_flags() {
        let d = extract_makefile_defines(
            "\t$(MAKE) EXTRA_CXXFLAGS=\"$(EXTRA_CXXFLAGS) -DDISABLE_OUTPUT\" CXX=\"$(CXX)\" all\n",
        );
        assert!(d.contains(&"DISABLE_OUTPUT".to_string()));
    }

    #[test]
    fn a_valueless_define_canonicalizes_to_one() {
        // `-DENABLE_PREFETCH` and `-DENABLE_PREFETCH=1` must compare equal, or
        // the mirrored-value check reports drift on an unchanged Makefile.
        assert_eq!(
            canonical_define("ENABLE_PREFETCH"),
            ("ENABLE_PREFETCH", "1")
        );
        assert_eq!(
            canonical_define("ENABLE_PREFETCH=1"),
            ("ENABLE_PREFETCH", "1")
        );
    }

    #[test]
    fn a_changed_value_is_visible_to_the_mirror_check() {
        // The case a name-only comparison misses: upstream flipping V17 or
        // MATE_SORT while build.rs keeps compiling with the old value.
        assert_eq!(canonical_define("V17=0"), ("V17", "0"));
        assert_eq!(canonical_define("MATE_SORT=1"), ("MATE_SORT", "1"));
    }

    #[test]
    fn only_the_first_equals_splits_name_from_value() {
        // `-DKERNEL_VARIANT=_avx2` and make-variable forwards like
        // `-DX=$(Y)` must keep their full value rather than being truncated.
        assert_eq!(
            canonical_define("KERNEL_VARIANT=_avx2"),
            ("KERNEL_VARIANT", "_avx2")
        );
        assert_eq!(
            canonical_define("DISABLE_BATCHED_MATESW=$(DISABLE_BATCHED_MATESW)"),
            ("DISABLE_BATCHED_MATESW", "$(DISABLE_BATCHED_MATESW)")
        );
    }
}

fn generate_bindings(manifest: &Path, _vendor_src: &Path, out: &Path) {
    let shim_dir = manifest.join("shim");
    let bindings = bindgen::Builder::default()
        .header(shim_dir.join("bwa_shim.h").to_string_lossy())
        .clang_arg(format!("-I{}", shim_dir.display()))
        .allowlist_type("BwaReadPair")
        .allowlist_type("BwaIndex")
        .allowlist_type("BwaSeeds")
        .allowlist_type("BwaBatch")
        .allowlist_type("mem_opt_t")
        .allowlist_type("mem_pestat_t")
        .allowlist_function("bwa_shim_.*")
        .allowlist_var("MEM_F_.*")
        .derive_default(true)
        .generate()
        .expect("bindgen failed");
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("failed to write bindings.rs");
}
