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
    println!("cargo:rerun-if-changed=shim/bwa_shim.h");
    println!("cargo:rerun-if-changed=shim/bwa_shim_types.h");
    println!("cargo:rerun-if-changed=vendor/COMMIT");
    println!("cargo:rerun-if-changed=patches");

    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let build_dir = out.join("build");

    // 1. Verify MATE_SORT=0 default in vendored Makefile (shim semantics depend on it).
    let makefile_path = manifest.join("vendor/bwa-mem3/Makefile");
    let makefile = fs::read_to_string(&makefile_path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", makefile_path.display(), e));
    assert!(
        makefile
            .lines()
            .any(|l| l.contains("CPPFLAGS") && l.contains("-DMATE_SORT=0")),
        "vendored bwa-mem3 Makefile must retain `-DMATE_SORT=0` default; shim pairing logic depends on it",
    );

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
    ];
    // fastmap.cpp used to be excluded (CLI-side batch driver) but is now
    // built to expose worker_alloc/worker_free. Its entry point is
    // `main_mem`, not `main`, so no collision with the Rust test harness.
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
