# CLAUDE.md

Project-specific guidance for AI assistants working in this repository.

## Project overview

`bwa-mem3-rs` is a Rust FFI crate for [bwa-mem3] that emits packed BAM records
(per the BAM spec) and leaves all parallelism to the caller. It ships as a
three-crate workspace:

- `bwa-mem3-sys` — unsafe FFI + vendored C/C++ + custom shim.
- `bwa-mem3-rs` — safe wrapper.
- `bwa-mem3-rs-cli` — thin `bwa-rs` CLI binary.

Hosted under [`fg-labs/bwa-mem3-rs`](https://github.com/fg-labs/bwa-mem3-rs).
Vendors [`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3) at the
`main` branch.

[bwa-mem3]: https://github.com/fg-labs/bwa-mem3

## Build and test

```bash
cargo ci-build     # cargo build --workspace --all-targets --locked
cargo ci-test      # cargo test --workspace --locked
cargo ci-fmt       # cargo fmt --all -- --check
cargo ci-lint      # cargo clippy --workspace --all-targets -- -D warnings

# Integration tests (require a prebuilt bwa-mem3 index; skipped otherwise):
BWA_MEM3_RS_TEST_REF=/Users/nhomer/work/references/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
    cargo test --workspace

# E2E (requires bwa-mem3 CLI + samtools on PATH or BWA_MEM3_BIN set):
cargo test -p bwa-mem3-rs-cli --test e2e
```

## Architecture

Ownership layers:

| Layer | What's there |
|---|---|
| `fg-labs/bwa-mem3` / `main` branch | Our integration fork of upstream bwa-mem3. Everything-we've-merged: Apple Silicon NEON, Linux-arm64 Makefile, drop-stat fix, future bwa-meth / XB / etc. |
| `bwa-mem3-sys/vendor/bwa-mem3/` | Pristine snapshot of `fg-labs/bwa-mem3@main` at `vendor/COMMIT`. Never edit in place. Refresh via `scripts/refresh-bwa-mem3.sh`. |
| `bwa-mem3-sys/patches/` | Numbered `.patch` files applied to the vendored source at build time. Goal state: empty (fixes landed in `main` instead). |
| `bwa-mem3-sys/shim/` | Our C/C++ shim. `bwa_shim.h` is the public header bindgen consumes. `bwa_shim.cpp` uses our POD copies of `mem_opt_t` / `mem_pestat_t` from `bwa_shim_types.h`. `bwa_shim_align.cpp` includes upstream's real `bwamem.h` / `FMI_search.h` and bridges to `bwa_shim.cpp` via opaque pointers (layouts match; the POD-vs-upstream match is guarded at build time by `bwa_shim_layout_assert.cpp`, and bindgen's own tests check the Rust-vs-POD match). |
| `bwa-mem3-rs/` | Safe Rust API: `BwaIndex`, `MemOpts`, `MemPeStat`, `ReadPair`, `Seeds`, `AlignmentBatch`, `Record`. Phase split: `seed_batch` → `extend_batch`. Convenience: `align_batch`. Pilot-only: `estimate_pestat`. |
| `bwa-mem3-rs-cli/` | Minimal `bwa-rs mem` CLI. Reads paired FASTQ (gzip-aware), writes proper BGZF-BAM. |

## Gotchas

### 1. Vendored Makefile must retain `MATE_SORT=0`

Our shim delegates the paired-end decision to upstream's `mem_pair_resolve` (exposed by fg-labs/bwa-mem3 PR #9) and then runs our own BAM emission. `mem_pair_resolve`'s internal branching is guarded on `#if MATE_SORT` vs. the default; only the `MATE_SORT=0` path is exercised (and audited) by our shim. If `main` ever flips the default, the pairing logic would swap to an untested branch. `build.rs` asserts `-DMATE_SORT=0` at build time.

### 2. `mem_opt_t` / `mem_pestat_t` layouts are mirrored in two places

They're in `shim/bwa_shim_types.h` (what bindgen reads) and in upstream's `bwamem.h` (what `bwa_shim_align.cpp` includes). Both must stay byte-identical — the shim allocates the real struct via upstream's `mem_opt_init()`, so any field-order/size drift silently corrupts the offsets Rust reads/writes through the bindgen view. **`shim/bwa_shim_layout_assert.cpp` guards this**: it `#include`s the POD under renamed tags (`#define mem_opt_t pod_mem_opt_t`) and then the real `bwamem.h`, and `static_assert`s `offsetof`/`sizeof` field by field, so drift fails `cargo build`. (The bindgen `bindgen_test_layout_mem_opt_t` test only checks the Rust struct against the C compiler's view of the *POD copy*, **not** the POD against upstream — that gap is what the guard TU closes.) On `refresh-bwa-mem3.sh`, diff `vendor/bwa-mem3/src/bwamem.h` around lines 95–156 (`mem_opt_t`) and 239–243 (`mem_pestat_t`); if either changed, update `shim/bwa_shim_types.h` **and** the field list in `bwa_shim_layout_assert.cpp` to match.

### 3. macOS deployment target mismatch → SIGBUS at test-binary startup

`build.rs` sets `MACOSX_DEPLOYMENT_TARGET=11.0` explicitly when building on macOS to keep `cc`'s emitted objects aligned with what `cargo`/`rustc` links. Without this, linked binaries can fault at startup on macOS 26+.

### 4. Shadowing libc

Upstream `bwamem.cpp` had an unused file-scope `int stat;` that shadowed libc's `stat()` syscall wrapper → SIGBUS on test-harness startup. The fix is carried on `main`; should be fixed-forward there rather than as a patch in our crate.

### 5. Some libbwa-mem3 symbols have C++ linkage

`mem_matesw` (and various internal helpers) are not `extern "C"` in `bwamem.h`. If you add forward declarations in the shim for any such symbol, put them outside `extern "C"` blocks so the mangled names match. As of the `mem_pair_resolve` adoption, the shim no longer forward-declares any of these directly.

### 6. bwa-mem3 CIGAR opcodes use a 5-char table, not BAM spec

bwa-mem3's internal `mem_aln_t.cigar` uses opcode table `MIDSH` (M=0 I=1 D=2 S=3 H=4). The BAM spec uses `MIDNSHP=X` (M=0 I=1 D=2 N=3 S=4 H=5). Our emitter (`bwa_cigar_to_bam` in `bwa_shim_align.cpp`) remaps before writing packed BAM. Don't copy bwa-mem3's opcodes verbatim into BAM output.

### 7. `s->seq` after `mem_kernel1_core` is 2-bit encoded, not ASCII

`mem_kernel1_core` rewrites each read's `seq` in place via `nst_nt4_table`: bytes become 0=A, 1=C, 2=G, 3=T, 4=N. Our emitter uses the `bwa2_to_bam4` table to map these to BAM 4-bit nibbles (1/2/4/8/15), plus `bwa2_complement` for the reverse-strand path.

### 8. Files excluded from the C++ build (see `build.rs`)

- `main.cpp` — CLI entry point
- `bwtindex.cpp` — index builder (not our concern)
- `bam_writer.cpp`, `meth_bam.cpp` — htslib-dependent; shim emits BAM directly
- `fm_index_writer.cpp`, `index_prelude.cpp`, `libsais_build.cpp` — index builder
- `fastmap.cpp` — CLI batch driver (see below)
- `fast_reader.c`, `fast_reader_bseq.c`, `fr_fastq.c` — the 0.6.0 fast FASTQ
  reader (libdeflate + zlib-ng). These are `.c`, and `build.rs` only globs
  `src/*.cpp`, so they are never picked up; they sit unused in the vendor tree.

The htslib- and libsais-dependent TUs would also fail to link: the
refresh script prunes `ext/htslib`, `ext/libsais`, `ext/mimalloc`,
`ext/doctest`, and (as of 0.6.0) `ext/zlib-ng` from the vendor tree. If a
future refresh wants to re-enable any of them, restore the relevant
submodule first (or `scripts/refresh-bwa-mem3.sh`'s `DROP_SUBTREES`).

`fastmap.cpp` was compiled at 0.2.x purely to export `worker_alloc` /
`worker_free`. As of bwa-mem3 0.6.0 it is transitively coupled to the new
`fast_reader` FASTQ path (libdeflate + zlib-ng), `numa`, and htslib — none of
which this crate vendors — so it is excluded again. We instead carry verbatim
file-local copies of just `worker_alloc` / `worker_free` in
`shim/bwa_shim_align.cpp`; **keep them in sync** with upstream's `fastmap.cpp`
on every vendor refresh (the buffer set must match what `mem_kernel1_core` /
`mem_kernel2_core` expect).

`runsimd.cpp` is gone in bwa-mem3 v0.2.0 — the multi-binary launcher (`bwa-mem3.<tier>` companions) was replaced by single-binary SIMD dispatch in `simd_dispatch.cpp`. The build no longer needs to skip an unguarded `main()`.

### 9. Per-tier kernel build on x86_64 (v0.2.0+)

bwa-mem3 v0.2.0 splits four SIMD-bearing TUs — `bandedSWA.cpp`, `kswv.cpp`,
`ksw.cpp`, `sam_encode.cpp` — into per-tier kernel compilations. On x86_64
each is compiled once per tier (`sse41`, `sse42`, `avx`, `avx2`, `avx512bw`)
with `-DKERNEL_VARIANT=_<tier>` plus the tier's `-m...` flag set; `kernel_dispatch.h`
mangles every exported symbol to `<name><tier>`. `simd_dispatch.cpp` (compiled
once at the baseline) provides unmangled wrappers (`make_banded_pair_wise_sw`,
`make_kswv`, `ksw_extend2`, `sam_encode_seq_fwd`, …) that pick a tier at
runtime via `__builtin_cpu_supports`. `build.rs` walks
`KERNEL_TIERS_X86` to emit the five `bwa-mem3-kernel-<tier>` archives and
excludes the kernel TUs from the baseline build. On aarch64, kernel TUs are
compiled once with `KERNEL_VARIANT` unset and the dispatcher's `#else` branch
calls them directly.

### 10. The Python stub mirrors `bwa-mem3-py/src/lib.rs` by hand

`bwa-mem3-py/python/bwa_mem3/_bwa_mem3.pyi` is a hand-written type stub for
the compiled extension; it ships in the wheel (next to `py.typed`) and is
what mypy type-checks against. mypy only verifies the stub is *self*-
consistent — it never imports the `.so` — so a stub that drifts from the
PyO3 bindings passes silently. When you add or change a `#[pyclass]`,
`#[pymethods]` member, getter/setter, `#[pyfunction]`, or the `shm`
submodule in `lib.rs`, update the `.pyi` in the same commit. (`mypy.stubtest`
isn't wired up: the stub deliberately types `shm` as a `_Shm` instance while
the runtime `shm` is a submodule, which stubtest would flag without an
allowlist.)

### 11. Bisulfite (`--meth`) is dual-coordinate

`--meth` (D3) alignment needs a **dual index** built by `bwa-mem3 index --meth`:
the converted, f/r-doubled seed FM-index (`<ref>.meth.*`) plus the original
un-converted reference (`<ref>.*`). Load both via `BwaIndex::load_meth(seed,
orig)` → `shim_align_idx_load_meth`; the shim keeps the original `bns`/`pac`
and a second unpacked `ref_string` resident on `BwaShimIndex`. Seeding runs
against the converted index; **everything after the seed→original remap in
`mem_kernel1_core`/`mem_kernel2_core` runs in original coordinates** — insert
size, pairing, `mem_reg2aln`, output rids, and the reported contigs
(`shim_header_bns` returns the original `bns` so the BAM header matches the
emitted rids). Per read the shim retains the unconverted bases in
`bseq1_t.meth_orig_seq` and projects `seq` in place (R1 C→T, R2 G→A) before
seeding; `mem_reg2aln` takes `meth_orig_seq` so NM/MD/CIGAR reflect the
original read, and `append_bam_record` emits Bismark `XR:Z` (read conversion,
from R1/R2), `XG:Z` (genome strand, from `mem_aln_t.meth_hypothesis`), and
`XM:Z` (via upstream `meth_build_xm`, which is compiled — only `meth_bam.cpp`,
the htslib writer, is excluded). All meth code is gated on `opt->meth_mode` /
non-NULL `meth_orig_*`, so the non-meth path is unchanged. Output matches the
CLI byte-for-byte on every record including secondaries/`XA:Z`
(`bwa-mem3-rs-cli/tests/meth_e2e.rs` pins this), because `pair_and_emit`
replicates `mem_reg2sam`'s XA folding (see gotcha #12).

### 12. `pair_and_emit` folds secondaries into `XA:Z` like `mem_reg2sam`

The shim emits records from the per-read alnreg list itself rather than calling
upstream's `mem_reg2sam`, so it must reproduce that function's output policy:
after `mem_pair_resolve` (which runs `mem_mark_primary_se`), `pair_and_emit`
calls `mem_gen_alt` to build each read's `XA:Z` string and then, per alnreg,
**skips** any region that is secondary (`ar->secondary >= 0`, folded into the
primary's `XA:Z`), below `opt->T`, or below `drop_ratio` — emitting only
primaries + supplementaries (2nd+ emitted region gets `0x800` and its MAPQ
lowered to the primary's). Without this the shim emitted every surviving
region, so multi-mapping reads got a record per hit instead of one record with
an `XA:Z` tag — harmless for unique mappers but a large divergence on
repetitive reads and on `--meth` (whose collapsed C/T scoring surfaces extra
weak hits). `lists[k]` stays 1:1 with `a[k]` so the pairing indices (`z[k]`)
and mate/SA logic are unaffected; a parallel `emit[k]` mask gates the append.

## Commit / PR conventions

- Conventional Commits; sign with `-S`; see `CONTRIBUTING.md`.
- Never mention AI/Claude/Anthropic in commits or GitHub-visible content
  (enforced globally in `~/.claude/CLAUDE.md`).
- 1Password signing flakes when the user is away from the machine. After
  three failed attempts, commit with `--no-gpg-sign` and note it in the
  commit message for re-signing later.

## Updating upstream

See `CONTRIBUTING.md` → "Updating the vendored `bwa-mem3` source."
