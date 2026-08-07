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

They're in `shim/bwa_shim_types.h` (what bindgen reads) and in upstream's `bwamem.h` (what `bwa_shim_align.cpp` includes). Both must stay byte-identical — the shim allocates the real struct via upstream's `mem_opt_init()`, so any field-order/size drift silently corrupts the offsets Rust reads/writes through the bindgen view. **`shim/bwa_shim_layout_assert.cpp` guards this**: it `#include`s the POD under renamed tags (`#define mem_opt_t pod_mem_opt_t`) and then the real `bwamem.h`, and `static_assert`s `offsetof`/`sizeof` field by field, so drift fails `cargo build`. (The bindgen `bindgen_test_layout_mem_opt_t` test only checks the Rust struct against the C compiler's view of the *POD copy*, **not** the POD against upstream — that gap is what the guard TU closes.) The same TU also `static_assert`s the `MEM_F_*` flag *values* the POD hardcodes against upstream's, since bindgen's `allowlist_var("MEM_F_.*")` reads only the POD and would never notice a renumbered flag on its own. On `refresh-bwa-mem3.sh`, diff `vendor/bwa-mem3/src/bwamem.h` around lines 182–295 (`mem_opt_t`) and 297–301 (`mem_pestat_t`) — the ranges move on most refreshes, so locate them by the `} mem_opt_t;` / `} mem_pestat_t;` closers rather than trusting these numbers; if either changed, update `shim/bwa_shim_types.h` **and** the field list in `bwa_shim_layout_assert.cpp` to match.

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

The htslib- and libsais-dependent TUs would also fail to link:
`scripts/vendor-drop-subtrees.txt` (the shared pruning list read by both
`refresh-bwa-mem3.sh` and `bwa-mem3-drift-report.sh`) prunes `ext/htslib`,
`ext/libsais`, `ext/mimalloc`, `ext/doctest`, and (as of 0.6.0) `ext/zlib-ng`
from the vendor tree — along with upstream's own `test/` suite, which is
vendored-in-principle but never compiled. If a future refresh wants to
re-enable any of them, restore the relevant submodule first (or edit the
shared list).

`fastmap.cpp` was compiled at 0.2.x purely to export `worker_alloc` /
`worker_free`. As of bwa-mem3 0.6.0 it is transitively coupled to the new
`fast_reader` FASTQ path (libdeflate + zlib-ng), `numa`, and htslib — none of
which this crate vendors — so it is excluded again. We instead carry verbatim
file-local copies of just `worker_alloc` / `worker_free` in
`shim/bwa_shim_align.cpp`; **keep them in sync** with upstream's `fastmap.cpp`
on every vendor refresh (the buffer set must match what `mem_kernel1_core` /
`mem_kernel2_core` expect).

**They are no longer verbatim.** As of v0.9.0 upstream sizes `chain_scratch` /
`seed_scratch` (renamed from `chain_ar` / `seedBuf`) as `nthreads * BATCH_SIZE`
and indexes them by `tid`, because it fused its two `kt_for` passes into
`worker_bwt_aln` so chains no longer outlive a work item. This crate's public
API **is** that barrier — `seed_batch` → `extend_batch` — so our copy keeps the
pre-0.9.0 `nreads` sizing and `seq_id` indexing. The kernels take both arrays as
parameters and index them `[0, nseq)`, and `tid` selects only `mmc`, so a
caller-owned per-read array stays valid. Upstream renamed the fields precisely
so out-of-tree callers would fail to compile here rather than silently run off
the end; do **not** take the compiler's `auxSeedBuf` did-you-mean, which is a
vestigial member nothing allocates.

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
seeding. It must also set **`bseq1_t.meth_base_ot`** (R1 = OT = 1, R2 = OB = 0):
that field feeds the seed-chemistry filter in `meth_seed_to_orig`, which drops
any seed whose genomic strand disagrees with the read's. Leaving it zero is not
inert — `0` is the *valid* encoding for OB, so the filter runs and silently
discards every R1 seed (it only self-disables at `< 0`). `mem_reg2aln` and
`mem_gen_alt` both take `meth_orig_seq` so NM/MD/CIGAR — and every `XA:Z`
sub-entry's NM — reflect the original read rather than the projected one, and `append_bam_record` emits Bismark `XR:Z` (read conversion,
from R1/R2), `XG:Z` (genome strand, from `mem_aln_t.meth_hypothesis`), and
`XM:Z` (via upstream `meth_build_xm`, which is compiled — only `meth_bam.cpp`,
the htslib writer, is excluded; it takes a `meth_chem_t` chemistry argument as
of v0.9.0, sourced from `opt->meth_chem`).

`HN:i` **is** emitted under `--meth` as of v0.9.0. It previously was not, because
upstream's meth writer never emitted it; v0.9.0 put both writers behind a
compat-target switch (`p.HN >= 0 && opt->compat->emit_hn`, `meth_bam.cpp:663` /
`bam_writer.cpp:458`) and `COMPAT_TARGET_OFF` — the default native-output
target — sets `emit_hn = 1`. The shim gates on the same expression. Note the
same struct carries `emit_mq = 1`, which bears on the `MQ:i` tag this crate
deliberately never emits (`DELIBERATELY_ASYMMETRIC_TAG_KEYS`). All meth code is gated on `opt->meth_mode` /
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

One paired-branch subtlety the unified emission must reproduce: when
`mem_pair` selects a non-top region (`z[k] != 0`), `mem_pair_resolve` promotes
`a[k].a[z[k]]` (sets its `secondary` to `-2`) and runs the `secondary_all`
switch, which reassigns the old SE-primary (region 0) into `z[k]`'s group —
leaving it with `secondary < 0` but `secondary_all >= 0`. `mem_gen_alt` folds
that region into `z[k]`'s `XA:Z`, and upstream `mem_sam_pe`'s paired block emits
**only** `z[k]` as primary, never the switched-away region. Because our emit
filter keys off `secondary` alone, the filter also drops any region with
`secondary < 0 && secondary_all >= 0` on the paired branch; without it the shim
would surface the old primary as an extra record and demote the `z[k]`
pair-primary to a `0x800` supplementary. `cli_parity_pair_select.rs` pins this
(two near-identical motif copies + a mate that anchors R1 to the lower-scoring
copy, so `z[0] != 0`).

### 13. Upstream's defaulted parameters are a silent-breakage surface

Several `bwamem.h` functions the shim calls take C++ **default arguments**, and
upstream keeps adding them:

| function | defaulted param(s) |
|---|---|
| `mem_gen_alt` | `meth_orig_query` |
| `mem_reg2aln` | `meth_orig_query` |
| `mem_kernel1_core` / `mem_kernel2_core` | `meth_orig_bns`, `meth_orig_pac` |
| `mem_matesw` / `mem_matesw_batch_post` | `ms_orig`, `mat`, `mate_meth_ot` |

A new one **compiles clean at every existing call site** and silently selects
the legacy path. The same hazard applies to new `bseq1_t` / `mem_opt_t` fields
the caller is expected to populate: they come up zero from a `calloc`, and zero
is often a meaningful value rather than an "unset" sentinel.

The v0.9.0 bump hit this three times — `bseq1_t::meth_base_ot` (unset ⇒ every R1
seed discarded), `mem_gen_alt`'s `meth_orig_query` (unset ⇒ every `XA:Z`
sub-entry's NM counted bisulfite conversions), and the `HN` compat switch. None
produced a warning, and the non-meth suite stayed green throughout.

**On every refresh**, grep `bwamem.h` for `= NULL)` / `= -1)` / `= NULL,` /
`= -1,` and check each against the shim's call sites, then diff `bwa.h`'s
`bseq1_t` for new fields. Where upstream exposes a helper that encapsulates the
policy (`mem_proper_pair_extra_flag`, `mem_opt_apply_meth_defaults`,
`mem_aln_ref_string`), prefer calling it over replicating it — a replicated
policy is a twin that drifts, and upstream factors these out precisely because
its own twins drifted.

### 14. `XA:Z` sub-entry ORDER can differ from the CLI on tied hits

Byte-parity against the CLI holds for the `XA:Z` *set*, not always its order. On
a 2,000-pair hg38 fixture: 26 of 27 `XA`-bearing records byte-identical, 1 with
the same five entries in a different order, 0 with a different set. The two
swapped entries had equal NM.

This is upstream's, not ours. With the default `alnreg_sort_fast == 0` the
surviving alnreg set "is defined by the *permutation* klib's unstable introsort
happens to produce" (`bwamem.cpp:560-575`), and upstream measured that sort
tying on ~0.98% of calls over 186M regions. When the comparator sees a tie the
resulting order is not determined by the input alone, so two processes linking
the same `mem_gen_alt` can legitimately disagree.

Practical consequences:

- A parity assertion over `XA:Z` must compare the entry **set**, not the string,
  unless the fixture is known tie-free. The PhiX-scale suites happen to be, which
  is why they compare strings and pass.
- Do not "fix" this in the shim. The regions, coordinates, CIGARs and NMs all
  match; only the order of equal-NM entries moves.
- `opt->alnreg_sort_fast` (v0.9.0, `--fast`) swaps in a strict-total-order
  comparator plus pdqsort, which makes the order deterministic — but it is a
  different, non-bwa-mem2-compatible surviving set, so it is not a parity fix.

## Commit / PR conventions

- Conventional Commits; sign with `-S`; see `CONTRIBUTING.md`.
- Never mention AI/Claude/Anthropic in commits or GitHub-visible content
  (enforced globally in `~/.claude/CLAUDE.md`).
- 1Password signing flakes when the user is away from the machine. After
  three failed attempts, commit with `--no-gpg-sign` and note it in the
  commit message for re-signing later.

## Updating upstream

See `CONTRIBUTING.md` → "Updating the vendored `bwa-mem3` source."
