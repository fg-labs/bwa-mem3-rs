# Contributing to bwa-mem3-rs

## Development setup

### Prerequisites

- Rust (stable toolchain — see `rust-toolchain.toml` for the pinned version).
- A C++17 compiler (clang on macOS, gcc on Linux).
- `zlib` development headers (`zlib1g-dev` on Debian/Ubuntu; included on macOS).
- For integration tests: a prebuilt bwa-mem3 index — set `BWA_MEM3_RS_TEST_REF`
  to its prefix (the path such that `<prefix>.bwt.2bit.64` exists). Build one
  with `bwa-mem3 index <ref.fa>`.

### Install git hooks

Pre-commit hooks run `cargo ci-fmt` and `cargo ci-lint` before each commit:

```bash
./scripts/install-hooks.sh
```

The hook is skippable with `git commit --no-verify` when you know what you're
doing.

### Run checks manually

```bash
cargo ci-fmt      # format check
cargo ci-lint     # clippy with -D warnings
cargo ci-build    # build everything
cargo ci-test     # unit + integration tests

# or everything at once:
cargo ci-fmt && cargo ci-lint && cargo ci-build && cargo ci-test
```

### Run integration tests against a real index

```bash
BWA_MEM3_RS_TEST_REF=/path/to/hg38.fa cargo test --workspace
```

Without `BWA_MEM3_RS_TEST_REF` the integration tests skip gracefully.

### End-to-end regression

`bwa-mem3-rs-cli` ships a `tests/e2e.rs` test that embeds PhiX174, builds an
index via `bwa-mem3 index`, simulates 1,000 paired reads, runs our `bwa-rs`
CLI, and verifies the output via `samtools quickcheck` + `samtools view -c`.
Runs in CI when `bwa-mem3` + `samtools` are available; locally set
`BWA_MEM3_BIN=/path/to/bwa-mem3` if the CLI is not on `PATH`.

### CLI-parity test (optional)

The `cli_parity` test compares our packed-BAM output against `bwa-mem3 mem`
CLI output. It requires:

- `BWA_MEM3_BIN=/path/to/bwa-mem3`
- `BWA_MEM3_RS_TEST_REF=/path/to/index/prefix`
- `samtools` on `PATH`

## Code style

- Follow the Rust API Guidelines; prefer `impl Trait` over `Box<dyn Trait>` on
  returns; `&[T]` over `&Vec<T>` on borrows.
- `#![forbid(unsafe_code)]` would be ideal but is infeasible for an FFI crate;
  every `unsafe impl Send/Sync` and unsafe block must carry a safety comment
  naming the invariant it relies on (e.g. "BwaIndex has no mutable state after
  load" or "Seeds owns all memory; no aliasing back into shared state").
- Errors return `Result<_, bwa_mem3_rs::Error>`; no panics across the FFI
  boundary.
- Clippy runs with `-D warnings`.

## Commit conventions

- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, etc.
  per <https://conventionalcommits.org>).
- Sign commits (`-S`). Session commits made without signing should be re-signed
  via `git rebase --exec 'git commit --amend --no-edit -S'` before merge.
- Branch names: `JIRA-1234/initials_short-description` or similar.
- Never commit directly to `main`.

## Updating the vendored `bwa-mem3` source

The crate vendors a snapshot of [`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3)
at the `main` branch under `bwa-mem3-sys/vendor/bwa-mem3/`. `main` is
our integration branch tracking upstream `bwa-mem2/bwa-mem2` with:

- Apple Silicon / NEON hot-path kernels (PR #288 equivalent).
- `sse2neon` bridge + Linux/macOS aarch64 Makefile targets (PR #271
  equivalent).
- `drop-unused-global-stat` fix (the `int stat;` in `bwamem.cpp` that
  collides with libc on macOS).
- Future: bwa-meth support, XB tag, etc.

### Automated bump PRs

A nightly workflow (`.github/workflows/vendor-bump.yml`) watches
[`fg-labs/bwa-mem3`](https://github.com/fg-labs/bwa-mem3) for new releases. On
finding one it refreshes the snapshot, builds it, and opens a tracking issue
plus a **draft PR**. That PR is expected to fail CI — a bump needs shim and
`build.rs` adaptation. Its drift report is the worklist: compiler errors,
`mem_opt_t` field changes, flag and enum set changes, contract changes in the
upstream functions the shim calls, new TUs, dependency changes, and a manual
checklist of what CI cannot verify.

Normal flow: push your adaptation commits onto the bot's branch, squash into a
logical history, sign, and mark ready for review.

To bump off-schedule, or to a specific release, dispatch it manually:

```bash
gh workflow run vendor-bump.yml -f tag=v0.8.0
gh workflow run vendor-bump.yml -f tag=v0.8.0 -f dry_run=true   # report only
```

Passing `-f tag=` implies **force**: it is the deliberate bypass of the
tracking-issue dedup, so an already-open issue for that release no longer stops
the run. It does *not* bypass an existing PR or branch — if
`vendor-bump/bwa-mem3-<version>` already exists, the run fails immediately and
creates nothing, because reusing that branch would either fail as a
non-fast-forward push or require force-pushing over commits that may already
carry review comments. Close the PR and delete the branch first, then
re-dispatch. A cron tick never sets `tag`, so only a human can reach this path.

Two cases produce an issue but no PR: a refresh that fails outright (e.g.
upstream dropped `-DMATE_SORT=0`), and upstream adding an `ext/` submodule the
refresh could not know to prune — committing then would vendor the whole
subtree. For the latter, add it to `scripts/vendor-drop-subtrees.txt` and
re-dispatch. A third case produces neither: if the tracking issue itself fails
to open (e.g. a permissions gap on the bot's GitHub App), the run ends there —
no branch is pushed and no PR is opened, so a broken issue never leaves an
orphaned PR behind.

The two scripts run locally too, and need authenticated `gh` plus `jq`:

```bash
scripts/bwa-mem3-latest-release.sh              # what is new upstream?
scripts/refresh-bwa-mem3.sh <40-hex-sha>        # then, in a dirty tree:
scripts/bwa-mem3-drift-report.sh                # ...the worklist
```

To refresh the vendored snapshot to a new `main` tip:

```bash
scripts/refresh-bwa-mem3.sh <commit-hash> [local-source-path]
```

The script:

1. Clones (or rsyncs from a local tree) the target commit into `vendor/`.
2. Writes the commit hash into `vendor/COMMIT`.
3. Verifies the Makefile still has `-DMATE_SORT=0` (the shim's pairing logic
   depends on that default).

After refreshing, run the full test suite against a real index to catch any
upstream-API drift (`mem_kernel1_core`, `mem_kernel2_core`, `mem_matesw`,
etc.).

## Releases

Releases are managed by [release-plz](https://release-plz.dev). Merging to
`main` triggers a release PR that bumps versions and generates `CHANGELOG.md`
entries from conventional commits. Merging the release PR triggers the
`publish.yml` workflow, which:

1. `cargo publish`es `bwa-mem3-sys`, `bwa-mem3-rs`, and `bwa-mem3-rs-cli` to
   crates.io via Trusted Publishing.
2. Creates a GitHub release tagged `v<VERSION>`.
3. The GitHub release fires `pypi.yml`, which builds the `bwa-mem3` Python
   wheels (Linux x86_64 / aarch64, macOS arm64) and uploads them to PyPI via
   Trusted Publishing.

`bwa-mem3-py` is a workspace member that inherits `[workspace.package].version`,
so release-plz bumps it with the rest and no manual step is needed at release
time. It is deliberately excluded from `default-members` and from the `ci-*`
aliases in `.cargo/config.toml`, so the pure-Rust checks never need a Python
interpreter; maturin builds it from `bwa-mem3-py/` in the Python workflow. Its
dependency on `bwa-mem3-rs` is path-only for the same reason -- a `version`
requirement there is used only when publishing to crates.io, and the crate is
`publish = false`. Because it is now in the workspace, the publish workflow's
"Verify lockstep workspace versions" gate covers it too, which it could not do
while the crate sat outside.

### PyPI Trusted Publishing (one-time setup)

Before the first PyPI release, the project maintainer must register a
[pending publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
on PyPI for the `bwa-mem3` project:

- PyPI Project Name: `bwa-mem3`
- Owner: `fg-labs`
- Repository name: `bwa-mem3-rs`
- Workflow filename: `pypi.yml`
- Environment name: `pypi`

The `pypi` GitHub Environment should also exist on the repo (`Settings →
Environments → New environment`) so the `environment:` clause in `pypi.yml`
resolves; protection rules (required reviewers, branch policy) are optional
but recommended for a project that auto-publishes on release.

No `PYPI_API_TOKEN` secret is needed — the workflow exchanges its OIDC token
for a short-lived upload credential at run time.
