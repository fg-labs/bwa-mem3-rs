#!/usr/bin/env bash
# Replay the v0.2.2 -> v0.6.0 vendor bump and assert the drift report surfaces
# every finding that bump actually required. This is the acceptance test for
# scripts/bwa-mem3-drift-report.sh: the correct answer is already in git.
#
# Slow (full upstream clone + cargo ci-build) and needs network + gh, so it is
# deliberately NOT part of the `scripts` CI job. Run it when you change the
# report.
#
# Usage: bash scripts/tests/drift-report-replay.sh
set -euo pipefail

# This repo's commit immediately BEFORE the v0.6.0 bump.
PRE_BUMP_REF="69f81e6"
# Upstream's sha for v0.6.0. NOT this repo's bump commit — the report and the
# refresh script both query fg-labs/bwa-mem3, where our shas do not exist.
UPSTREAM_V060="a887e36cb8fbdc54bd5a3543cdfd1850bf7e2f52"

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"

# A shallow clone would not have $PRE_BUMP_REF's history available to
# `git worktree add`; fill it in rather than failing partway through.
if ! git -C "$repo_root" rev-parse --quiet --verify "${PRE_BUMP_REF}^{commit}" >/dev/null; then
    echo "=== $PRE_BUMP_REF not found locally; unshallowing ===" >&2
    git -C "$repo_root" fetch --unshallow
fi

work="$(mktemp -d)"
# Single trap covering both "worktree was never created" (the `|| true` makes
# the no-op removal harmless) and "worktree exists and must go" — two
# separate traps that overwrite each other (the pattern the brief used) is
# one too many moving parts for something that runs unconditionally on every
# exit path, success or failure.
#
# INT and TERM are trapped alongside EXIT because bash runs an EXIT trap on
# death by signal only when that signal is trapped too. This script drives the
# report through `cargo ci-build` twice, so there is a multi-minute window in
# which Ctrl-C is likely — and an untrapped interrupt would skip `git worktree
# remove`, leaving a worktree registered in the developer's real repository at
# a path under /tmp that later `git worktree list` calls keep reporting.
# bwa-mem3-drift-report.sh traps the same three signals for the same reason.
cleanup() {
    git -C "$repo_root" worktree remove --force "$work/tree" >/dev/null 2>&1 || true
    rm -rf "$work"
}
trap cleanup EXIT INT TERM

echo "=== Creating a worktree at $PRE_BUMP_REF ==="
git -C "$repo_root" worktree add --detach "$work/tree" "$PRE_BUMP_REF" >/dev/null

# Copy the CURRENT scripts in. The worktree predates them, and the report must
# be the version under test — but note it will read the worktree's build.rs,
# where fastmap.cpp is not yet in skip_common. That is deliberate: it is what
# the report would have seen at the time.
cp "$repo_root/scripts/bwa-mem3-drift-report.sh" \
   "$repo_root/scripts/refresh-bwa-mem3.sh" \
   "$repo_root/scripts/vendor-drop-subtrees.txt" \
   "$work/tree/scripts/"

# The pre-bump prune list must not include entries the bump introduced, or the
# replay cannot rediscover them. ext/zlib-ng is the finding we are testing for.
grep -v 'ext/zlib-ng' "$repo_root/scripts/vendor-drop-subtrees.txt" \
    > "$work/tree/scripts/vendor-drop-subtrees.txt"

# Commit the substituted scripts into the (throwaway, detached-HEAD, never
# pushed) worktree. Without this, `bwa-mem3-drift-report.sh`'s own
# `assert_refreshed_tree` guard sees `scripts/refresh-bwa-mem3.sh` as modified
# and the other two as untracked — both outside `bwa-mem3-sys/vendor/` — and
# refuses to run at all ("files outside bwa-mem3-sys/vendor/ are dirty"),
# before the vendor refresh below even happens. Committing first means the
# refresh is the only dirty thing left for that guard to see, matching what a
# real refresh-on-a-clean-checkout looks like. This commit never leaves the
# worktree and is unsigned on purpose: it is not repo history, it is a fixture
# for this replay, and the worktree (and this commit with it) is deleted on
# every exit path via the trap above.
git -C "$work/tree" add scripts/bwa-mem3-drift-report.sh scripts/refresh-bwa-mem3.sh scripts/vendor-drop-subtrees.txt
git -C "$work/tree" \
    -c user.name="drift-report-replay" -c user.email="drift-report-replay@localhost" \
    commit --no-gpg-sign -q -m "test fixture: substitute current drift-report scripts"

echo "=== Refreshing to upstream v0.6.0 ($UPSTREAM_V060) ==="
(cd "$work/tree" && scripts/refresh-bwa-mem3.sh "$UPSTREAM_V060")

echo "=== Running the drift report ==="
report="$work/report.md"
# The status is captured rather than discarded: the report exits 0 whether or
# not it finds drift, so a non-zero status is always misuse (a rejected tree,
# a missing dependency) rather than a finding. Without this, such a run shows
# up only as every assertion below failing at once, which reads like eight
# regressions instead of one harness problem.
#
# It is fatal, not a warning: a report that died partway can still have emitted
# every pattern the assertions look for before it stopped, so continuing would
# let a broken run finish with "9 passed, 0 failed" — precisely the false pass
# an acceptance test exists to prevent.
report_status=0
(cd "$work/tree" && MISSED_TAGS="v0.3.0 v0.4.0 v0.5.0 v0.6.0" \
    scripts/bwa-mem3-drift-report.sh) > "$report" 2>&1 || report_status=$?
if [ "$report_status" -ne 0 ]; then
    echo "ERROR: the drift report exited $report_status — that is misuse, not" >&2
    echo "       drift, so nothing below would be meaningful. Last 20 lines:" >&2
    tail -20 "$report" >&2
    exit 1
fi

echo "=== Asserting expected findings ==="
pass=0
fail=0
# -E (ERE) so alternation is plain `|`. The BRE spelling `\|` is a GNU
# extension, and this script is documented as a by-hand developer tool rather
# than a CI job, so it has to hold up on a grep that does not implement it.
expect() { # expect <description> <pattern>
    if grep -qiE -- "$2" "$report"; then
        pass=$((pass + 1))
        echo "  ok: $1"
    else
        fail=$((fail + 1))
        echo "  FAIL: $1 (no match for '$2')" >&2
    fi
}

# Every pattern below except the last two must match text the report emits ONLY
# when it has the corresponding finding. That is easy to get wrong: `expect`
# greps the whole report unanchored, and the report also prints the gate's
# status line and a manual checklist on EVERY run, drift or not — so bare
# `NEW_SUBMODULE` and bare `worker_alloc` both matched boilerplate and could
# not fail. They are anchored on the finding's own form instead: the gate's
# bare marker line, and check 5's `### <fn> changed|NOT FOUND` heading.
expect "flags the new ext/zlib-ng submodule"      "zlib-ng"
expect "emits the hard-gate marker"               '^NEW_SUBMODULE$'
expect "reports mem_opt_t layout drift"           "mem_opt_t\` changed"
expect "names a new mem_opt_t field"               "smem_dedup|min_ext_len|band_start"
expect "reports the new .c reader TUs"            "fast_reader"
# mem_kernel1_core, NOT worker_alloc: this assertion exists to prove check 5
# emits a real contract finding, and worker_alloc cannot serve that purpose
# here because its body is byte-identical across this replay's two endpoints
# (69f81e6 -> a887e36) -- the report is correct to stay silent about it. The
# assertion only ever passed because an unanchored, case-insensitive grep for
# `worker_alloc` also matched the manual checklist, which is printed on every
# run. mem_kernel1_core is the shim's main entry point and did change.
# shellcheck disable=SC2016 # backticks are literal markdown in the report,
# not command substitution -- single quotes are what keeps them that way.
expect "reports an upstream contract change"      '^### `mem_kernel1_core` (changed|NOT FOUND)'
# Unconditional by design: this one asserts the checklist section is present at
# all, not that a particular finding fired.
expect "includes the manual checklist"            "Manual checklist"
expect "includes release notes for a missed tag"  "v0.5.0"

echo "=== Synthetic prune: removing a needed header must fail the build check ==="
rm -rf "$work/tree/bwa-mem3-sys/vendor/bwa-mem3/ext/pdqsort"
synth="$work/synthetic.md"
(cd "$work/tree" && scripts/bwa-mem3-drift-report.sh) > "$synth" 2>&1 || true
if grep -q "Build FAILED" "$synth"; then
    pass=$((pass + 1)); echo "  ok: pruning ext/pdqsort fails the build check"
else
    fail=$((fail + 1)); echo "  FAIL: build check did not notice a missing needed header" >&2
fi

echo
echo "$pass passed, $fail failed"
echo "Full report: $report (removed on exit; copy it out if you need it)"
[ "$fail" -eq 0 ]
