#!/usr/bin/env bash
# Unit-tests the pure logic in bwa-mem3-latest-release.sh. Sourcing the script
# must not run main(), which is what the BASH_SOURCE guard at its foot is for.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "$here/../bwa-mem3-latest-release.sh"

pass=0
fail=0
check() { # check <description> <expected> <actual>
    if [ "$2" = "$3" ]; then
        pass=$((pass + 1))
    else
        fail=$((fail + 1))
        echo "FAIL: $1 — expected '$2', got '$3'" >&2
    fi
}
check_true()  { if "$@" >/dev/null 2>&1; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL: expected success: $*" >&2; fi; }
check_false() { if "$@" >/dev/null 2>&1; then fail=$((fail+1)); echo "FAIL: expected failure: $*" >&2; else pass=$((pass+1)); fi; }

# --- version_gt: numeric field-by-field, not lexical, not sort -V ---
check_true  version_gt 0.8.0 0.6.0
check_false version_gt 0.6.0 0.8.0
check_false version_gt 0.6.0 0.6.0          # equal is not greater
check_true  version_gt 0.10.0 0.9.0         # lexical compare would fail this
check_true  version_gt 0.3.0 0.2.2
check_true  version_gt 1.0.0 0.99.99
check_true  version_gt 0.6.1 0.6.0
# Tolerate a leading v on either side.
check_true  version_gt v0.8.0 v0.6.0
check_false version_gt v0.6.0 v0.8.0

# --- REPO_ROOT / COMMIT_FILE: regression test for the $0-vs-BASH_SOURCE bug ---
# $0 is not rewritten by `source`, so a REPO_ROOT computed from $0 resolves
# relative to *this* test file's directory (scripts/tests/) rather than the
# sourced script's own directory (scripts/) — one level too shallow, landing
# COMMIT_FILE on a path that doesn't exist. This must hold on every sourced
# invocation, not just when the script is executed directly.
check_true test -f "$COMMIT_FILE"

# --- vendored_version: offline, via a stubbed resolve_tag_sha ---
# resolve_tag_sha is redefined below purely for this test process (sourcing
# makes it a plain shell function, so redefining it here does not touch the
# real script on disk) so these assertions never depend on GitHub being
# reachable or on upstream's current release list.
fixture_commit_file="$(mktemp)"
# EXIT INT TERM, like every other script in this directory: bash runs an EXIT
# trap on death by signal only when that signal is trapped too. (Nothing is
# clobbered here -- bwa-mem3-latest-release.sh installs no trap of its own,
# unlike bwa-mem3-drift-report.sh; see the note in test-drift-report.sh.)
trap 'rm -f "$fixture_commit_file"' EXIT INT TERM
COMMIT_FILE="$fixture_commit_file"
fixture_tags=$'v0.3.0\nv0.2.0\nv0.1.0'

# Invoked indirectly by vendored_version, which lives in the sourced script —
# invisible to shellcheck here because of the `source=/dev/null` directive above.
# shellcheck disable=SC2329
resolve_tag_sha() {
    case "$1" in
        v0.3.0) printf '%s' "3333333333333333333333333333333333333333" ;;
        v0.2.0) printf '%s' "2222222222222222222222222222222222222222" ;;
        v0.1.0) printf '%s' "1111111111111111111111111111111111111111" ;;
        *) return 1 ;;
    esac
}

printf '%s' "2222222222222222222222222222222222222222" > "$fixture_commit_file"
check "vendored_version: finds a mid-list match" "v0.2.0" "$(vendored_version "$fixture_tags")"

# Genuine no-match: every fixture tag resolves successfully, none matches ->
# exit 1 (retag/absence — NOT the same as an API failure, exit 2, below).
printf '%s' "deadbeefdeadbeefdeadbeefdeadbeefdeadbeef" > "$fixture_commit_file"
rc=0
out="$(vendored_version "$fixture_tags" 2>/dev/null)" || rc=$?
check "vendored_version: no match -> empty output" "" "$out"
check "vendored_version: no match -> exit 1" "1" "$rc"

# --- Regression test for the API-failure-misreported-as-retag bug ---
# Stub resolve_tag_sha to fail for every tag, as a rate-limited or
# network-blipped `gh api` call would. Before the fix, a failed lookup for
# the matching tag was indistinguishable from "this tag doesn't match," so
# the loop ran to completion and returned the same exit code as a genuine
# no-match (1) — which main() reports as a possible upstream retag. The fix
# must return a distinct code (2) so main() can report an API failure and
# say "retry" instead of accusing upstream of retagging.
resolve_tag_sha() { return 1; }
printf '%s' "2222222222222222222222222222222222222222" > "$fixture_commit_file"
rc=0
out="$(vendored_version "$fixture_tags" 2>/dev/null)" || rc=$?
check "vendored_version: resolve failure -> empty output" "" "$out"
check "vendored_version: resolve failure -> exit 2, not 1 (not a retag)" "2" "$rc"

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
