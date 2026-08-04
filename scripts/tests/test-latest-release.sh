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

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
