#!/usr/bin/env bash
# Run every script unit-test. Wired into CI as the `scripts` job so the bash
# that drives the vendor-bump automation is not the one untested thing in the
# repo.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
status=0
for t in "$here"/test-*.sh; do
    echo "=== $(basename "$t") ==="
    if ! bash "$t"; then
        status=1
    fi
done
exit "$status"
