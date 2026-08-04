#!/usr/bin/env bash
# Unit-tests the pure helpers in bwa-mem3-drift-report.sh. The git- and
# cargo-dependent checks are covered by the replay harness (Task 13).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "$here/../bwa-mem3-drift-report.sh"

pass=0
fail=0
check() { if [ "$2" = "$3" ]; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL: $1 — expected '$2', got '$3'" >&2; fi; }
contains() { if printf '%s' "$2" | grep -q -- "$3"; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL: $1 — '$3' not found in output" >&2; fi; }
absent() { if printf '%s' "$2" | grep -q -- "$3"; then fail=$((fail+1)); echo "FAIL: $1 — '$3' unexpectedly present" >&2; else pass=$((pass+1)); fi; }
check_true()  { if "$@" >/dev/null 2>&1; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL: expected success: $*" >&2; fi; }
check_false() { if "$@" >/dev/null 2>&1; then fail=$((fail+1)); echo "FAIL: expected failure: $*" >&2; else pass=$((pass+1)); fi; }

tmp="$(mktemp -d)"
# Both directories, and the same three signals the sourced script uses: bash
# keeps only ONE handler per signal, so a bare `trap ... EXIT` here would
# REPLACE the `trap 'rm -rf "$TMPDIR_SELF"' EXIT INT TERM` that
# bwa-mem3-drift-report.sh installs when sourced above. Every helper this file
# exercises (check_gitmodules, check_call_site_contracts,
# check_source_inventory, check_new_opts) allocates through mktemp_tracked
# inside $TMPDIR_SELF, so replacing it leaked one temp directory per run and
# dropped INT/TERM handling with it.
trap 'rm -rf "$tmp" "$TMPDIR_SELF"' EXIT INT TERM

# --- gitmodules_new_submodules: names only what the new file adds ---
printf '[submodule "ext/a"]\n\tpath = ext/a\n' > "$tmp/old"
printf '[submodule "ext/a"]\n\tpath = ext/a\n[submodule "ext/b"]\n\tpath = ext/b\n' > "$tmp/new"
out="$(gitmodules_new_submodules "$tmp/old" "$tmp/new")"
check "one new submodule detected" "ext/b" "$out"

# --- unchanged .gitmodules yields nothing ---
out="$(gitmodules_new_submodules "$tmp/old" "$tmp/old")"
check "no drift on identical files" "" "$out"

# --- a REMOVED submodule is not reported as new ---
out="$(gitmodules_new_submodules "$tmp/new" "$tmp/old")"
check "removal is not a new submodule" "" "$out"

# --- TWO new submodules at once: the exact case the tr-newline bug corrupted ---
# (that bug joined multi-line grep output into one token, e.g. "ext/cext/d",
# via `tr -d '[:space:]'` stripping the newline between them — a single-new
# fixture can't catch that because there is no second line to glue on to.)
printf '[submodule "ext/a"]\n\tpath = ext/a\n' > "$tmp/old2"
printf '[submodule "ext/a"]\n\tpath = ext/a\n[submodule "ext/c"]\n\tpath = ext/c\n[submodule "ext/d"]\n\tpath = ext/d\n' > "$tmp/new2"
out="$(gitmodules_new_submodules "$tmp/old2" "$tmp/new2")"
check "two new submodules detected, not glued into one token" "$(printf 'ext/c\next/d')" "$out"

# --- struct_body: anonymous typedef must extract (mem_pestat_t's shape) ---
cat > "$tmp/hdr.h" <<'EOF'
typedef struct mem_opt_t {
    int a, b;
    int w;
} mem_opt_t;

typedef struct {
    int low, high;
    double avg, std;
} mem_pestat_t;
EOF
out="$(struct_body "$tmp/hdr.h" mem_pestat_t)"
contains "anonymous typedef extracted" "$out" "low, high"
absent "anonymous extraction did not bleed into mem_opt_t" "$out" "int w"
out="$(struct_body "$tmp/hdr.h" mem_opt_t)"
contains "named typedef extracted" "$out" "int w"
absent "named extraction did not bleed into mem_pestat_t" "$out" "avg, std"

# --- assert_refreshed_tree: all three branches, via a stubbed `git` ---
# Overriding `git` as a shell function (rather than touching real repo state)
# keeps this offline and independent of whatever is actually dirty in the
# working tree the test happens to run in. The first call inside
# assert_refreshed_tree passes the literal pathspec "bwa-mem3-sys/vendor/";
# the second does not — that's how the stub tells the two calls apart.
#
# Case (b) — vendor dirty, nothing else dirty — is the regression test for
# the missing `|| true` bug: `git status --porcelain | grep -v
# 'bwa-mem3-sys/vendor/'` exits 1 when grep filters out every line (the
# success case), which under `pipefail`/`set -e` aborted the whole script
# before this branch's own `[ -n "$other_dirty" ]` check ever ran. Removing
# the `|| true` at the `other_dirty=` assignment reproduces that: this test
# then fails by never completing (the sourced script's `set -e` kills the
# test process), rather than by a wrong exit code — confirmed by hand while
# fixing the bug, not asserted here since a killed process can't run further
# assertions in the same test file.
# Invoked indirectly by assert_refreshed_tree, which lives in the sourced
# script — invisible to shellcheck here because of the `source=/dev/null`
# directive above.
# shellcheck disable=SC2329
git() {
    local arg
    for arg in "$@"; do
        if [ "$arg" = "bwa-mem3-sys/vendor/" ]; then
            printf '%s\n' "$STUB_VENDOR_STATUS"
            return 0
        fi
    done
    printf '%s\n' "$STUB_ALL_STATUS"
    return 0
}

# (a) vendor clean -> nothing was refreshed.
STUB_VENDOR_STATUS=""
STUB_ALL_STATUS=""
err="$(assert_refreshed_tree 2>&1 1>/dev/null)" || true
check_false assert_refreshed_tree
contains "clean tree: nothing-refreshed message" "$err" "nothing was refreshed"

# (b) vendor dirty, nothing else dirty -> success. Regression test for the
# missing `|| true`.
STUB_VENDOR_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
STUB_ALL_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
check_true assert_refreshed_tree

# (c) vendor dirty, plus a non-vendor file dirty -> files-outside message.
STUB_VENDOR_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
STUB_ALL_STATUS=$' M bwa-mem3-sys/vendor/bwa-mem3/foo.c\n M Cargo.toml'
err="$(assert_refreshed_tree 2>&1 1>/dev/null)" || true
check_false assert_refreshed_tree
contains "extra dirty file: files-outside message" "$err" "files outside"
contains "extra dirty file: names the offending path" "$err" "Cargo.toml"

unset -f git

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
