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
# the missing `|| true` bug, and it CANNOT go through check_true/check_false:
# bash suspends errexit for the ENTIRE call stack of a command used as an
# if/while/until condition (or as either side of &&/||) — check_true's `if
# "$@"; then` included. That makes a pipefail abort inside
# assert_refreshed_tree's `other_dirty=` assignment unobservable through
# check_true regardless of whether `|| true` is present, so a version of
# this test that called `check_true assert_refreshed_tree` would pass on a
# mutated script that dropped the `|| true` just as readily as on the fixed
# one — verified by mutation, see the task report. Case (b) below instead
# runs assert_refreshed_tree as a BARE statement inside a child `bash -c`
# with errexit active, where an abort actually terminates the child and
# shows up in its exit status.
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

# (b) vendor dirty, nothing else dirty -> success, with the abort this test
# guards against actually observable. `export -f` (function) and `export`
# (variables) make the stub visible to the separate `bash -c` process below —
# a plain subshell wouldn't need this, but a freshly spawned bash does.
STUB_VENDOR_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
STUB_ALL_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
export -f git
export STUB_VENDOR_STATUS STUB_ALL_STATUS here
rc=0
bash -c '
    source "$here/../bwa-mem3-drift-report.sh"
    assert_refreshed_tree
' >/dev/null 2>&1 || rc=$?
check "vendor-only-dirty tree: bare statement under errexit does not abort" "0" "$rc"

# (c) vendor dirty, plus a non-vendor file dirty -> files-outside message.
STUB_VENDOR_STATUS=" M bwa-mem3-sys/vendor/bwa-mem3/foo.c"
STUB_ALL_STATUS=$' M bwa-mem3-sys/vendor/bwa-mem3/foo.c\n M Cargo.toml'
err="$(assert_refreshed_tree 2>&1 1>/dev/null)" || true
check_false assert_refreshed_tree
contains "extra dirty file: files-outside message" "$err" "files outside"
contains "extra dirty file: names the offending path" "$err" "Cargo.toml"

unset -f git

# --- flag_defines: extracts MEM_F_* name/value pairs ---
cat > "$tmp/flags.h" <<'EOF'
/* MEM_F_* flag bits (from bwamem.h). */
#define MEM_F_PE             0x2
#define MEM_F_NOPAIRING      0x4
#define SOMETHING_ELSE       0x8
EOF
out="$(flag_defines "$tmp/flags.h")"
check "flag_defines extracts pairs" "MEM_F_NOPAIRING 0x4
MEM_F_PE 0x2" "$out"
absent "flag_defines ignores unrelated defines" "$out" "SOMETHING_ELSE"

# --- enum_bodies: extracts full enum declarations, single- and multi-line ---
# Modeled on bwamem.h's actual two enum shapes: a single-line named `enum
# mem_meth_scoring { ... };` and a multi-line anonymous `typedef enum { ... }
# seed_order_t;`. A grep keyed on a `MEM_` prefix (an earlier draft of this
# check) would silently miss seed_order_t entirely — none of its members
# (SEED_ORDER_*) start with MEM_ — which is exactly the enum opts.rs's
# SeedOrder mirrors, so that gap would defeat the check's own stated purpose.
cat > "$tmp/mixed.h" <<'EOF'
typedef struct mem_opt_t {
    int a;
} mem_opt_t;

enum mem_meth_scoring { MEM_METH_SCORING_COLLAPSED = 0, MEM_METH_SCORING_GENOMIC = 1 };

typedef enum {
    SEED_ORDER_OFF = 0,
    SEED_ORDER_GLOBAL_LONGEST
} seed_order_t;

typedef struct {
    int low;
} mem_pestat_t;
EOF
out="$(enum_bodies "$tmp/mixed.h")"
contains "enum_bodies finds the single-line named enum" "$out" "MEM_METH_SCORING_GENOMIC"
contains "enum_bodies finds the multi-line anonymous typedef enum" "$out" "SEED_ORDER_GLOBAL_LONGEST"
contains "enum_bodies finds the multi-line enum's typedef name" "$out" "seed_order_t"
absent "enum_bodies does not leak a preceding struct's fields" "$out" "int a;"
absent "enum_bodies does not leak a following struct's fields" "$out" "int low;"

# --- enum_bodies: member ORDER is preserved, not sorted ---
# Enum members without an explicit `= N` take their value from position, so
# unlike flag_defines (independent #define macros — order-insensitive by
# construction), sorting an enum's member lines would hide a reorder that
# silently renumbers every unlabeled member after it. Two files with the same
# two members in opposite order must come back as two DIFFERENT exact strings.
cat > "$tmp/enumA.h" <<'EOF'
typedef enum {
    SEED_ORDER_OFF = 0,
    SEED_ORDER_GLOBAL_LONGEST = 1
} seed_order_t;
EOF
cat > "$tmp/enumB.h" <<'EOF'
typedef enum {
    SEED_ORDER_GLOBAL_LONGEST = 1,
    SEED_ORDER_OFF = 0
} seed_order_t;
EOF
out="$(enum_bodies "$tmp/enumA.h")"
check "enum_bodies output A keeps source order" "$(printf 'typedef enum {\n    SEED_ORDER_OFF = 0,\n    SEED_ORDER_GLOBAL_LONGEST = 1\n} seed_order_t;')" "$out"
out="$(enum_bodies "$tmp/enumB.h")"
check "enum_bodies output B keeps source order (reordered vs A)" "$(printf 'typedef enum {\n    SEED_ORDER_GLOBAL_LONGEST = 1,\n    SEED_ORDER_OFF = 0\n} seed_order_t;')" "$out"

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
