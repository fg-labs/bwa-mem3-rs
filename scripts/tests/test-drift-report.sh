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

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
