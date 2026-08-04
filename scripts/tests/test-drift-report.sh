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
# -F (fixed-string): callers pass literal markdown snippets (backticks,
# `**bold**`, `[i]`-shaped code fragments), and BRE would otherwise
# misinterpret `[...]` as a bracket expression or a leading `**` as an
# invalid repetition operator with no preceding element.
contains() { if printf '%s' "$2" | grep -qF -- "$3"; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL: $1 — '$3' not found in output" >&2; fi; }
absent() { if printf '%s' "$2" | grep -qF -- "$3"; then fail=$((fail+1)); echo "FAIL: $1 — '$3' unexpectedly present" >&2; else pass=$((pass+1)); fi; }
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

# --- function_body: single-line signature, straight-line extraction ---
cat > "$tmp/fb1.cpp" <<'EOF'
void unrelated_before() {
    int x = 1;
}

void worker_alloc(int a, int b)
{
    int c = a + b;
    return;
}

void unrelated_after() {
    int y = 2;
}
EOF
out="$(function_body "$tmp/fb1.cpp" worker_alloc)"
contains "function_body: extracts the signature line" "$out" "void worker_alloc(int a, int b)"
contains "function_body: extracts the body" "$out" "int c = a + b;"
absent "function_body: does not bleed into the following function" "$out" "int y = 2;"
absent "function_body: does not bleed in the preceding function" "$out" "int x = 1;"

# --- function_body: signature spanning two lines (worker_alloc's real shape) ---
cat > "$tmp/fb2.cpp" <<'EOF'
void worker_alloc(const mem_opt_t *opt, worker_t &w, int32_t nreads,
                   int32_t nthreads)
{
    w.nthreads = nthreads;
}
EOF
out="$(function_body "$tmp/fb2.cpp" worker_alloc)"
contains "function_body: multi-line signature captured" "$out" "int32_t nthreads)"
contains "function_body: multi-line signature body captured" "$out" "w.nthreads = nthreads;"

# --- function_body: a forward declaration (ends in ';') must not trigger ---
cat > "$tmp/fb3.cpp" <<'EOF'
void mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns);

void mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns)
{
    int real_body_marker = 1;
}
EOF
out="$(function_body "$tmp/fb3.cpp" mem_sam_pe)"
contains "function_body: skips a semicolon-terminated forward decl, finds the real definition" "$out" "real_body_marker"

# --- function_body: name absent entirely -> empty output (the NOT-FOUND signal check 5 relies on) ---
out="$(function_body "$tmp/fb1.cpp" totally_absent_name)"
check "function_body: absent name yields empty output" "" "$out"

# --- function_body: a macro invocation of the same name at column 0 must not trigger ---
cat > "$tmp/fb4.cpp" <<'EOF'
#define WRAP_IT(x) worker_alloc(x)

void worker_alloc(int a)
{
    int real_body_marker = 1;
}
EOF
out="$(function_body "$tmp/fb4.cpp" worker_alloc)"
contains "function_body: does not trigger on a #define line, finds the real definition" "$out" "real_body_marker"
absent "function_body: does not capture the macro line itself" "$out" "WRAP_IT"

# --- function_body: REGRESSION — an indented, multi-line CALL SITE that
# textually precedes the real column-0 definition must not hijack the scan.
# This is the exact shape found in bwa-mem3's bwamem.cpp for mem_reg2sam and
# mem_mark_primary_se: both are called, indented, inside an earlier function,
# with the call's first physical line ending in a comma (not a semicolon) —
# which satisfies "contains name( and isn't semicolon-terminated" just as
# well as a real definition unless the scan also requires column 0.
cat > "$tmp/fb5.cpp" <<'EOF'
void some_earlier_function(worker_t *w)
{
    for (int i = 0; i < 10; i++)
    {
        mem_reg2sam(w->opt, mem_aln_bns(w), mem_aln_pac(w), &w->seqs[i],
                    &w->regs[i], 0, 0);
        free(w->regs[i].a);
    }
}

void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                 bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
    int real_definition_marker = 1;
}
EOF
out="$(function_body "$tmp/fb5.cpp" mem_reg2sam)"
contains "function_body: call-before-def fixture finds the real definition" "$out" "real_definition_marker"
absent "function_body: call-before-def fixture does not capture the call site's own scope" "$out" "free(w->regs[i].a)"

# Mutation check for the above: reproduce the brief's ORIGINAL function_body
# (no column-0 requirement on the trigger) and confirm it fails exactly this
# assertion — i.e. it captures the call site's enclosing scope instead of the
# definition. This is the bug found while implementing check 5 against the
# real vendored tree (mem_reg2sam and mem_mark_primary_se both reproduce it
# there today); the fixture above is the minimal repro.
function_body_unanchored_buggy() {
    local file="$1" name="$2"
    mawk -v name="$name" '
        index($0, name "(") && !started && $0 !~ /;[[:space:]]*$/ { started = 1 }
        started { print }
        started && /^\}/ { exit }
    ' "$file"
}
out="$(function_body_unanchored_buggy "$tmp/fb5.cpp" mem_reg2sam)"
absent "MUTATION (expected red): unanchored function_body wrongly captures the call site, missing the real definition marker" "$out" "real_definition_marker"
contains "MUTATION (expected red): unanchored function_body's wrong capture includes the call site's own scope" "$out" "free(w->regs[i].a)"

# --------------------------------------------------------------------------
# Fixture-repo integration tests for checks 5, 6, and 8.
#
# These are the checks Task 12's brief flagged as able to silently do
# nothing: check 5 (call-site contracts) depends on UPSTREAM_CONTRACTS
# pointing at the right file, check 6 (source inventory) depends on a regex
# that must match every skip_common entry, and check 8 (new mem_opt_t
# fields) has a bare `diff | grep | sed` command substitution that aborts
# the whole script under `set -o pipefail` in the common case (diff exits 1
# because the two struct bodies actually differ).
#
# A local `git init` is not a clone and touches no network, so this stays
# within "offline and fast." REPO_ROOT/VENDOR/VENDOR_ABS are reassigned for
# the duration of this block (they are ordinary globals in the sourced
# script, not readonly) and restored afterward.
# --------------------------------------------------------------------------
orig_repo_root="$REPO_ROOT"
orig_vendor="$VENDOR"
orig_vendor_abs="$VENDOR_ABS"

fixture="$tmp/fixture"
# Exported immediately: several tests below spawn a child `bash -c` that
# reads $fixture from its environment (the script it sources runs under
# `set -u`, so an unexported reference there is an "unbound variable" abort,
# not a silent empty string).
export fixture
mkdir -p "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src" \
         "$fixture/bwa-mem3-sys/shim" \
         "$fixture/bwa-mem3-rs/src"

cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/foo.cpp" <<'EOF'
void foo(int a, int b)
{
    int original_marker = 1;
}
EOF
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/keep.cpp" <<'EOF'
void keep_fn() {}
EOF
cat > "$fixture/bwa-mem3-sys/shim/dummy.cpp" <<'EOF'
// See gotcha #12 for how foo's return value is used downstream.
void call_it()
{
    foo(3,
        4);
}
EOF
cat > "$fixture/bwa-mem3-sys/build.rs" <<'EOF'
fn main() {
    let skip_common: &[&str] = &[
        "foo.cpp",
        "stale_tool.cpp",
        "FMI_search2.cpp",
        "reader_v2.cpp",
    ];
}
EOF
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h" <<'EOF'
typedef struct mem_opt_t {
    int w;
} mem_opt_t;
EOF
cat > "$fixture/bwa-mem3-rs/src/opts.rs" <<'EOF'
// existing_field is already wired up.
pub struct MemOpts { existing_field: i32 }
EOF

# commit.gpgsign is forced off, not just identity-overridden: a developer with
# global signing on (this repo's CONTRIBUTING asks for it) would otherwise have
# this fixture commit signed, and it dies with "fatal: failed to write commit
# object" whenever the signing agent is locked or unavailable. CI has no signing
# config, so the failure reproduces only on a developer machine.
( cd "$fixture" && git init -q \
    && git -c user.email=t@t -c user.name=t add \
        bwa-mem3-sys/vendor/bwa-mem3/src/foo.cpp \
        bwa-mem3-sys/vendor/bwa-mem3/src/keep.cpp \
        bwa-mem3-sys/shim/dummy.cpp \
        bwa-mem3-sys/build.rs \
        bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h \
        bwa-mem3-rs/src/opts.rs \
    && git -c user.email=t@t -c user.name=t -c commit.gpgsign=false commit -q -m init )

REPO_ROOT="$fixture"
VENDOR="bwa-mem3-sys/vendor/bwa-mem3"
VENDOR_ABS="$REPO_ROOT/$VENDOR"

# --- check 5: unchanged tracked contract ---
UPSTREAM_CONTRACTS=("foo:src/foo.cpp")
out="$(check_call_site_contracts)"
contains "check 5: unchanged contract reported as such" "$out" "No change to any tracked contract"

# --- check 5: changed body surfaces the diff AND the shim call site ---
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/foo.cpp" <<'EOF'
void foo(int a, int b)
{
    int original_marker = 1;
    int added_marker = 2;
}
EOF
out="$(check_call_site_contracts)"
contains "check 5: changed contract is flagged" "$out" "\`foo\` changed"
contains "check 5: changed contract shows the diff" "$out" "added_marker"
contains "check 5: changed contract lists the shim call site" "$out" "dummy.cpp"

# --- check 5: call-site listing keeps a multi-line call's full argument
# list, shows repo-relative paths, and drops a paren-less prose mention ---
# This is the review-flagged gap: a bare first-line grep match previously
# showed only "foo(3," and hid the "4);" continuation line entirely, which
# is exactly what happened on the real tree for mem_kernel1_core /
# mem_kernel2_core / mem_pair_resolve / mem_reg2aln (all called across
# several lines in bwa_shim_align.cpp -- verified separately, see the task
# report). The fixture's dummy.cpp above is the minimal repro: a two-line
# call plus a comment that mentions "foo" without an adjacent "(".
contains "check 5: multi-line call site keeps the second argument line" "$out" "4);"
absent "check 5: paren-less prose mention of the function name is excluded" "$out" "foo's return value"
absent "check 5: call-site paths are repo-relative, not $fixture's absolute prefix" "$out" "$fixture"

# MUTATION (expected red): reproduce the ORIGINAL bare `grep -rn "$fn"` (no
# -A, no paren anchor, no path-stripping) against this exact fixture and
# confirm it fails both of the properties above -- it truncates the
# multi-line call to its first line only, AND it surfaces the paren-less
# comment as a false "call site".
out_unfixed="$(grep -rn "foo" "$fixture/bwa-mem3-sys/shim/" || printf '(none found)\n')"
absent "MUTATION (expected red): unfixed grep truncates the multi-line call, missing the second argument line" "$out_unfixed" "4);"
contains "MUTATION (expected red): unfixed grep surfaces the paren-less prose mention as a call site" "$out_unfixed" "foo's return value"

# --- check 5: a real change must not abort main()'s remaining checks ---
# check_call_site_contracts's last statement is `[ "$changed" -eq 0 ] &&
# printf ...`, which is itself false (and so returns non-zero) whenever a
# contract DID change — exactly the case this check exists to surface. Since
# main() calls it as a bare statement, that would silently abort the rest of
# the report (checks 6-11 never run) precisely when there is real drift to
# report. The fix is a trailing `return 0`; this asserts a command AFTER the
# call is actually reached, using the bare-statement-in-a-child-`bash -c`
# pattern (command substitution would swallow the abort the same way
# check_true/check_false do — see the note on the check_new_opts assertion
# below). Confirmed by mutation: deleting the trailing `return 0` makes this
# assertion fail (and, more drastically, aborts the whole test file with no
# output at all, since the existing `out="$(check_call_site_contracts)"`
# assignments above are themselves bare top-level statements too).
rc=0
bash -c '
    source "$here/../bwa-mem3-drift-report.sh"
    REPO_ROOT="$fixture"
    VENDOR="bwa-mem3-sys/vendor/bwa-mem3"
    VENDOR_ABS="$REPO_ROOT/$VENDOR"
    UPSTREAM_CONTRACTS=("foo:src/foo.cpp")
    check_call_site_contracts >/dev/null
    echo reached
' >/dev/null 2>&1 || rc=$?
check "check 5: a changed contract does not abort the caller" "0" "$rc"

# --- check 5: NOT FOUND when UPSTREAM_CONTRACTS names the wrong file ---
# This is exactly the mem_gen_alt path bug found against the real vendor tree
# (the brief's UPSTREAM_CONTRACTS pointed mem_gen_alt at src/bwamem.cpp, but
# it is defined in src/bwamem_extra.cpp) — pointing at a file that doesn't
# define the function must report NOT FOUND, not silently pass.
# Read by check_call_site_contracts, defined in the sourced (and
# source=/dev/null-annotated) script — invisible to shellcheck here, same as
# the git() stub function above.
# shellcheck disable=SC2034
UPSTREAM_CONTRACTS=("foo:src/keep.cpp")
out="$(check_call_site_contracts)"
contains "check 5: wrong-file contract path reports NOT FOUND" "$out" "NOT FOUND"

# --- check 6: skip_common triage, including the exact hazard the brief warns about ---
# stale_tool.cpp: absent from src/ -> GONE. FMI_search2.cpp / reader_v2.cpp:
# absent too, but their names carry an uppercase letter and a digit
# respectively — deliberately chosen to probe whether the extraction regex
# actually matches every entry, per the task's explicit instruction to verify
# this rather than transcribe it.
out="$(check_source_inventory)"
contains "check 6: present skip_common entry" "$out" "\`foo.cpp\` — present"
contains "check 6: stale lowercase skip_common entry reported GONE" "$out" "\`stale_tool.cpp\` — **GONE**"
contains "check 6: stale skip_common entry with an uppercase letter is still matched" "$out" "\`FMI_search2.cpp\` — **GONE**"
contains "check 6: stale skip_common entry with a digit is still matched" "$out" "\`reader_v2.cpp\` — **GONE**"

# MUTATION (expected red): the brief's original regex character class
# (`[a-z_]+`) silently drops any skip_common entry containing an uppercase
# letter or a digit — the exact vacuous-check failure mode the task flagged.
# All of the repo's REAL skip_common entries happen to be all-lowercase
# today (verified separately), so this only bites a future entry; this
# fixture is what makes that failure mode visible now.
out_narrow="$(mawk '/let skip_common/,/\];/' "$fixture/bwa-mem3-sys/build.rs" \
    | grep -oE '"[a-z_]+\.cpp"' | tr -d '"')"
absent "MUTATION (expected red): narrow regex silently misses the uppercase entry" "$out_narrow" "FMI_search2.cpp"
absent "MUTATION (expected red): narrow regex silently misses the digit-bearing entry" "$out_narrow" "reader_v2.cpp"
contains "MUTATION baseline: narrow regex still matches the plain lowercase entry" "$out_narrow" "stale_tool.cpp"

# --- check 6: a stale extraction pattern must be loud, not silent ---
# Review-flagged gap: build.rs already declares a sibling list in a
# different style (`const KERNEL_SRCS: &[&str] = &[...]`), so a plausible
# future rename of skip_common to a `const` binding is not hypothetical.
# Without a guard, a pattern that stops matching prints the section header
# and then NOTHING AT ALL -- exit 0, no warning -- reading as "no stale
# exclusions" when the truth is "the extraction itself is broken." That is
# the exact failure mode this check exists to prevent.
cp "$fixture/bwa-mem3-sys/build.rs" "$tmp/build.rs.let-style"

# First: confirm the WIDENED anchor (shipped fix) copes with a rename to
# `const SKIP_COMMON` gracefully, with no warning needed.
cat > "$fixture/bwa-mem3-sys/build.rs" <<'EOF'
fn main() {
    const SKIP_COMMON: &[&str] = &[
        "foo.cpp",
        "stale_tool.cpp",
    ];
}
EOF
out="$(check_source_inventory)"
contains "check 6: a renamed skip_common binding (const SKIP_COMMON) is still matched" "$out" "\`foo.cpp\` — present"
contains "check 6: a renamed binding's stale entry is still matched too" "$out" "\`stale_tool.cpp\` — **GONE**"
absent "check 6: a renamed binding the widened anchor DOES match does not warn" "$out" "Extraction matched nothing"

# Second: a rename the widened anchor can't anticipate either (e.g. neither
# `let`/`const skip_common` nor `SKIP_COMMON`) must trigger the loud warning
# instead of silently printing nothing.
cat > "$fixture/bwa-mem3-sys/build.rs" <<'EOF'
fn main() {
    static UNANTICIPATED_RENAME: &[&str] = &[
        "foo.cpp",
    ];
}
EOF
out="$(check_source_inventory)"
contains "check 6: an anchor-defeating rename triggers the loud warning, not silence" "$out" "Extraction matched nothing"

# MUTATION (expected red): reproduce the ORIGINAL, unwidened `/let
# skip_common/` anchor (no guard) against the const-SKIP_COMMON fixture from
# above and confirm it goes silent -- exactly the reviewer's demonstration.
cat > "$fixture/bwa-mem3-sys/build.rs" <<'EOF'
fn main() {
    const SKIP_COMMON: &[&str] = &[
        "foo.cpp",
        "stale_tool.cpp",
    ];
}
EOF
check_source_inventory_original_anchor() {
    local entries
    entries="$(mawk '/let skip_common/,/\];/' "$fixture/bwa-mem3-sys/build.rs" \
                 | grep -oE '"[a-zA-Z_0-9]+\.cpp"' | tr -d '"' || true)"
    if [ -z "$entries" ]; then
        printf 'nothing printed after the header (exit 0, no warning)\n'
    else
        printf 'entries=[%s]\n' "$entries"
    fi
}
out_stale="$(check_source_inventory_original_anchor)"
contains "MUTATION (expected red): original /let skip_common/ anchor goes silent on a renamed const binding" "$out_stale" "nothing printed after the header"

# Restore the fixture's build.rs to its normal let-style form for the tests below.
cp "$tmp/build.rs.let-style" "$fixture/bwa-mem3-sys/build.rs"

# --- check 6: a newly-added .cpp not yet tracked at HEAD is reported ---
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/new_dep.cpp" <<'EOF'
void new_dep_fn() {}
EOF
out="$(check_source_inventory)"
contains "check 6: a new untracked .cpp is reported in the inventory diff" "$out" "new_dep.cpp"

# --- check 8: new field triage (skip / already-referenced / consider-exposing) ---
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h" <<'EOF'
typedef struct mem_opt_t {
    int w;
    int existing_field; // consider
    int internal_state; // NOT a user option
    double new_score; // consider
} mem_opt_t;
EOF
out="$(check_new_opts)"
contains "check 8: internal-marked field is skipped" "$out" "internal_state"
contains "check 8: internal-marked field says skip" "$out" "skip (upstream marks it internal)"
contains "check 8: field already referenced in opts.rs is flagged as such" "$out" "existing_field"
contains "check 8: field already referenced in opts.rs says already referenced" "$out" "already referenced in \`opts.rs\`"
contains "check 8: unreferenced new field is flagged for exposure" "$out" "new_score"
contains "check 8: unreferenced new field says consider exposing" "$out" "consider exposing"

# --- check 8: no new fields -> trivial "no new fields" message ---
cp "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h" "$tmp/bwamem_with_new_fields.h"
cat > "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h" <<'EOF'
typedef struct mem_opt_t {
    int w;
} mem_opt_t;
EOF
out="$(check_new_opts)"
contains "check 8: unchanged struct reports no new fields" "$out" "No new fields."

# MUTATION (expected red): reproduce check_new_opts's diff|grep|sed line
# WITHOUT the `|| true` fallback, restore the struct that actually differs
# from HEAD, and confirm it aborts the process instead of reporting — this
# is the exact hazard flagged repeatedly across this plan (a diff that finds
# real differences exits 1, which is the ORDINARY case here, not an error).
# Run as a bare statement in a child `bash -c`, per the same reasoning as
# assert_refreshed_tree's case (b): a `set -e` abort inside a function called
# through check_true/check_false is unobservable, because bash suspends
# errexit for the entire call stack of a command used as an if-condition.
cp "$tmp/bwamem_with_new_fields.h" "$fixture/bwa-mem3-sys/vendor/bwa-mem3/src/bwamem.h"
export REPO_ROOT VENDOR VENDOR_ABS here fixture
rc=0
bash -c '
    source "$here/../bwa-mem3-drift-report.sh"
    REPO_ROOT="$fixture"
    VENDOR="bwa-mem3-sys/vendor/bwa-mem3"
    VENDOR_ABS="$REPO_ROOT/$VENDOR"
    check_new_opts_unfixed() {
        report_section "8. New mem_opt_t fields (unfixed)"
        local old new o n field
        old="$(mktemp)"; new="$(mktemp)"
        old_file "$VENDOR/src/bwamem.h" > "$old"
        new_file "$VENDOR/src/bwamem.h" > "$new"
        o="$(mktemp)"; n="$(mktemp)"
        struct_body "$old" mem_opt_t > "$o"
        struct_body "$new" mem_opt_t > "$n"
        local added
        added="$(diff "$o" "$n" | grep "^>" | sed "s/^> //")"
        printf "unreachable: %s\n" "$added"
    }
    check_new_opts_unfixed
' >/dev/null 2>&1 || rc=$?
check "MUTATION (expected red): check_new_opts without || true aborts instead of reporting" "1" "$rc"

# The FIXED check_new_opts, same fixture, must not abort and must report.
#
# Deliberately NOT written as `out="$(check_new_opts)" || rc=$?`: wrapping the
# call in a command substitution spawns a subshell that, by default (without
# `shopt -s inherit_errexit`), does NOT inherit errexit from the caller — so
# an internal abort inside check_new_opts would be swallowed there too, blind
# in the same way check_true/check_false are, just via a different mechanism.
# Verified directly: a version of this assertion using `out="$(check_new_opts)"`
# reported rc=0 even against the real script with `|| true` deleted from the
# diff|grep|sed line. Calling check_new_opts as the bare, final statement of a
# freshly spawned `bash -c` (same pattern as assert_refreshed_tree's case (b))
# keeps errexit live for the duration of the call, so a real abort is
# observable; stdout is captured via redirection instead of substitution.
outfile="$tmp/check_new_opts_fixed.out"
rc=0
bash -c '
    source "$here/../bwa-mem3-drift-report.sh"
    REPO_ROOT="$fixture"
    VENDOR="bwa-mem3-sys/vendor/bwa-mem3"
    VENDOR_ABS="$REPO_ROOT/$VENDOR"
    check_new_opts
' > "$outfile" 2>&1 || rc=$?
check "check_new_opts (fixed): does not abort on a real struct diff" "0" "$rc"
contains "check_new_opts (fixed): reports the new field despite diff exiting 1 internally" "$(cat "$outfile")" "new_score"

REPO_ROOT="$orig_repo_root"
VENDOR="$orig_vendor"
VENDOR_ABS="$orig_vendor_abs"

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
