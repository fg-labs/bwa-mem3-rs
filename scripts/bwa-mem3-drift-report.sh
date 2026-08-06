#!/usr/bin/env bash
# Report what a vendor refresh broke, as markdown on stdout.
#
# Run AFTER scripts/refresh-bwa-mem3.sh has rewritten the working tree and
# BEFORE anything is committed: the "before" side of every comparison comes
# from `git show HEAD:<path>`.
#
# Usage: scripts/bwa-mem3-drift-report.sh
#
# Exit status is 0 whatever it finds — a refresh that breaks the build is the
# expected case, and the report is the deliverable. Non-zero means misuse
# (nothing refreshed, or unrelated files dirty).
#
# The primary check is `cargo ci-build`. The compiler subsumes pruned
# submodules, Makefile-generated headers, new external dependencies, changed
# function arity, and renamed struct fields as a single class, with none of
# the false positives a hand-rolled include scan produces (~20 today on an
# unchanged tree, all behind #ifdef USE_MALLOC_WRAPPERS / VTUNE_ANALYSIS /
# arch guards) and without tripping over ext/sse2neon being a submodule.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENDOR="bwa-mem3-sys/vendor/bwa-mem3"
# Not read by this task's two checks; kept as the absolute-path counterpart to
# $VENDOR for the struct/flag/enum-diff checks (Task 11/12) built on this
# harness, which resolve vendored headers by absolute path.
# shellcheck disable=SC2034
VENDOR_ABS="$REPO_ROOT/$VENDOR"

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

# One private directory for every temp file this script allocates, released by
# the EXIT trap on ANY exit path — including the one that actually happens in
# practice: check_build shells out to a multi-minute `cargo ci-build`, and
# main() runs each check inside a command substitution, so a Ctrl-C lands
# mid-capture with several temp files allocated and the per-function `rm -f`
# at the tail of that check never reached. The eager `rm -f` calls are kept as
# well, so a long run releases files as it goes rather than only at exit.
#
# A directory, deliberately, rather than an array of file paths: every caller
# invokes mktemp_tracked as `x="$(mktemp_tracked)"`, and a command
# substitution runs in a SUBSHELL — so `ARRAY+=(...)` inside the function
# would mutate a copy that is discarded on return, leaving the parent's array
# empty and the trap removing nothing. TMPDIR_SELF is assigned once here at
# top level, so the trap closes over a path that is correct no matter how many
# subshells deep the allocation happened.
TMPDIR_SELF="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_SELF"' EXIT INT TERM

# mktemp, inside this script's auto-cleaned directory. Use instead of bare
# `mktemp` so an interrupted run leaks nothing.
mktemp_tracked() {
    mktemp "$TMPDIR_SELF/drift.XXXXXX"
}

report_section() { printf '\n## %s\n\n' "$1"; }

# The committed (pre-refresh) content of a vendored path, or empty if it is new.
old_file() {
    git -C "$REPO_ROOT" show "HEAD:$1" 2>/dev/null || true
}

# The working-tree (post-refresh) content of a path, or empty if it is gone.
new_file() {
    cat "$REPO_ROOT/$1" 2>/dev/null || true
}

# The pruned-subtree paths from scripts/vendor-drop-subtrees.txt, one per
# line, parsed exactly the way refresh-bwa-mem3.sh parses the same file
# (strip a trailing comment, strip all whitespace, drop empty lines) so the
# two scripts cannot disagree about what "already pruned" means. Deliberately
# not a second hardcoded copy of the list — see the file's own header comment
# on why.
#
# A missing list file writes an ERROR to stderr and the caller then sees an
# empty drop set, so every added path stays in `added` and the gate FIRES.
# That is the safe direction (a refresh is blocked, not silently un-gated),
# but note the `exit 1` below does NOT propagate: bash does not honour
# `errexit` for a failing bare assignment inside a command substitution, and
# every caller of this function runs it as `$(drop_subtrees ...)`. Do not
# rely on it to abort the script.
drop_subtrees() {
    local list="$REPO_ROOT/scripts/vendor-drop-subtrees.txt"
    [ -f "$list" ] || { echo "ERROR: missing $list" >&2; exit 1; }
    local line
    while IFS= read -r line; do
        line="${line%%#*}"
        line="$(printf '%s' "$line" | tr -d '[:space:]')"
        [ -n "$line" ] && printf '%s\n' "$line"
    done < "$list"
}

# Submodule paths present in the new .gitmodules but not the old one, MINUS
# anything scripts/vendor-drop-subtrees.txt already prunes.
#
# A new submodule means refresh-bwa-mem3.sh could not know to prune it, so the
# clone kept it and a commit would carry its whole source tree. The caller uses
# this to gate PR creation.
#
# The MINUS half exists because refresh-bwa-mem3.sh prunes the pruned
# submodule's *subtree* but leaves .gitmodules itself fully tracked (it never
# edits it) — so a submodule that is already on the drop list still shows up
# as "new" here on every single future refresh, forever. Without this filter
# the documented remedy ("add it to vendor-drop-subtrees.txt and re-dispatch")
# would have zero effect: the gate would still fire on the very next run.
#
# Extracts the path value directly with grep -oE rather than trimming
# whitespace afterward with `tr -d '[:space:]'`: that would delete the
# newlines between grep's output lines too, gluing multiple submodule paths
# into one bogus token (e.g. "ext/aext/b") and hiding real drift. Uses
# POSIX [[:space:]] instead of the GNU-only \s shorthand, which BSD sed
# (macOS) does not expand and silently leaves as a literal "s".
gitmodules_new_submodules() {
    local old="$1" new="$2"
    local old_paths new_paths added dropped
    old_paths="$(grep -oE '^[[:space:]]*path[[:space:]]*=[[:space:]]*[^[:space:]]+' "$old" 2>/dev/null | sed -E 's/^[[:space:]]*path[[:space:]]*=[[:space:]]*//' | sort -u || true)"
    new_paths="$(grep -oE '^[[:space:]]*path[[:space:]]*=[[:space:]]*[^[:space:]]+' "$new" 2>/dev/null | sed -E 's/^[[:space:]]*path[[:space:]]*=[[:space:]]*//' | sort -u || true)"
    added="$(comm -13 <(printf '%s\n' "$old_paths") <(printf '%s\n' "$new_paths") | grep -v '^$' || true)"
    [ -n "$added" ] || return 0
    dropped="$(drop_subtrees | sort -u)"
    comm -23 <(printf '%s\n' "$added" | sort -u) <(printf '%s\n' "$dropped") | grep -v '^$' || true
}

# The body of a C struct typedef, handling BOTH shapes upstream uses:
#   typedef struct mem_opt_t { ... } mem_opt_t;   (named)
#   typedef struct { ... } mem_pestat_t;          (anonymous)
#
# Scans BACKWARD from the closing `} <name>;` to the nearest preceding
# `typedef struct`, because a forward range keyed on `<name> {` finds nothing
# for the anonymous form — which is exactly how mem_pestat_t is declared
# (bwamem.h:239-243).
struct_body() {
    local file="$1" name="$2"
    mawk -v name="$name" '
        { lines[NR] = $0 }
        $0 ~ "^\\}[[:space:]]*" name "[[:space:]]*;" { close_at = NR }
        END {
            if (close_at == 0) exit 0
            for (i = close_at; i >= 1; i--) {
                if (lines[i] ~ /typedef[[:space:]]+struct/) { open_at = i; break }
            }
            if (open_at == 0) exit 0
            for (i = open_at; i <= close_at; i++) print lines[i]
        }
    ' "$file"
}

# `NAME VALUE` for every MEM_F_* define in a file, sorted by name so a diff is
# insensitive to declaration order. Sorting is safe here specifically because
# each #define is an independent macro — unlike enum members (see
# enum_bodies below), nothing about a flag's meaning depends on where among
# its siblings it is declared.
flag_defines() {
    grep -E '^[[:space:]]*#define[[:space:]]+MEM_F_[A-Z0-9_]+[[:space:]]+' "$1" 2>/dev/null \
        | mawk '{ print $2, $3 }' | sort || true
}

# Every enum declaration in a header, verbatim and in source order.
#
# bwamem.h uses BOTH shapes: a single-line named form
# (`enum mem_meth_scoring { A = 0, B = 1 };`) and a multi-line anonymous
# `typedef` form (`typedef enum {\n ... \n} seed_order_t;`). Triggers on the
# bare word `enum` (word-bounded by hand since mawk has no `\b`) and closes on
# the first `;` seen after that — which ends a single-line enum immediately
# and an anonymous typedef enum at its `} name;` line — without hardcoding
# either enum's name, so a brand-new enum is picked up for free.
#
# Deliberately does NOT sort (contrast flag_defines above): an enum member
# with no explicit `= N` takes its value from its position among its
# siblings, so reordering members silently renumbers every unlabeled one
# after the moved member — a real semantic change a sorted, order-blind diff
# would hide. An earlier draft of this check used a `MEM_[A-Z_]+\s*=` grep,
# which also missed seed_order_t entirely: none of its members
# (SEED_ORDER_*) carry a MEM_ prefix, even though that enum is exactly what
# opts.rs's SeedOrder mirrors.
#
# Known blind spot: the trigger is the bare word "enum" anywhere on a line,
# not just at a declaration, so a future comment containing that word (e.g.
# "// like a C enum") would open a spurious capture through to the next `;`.
# Harmless today — this check is advisory, and no such comment exists in
# bwamem.h now — but worth knowing if this section ever prints noise.
enum_bodies() {
    mawk '
        /(^|[^a-zA-Z_])enum([^a-zA-Z_]|$)/ { in_enum = 1 }
        in_enum { print }
        in_enum && /;/ { in_enum = 0 }
    ' "$1" 2>/dev/null || true
}

# The definition of a C/C++ function, from its signature line to the first
# column-0 closing brace. Good enough for upstream's style, which puts the
# opening brace on its own line at column 0 and closes at column 0.
#
# The trigger requires the signature line ITSELF to start at column 0 (not
# just be un-semicolon-terminated): upstream's call sites are indented, but a
# multi-line call — one whose first physical line ends in a comma rather than
# a semicolon, e.g.
#   mem_reg2sam(w->opt, mem_aln_bns(w), mem_aln_pac(w), &w->seqs[i],
#               &w->regs[i], 0, 0);
# — satisfies the "contains name( and doesn't end in ;" test just as well as
# a real definition. Without the column-0 requirement this hijacks the scan
# on any call site that textually precedes the definition in the same file,
# which is exactly what happens for mem_reg2sam and mem_mark_primary_se in
# bwamem.cpp today (both are called, indented, earlier in the file than they
# are defined) — verified against the vendored tree while writing this check.
# A leading `#` is also excluded so a macro invocation of the same name
# (`#define FOO(x) name(x)`) can't trigger it either.
function_body() {
    local file="$1" name="$2"
    mawk -v name="$name" '
        $0 ~ /^[^[:space:]#]/ && index($0, name "(") && !started && $0 !~ /;[[:space:]]*$/ { started = 1 }
        started { print }
        started && /^\}/ { exit }
    ' "$file"
}

# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

# Formats the "Build FAILED" body for a given log file and exit status. Split
# out from check_build so scripts/tests/test-drift-report.sh can exercise the
# filtering/tail logic against fixture logs without actually running
# `cargo ci-build` — the thing under test here is the report's OUTPUT SHAPE,
# not whether this checkout currently builds.
format_build_failure() {
    local log="$1" status="$2"
    printf 'Build FAILED (exit %s).\n\n' "$status"
    # Two sections below, not one filtered dump. A pattern list is never
    # exhaustive — this exact line was already widened twice (first for
    # `: error:`, then for `ld:`/`Undefined symbols`) and STILL missed
    # clang's `fatal error:` on a real historical bump (the v0.2.2 ->
    # v0.6.0 replay in scripts/tests/drift-report-replay.sh caught it: the
    # actual failure was `fastmap.cpp:47:10: fatal error: 'version.h' file
    # not found`). The real defect was never the pattern; it was treating
    # "grep found a hit" as "the report is complete" — cargo's own wrapper
    # line ("error: failed to run custom build command...") matches
    # `^error` on EVERY build-script failure, so the old single-shot
    # `grep ... || tail -40` fallback never reached the tail once the
    # filter had ANY hit, even when that hit was the wrapper line telling
    # the human nothing. So: print the filtered hits (if any) as a targeted
    # worklist, then UNCONDITIONALLY print a raw tail too, so a diagnostic
    # shaped like nothing on the list below is still visible. `fatal
    # error:` is added to the pattern as belt-and-braces, not as the fix —
    # the fix is that this section no longer depends on the pattern being
    # complete.
    local hits
    # `|| true` INSIDE the substitution, not chained after it: grep exits 1
    # when nothing matches, and under `set -o pipefail` that would make the
    # bare assignment abort the whole script via `set -e` before
    # `[ -n "$hits" ]` below ever runs.
    hits="$(grep -E '(^error|: (fatal )?error:|^warning: unused|  -->|ld:|^Undefined symbols)' "$log" | head -40 || true)"
    if [ -n "$hits" ]; then
        printf 'Lines matching known error patterns:\n\n```\n'
        printf '%s\n' "$hits"
        printf '```\n\n'
    else
        printf 'No line matched a known error pattern.\n\n'
    fi
    # Deliberately unconditional, and deliberately NOT deduped against the
    # filtered hits above: the two sections answer different questions
    # (targeted-by-pattern vs. "the actual end of the log, in case the
    # filter missed it"), and reliably detecting the overlap is more
    # machinery than the handful of possibly-repeated lines it would save.
    # Accept the overlap; the label says what this section is.
    printf 'Last 40 lines of the raw log (may repeat lines already shown above):\n\n```\n'
    tail -40 "$log"
    printf '```\n\n'
    printf 'Full log is in the workflow run artifacts.\n'
}

# Check 2 (primary): does the refreshed tree still build?
check_build() {
    report_section "2. Build (\`cargo ci-build\`)"
    local log status
    log="$(mktemp_tracked)"
    status=0
    (cd "$REPO_ROOT" && cargo ci-build) > "$log" 2>&1 || status=$?
    if [ "$status" -eq 0 ]; then
        printf 'Builds clean. Note this proves compilation only — it does not\n'
        printf 'prove byte-parity with the CLI. See the manual checklist below.\n'
    else
        format_build_failure "$log" "$status"
    fi
    rm -f "$log"
}

# Check 1 (hard gate): a new submodule makes a refresh commit wrong.
check_gitmodules() {
    local old new
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    old_file "$VENDOR/.gitmodules" > "$old"
    new_file "$VENDOR/.gitmodules" > "$new"
    local added
    added="$(gitmodules_new_submodules "$old" "$new")"
    rm -f "$old" "$new"

    report_section "1. Submodules (\`.gitmodules\`)"
    if [ -z "$added" ]; then
        printf 'No new submodules.\n'
        return 0
    fi
    # NEW_SUBMODULE is the machine-readable marker the workflow greps for to
    # gate PR creation. It is NOT the first line of the report or even of
    # this check's own output — report_section above already emitted the
    # "## 1. Submodules" heading. Consumers must match it anchored against
    # the whole report (`grep -qE '^NEW_SUBMODULE$'`), never by position
    # (e.g. `head -1`), since a bare-line anchor can't collide with prose or
    # pasted compiler output and doesn't break if a check is inserted earlier.
    printf 'NEW_SUBMODULE\n\n'
    printf 'Upstream added submodule(s) this refresh could not know to prune:\n\n'
    # Backticks below are literal markdown code-formatting, not command
    # substitution — single-quoting deliberately prevents expansion.
    # shellcheck disable=SC2016
    printf '%s\n' "$added" | sed 's/^/- `/; s/$/`/'
    # shellcheck disable=SC2016
    printf '\n**No PR was created.** `refresh-bwa-mem3.sh` clones with\n'
    # shellcheck disable=SC2016
    printf '`--recurse-submodules`, so committing now would vendor the whole\n'
    printf 'subtree. Decide whether each is needed at build time:\n\n'
    # `--` stops bash's printf builtin from parsing the leading "-" of the
    # format string itself as an option (it otherwise fails with "invalid
    # option" and silently drops the line).
    # shellcheck disable=SC2016
    printf -- '- not needed → add it to `scripts/vendor-drop-subtrees.txt`\n'
    printf -- '- needed → leave it, and confirm the build picks up its headers\n\n'
    printf 'Then re-run the workflow.\n'
}

# Check 3: mem_opt_t / mem_pestat_t body drift.
#
# Advisory only: shim/bwa_shim_layout_assert.cpp already static_asserts every
# field offset against upstream, so real drift fails `cargo build` with a
# precise message. This just names the fields up front so the human knows what
# they are walking into.
check_structs() {
    report_section "3. \`mem_opt_t\` / \`mem_pestat_t\` layout"
    local old new o n s
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    old_file "$VENDOR/src/bwamem.h" > "$old"
    new_file "$VENDOR/src/bwamem.h" > "$new"
    local any=0
    for s in mem_opt_t mem_pestat_t; do
        o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
        struct_body "$old" "$s" > "$o"
        struct_body "$new" "$s" > "$n"
        if ! diff -q "$o" "$n" >/dev/null; then
            any=1
            # Backticks below are literal markdown, not command substitution.
            # shellcheck disable=SC2016
            printf '### `%s` changed\n\n```diff\n' "$s"
            diff -u "$o" "$n" | tail -n +3 || true
            printf '```\n\n'
        fi
        rm -f "$o" "$n"
    done
    rm -f "$old" "$new"
    if [ "$any" -eq 0 ]; then
        printf 'No change.\n'
    else
        # shellcheck disable=SC2016
        printf 'Update `bwa-mem3-sys/shim/bwa_shim_types.h` and the field list\n'
        # shellcheck disable=SC2016
        printf 'in `shim/bwa_shim_layout_assert.cpp` to match. The layout\n'
        printf 'assertions will fail the build until you do — that is the real\n'
        printf 'guard; this section is just the heads-up (gotcha #2).\n'
    fi
}

# Check 4: flag and enum SET drift.
#
# Per-flag static_asserts cannot see a flag that does not yet exist, and
# bindgen's allowlist_var("MEM_F_.*") reads our POD copy rather than upstream,
# so a NEW upstream flag is invisible to both. This covers that direction, and
# the same argument applies to enum values: opts.rs mirrors them by hand.
#
# A NEW flag needs edits in TWO files: `bwa_shim_types.h` (the #define block
# bindgen reads), and, inside `bwa_shim_layout_assert.cpp`, THREE coordinated
# edits — the `pod_flags::MEM_F_X_v` capture, its paired `#undef`, and the
# `BWA_SHIM_CK_FLAG(MEM_F_X)` invocation. Miss the `#undef` and the assert
# silently becomes a tautology (both sides read upstream's post-include
# value) — see the capture-step comment in that file for why.
check_flags_and_enums() {
    report_section "4. \`MEM_F_*\` flags and \`bwamem.h\` enums"
    local old new o n
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    old_file "$VENDOR/src/bwamem.h" > "$old"
    new_file "$VENDOR/src/bwamem.h" > "$new"

    o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
    flag_defines "$old" > "$o"
    flag_defines "$new" > "$n"
    if diff -q "$o" "$n" >/dev/null; then
        printf 'Flag set unchanged.\n\n'
    else
        # Backticks below are literal markdown, not command substitution.
        # shellcheck disable=SC2016
        printf '### `MEM_F_*` set changed\n\n```diff\n'
        diff -u "$o" "$n" | tail -n +3 || true
        printf '```\n\n'
        # shellcheck disable=SC2016
        printf 'Update `shim/bwa_shim_types.h` (the #define block bindgen\n'
        # shellcheck disable=SC2016
        printf 'reads), then in `shim/bwa_shim_layout_assert.cpp` add all\n'
        # shellcheck disable=SC2016
        printf 'THREE of: the `pod_flags::MEM_F_X_v` capture, its paired\n'
        # shellcheck disable=SC2016
        printf '`#undef`, and the `BWA_SHIM_CK_FLAG(MEM_F_X)` invocation —\n'
        # shellcheck disable=SC2016
        printf 'miss the `#undef` and the assert silently becomes a\n'
        printf 'tautology instead of catching a renumbered flag. A NEW flag\n'
        printf 'never fails the build on its own, and also never reaches\n'
        printf 'Rust because bindgen reads the POD copy.\n\n'
    fi
    rm -f "$o" "$n"

    o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
    enum_bodies "$old" > "$o"
    enum_bodies "$new" > "$n"
    if diff -q "$o" "$n" >/dev/null; then
        printf 'Enum values unchanged.\n'
    else
        printf '### Enum values changed\n\n```diff\n'
        diff -u "$o" "$n" | tail -n +3 || true
        printf '```\n\n'
        # shellcheck disable=SC2016
        printf 'Check the hand-written mirrors in `bwa-mem3-rs/src/opts.rs`\n'
        # shellcheck disable=SC2016
        printf '(`SeedOrder`, `MethScoring`). Their `TryFrom<i32>` returns\n'
        # shellcheck disable=SC2016
        printf '`Error::UnrecognizedEnum` for an unknown value, so a new\n'
        printf 'variant surfaces at runtime rather than mis-mapping — but the\n'
        printf 'mirror still needs teaching, and any new enum needs one too.\n'
    fi
    rm -f "$o" "$n" "$old" "$new"
}

# Upstream functions the shim calls or carries verbatim copies of. Format:
# <function>:<upstream file>. Keep in sync when the shim adopts a new upstream
# entry point — check 5 is the only thing watching these contracts.
#
# mem_gen_alt's definition lives in bwamem_extra.cpp, not bwamem.cpp (it is
# only declared, across two lines, in bwamem.h) — verified against the
# vendored tree; pointing this at the wrong file makes the check silently
# report a real function as "NOT FOUND" on every run.
UPSTREAM_CONTRACTS=(
    "worker_alloc:src/fastmap.cpp"
    "worker_free:src/fastmap.cpp"
    "mem_kernel1_core:src/bwamem.cpp"
    "mem_kernel2_core:src/bwamem.cpp"
    "mem_pair_resolve:src/bwamem_pair.cpp"
    "mem_gen_alt:src/bwamem_extra.cpp"
    "mem_reg2aln:src/bwamem.cpp"
    "mem_mark_primary_se:src/bwamem.cpp"
    "mem_reg2sam:src/bwamem.cpp"
    "mem_sam_pe:src/bwamem_pair.cpp"
)

# Check 5: contracts of upstream functions the shim calls or copies.
#
# A body diff alone is not enough. The shim passes specific arguments and
# replicates specific output policy, so a change in what an argument MEANS is
# silent even when the signature is unchanged — and a rename that makes it fail
# to compile can hide a semantic change underneath.
check_call_site_contracts() {
    report_section "5. Upstream contracts the shim depends on"
    printf 'The shim calls or carries copies of these. For each that changed,\n'
    printf 'read the diff and then re-check the shim call sites listed under it\n'
    printf '— a rename that breaks the build can hide a semantic change, and a\n'
    printf 'changed argument meaning is silent even with an identical signature\n'
    printf '(gotchas #8 and #12).\n\n'
    local entry fn path old new o n changed=0
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    for entry in "${UPSTREAM_CONTRACTS[@]}"; do
        fn="${entry%%:*}"
        path="${entry#*:}"
        old_file "$VENDOR/$path" > "$old"
        new_file "$VENDOR/$path" > "$new"
        o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
        function_body "$old" "$fn" > "$o"
        function_body "$new" "$fn" > "$n"
        if [ ! -s "$n" ]; then
            changed=1
            # Backticks below are literal markdown, not command substitution.
            # shellcheck disable=SC2016
            printf '### `%s` NOT FOUND in `%s`\n\n' "$fn" "$path"
            printf 'It was renamed, moved, or removed upstream. Find where it\n'
            # shellcheck disable=SC2016
            printf 'went and update `UPSTREAM_CONTRACTS` in this script.\n\n'
        elif ! diff -q "$o" "$n" >/dev/null; then
            changed=1
            # shellcheck disable=SC2016
            printf '### `%s` changed (`%s`)\n\n```diff\n' "$fn" "$path"
            diff -u "$o" "$n" | tail -n +3 || true
            # shellcheck disable=SC2016
            printf '```\n\nShim call sites:\n\n```\n'
            # -A3: a bare first-line match hides the argument list for every
            # multi-line call — verified against the real tree, where
            # mem_kernel1_core/mem_kernel2_core/mem_pair_resolve/mem_reg2aln
            # are all called across several lines, so a human previously saw
            # only e.g. "mem_kernel1_core(s->w.fmi, s->w.opt," and had to open
            # the file to see the rest. The pattern requires "(" (optionally
            # preceded by whitespace) right after the name so prose mentions
            # in comments (e.g. "matching mem_reg2sam's...", "mirrors
            # mem_sam_pe's...", never followed by "(") drop out — every real
            # call site in the shim today is written as `name(` or `name (`,
            # so this does not risk missing one; if a future call site is
            # ever split across a name-then-newline-then-"(" shape, this
            # would miss it and this comment's claim would need revisiting.
            # Paths are shown repo-relative, not as $REPO_ROOT's absolute path.
            grep -rn -A3 -E "${fn}[[:space:]]*\\(" "$REPO_ROOT/bwa-mem3-sys/shim/" 2>/dev/null \
                | sed "s|^$REPO_ROOT/||" \
                || printf '(no direct call site in shim/ -- see gotcha #12 if this\nfunction is mirrored rather than called directly)\n'
            printf '```\n\n'
        fi
        rm -f "$o" "$n"
    done
    rm -f "$old" "$new"
    [ "$changed" -eq 0 ] && printf 'No change to any tracked contract.\n'
    return 0
}

# Check 6: source inventory. A new .cpp is auto-globbed by build.rs and may
# fail to link; a new .c is silently ignored (the glob is *.cpp only); a
# skip_common entry that no longer exists upstream is a stale exclusion.
check_source_inventory() {
    report_section "6. Source inventory"
    local o n
    o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
    # 2>/dev/null: a first-ever-vendored src/ (no HEAD:$VENDOR/src to show)
    # makes git print a raw "fatal: Not a valid object name" to stderr, which
    # would otherwise land in the middle of the markdown report.
    git -C "$REPO_ROOT" ls-tree --name-only "HEAD:$VENDOR/src" 2>/dev/null | grep -E '\.(c|cpp)$' | sort > "$o" || true
    (cd "$VENDOR_ABS/src" && ls) | grep -E '\.(c|cpp)$' | sort > "$n" || true
    if diff -q "$o" "$n" >/dev/null; then
        printf 'No TUs added or removed.\n\n'
    else
        printf '```diff\n'
        diff -u "$o" "$n" | tail -n +3 || true
        printf '```\n\n'
        # Backticks below are literal markdown, not command substitution.
        # shellcheck disable=SC2016
        printf 'A new `.cpp` is auto-globbed by `build.rs` and may fail to\n'
        # shellcheck disable=SC2016
        printf 'link — add it to `skip_common` if it is out of scope. A new\n'
        # shellcheck disable=SC2016
        printf '`.c` is silently ignored (the glob is `*.cpp` only), so it sits\n'
        printf 'unused; that is usually fine but worth knowing.\n\n'
    fi
    rm -f "$o" "$n"

    # shellcheck disable=SC2016
    printf '### `skip_common` entries still present upstream\n\n'
    # Extracted into a variable (rather than piped straight into the while
    # loop's process substitution) so an empty result is itself detectable.
    # Anchored on `(let|const)[[:space:]]+skip_common|SKIP_COMMON` rather
    # than just `let skip_common`: build.rs already declares a sibling list
    # in a different style (`const KERNEL_SRCS: &[&str] = &[...]`), so a
    # plausible future rename of skip_common to a `const` (or to
    # SCREAMING_CASE, matching that sibling) would make the narrower anchor
    # match nothing. Without this guard, a stale pattern here used to fail
    # silently: the section header would print, then nothing at all, exit 0
    # -- reading as "no stale exclusions" when it was actually "the
    # extraction itself is broken." That is the one failure mode this check
    # exists to prevent, so a stale pattern must be loud, not silent.
    local name entries
    entries="$(mawk '/(let|const)[[:space:]]+skip_common|SKIP_COMMON/,/\];/' "$REPO_ROOT/bwa-mem3-sys/build.rs" \
                 | grep -oE '"[a-zA-Z_0-9]+\.cpp"' | tr -d '"' || true)"
    if [ -z "$entries" ]; then
        printf 'Extraction matched nothing -- the skip_common pattern in\n'
        printf 'check_source_inventory (scripts/bwa-mem3-drift-report.sh) is\n'
        printf 'stale relative to build.rs. Update the mawk anchor to match\n'
        printf 'however skip_common is declared now, then re-run this report.\n'
        return 0
    fi
    while IFS= read -r name; do
        [ -n "$name" ] || continue
        if [ -f "$VENDOR_ABS/src/$name" ]; then
            # shellcheck disable=SC2016
            printf -- '- `%s` — present\n' "$name"
        else
            # shellcheck disable=SC2016
            printf -- '- `%s` — **GONE**: stale exclusion in `build.rs`\n' "$name"
        fi
    done <<< "$entries"
}

# Check 7: upstream build dependencies. check.yml's e2e job and python.yml's
# integration job build upstream from the FULL Makefile (htslib, fast_reader),
# so they inherit upstream's whole dependency set even though this crate
# compiles none of it. A new -lfoo breaks e2e with the cause buried in a log.
check_dependencies() {
    report_section "7. Upstream build dependencies"
    local old new o n
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    old_file "$VENDOR/Makefile" > "$old"
    new_file "$VENDOR/Makefile" > "$new"
    o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
    # The anchor covers indirect dependency variables, not just LIBS/LDLIBS/
    # LDFLAGS themselves: LIBS_EXTRA, LIBSAIS_OPENMP_LIBS, and
    # MIMALLOC_LDFLAGS are all real Makefile variables (verified against the
    # vendored tree) that feed into LIBS/LDFLAGS via `LIBS += $(LIBS_EXTRA)`
    # -- an aggregate line like that does not change when the variable IT
    # REFERENCES gains a new -lfoo, so anchoring on the exact names LIBS/
    # LDLIBS/LDFLAGS alone would miss exactly the kind of new dependency this
    # check exists to catch.
    grep -E '^[[:space:]]*[A-Za-z_]*(LIBS|LDLIBS|LDFLAGS)[A-Za-z_]*[[:space:]]*[+:]?=|pkg-config|\./configure' "$old" | sort -u > "$o" || true
    grep -E '^[[:space:]]*[A-Za-z_]*(LIBS|LDLIBS|LDFLAGS)[A-Za-z_]*[[:space:]]*[+:]?=|pkg-config|\./configure' "$new" | sort -u > "$n" || true
    if diff -q "$o" "$n" >/dev/null; then
        # Backticks below are literal markdown, not command substitution.
        # shellcheck disable=SC2016
        printf 'No change to `LIBS` / `pkg-config` / `./configure`.\n'
    else
        printf '```diff\n'
        diff -u "$o" "$n" | tail -n +3 || true
        printf '```\n\n'
        printf 'A new external library may need adding to the apt list in\n'
        # shellcheck disable=SC2016
        printf '`.github/workflows/check.yml` (e2e job) and\n'
        # shellcheck disable=SC2016
        printf '`.github/workflows/python.yml` (integration job). Not a\n'
        # shellcheck disable=SC2016
        printf 'mechanical mapping: `-lhts`/`-lz-ng` come from vendored\n'
        # shellcheck disable=SC2016
        printf 'submodules and `-lbwa` is upstream'"'"'s own archive.\n'
    fi
    rm -f "$o" "$n" "$old" "$new"
}

# Check 8: new mem_opt_t fields worth exposing on MemOpts.
#
# Filtered by trailing-comment triage, because most new fields are not user
# options: 0.6.0 added est_insert_high (marked "runtime state (NOT a user
# option)"), the derived mat_ot/mat_ob matrices, and bam_mode/bam_level (moot —
# the shim emits BAM itself); v0.8.0 adds a non-settable `compat` pointer.
check_new_opts() {
    report_section "8. New \`mem_opt_t\` fields vs. the Rust API"
    local old new o n field
    old="$(mktemp_tracked)"; new="$(mktemp_tracked)"
    old_file "$VENDOR/src/bwamem.h" > "$old"
    new_file "$VENDOR/src/bwamem.h" > "$new"
    o="$(mktemp_tracked)"; n="$(mktemp_tracked)"
    struct_body "$old" mem_opt_t > "$o"
    struct_body "$new" mem_opt_t > "$n"
    local added
    # `|| true`: diff exits 1 whenever the two struct bodies differ, which is
    # the ordinary case this check exists to report — without the fallback,
    # `set -o pipefail` would make that expected diff abort the whole script
    # before `added` is ever inspected, even though grep/sed both succeed.
    added="$(diff "$o" "$n" | grep '^>' | sed 's/^> //' || true)"
    rm -f "$o" "$n" "$old" "$new"
    if [ -z "$added" ]; then
        printf 'No new fields.\n'
        return 0
    fi
    printf 'New field declarations:\n\n'
    while IFS= read -r field; do
        [ -n "$field" ] || continue
        case "$field" in
            *"NOT a user option"*|*"runtime state"*|*derived*)
                # Backticks below are literal markdown, not command substitution.
                # shellcheck disable=SC2016
                printf -- '- `%s` — skip (upstream marks it internal)\n' "$field" ;;
            *)
                local name
                name="$(printf '%s' "$field" | sed 's/;.*//' | mawk '{print $NF}' | tr -d '*')"
                case "$name" in
                    ''|/*)
                        # Nothing usable to check: an empty extraction (a
                        # whitespace-only or comment-only field line) or a
                        # name starting with "/" (a stray comment fragment
                        # picked up as the last token). Either way, don't
                        # feed it to grep -- flag for a human instead of
                        # guessing.
                        # shellcheck disable=SC2016
                        printf -- '- `%s` — could not extract a field name; check manually\n' "$field" ;;
                    *)
                        # -F (fixed-string): $name is used as a literal
                        # identifier, not a regex. An added array field like
                        # `int mat_ot[25];` extracts as `mat_ot[25]`, which
                        # plain grep parses as "mat_ot" followed by one
                        # character from the class {2,5} -- silently matching
                        # an unrelated `mat_ot2` in opts.rs and misreporting
                        # "already referenced". -F treats it literally.
                        if grep -qF "$name" "$REPO_ROOT/bwa-mem3-rs/src/opts.rs" 2>/dev/null; then
                            # shellcheck disable=SC2016
                            printf -- '- `%s` — already referenced in `opts.rs`\n' "$field"
                        else
                            # shellcheck disable=SC2016
                            printf -- '- `%s` — **consider exposing** on `MemOpts`\n' "$field"
                        fi ;;
                esac ;;
        esac
    done <<< "$added"
    printf '\nExposing a field is optional and separable from the refresh —\n'
    printf 'file it as follow-up rather than growing this PR.\n'
}

# Check 9: an auto-generated re-vendoring commit must not drop attribution.
#
# ext/pdqsort carries no standalone license file — its MIT-style notice is
# embedded in the header comment of pdqsort.h itself (verified against the
# vendored tree) — so that is the path checked, not a license.txt that has
# never existed and would report a false "MISSING" on every single run.
check_licences() {
    report_section "9. Vendored licences"
    local p
    for p in "LICENSE" "ext/sse2neon/LICENSE" "ext/pdqsort/pdqsort.h"; do
        if [ -e "$VENDOR_ABS/$p" ]; then
            if [ "$(old_file "$VENDOR/$p")" = "$(new_file "$VENDOR/$p")" ]; then
                # Backticks below are literal markdown, not command substitution.
                # shellcheck disable=SC2016
                printf -- '- `%s` — present, unchanged\n' "$p"
            else
                # shellcheck disable=SC2016
                printf -- '- `%s` — present, **CHANGED**: read the new terms\n' "$p"
            fi
        else
            # shellcheck disable=SC2016
            printf -- '- `%s` — **MISSING**: attribution lost, do not merge\n' "$p"
        fi
    done
}

# Check 10: what actually changed upstream, per release.
check_release_notes() {
    report_section "10. Upstream release notes"
    local tags t
    tags="${MISSED_TAGS:-}"
    if [ -z "$tags" ]; then
        printf '(no tag list supplied via MISSED_TAGS)\n'
        return 0
    fi
    # shellcheck disable=SC2086 # deliberate word-splitting: MISSED_TAGS is a
    # space-separated list of tag names, not a single path or token.
    for t in $tags; do
        printf '### %s\n\n' "$t"
        gh api "repos/fg-labs/bwa-mem3/releases/tags/$t" --jq '.body' 2>/dev/null \
            | sed 's/^/> /' || printf '> (release notes unavailable)\n'
        printf '\n'
    done
}

# Check 11: what CI structurally cannot verify.
print_manual_checklist() {
    report_section "11. Manual checklist (CI cannot do these)"
    cat <<'EOF'
Green CI is necessary but **not sufficient**. Before merging:

- [ ] Byte-parity against a real index. The PhiX-scale parity targets run in
      CI, but the library-level ones need a large reference:
      `BWA_MEM3_RS_TEST_REF=/path/to/hg38.fa cargo test --workspace`
      covers `align_smoke`, `cli_parity`, `concurrency`, `phase_split`, and
      `shm_lifecycle`.
- [ ] `bwa-mem3-py` (gotcha #10). It is outside the cargo workspace, so
      nothing in CI compiles it. If any `mem_opt_t` field or enum changed,
      check `bwa-mem3-py/src/lib.rs` **and** the hand-written stub
      `bwa-mem3-py/python/bwa_mem3/_bwa_mem3.pyi` — mypy only checks the stub
      for self-consistency and never imports the `.so`, so drift is silent.
- [ ] `shim/bwa_shim_types.h`'s `MEM_F_*` `#define` block, and the three
      coordinated edits in `shim/bwa_shim_layout_assert.cpp` (the
      `pod_flags::MEM_F_X_v` capture, its paired `#undef`, and the
      `BWA_SHIM_CK_FLAG(MEM_F_X)` invocation), match upstream (check 4 above
      diffs it, but confirm you acted on it).
- [ ] `worker_alloc` / `worker_free` copies in `bwa_shim_align.cpp` still
      match upstream's, **including the buffer sizing**, not just field names.
- [ ] Commits are signed (`-S`). The bot's commit is not.
- [ ] Squash the bot commit together with your adaptation commits into a
      logical history before merging.
EOF
}

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

assert_refreshed_tree() {
    local vendor_dirty other_dirty
    # `|| true` for the same reason as the other_dirty pipeline below, but a
    # different failure: `head -1` exits after one line and closes the pipe, so
    # `git status` dies on SIGPIPE and the pipeline reports 141. A refresh
    # rewrites the whole vendor subtree -- thousands of lines, far past the
    # 64 KiB pipe buffer -- so git is guaranteed to still be writing when head
    # goes away. Unguarded, `set -e` then aborts the script at exit 141 before
    # a single line of report is printed, and the workflow reads that as "the
    # refresh failed" rather than "here is the drift".
    vendor_dirty="$(git -C "$REPO_ROOT" status --porcelain -- "bwa-mem3-sys/vendor/" | head -1 || true)"
    if [ -z "$vendor_dirty" ]; then
        echo "ERROR: bwa-mem3-sys/vendor/ is clean — nothing was refreshed." >&2
        echo "       Run scripts/refresh-bwa-mem3.sh <sha> first, and run this" >&2
        echo "       BEFORE committing (the 'before' side is git HEAD)." >&2
        return 1
    fi
    # `|| true`: grep -v exits 1 when EVERY line is filtered out, which is
    # exactly the success case (only bwa-mem3-sys/vendor/ is dirty) — under
    # `pipefail` that would otherwise abort the whole script via `set -e`
    # before the `[ -n "$other_dirty" ]` check below ever runs.
    other_dirty="$(git -C "$REPO_ROOT" status --porcelain | grep -v 'bwa-mem3-sys/vendor/' | head -1 || true)"
    if [ -n "$other_dirty" ]; then
        echo "ERROR: files outside bwa-mem3-sys/vendor/ are dirty: $other_dirty" >&2
        echo "       Commit or stash them so the report reflects the refresh alone." >&2
        return 1
    fi
    return 0
}

# A ~5-line lede printed before the numbered sections, so a human opening a
# two-release bump (potentially thousands of lines once #5's function diffs
# and #10's release notes stack up) gets the verdict first instead of having
# to read linearly to find it. Takes each check's ALREADY-CAPTURED output as
# an argument (main buffers every check before printing anything, precisely
# so this can exist) rather than re-running any check a second time.
#
# Detection is by grepping for the literal "nothing changed" sentence each
# check prints on its own no-drift path — e.g. check_structs's "No change."
# — rather than adding a second reporting channel (return codes, globals) to
# every check function just for this. That keeps the checks themselves as the
# single source of truth for their own text; this function only reads it.
print_summary() {
    local gitmodules_out="$1" build_out="$2" structs_out="$3" flags_out="$4" \
        contracts_out="$5" inventory_out="$6" deps_out="$7" newopts_out="$8" \
        licences_out="$9"

    local gate="not tripped" build_status="passed"
    grep -q '^NEW_SUBMODULE$' <<< "$gitmodules_out" && gate="**TRIPPED** — no PR will be opened"
    grep -q '^Builds clean\.' <<< "$build_out" || build_status="**FAILED**"

    local drifted=()
    grep -qF 'No change.' <<< "$structs_out" || drifted+=("struct layout (§3)")
    grep -q '^### ' <<< "$flags_out" && drifted+=("flags/enums (§4)")
    grep -qF 'No change to any tracked contract.' <<< "$contracts_out" || drifted+=("upstream contracts (§5)")
    { grep -qF 'GONE' <<< "$inventory_out" || ! grep -qF 'No TUs added or removed.' <<< "$inventory_out"; } \
        && drifted+=("source inventory (§6)")
    # Backticks below are literal markdown, not command substitution.
    # shellcheck disable=SC2016
    grep -qF 'No change to `LIBS`' <<< "$deps_out" || drifted+=("build deps (§7)")
    grep -qF 'No new fields.' <<< "$newopts_out" || drifted+=("new mem_opt_t fields (§8)")
    grep -qE 'CHANGED|MISSING' <<< "$licences_out" && drifted+=("licences (§9)")

    report_section "Summary"
    # shellcheck disable=SC2016
    printf -- '- Gate (`NEW_SUBMODULE`): %s\n' "$gate"
    # shellcheck disable=SC2016
    printf -- '- Build (`cargo ci-build`): %s\n' "$build_status"
    if [ "${#drifted[@]}" -eq 0 ]; then
        printf -- '- Drift: none of §3–§9 report a change.\n'
    else
        local joined
        joined="$(printf '%s, ' "${drifted[@]}")"
        joined="${joined%, }"
        printf -- '- Drift (%d/7): %s.\n' "${#drifted[@]}" "$joined"
    fi
    printf -- "- Green CI is **not sufficient** — §11's manual checklist is not run by CI\n"
    printf '  and sits last, after potentially thousands of diff lines below; do not skip it.\n'
}

main() {
    assert_refreshed_tree

    # Every check runs and is captured BEFORE anything is printed, so
    # print_summary (above) can report on all of them without re-running any
    # check a second time. Capturing via command substitution is safe here:
    # every check function's last statement exits 0 on its own (verified
    # against each one), so a captured assignment cannot trip `set -e`
    # spuriously the way an un-guarded `diff`/`grep` mid-function could.
    local gitmodules_out build_out structs_out flags_out contracts_out
    local inventory_out deps_out newopts_out licences_out release_out
    gitmodules_out="$(check_gitmodules)"
    build_out="$(check_build)"
    structs_out="$(check_structs)"
    flags_out="$(check_flags_and_enums)"
    contracts_out="$(check_call_site_contracts)"
    inventory_out="$(check_source_inventory)"
    deps_out="$(check_dependencies)"
    newopts_out="$(check_new_opts)"
    licences_out="$(check_licences)"
    release_out="$(check_release_notes)"

    printf '# bwa-mem3 vendor refresh — drift report\n'
    print_summary "$gitmodules_out" "$build_out" "$structs_out" "$flags_out" \
        "$contracts_out" "$inventory_out" "$deps_out" "$newopts_out" "$licences_out"
    printf '%s\n' "$gitmodules_out"
    printf '%s\n' "$build_out"
    printf '%s\n' "$structs_out"
    printf '%s\n' "$flags_out"
    printf '%s\n' "$contracts_out"
    printf '%s\n' "$inventory_out"
    printf '%s\n' "$deps_out"
    printf '%s\n' "$newopts_out"
    printf '%s\n' "$licences_out"
    printf '%s\n' "$release_out"
    print_manual_checklist
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    main "$@"
fi
