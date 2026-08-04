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

report_section() { printf '\n## %s\n\n' "$1"; }

# The committed (pre-refresh) content of a vendored path, or empty if it is new.
old_file() {
    git -C "$REPO_ROOT" show "HEAD:$1" 2>/dev/null || true
}

# The working-tree (post-refresh) content of a path, or empty if it is gone.
new_file() {
    cat "$REPO_ROOT/$1" 2>/dev/null || true
}

# Submodule paths present in the new .gitmodules but not the old one.
#
# A new submodule means refresh-bwa-mem3.sh could not know to prune it, so the
# clone kept it and a commit would carry its whole source tree. The caller uses
# this to gate PR creation.
#
# Extracts the path value directly with grep -oE rather than trimming
# whitespace afterward with `tr -d '[:space:]'`: that would delete the
# newlines between grep's output lines too, gluing multiple submodule paths
# into one bogus token (e.g. "ext/aext/b") and hiding real drift. Uses
# POSIX [[:space:]] instead of the GNU-only \s shorthand, which BSD sed
# (macOS) does not expand and silently leaves as a literal "s".
gitmodules_new_submodules() {
    local old="$1" new="$2"
    local old_paths new_paths
    old_paths="$(grep -oE '^[[:space:]]*path[[:space:]]*=[[:space:]]*[^[:space:]]+' "$old" 2>/dev/null | sed -E 's/^[[:space:]]*path[[:space:]]*=[[:space:]]*//' | sort -u || true)"
    new_paths="$(grep -oE '^[[:space:]]*path[[:space:]]*=[[:space:]]*[^[:space:]]+' "$new" 2>/dev/null | sed -E 's/^[[:space:]]*path[[:space:]]*=[[:space:]]*//' | sort -u || true)"
    comm -13 <(printf '%s\n' "$old_paths") <(printf '%s\n' "$new_paths") | grep -v '^$' || true
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
    grep -E '^\s*#define\s+MEM_F_[A-Z0-9_]+\s+' "$1" 2>/dev/null \
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

# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------

# Check 1 (primary): does the refreshed tree still build?
check_build() {
    report_section "1. Build (\`cargo ci-build\`)"
    local log status
    log="$(mktemp)"
    status=0
    (cd "$REPO_ROOT" && cargo ci-build) > "$log" 2>&1 || status=$?
    if [ "$status" -eq 0 ]; then
        printf 'Builds clean. Note this proves compilation only — it does not\n'
        printf 'prove byte-parity with the CLI. See the manual checklist below.\n'
    else
        printf 'Build FAILED (exit %s). These errors are the worklist:\n\n' "$status"
        printf '```\n'
        # `: error:` (unanchored) alongside the anchored `^error`: cargo's own
        # "failed to run custom build command" wrapper matches ^error, but the
        # actual C/C++ diagnostic that matters after a vendor refresh is a
        # couple of lines further down, indented under "--- stderr" as
        # `file:line:col: error: ...` — clang/g++'s mid-line form, which
        # ^error alone would never match, silently hiding the one line a
        # human needs behind the useless wrapper. `ld:`/`Undefined symbols`
        # catch the equivalent linker-error case (e.g. a renamed/removed
        # exported symbol after a refresh); still not exhaustive of every
        # compiler/linker's diagnostic format — that gap is tracked
        # separately, not solved by this alternation.
        grep -E '(^error|: error:|^warning: unused|  -->|ld:|^Undefined symbols)' "$log" | head -40 || tail -40 "$log"
        printf '```\n\n'
        printf 'Full log is in the workflow run artifacts.\n'
    fi
    rm -f "$log"
}

# Check 2 (hard gate): a new submodule makes a refresh commit wrong.
check_gitmodules() {
    local old new
    old="$(mktemp)"; new="$(mktemp)"
    old_file "$VENDOR/.gitmodules" > "$old"
    new_file "$VENDOR/.gitmodules" > "$new"
    local added
    added="$(gitmodules_new_submodules "$old" "$new")"
    rm -f "$old" "$new"

    report_section "2. Submodules (\`.gitmodules\`)"
    if [ -z "$added" ]; then
        printf 'No new submodules.\n'
        return 0
    fi
    # NEW_SUBMODULE is the machine-readable marker the workflow greps for to
    # gate PR creation. It is NOT the first line of the report or even of
    # this check's own output — report_section above already emitted the
    # "## 2. Submodules" heading. Consumers must match it anchored against
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
    old="$(mktemp)"; new="$(mktemp)"
    old_file "$VENDOR/src/bwamem.h" > "$old"
    new_file "$VENDOR/src/bwamem.h" > "$new"
    local any=0
    for s in mem_opt_t mem_pestat_t; do
        o="$(mktemp)"; n="$(mktemp)"
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
    old="$(mktemp)"; new="$(mktemp)"
    old_file "$VENDOR/src/bwamem.h" > "$old"
    new_file "$VENDOR/src/bwamem.h" > "$new"

    o="$(mktemp)"; n="$(mktemp)"
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

    o="$(mktemp)"; n="$(mktemp)"
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

main() {
    assert_refreshed_tree
    printf '# bwa-mem3 vendor refresh — drift report\n'
    check_gitmodules
    check_build
    check_structs
    check_flags_and_enums
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    main "$@"
fi
