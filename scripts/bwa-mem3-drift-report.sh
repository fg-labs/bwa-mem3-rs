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
        # human needs behind the useless wrapper.
        grep -E '(^error|: error:|^warning: unused|  -->)' "$log" | head -40 || tail -40 "$log"
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
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    main "$@"
fi
