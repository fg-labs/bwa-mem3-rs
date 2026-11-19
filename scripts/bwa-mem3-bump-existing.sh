#!/usr/bin/env bash
# Source-only: vendor-bump.yml's "Check for an existing PR or issue" step
# sources this file; running it does nothing. It holds that step's decisions
# about what already exists for a release, as pure functions with no network
# access, so scripts/tests/test-bump-existing.sh can unit-test them. The step
# does the lookups (gh, git ls-remote) and passes the results in.

# forced_verdict <ls-remote rc> <issue lookup ok: 1|0> <vendored sha> <target sha>
#
# Whether a forced (`-f tag=`) run may refresh to the target release and open
# its branch and PR. Prints one of:
#   proceed         go ahead, even when a closed or merged PR exists for the
#                   branch: closing the PR and deleting the branch is the
#                   documented way to retry a release (CONTRIBUTING.md), and a
#                   merged release whose bump was reverted can be re-forced.
#   vendored        the target is already the vendored snapshot, so a refresh
#                   would change nothing and the drift report would fail.
#   branch          the branch still exists. Pushing fresh history to it would
#                   fail as a non-fast-forward, or need a force push that
#                   orphans review comments on its commits. (An open PR always
#                   has a live branch: deleting the branch closes the PR.)
#   branch-unknown  the branch lookup failed. Fails closed: a run that cannot
#                   tell whether the branch exists must not try to create it.
#   issue-unknown   the open-issue lookup failed, so the run cannot tell which
#                   tracking issue to reuse and would file a duplicate.
# `git ls-remote --exit-code` returns 0 when the ref exists, 2 when it does
# not, and anything else on error.
forced_verdict() {
    local ls_rc="$1" issue_ok="$2" vendored="$3" target="$4"
    if [ -n "$vendored" ] && [ "$vendored" = "$target" ]; then
        echo vendored
    elif [ "$ls_rc" -eq 0 ]; then
        echo branch
    elif [ "$ls_rc" -ne 2 ]; then
        echo branch-unknown
    elif [ "$issue_ok" != 1 ]; then
        echo issue-unknown
    else
        echo proceed
    fi
}

# forced_message <verdict> <tag> <branch> <PRs found, e.g. "98:CLOSED", or ""> <ls-remote rc>
#
# The step-summary text for a forced run forced_verdict refused, with the
# remedy that fits that verdict. Returns 1 for a verdict it does not know.
forced_message() {
    local verdict="$1" tag="$2" branch="$3" found="$4" ls_rc="$5"
    case "$verdict" in
        vendored)
            echo "$tag is already the vendored snapshot (bwa-mem3-sys/vendor/COMMIT),"
            echo "so a forced refresh to it would change nothing."
            ;;
        branch)
            echo "Branch $branch still exists (PRs for it: ${found:-none})."
            echo "Force-dispatching --tag=$tag would either fail to push"
            echo "(non-fast-forward) or require force-pushing over the branch,"
            echo "orphaning any review comments on its commits. Close its PR if"
            echo "one is open, delete branch $branch, then re-dispatch."
            ;;
        branch-unknown)
            echo "Could not tell whether branch $branch exists"
            echo "(git ls-remote exited $ls_rc). See the step log, then"
            echo "re-dispatch once the lookup works."
            ;;
        issue-unknown)
            echo "The open tracking-issue lookup failed, so this run could not"
            echo "tell which issue to reuse and might file a duplicate. See the"
            echo "step log, then re-dispatch once the lookup works."
            ;;
        *)
            echo "forced_message: unknown verdict '$verdict'" >&2
            return 1
            ;;
    esac
    echo "Nothing was created or pushed by this run."
}

# The jq filter the step runs over the open-issues API response to produce the
# lines issue_number_for_title reads. `@tsv` escapes a line break, tab or
# backslash inside a field, so every issue is exactly one line: GitHub allows a
# line break in a title, and an unescaped one would let anyone who can open an
# issue forge a second "<number>\t<login>\t<title>" record naming the bot, and
# steer a forced run's report onto an issue of their choosing.
# shellcheck disable=SC2034 # read by vendor-bump.yml and the tests
ISSUE_LINES_JQ='.[] | select(.pull_request == null) | [.number, .user.login, .title] | @tsv'

# issue_number_for_title <title> [author]
#
# Reads "<number>\t<login>\t<title>" lines (one open issue each, as
# ISSUE_LINES_JQ emits them) on stdin and
# prints the number of the first whose title equals <title> exactly, and whose
# author is <author> when one is given, or nothing. Exact: not a prefix, a
# pattern or a whitespace-trimmed match, so "...v0.18.0" never matches
# "...v0.18.0-rc1". Each line is split on its first two tabs only, so a title
# keeps any tab or surrounding space it has, and a last line without a
# trailing newline is still read.
issue_number_for_title() {
    local want="$1" author="${2:-}" line number login title
    while IFS= read -r line || [ -n "$line" ]; do
        case "$line" in
            *$'\t'*$'\t'*) ;;
            *) continue ;; # not "<number>\t<login>\t<title>": skip it
        esac
        number="${line%%$'\t'*}"
        line="${line#*$'\t'}"
        login="${line%%$'\t'*}"
        title="${line#*$'\t'}"
        if [ "$title" = "$want" ] && { [ -z "$author" ] || [ "$login" = "$author" ]; }; then
            printf '%s\n' "$number"
            return 0
        fi
    done
    return 0
}
