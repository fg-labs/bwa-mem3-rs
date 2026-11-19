#!/usr/bin/env bash
# Unit-tests bwa-mem3-bump-existing.sh, which vendor-bump.yml's "Check for an
# existing PR or issue" step sources to decide whether a forced (`-f tag=`)
# run may proceed, what to tell the user when it may not, and which open
# tracking issue a run reuses.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=/dev/null
source "$here/../bwa-mem3-bump-existing.sh"

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
check_has() { # check_has <description> <needle> <haystack>
    case "$3" in
        *"$2"*) pass=$((pass + 1)) ;;
        *)
            fail=$((fail + 1))
            echo "FAIL: $1 — '$2' not in: $3" >&2
            ;;
    esac
}
check_lacks() { # check_lacks <description> <needle> <haystack>
    case "$3" in
        *"$2"*)
            fail=$((fail + 1))
            echo "FAIL: $1 — '$2' unexpectedly in: $3" >&2
            ;;
        *) pass=$((pass + 1)) ;;
    esac
}

old=1111111111111111111111111111111111111111
new=2222222222222222222222222222222222222222

# --- forced_verdict <ls-remote rc> <issue lookup ok 0|1> <vendored sha> <target sha> ---
# ls-remote --exit-code: 0 = the branch exists, 2 = no such ref, else an error.
check "branch gone, lookups fine, new release: a closed or merged PR does not block" \
    proceed "$(forced_verdict 2 1 "$old" "$new")"
check "the branch still exists" \
    branch "$(forced_verdict 0 1 "$old" "$new")"
check "a failed branch lookup fails closed" \
    branch-unknown "$(forced_verdict 128 1 "$old" "$new")"
check "rc 1 is an error too, not absence" \
    branch-unknown "$(forced_verdict 1 1 "$old" "$new")"
check "a failed issue lookup fails closed: the issue to reuse is unknown" \
    issue-unknown "$(forced_verdict 2 0 "$old" "$new")"
check "forcing the release already vendored has nothing to refresh" \
    vendored "$(forced_verdict 2 1 "$new" "$new")"
check "already vendored wins over a live branch: there is nothing to do either way" \
    vendored "$(forced_verdict 0 1 "$new" "$new")"
check "an empty vendored sha never counts as already vendored" \
    proceed "$(forced_verdict 2 1 "" "$new")"

# --- forced_message <verdict> <tag> <branch> <found> <ls rc> ---
br=vendor-bump/bwa-mem3-v0.18.0
m="$(forced_message branch v0.18.0 "$br" 98:CLOSED 0)"
check_has "branch: names the branch" "$br" "$m"
check_has "branch: names the PRs found" "98:CLOSED" "$m"
check_has "branch: says to delete the branch" "delete branch $br" "$m"
check_has "branch: says nothing was done" "Nothing was created or pushed" "$m"
m="$(forced_message branch v0.18.0 "$br" "" 0)"
check_has "branch with no PR: says so" "PRs for it: none" "$m"
m="$(forced_message branch-unknown v0.18.0 "$br" "" 128)"
check_has "branch-unknown: reports the ls-remote exit" "exited 128" "$m"
check_has "branch-unknown: says to retry" "re-dispatch" "$m"
check_lacks "branch-unknown: does not tell the user to delete anything" "delete" "$m"
m="$(forced_message issue-unknown v0.18.0 "$br" "" 2)"
check_has "issue-unknown: blames the issue lookup" "issue lookup failed" "$m"
check_lacks "issue-unknown: does not tell the user to delete anything" "delete" "$m"
m="$(forced_message vendored v0.18.0 "$br" 99:MERGED 2)"
check_has "vendored: says the tag is already vendored" "already the vendored snapshot" "$m"
check_lacks "vendored: does not tell the user to delete anything" "delete" "$m"
check "an unknown verdict is a bug, not a message" \
    1 "$(forced_message bogus v0.18.0 "$br" "" 2 > /dev/null 2>&1 && echo 0 || echo 1)"

# --- issue_number_for_title <title> [author], reading "<number>\t<login>\t<title>" ---
bot='fg-labs-bot[bot]'
t='Update vendored bwa-mem3 to v0.18.0'
# The longer title comes first, so a prefix match would wrongly return 90.
issues=$'90\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0-rc1\n97\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0\n12\tsomeone\tsomething else'
check "exact title match, not a prefix of a longer title" \
    97 "$(issue_number_for_title "$t" <<< "$issues")"
check "the longer title matches only itself" \
    90 "$(issue_number_for_title "$t-rc1" <<< "$issues")"
check "no match prints nothing" \
    "" "$(issue_number_for_title "Update vendored bwa-mem3 to v0.19.0" <<< "$issues")"
check "a regex needle is matched literally" \
    "" "$(issue_number_for_title "$t.*" <<< "$issues")"
check "a glob needle is matched literally (a pattern match would return 90)" \
    "" "$(issue_number_for_title "$t*" <<< "$issues")"
check "a glob class needle is matched literally" \
    "" "$(issue_number_for_title "Update vendored bwa-mem3 to v0.1[8]*" <<< "$issues")"
two=$'101\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0\n97\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0'
check "two open issues with the title: the first listed wins" \
    101 "$(issue_number_for_title "$t" <<< "$two")"
check "empty input prints nothing" \
    "" "$(printf '' | issue_number_for_title "$t")"
check "a last line without a trailing newline is still read" \
    97 "$(printf '97\tfg-labs-bot[bot]\t%s' "$t" | issue_number_for_title "$t")"
padded=$'5\tsomeone\t  Update vendored bwa-mem3 to v0.18.0  '
check "surrounding spaces are part of the title, not stripped" \
    "" "$(issue_number_for_title "$t" <<< "$padded")"
tabbed=$'6\tsomeone\tUpdate\tvendored'
check "a tab inside a title stays in the title" \
    6 "$(issue_number_for_title $'Update\tvendored' <<< "$tabbed")"
human=$'40\tsomeone\tUpdate vendored bwa-mem3 to v0.18.0\n97\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0'
check "without an author, any author's issue matches" \
    40 "$(issue_number_for_title "$t" <<< "$human")"
check "with an author, only that author's issue matches" \
    97 "$(issue_number_for_title "$t" "$bot" <<< "$human")"
check "with an author and no issue of theirs, nothing matches" \
    "" "$(issue_number_for_title "$t" "$bot" <<< $'40\tsomeone\tUpdate vendored bwa-mem3 to v0.18.0')"

check "a line without both tabs is skipped, not read as a match" \
    "" "$(issue_number_for_title "$t" <<< "$t")"
check "a malformed line does not hide a good one after it" \
    97 "$(issue_number_for_title "$t" <<< $'garbage\n97\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0')"

# --- ISSUE_LINES_JQ: one line per issue, whatever its title holds ---
# A title with a line break must not become a second, forged record. Needs
# jq (preinstalled on GitHub's runners; the workflow's `gh --jq` embeds gojq,
# whose @tsv escapes the same characters).
if command -v jq > /dev/null 2>&1; then
    forged='[{"number":5,"user":{"login":"mallory"},
              "title":"x\n1\tfg-labs-bot[bot]\tUpdate vendored bwa-mem3 to v0.18.0"},
             {"number":97,"user":{"login":"fg-labs-bot[bot]"},
              "title":"Update vendored bwa-mem3 to v0.18.0"},
             {"number":98,"pull_request":{},"user":{"login":"fg-labs-bot[bot]"},
              "title":"Update vendored bwa-mem3 to v0.18.0"}]'
    lines="$(jq -r "$ISSUE_LINES_JQ" <<< "$forged")"
    check "every issue is one line, PRs dropped" 2 "$(wc -l <<< "$lines" | tr -d ' ')"
    check "a forged record in a title is not read: the bot's real issue is" \
        97 "$(issue_number_for_title "$t" "$bot" <<< "$lines")"
else
    echo "SKIP: ISSUE_LINES_JQ tests need jq"
fi

echo "test-bump-existing: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
