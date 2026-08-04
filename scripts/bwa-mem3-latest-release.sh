#!/usr/bin/env bash
# Report whether fg-labs/bwa-mem3 has a release newer than the vendored
# snapshot, and resolve it to an immutable commit sha.
#
# Usage:
#   scripts/bwa-mem3-latest-release.sh [--tag vX.Y.Z]
#
# --tag forces a specific release and implies force: the script reports on it
# and sets forced=true even when it is older than or equal to what we vendor,
# so an explicit human dispatch is never a silent no-op.
#
# Emits key=value lines on stdout, and appends the same to $GITHUB_OUTPUT when
# that variable is set:
#   current_version  vendored release tag, e.g. v0.6.0
#   current_sha      contents of vendor/COMMIT
#   latest_version    newest non-draft, non-prerelease release tag
#   latest_sha       that tag's commit sha (40 hex)
#   missed_tags      space-separated tags strictly newer than current, oldest first
#   needs_bump       true|false
#   forced           true|false
#
# Requires authenticated `gh` and `jq`.
set -euo pipefail

UPSTREAM_REPO="fg-labs/bwa-mem3"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
COMMIT_FILE="$REPO_ROOT/bwa-mem3-sys/vendor/COMMIT"
VERSION_TXT="$REPO_ROOT/bwa-mem3-sys/vendor/bwa-mem3/version.txt"

# Strip a leading v and any pre-release suffix, leaving X.Y.Z.
normalize_version() {
    printf '%s' "${1#v}" | sed 's/-.*$//'
}

# Exit 0 iff $1 > $2, comparing X.Y.Z numerically field by field.
#
# Not `sort -V`: it is not a semver comparator, BSD and GNU disagree, and on
# macOS `printf 'v0.1.0-pre\nv0.1.0\n' | sort -V | tail -1` yields the
# prerelease. Equal versions are NOT greater — callers depend on being able to
# distinguish `<` from `==`.
version_gt() {
    local a b
    a="$(normalize_version "$1")"
    b="$(normalize_version "$2")"
    local ia ib
    for i in 1 2 3; do
        ia="$(printf '%s' "$a" | cut -d. -f"$i")"
        ib="$(printf '%s' "$b" | cut -d. -f"$i")"
        ia="${ia:-0}"
        ib="${ib:-0}"
        if [ "$ia" -gt "$ib" ]; then return 0; fi
        if [ "$ia" -lt "$ib" ]; then return 1; fi
    done
    return 1
}

# Print the 40-hex commit sha a tag points at, dereferencing annotated tags.
#
# Never uses a release's target_commitish: upstream's v0.2.0 reports "main",
# and a mutable ref in vendor/COMMIT would make the snapshot
# non-reproducible (check.yml uses that value as a fetch refspec).
resolve_tag_sha() {
    local tag="$1" json type sha
    json="$(gh api "repos/$UPSTREAM_REPO/git/ref/tags/$tag")"
    type="$(printf '%s' "$json" | jq -r '.object.type')"
    sha="$(printf '%s' "$json" | jq -r '.object.sha')"
    if [ "$type" = "tag" ]; then
        sha="$(gh api "repos/$UPSTREAM_REPO/git/tags/$sha" --jq '.object.sha')"
    fi
    if ! printf '%s' "$sha" | grep -Eq '^[0-9a-f]{40}$'; then
        echo "ERROR: tag $tag did not resolve to a 40-hex sha (got '$sha')" >&2
        return 1
    fi
    printf '%s' "$sha"
}

# All non-draft, non-prerelease release tags, newest first.
release_tags() {
    gh api --paginate "repos/$UPSTREAM_REPO/releases" \
        --jq '.[] | select(.draft == false and .prerelease == false) | .tag_name'
}

# Look up which release tag (if any) vendor/COMMIT currently points at.
# Optional $1 is a pre-fetched `release_tags` list (newest first); omit it to
# have this function fetch its own, which is what makes it usable standalone.
#
# Reverse-looks-up vendor/COMMIT in the release list rather than trusting a
# marker file (vendor/bwa-mem3/version.txt): that also catches an upstream
# RETAG, where version.txt still says 0.6.0 but v0.6.0 now points somewhere
# else. Prints nothing (and returns non-zero) if no release tag matches.
vendored_version() {
    local tags="${1:-}" current_sha t
    [ -f "$COMMIT_FILE" ] || { echo "ERROR: missing $COMMIT_FILE" >&2; return 1; }
    current_sha="$(tr -d '[:space:]' < "$COMMIT_FILE")"

    [ -n "$tags" ] || tags="$(release_tags)"
    [ -n "$tags" ] || { echo "ERROR: no releases found for $UPSTREAM_REPO" >&2; return 1; }

    while IFS= read -r t; do
        if [ "$(resolve_tag_sha "$t")" = "$current_sha" ]; then
            printf '%s' "$t"
            return 0
        fi
    done <<< "$tags"
    return 1
}

main() {
    local forced_tag=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --tag) forced_tag="${2:?--tag needs a value}"; shift 2 ;;
            *) echo "ERROR: unknown argument: $1" >&2; exit 2 ;;
        esac
    done

    [ -f "$COMMIT_FILE" ] || { echo "ERROR: missing $COMMIT_FILE" >&2; exit 1; }
    local current_sha
    current_sha="$(tr -d '[:space:]' < "$COMMIT_FILE")"

    local tags
    tags="$(release_tags)"
    [ -n "$tags" ] || { echo "ERROR: no releases found for $UPSTREAM_REPO" >&2; exit 1; }

    # Reverse-look-up vendor/COMMIT in the release list. Stronger than reading
    # a version marker out of the tree: it also catches an upstream RETAG,
    # where version.txt still says 0.6.0 but v0.6.0 now points somewhere else.
    local current_version=""
    if ! current_version="$(vendored_version "$tags")"; then
        local marker="unknown"
        [ -f "$VERSION_TXT" ] && marker="$(tr -d '[:space:]' < "$VERSION_TXT")"
        echo "ERROR: vendor/COMMIT ($current_sha) matches no release tag." >&2
        echo "       The vendored tree's version.txt says '$marker'." >&2
        echo "       Either a main tip was vendored (documented in CONTRIBUTING)," >&2
        echo "       or upstream retagged/deleted a release. Pass --tag to" >&2
        echo "       proceed against a specific release deliberately." >&2
        [ -n "$forced_tag" ] || exit 1
    fi

    local latest_version latest_sha needs_bump forced missed
    # `|| true` for uniformity with the two `| head -1` sites in
    # bwa-mem3-drift-report.sh, where the producer outruns the pipe buffer and
    # dies on SIGPIPE when head exits. It cannot fire here -- $tags is at most
    # a few KB of release tags, so printf finishes writing long before head
    # goes away -- but the shape is identical and a future caller feeding this
    # a larger list should not have to rediscover the trap.
    latest_version="$(printf '%s\n' "$tags" | head -1 || true)"
    forced=false
    if [ -n "$forced_tag" ]; then
        latest_version="$forced_tag"
        forced=true
    fi
    latest_sha="$(resolve_tag_sha "$latest_version")"

    if [ "$forced" = true ]; then
        needs_bump=true
    elif [ -n "$current_version" ] && version_gt "$latest_version" "$current_version"; then
        needs_bump=true
    else
        needs_bump=false
    fi

    # Tags strictly newer than current_version, oldest first. $tags is
    # newest-first (as returned by release_tags/the GitHub API), so walk a
    # bash array back-to-front rather than shelling out to a reverse utility.
    # An earlier draft picked between BSD `tail -r` (macOS) and GNU `tac`
    # (Linux CI) with `cmd1 2>/dev/null || cmd2`: fragile, because it silently
    # swallows tail's stderr and falls through on ANY failure, not just
    # "flag unsupported" — a transient error would reverse-sort garbage
    # without complaint. Iterating the array we already have in bash needs
    # neither external tool nor a guess about which one is present.
    missed=""
    if [ -n "$current_version" ]; then
        local -a tag_array
        mapfile -t tag_array <<< "$tags"
        local idx
        for ((idx = ${#tag_array[@]} - 1; idx >= 0; idx--)); do
            t="${tag_array[idx]}"
            if version_gt "$t" "$current_version"; then
                missed="${missed:+$missed }$t"
            fi
        done
    fi

    local out
    out=$(cat <<EOF
current_version=$current_version
current_sha=$current_sha
latest_version=$latest_version
latest_sha=$latest_sha
missed_tags=$missed
needs_bump=$needs_bump
forced=$forced
EOF
)
    printf '%s\n' "$out"
    [ -n "${GITHUB_OUTPUT:-}" ] && printf '%s\n' "$out" >> "$GITHUB_OUTPUT"
    return 0
}

# Only run main when executed, not when sourced by the unit tests.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    main "$@"
fi
