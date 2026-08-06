#!/usr/bin/env bash
# Refresh the vendored bwa-mem3 snapshot to a given commit on fg-labs/bwa-mem3.
# Default branch vendored is `main` (our integration branch).
#
# Usage: scripts/refresh-bwa-mem3.sh <commit-hash> [source-path]
#   - If source-path is provided, the vendor is copied from a local working tree
#     (useful when the commit is not yet pushed). Otherwise, a fresh clone is made.
set -euo pipefail
HASH="${1:?usage: refresh-bwa-mem3.sh <commit-hash> [source-path]}"
SRC="${2:-}"

# A release's target_commitish is not necessarily a commit sha — upstream's
# v0.2.0 release reports "main". Writing that into vendor/COMMIT would give a
# mutable, non-reproducible marker, and check.yml uses that value as a fetch
# refspec, so the e2e job would build a different upstream than we vendor.
case "$HASH" in
    [0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]*)
        if [ "${#HASH}" -ne 40 ] || [ -n "$(printf '%s' "$HASH" | tr -d '0-9a-f')" ]; then
            echo "ERROR: HASH must be a full 40-hex commit sha, got: $HASH" >&2
            exit 1
        fi
        ;;
    *)
        echo "ERROR: HASH must be a full 40-hex commit sha, got: $HASH" >&2
        exit 1
        ;;
esac

REPO_URL="https://github.com/fg-labs/bwa-mem3.git"
DST="bwa-mem3-sys/vendor/bwa-mem3"

rm -rf "$DST"
mkdir -p "$(dirname "$DST")"

# Subtrees to prune, read from the shared list so the drift report and this
# script cannot disagree about what a pruned tree looks like.
DROP_LIST="$(dirname "$0")/vendor-drop-subtrees.txt"
[ -f "$DROP_LIST" ] || { echo "ERROR: missing $DROP_LIST" >&2; exit 1; }
DROP_SUBTREES=()
while IFS= read -r line; do
    line="${line%%#*}"                       # strip comments
    line="$(printf '%s' "$line" | tr -d '[:space:]')"
    [ -n "$line" ] && DROP_SUBTREES+=("$line")
done < "$DROP_LIST"
[ "${#DROP_SUBTREES[@]}" -gt 0 ] || { echo "ERROR: $DROP_LIST has no entries" >&2; exit 1; }

if [ -n "$SRC" ]; then
    # Copy from a local working tree at the given hash (submodules initialized).
    actual=$(git -C "$SRC" rev-parse HEAD)
    if [ "$actual" != "$HASH" ] && [ "${actual:0:${#HASH}}" != "$HASH" ]; then
        echo "ERROR: $SRC HEAD is $actual, not $HASH" >&2
        exit 1
    fi
    # rsync preserves submodule contents; exclude vcs, build artifacts, large tests.
    mkdir -p "$DST"
    rsync_excludes=(
        --exclude='.git' --exclude='.github'
        --exclude='*.o' --exclude='*.a' --exclude='obj/' --exclude='Debug/' --exclude='objtest/'
        --exclude='/bwa-mem3'
        --exclude='/bwa-mem3.sse41' --exclude='/bwa-mem3.sse42'
        --exclude='/bwa-mem3.avx' --exclude='/bwa-mem3.avx2' --exclude='/bwa-mem3.avx512bw'
        --exclude='/bwa-mem3.arm64' --exclude='/bwa-mem3.native'
    )
    for sub in "${DROP_SUBTREES[@]}"; do
        rsync_excludes+=(--exclude="/${sub}")
    done
    rsync -a "${rsync_excludes[@]}" "$SRC/" "$DST/"
else
    git clone --recurse-submodules "$REPO_URL" "$DST"
    git -C "$DST" checkout "$HASH"
    actual=$(git -C "$DST" rev-parse HEAD)
    if [ "$actual" != "$HASH" ]; then
        echo "ERROR: clone HEAD is $actual, not $HASH" >&2
        exit 1
    fi
    git -C "$DST" submodule update --init --recursive
    find "$DST" -name '.git' -print0 | xargs -0 rm -rf
    # Strip nested .github dirs (e.g. ext/sse2neon/.github); the loop
    # below only removes the top-level entry. The rsync path handles
    # both via the un-anchored `--exclude='.github'` above.
    find "$DST" -name '.github' -type d -print0 | xargs -0 rm -rf
    for sub in "${DROP_SUBTREES[@]}"; do
        rm -rf "${DST:?}/$sub"
    done
fi

# Upstream's own .gitattributes would let a future upstream linguist-* rule
# override the root .gitattributes that collapses this tree in PR diffs.
rm -f "$DST/.gitattributes"

# Record the commit; used by build.rs as a sanity marker.
echo "$HASH" > bwa-mem3-sys/vendor/COMMIT

# Verify MATE_SORT=0 default (see spec "Pairing without patching upstream").
if ! grep -E '^CPPFLAGS\+?=.*-DMATE_SORT=0' "$DST/Makefile" >/dev/null; then
    echo "ERROR: vendored Makefile does not have -DMATE_SORT=0; shim semantics would diverge." >&2
    exit 1
fi

echo "Vendored $HASH at $DST"
