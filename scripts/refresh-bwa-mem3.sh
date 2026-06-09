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
REPO_URL="https://github.com/fg-labs/bwa-mem3.git"
DST="bwa-mem3-sys/vendor/bwa-mem3"

rm -rf "$DST"
mkdir -p "$(dirname "$DST")"

# Vendor-only subtrees: dropped to keep the snapshot small. We don't compile
# the upstream CLI's BAM writer (htslib), the index builder (libsais), the
# allocator override (mimalloc), the unit-test harness (doctest), benches,
# CI configs, or upstream's docs/scripts. Keep `ext/sse2neon` (headers we
# compile against) and `ext/pdqsort` (the header-only sort dropped in at
# bwamem.cpp sort sites).
DROP_SUBTREES=(
    .github bench docs scripts
    ext/htslib ext/libsais ext/mimalloc ext/doctest
)

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
    git -C "$DST" submodule update --init --recursive
    find "$DST" -name '.git' -print0 | xargs -0 rm -rf
    # Strip nested .github dirs (e.g. ext/sse2neon/.github); the loop
    # below only removes the top-level entry. The rsync path handles
    # both via the un-anchored `--exclude='.github'` above.
    find "$DST" -name '.github' -type d -print0 | xargs -0 rm -rf
    for sub in "${DROP_SUBTREES[@]}"; do
        rm -rf "$DST/$sub"
    done
fi

# Record the commit; used by build.rs as a sanity marker.
echo "$HASH" > bwa-mem3-sys/vendor/COMMIT

# Verify MATE_SORT=0 default (see spec "Pairing without patching upstream").
if ! grep -E '^CPPFLAGS\+?=.*-DMATE_SORT=0' "$DST/Makefile" >/dev/null; then
    echo "ERROR: vendored Makefile does not have -DMATE_SORT=0; shim semantics would diverge." >&2
    exit 1
fi

echo "Vendored $HASH at $DST"
