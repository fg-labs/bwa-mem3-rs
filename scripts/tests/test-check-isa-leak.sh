#!/usr/bin/env bash
# Unit-tests check-isa-leak.sh against tiny x86_64 binaries built here. Each fixture puts one
# AVX-512 (or AVX2) instruction in a function via inline asm, so the result does not depend on
# what the compiler chooses to vectorize. Needs an x86_64 Linux host with cc and binutils;
# skipped elsewhere.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
check="$here/../check-isa-leak.sh"

if [ "$(uname -s)" != Linux ] || [ "$(uname -m)" != x86_64 ] \
    || ! command -v cc objdump nm strip mawk > /dev/null \
    || ! objdump --version | head -1 | grep -q 'GNU objdump'; then
    echo "SKIP: check-isa-leak tests need an x86_64 Linux host with cc, GNU binutils and mawk"
    exit 0
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
pass=0
fail=0

# build <name> <C source>: compile a fixture binary into $work/<name>.
build() {
    printf '%s\n' "$2" > "$work/$1.c"
    cc -O1 -o "$work/$1" "$work/$1.c"
}

# expect <description> <status> <binary> [text the output must contain]
expect() {
    local out status=0
    out=$("$check" "$3" 2>&1) || status=$?
    if [ "$status" -ne "$2" ]; then
        fail=$((fail + 1))
        echo "FAIL: $1 — expected exit $2, got $status: $out" >&2
    elif [ $# -ge 4 ] && [[ "$out" != *"$4"* ]]; then
        fail=$((fail + 1))
        echo "FAIL: $1 — output lacks '$4': $out" >&2
    else
        pass=$((pass + 1))
    fi
}

# A weak (COMDAT-like) helper with AVX-512 code is exactly the leak: every tier can bind to it.
build weak_leak '
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { shared_helper(); return 0; }'
expect "weak symbol with AVX-512 code is a leak" 1 "$work/weak_leak" "shared_helper"

# So is a global one outside the tier.
build global_leak '
__attribute__((noinline)) void fill_assign(void) { __asm__ volatile("vmovdqu16 %%ymm1, %%ymm2" ::: "xmm1", "xmm2"); }
int main(void) { fill_assign(); return 0; }'
expect "global non-tier symbol with an EVEX-only mnemonic is a leak" 1 "$work/global_leak" "fill_assign"

# The avx512bw tier itself may use AVX-512.
build tier_ok '
__attribute__((noinline)) void ksw_extend2_avx512bw(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { ksw_extend2_avx512bw(); return 0; }'
expect "avx512bw tier symbol is allowed" 0 "$work/tier_ok"

# A local (static) function cannot be bound to from outside its TU.
build static_ok '
static __attribute__((noinline)) void pack(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { pack(); return 0; }'
expect "static function is allowed" 0 "$work/static_ok"

# AVX2 code is not AVX-512.
build avx2_ok '
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile("vpxor %%ymm1, %%ymm1, %%ymm1" ::: "xmm1"); }
int main(void) { shared_helper(); return 0; }'
expect "AVX2 code is not a leak" 0 "$work/avx2_ok"

# AVX-512VL on low-numbered registers: no zmm, mask or xmm/ymm16-31, but still EVEX-encoded.
build vl_leak '
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile("vpxord %%ymm1, %%ymm1, %%ymm1" ::: "xmm1"); }
int main(void) { shared_helper(); return 0; }'
expect "EVEX-encoded AVX-512VL on low registers is a leak" 1 "$work/vl_leak" "shared_helper"

# Only the avx512bw tier's naming convention is exempt, not any name containing "avx512".
build name_leak '
__attribute__((weak, noinline)) void shared_avx512_helper(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { shared_avx512_helper(); return 0; }'
expect "a helper merely named avx512 is not exempt" 1 "$work/name_leak" "shared_avx512_helper"

# A stripped binary has no symbols to check: that is misuse, not a pass.
cp "$work/weak_leak" "$work/stripped"
strip "$work/stripped"
expect "stripped binary is misuse" 2 "$work/stripped" "no function symbols"

# The VEX-encoded opmask instructions are AVX-512 too, though no 0x62 byte marks them.
build kmask_leak '
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile("kmovw %%k1, %%k2" ::: ); }
int main(void) { shared_helper(); return 0; }'
expect "opmask (k-register) instruction is a leak" 1 "$work/kmask_leak" "shared_helper"

# An objdump whose output the parser cannot read (llvm-objdump's layout differs) must not pass.
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
printf '#!/bin/sh\n[ "$1" = --version ] && echo "LLVM version 21" ; exit 0\n' > "$work/llvm-objdump"
chmod +x "$work/llvm-objdump"
out=$(OBJDUMP="$work/llvm-objdump" "$check" "$work/weak_leak" 2>&1) && status=0 || status=$?
if [ "$status" -eq 2 ] && [[ "$out" == *"GNU objdump"* ]]; then pass=$((pass + 1)); else
    fail=$((fail + 1))
    echo "FAIL: a non-GNU objdump is misuse — expected exit 2, got $status: $out" >&2
fi

# Nor may a GNU objdump that emits no parsable instructions.
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
printf '#!/bin/sh\n[ "$1" = --version ] && echo "GNU objdump (GNU Binutils) 2.40"; exit 0\n' > "$work/empty-objdump"
chmod +x "$work/empty-objdump"
out=$(OBJDUMP="$work/empty-objdump" "$check" "$work/weak_leak" 2>&1) && status=0 || status=$?
if [ "$status" -eq 2 ] && [[ "$out" == *"no instructions"* ]]; then pass=$((pass + 1)); else
    fail=$((fail + 1))
    echo "FAIL: no parsed instructions is misuse — expected exit 2, got $status: $out" >&2
fi

# An objdump that fails must be a tool problem (exit 2), never read as "leaks" (exit 1).
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
printf '#!/bin/sh\n[ "$1" = --version ] && { echo "GNU objdump (GNU Binutils) 2.40"; exit 0; }\necho "cannot disassemble" >&2; exit 1\n' > "$work/failing-objdump"
chmod +x "$work/failing-objdump"
out=$(OBJDUMP="$work/failing-objdump" "$check" "$work/weak_leak" 2>&1) && status=0 || status=$?
if [ "$status" -eq 2 ] && [[ "$out" == *"cannot disassemble"* ]]; then pass=$((pass + 1)); else
    fail=$((fail + 1))
    echo "FAIL: a failing objdump is a tool problem — expected exit 2, got $status: $out" >&2
fi

# The stripped-binary diagnosis must not depend on nm's message language.
cp "$work/weak_leak" "$work/stripped_de"
strip "$work/stripped_de"
out=$(LANG=de_DE.UTF-8 LC_ALL=de_DE.UTF-8 "$check" "$work/stripped_de" 2>&1) && status=0 || status=$?
if [ "$status" -eq 2 ] && [[ "$out" == *"no function symbols"* ]]; then pass=$((pass + 1)); else
    fail=$((fail + 1))
    echo "FAIL: stripped binary under a German locale — expected 'no function symbols', got $status: $out" >&2
fi

expect "missing binary is misuse" 2 "$work/does-not-exist"

echo "check-isa-leak: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
