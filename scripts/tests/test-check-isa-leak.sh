#!/usr/bin/env bash
# Unit-tests check-isa-leak.sh against tiny x86_64 binaries built here. Each fixture puts one
# AVX-512 (or AVX2) instruction in a function via inline asm, so the result does not depend on
# what the compiler chooses to vectorize. Needs an x86_64 Linux host with cc, c++, GNU binutils
# and mawk; skipped elsewhere.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
check="$here/../check-isa-leak.sh"

skip() {
    echo "SKIP: check-isa-leak tests need an x86_64 Linux host with cc, c++, GNU binutils and mawk ($1)"
    exit 0
}
if [ "$(uname -s)" != Linux ] || [ "$(uname -m)" != x86_64 ]; then
    skip "not x86_64 Linux"
fi
# One tool per `command -v`: given several names it succeeds if any one of them exists.
for tool in cc c++ c++filt objdump nm strip mawk; do
    command -v "$tool" > /dev/null || skip "no $tool"
done
objdump --version | head -1 | grep -q 'GNU objdump' || skip "objdump is not GNU objdump"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
pass=0
fail=0

# build <name> <C source> [cc flags ...]: compile a C fixture binary into $work/<name>.
build() {
    local name=$1 src=$2
    shift 2
    printf '%s\n' "$src" > "$work/$name.c"
    cc -O1 "$@" -o "$work/$name" "$work/$name.c"
}

# build_cxx <name> <C++ source>: compile a C++ fixture binary into $work/<name>.
build_cxx() {
    printf '%s\n' "$2" > "$work/$1.cpp"
    c++ -O1 -o "$work/$1" "$work/$1.cpp"
}

# stub <name> <sh body>: write an executable /bin/sh script to $work/<name>.
stub() {
    printf '#!/bin/sh\n%s\n' "$2" > "$work/$1"
    chmod +x "$work/$1"
}

# expect <description> <status> <binary> [text the output must contain] [VAR=value ...]
# Runs the check with the optional environment and compares its exit status and output.
expect() {
    local desc=$1 want=$2 bin=$3 text=${4:-} out status=0
    shift $(($# < 4 ? $# : 4))
    out=$(env "$@" "$check" "$bin" 2>&1) || status=$?
    if [ "$status" -ne "$want" ]; then
        fail=$((fail + 1))
        echo "FAIL: $desc — expected exit $want, got $status: $out" >&2
    elif [ -n "$text" ] && [[ "$out" != *"$text"* ]]; then
        fail=$((fail + 1))
        echo "FAIL: $desc — output lacks '$text': $out" >&2
    else
        pass=$((pass + 1))
    fi
}

# --- What counts as a leak: AVX-512 code reachable from outside the avx512bw tier ---

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

# A leaked entry point with no AVX-512 of its own that calls a local AVX-512 helper still runs
# AVX-512 for every caller.
build local_leak '
static __attribute__((noinline)) void body(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
__attribute__((weak, noinline)) void shared_helper(void) { body(); }
int main(void) { shared_helper(); return 0; }'
expect "local AVX-512 helper behind a leaked entry is a leak" 1 "$work/local_leak" "body"

# The linker turns a hidden-visibility copy local in the output, but every caller still uses it.
build hidden_leak '
__attribute__((weak, noinline, visibility("hidden"))) void shared_helper(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { shared_helper(); return 0; }'
expect "hidden-visibility (local in the output) copy is a leak" 1 "$work/hidden_leak" "shared_helper"

# The avx512bw tier itself may use AVX-512.
build tier_ok '
__attribute__((noinline)) void ksw_extend2_avx512bw(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { ksw_extend2_avx512bw(); return 0; }'
expect "avx512bw tier symbol is allowed" 0 "$work/tier_ok"

# So may a local helper only the tier calls, even through another local helper.
build tier_local_ok '
static __attribute__((noinline)) void pack(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
static __attribute__((noinline)) void pack_outer(void) { pack(); __asm__ volatile("vpxord %%zmm2, %%zmm2, %%zmm2" ::: "xmm2"); }
__attribute__((noinline)) void ksw_extend2_avx512bw(void) { pack_outer(); }
int main(void) { ksw_extend2_avx512bw(); return 0; }'
expect "local helpers reached only from the tier are allowed" 0 "$work/tier_local_ok"

# A local AVX-512 function that code outside the tier calls is not.
build static_leak '
static __attribute__((noinline)) void pack(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { pack(); return 0; }'
expect "local AVX-512 function called from outside the tier is a leak" 1 "$work/static_leak" "pack"

# Nor is one the tier calls directly but whose address is also taken: a pointer to it can be called
# from anywhere. The pointer may live in data (a dynamic relocation's addend in a PIE, a plain word
# otherwise) or be formed in code (a rip-relative lea in a PIE, an immediate otherwise).
for pie in -pie -no-pie; do
    if [ "$pie" = -pie ]; then pic=-fPIE; else pic=-fno-PIE; fi
    build "ptr_data_leak$pie" '
static __attribute__((noinline)) void pack(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
__attribute__((noinline)) void ksw_extend2_avx512bw(void) { pack(); }
static void (*volatile fp)(void) = pack;
int main(void) { ksw_extend2_avx512bw(); fp(); return 0; }' "$pic" "$pie"
    expect "tier-called helper whose address is in data is a leak ($pie)" 1 "$work/ptr_data_leak$pie" "pack"
    build "ptr_code_leak$pie" '
static __attribute__((noinline)) void pack(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
__attribute__((noinline)) void ksw_extend2_avx512bw(void) { pack(); }
void (*volatile fp)(void);
int main(void) { fp = pack; ksw_extend2_avx512bw(); fp(); return 0; }' "$pic" "$pie"
    expect "tier-called helper whose address is taken in code is a leak ($pie)" 1 "$work/ptr_code_leak$pie" "pack"
done

# Only the avx512bw tier's naming convention is exempt, not any name containing "avx512".
build name_leak '
__attribute__((weak, noinline)) void shared_avx512_helper(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { shared_avx512_helper(); return 0; }'
expect "a helper merely named avx512 is not exempt" 1 "$work/name_leak" "shared_avx512_helper"

# The marker is "_avx512bw" ending an identifier or followed by "_", not a prefix of a longer word.
build near_miss_leak '
__attribute__((weak, noinline)) void ksw_extend2_avx512bw2(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { ksw_extend2_avx512bw2(); return 0; }'
expect "a name merely starting with the tier marker is not exempt" 1 "$work/near_miss_leak" "ksw_extend2_avx512bw2"

# bwa-mem3's own factories carry the marker mid-name (make_kswv_kernel_avx512bw_mat), and the
# compiler's clones keep the name with a dot suffix (.cold, .part.0, .isra.0).
build mid_marker_ok '
__attribute__((noinline)) void make_kswv_kernel_avx512bw_mat(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
__attribute__((noinline)) void clone(void) __asm__("ksw_extend2_avx512bw.part.0");
void clone(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main(void) { make_kswv_kernel_avx512bw_mat(); clone(); return 0; }'
expect "tier marker mid-name or before a clone suffix is allowed" 0 "$work/mid_marker_ok"

# Real tier symbols are C++ mangled, with the marker mid-name.
build_cxx cxx_tier_ok '
namespace kswv { __attribute__((noinline)) void ksw_extend2_avx512bw(int) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); } }
int main() { kswv::ksw_extend2_avx512bw(1); return 0; }'
expect "C++ mangled avx512bw tier symbol is allowed" 0 "$work/cxx_tier_ok"

# So are members of the tier's classes (BandedPairWiseSW_avx512bw::getScores8, ...).
build_cxx cxx_class_ok '
struct kswv_avx512bw { __attribute__((noinline)) void run(int); };
void kswv_avx512bw::run(int) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main() { kswv_avx512bw k; k.run(1); return 0; }'
expect "member of an avx512bw tier class is allowed" 0 "$work/cxx_class_ok"

# A mangled identifier that only starts with the marker is a different function.
build_cxx cxx_near_miss_leak '
namespace kswv { __attribute__((noinline)) void ksw_extend2_avx512bwExtra(int) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); } }
int main() { kswv::ksw_extend2_avx512bwExtra(1); return 0; }'
expect "C++ identifier merely starting with the tier marker is not exempt" 1 "$work/cxx_near_miss_leak" "ksw_extend2_avx512bwExtra"

# A leaked C++ template instantiation is reported by its demangled name.
build_cxx cxx_leak '
template <class T> __attribute__((noinline)) void fill_assign(T) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
int main() { fill_assign<short>(1); return 0; }'
expect "leaked C++ template is reported demangled" 1 "$work/cxx_leak" "void fill_assign<short>(short)"

# Other libraries' own runtime-dispatched AVX-512 code (crc32fast, zlib-ng, ...) is a tier root of
# its own; ISA_LEAK_ALLOW names it explicitly, and nothing else is exempt.
build dispatch_ok '
__attribute__((noinline)) void crc32_x86_vpclmulqdq_avx512_vl512(void) { __asm__ volatile("vpxord %%zmm1, %%zmm1, %%zmm1" ::: "xmm1"); }
static __attribute__((noinline)) void fold(void) { __asm__ volatile("vpxord %%zmm2, %%zmm2, %%zmm2" ::: "xmm2"); }
__attribute__((noinline)) void calculate_avx512(void) { fold(); }
int main(void) { crc32_x86_vpclmulqdq_avx512_vl512(); calculate_avx512(); return 0; }'
expect "self-dispatching AVX-512 code is a leak by default" 1 "$work/dispatch_ok" "crc32_x86_vpclmulqdq_avx512_vl512"
expect "ISA_LEAK_ALLOW makes self-dispatching code (and its helpers) a tier root" 0 "$work/dispatch_ok" "" \
    ISA_LEAK_ALLOW='crc32_x86_vpclmulqdq_avx512_vl(256|512)$|^calculate_avx512$'
expect "ISA_LEAK_ALLOW exempts only what it names" 1 "$work/weak_leak" "shared_helper" ISA_LEAK_ALLOW='^calculate_avx512$'

# --- What counts as AVX-512 ---

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

# A legacy (segment) prefix in front of the EVEX prefix: 64 62 ...
build seg_leak '
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile("vmovdqu64 %%fs:0, %%zmm1" ::: "xmm1"); }
int main(void) { shared_helper(); return 0; }'
expect "EVEX after a legacy prefix is a leak" 1 "$work/seg_leak" "shared_helper"

# A 0x62 byte that is not the opcode byte (here an immediate: b8 62 00 00 00) is not EVEX.
# shellcheck disable=SC2016 # $0x62 is assembler syntax, not a shell expansion
build imm_ok '
__attribute__((weak, noinline)) int shared_helper(void) { int x; __asm__ volatile("movl $0x62, %0" : "=r"(x)); return x; }
int main(void) { return shared_helper() == 0x62 ? 0 : 1; }'
expect "0x62 outside the opcode byte is not EVEX" 0 "$work/imm_ok"

# The VEX-encoded opmask instructions are AVX-512 too, though no 0x62 byte marks them.
for insn in "kmovw %%k1, %%k2" "kandw %%k1, %%k2, %%k3" "kandnw %%k1, %%k2, %%k3" "korw %%k1, %%k2, %%k3" \
    "kxorw %%k1, %%k2, %%k3" "kxnorw %%k1, %%k2, %%k3" "knotw %%k1, %%k2" "ktestw %%k1, %%k2" \
    "kortestw %%k1, %%k2" "kshiftlw \$1, %%k1, %%k2" "kunpckbw %%k1, %%k2, %%k3" "kaddw %%k1, %%k2, %%k3"; do
    name=kmask_$(echo "$insn" | mawk '{ print $1 }')
    build "$name" "
__attribute__((weak, noinline)) void shared_helper(void) { __asm__ volatile(\"$insn\" ::: ); }
int main(void) { shared_helper(); return 0; }"
    expect "opmask instruction ${insn%% *} is a leak" 1 "$work/$name" "shared_helper"
done

# --- Misuse and tool problems: exit 2, never a pass and never "leaks" ---

cp "$work/weak_leak" "$work/stripped"
strip "$work/stripped"
expect "stripped binary is misuse" 2 "$work/stripped" "no function symbols"

printf 'not an object file\n' > "$work/not_elf"
expect "a file nm cannot read is a tool problem" 2 "$work/not_elf" "nm failed"

expect "missing binary is misuse" 2 "$work/does-not-exist"

expect "a missing tool is a tool problem" 2 "$work/weak_leak" "not found" OBJDUMP="$work/no-such-objdump"

# An objdump whose output the parser cannot read (llvm-objdump's layout differs) must not pass.
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
stub llvm-objdump '[ "$1" = --version ] && echo "LLVM version 21"; exit 0'
expect "a non-GNU objdump is misuse" 2 "$work/weak_leak" "GNU objdump" OBJDUMP="$work/llvm-objdump"

# Nor may a GNU objdump that emits no parsable instructions.
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
stub empty-objdump '[ "$1" = --version ] && echo "GNU objdump (GNU Binutils) 2.40"; exit 0'
expect "no parsed instructions is misuse" 2 "$work/weak_leak" "no instructions" OBJDUMP="$work/empty-objdump"

# A failing objdump is a tool problem (exit 2), never read as "leaks" (exit 1).
# shellcheck disable=SC2016 # $1 belongs to the stub script, not to this one
stub failing-objdump '[ "$1" = --version ] && { echo "GNU objdump (GNU Binutils) 2.40"; exit 0; }
echo "cannot disassemble" >&2; exit 1'
expect "a failing objdump is a tool problem" 2 "$work/weak_leak" "cannot disassemble" OBJDUMP="$work/failing-objdump"

echo "check-isa-leak: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
