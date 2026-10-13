#!/usr/bin/env bash
# Fail if an x86_64 binary carries AVX-512 code outside the avx512bw kernel tier.
#
# Usage: scripts/check-isa-leak.sh <binary>
#
# bwa-mem3's kernel TUs are compiled once per SIMD tier (gotcha #9) and picked at runtime by
# simd_dispatch.cpp, so only `*_avx512bw` code may use AVX-512. Inline and template functions a
# kernel TU instantiates (std::vector<short>::_M_fill_assign, ...) are weak symbols, though: when
# the linker keeps the copy from the avx512bw TU, every tier calls it, and the program dies with
# SIGILL on an x86_64 CPU without AVX-512 the first time the AVX2 path reaches it.
#
# A function is reported when it contains an AVX-512 instruction, its symbol is global or weak (so
# other code can bind to it), and its mangled name does not carry the avx512bw tier suffix
# (`_avx512bw`, from kernel_dispatch.h). An instruction is AVX-512 when it is EVEX-encoded (in
# 64-bit mode an opcode byte of 0x62 after any legacy prefixes is always the EVEX prefix, whatever
# the registers: `vpxord %ymm1,%ymm1,%ymm1` is AVX-512VL) or is one of the VEX-encoded opmask
# (k-register) instructions. Local symbols (static functions and the compiler's .constprop/.isra
# clones) stay inside their TU and are not reported. A binary with no function symbols (stripped)
# cannot be checked and is rejected.
#
# Needs GNU binutils' objdump and nm for x86_64 (on a non-x86_64 host, set OBJDUMP and NM to the
# cross tools, e.g. x86_64-linux-gnu-objdump), and mawk. The parser reads GNU objdump's layout, so
# another objdump (macOS's default is llvm-objdump), or output with no instructions in it, is
# rejected rather than read as "no leak". Exit 0: no leak; 1: leak(s), listed on stdout;
# 2: misuse or a tool problem.
set -euo pipefail

if [ $# -ne 1 ] || [ ! -f "$1" ]; then
    echo "usage: $0 <x86_64 binary>" >&2
    exit 2
fi
BIN=$1
OBJDUMP=${OBJDUMP:-objdump}
NM=${NM:-nm}

die() { # die <message>: report a problem that prevents a check, exit 2
    echo "check-isa-leak: $1" >&2
    exit 2
}
for tool in "$OBJDUMP" "$NM" mawk; do
    command -v "$tool" > /dev/null || die "$tool not found"
done
"$OBJDUMP" --version 2> /dev/null | head -1 | grep -q 'GNU objdump' \
    || die "$OBJDUMP is not GNU objdump; set OBJDUMP to GNU binutils' objdump for x86_64"

# Defined function symbols: address -> binding (nm type letter) and mangled name.
SYMS=$(mktemp)
# Functions containing AVX-512 instructions: address -> count, from the disassembly.
EVEX=$(mktemp)

# Raw output of nm, and objdump's disassembly and errors.
NM_OUT=$(mktemp)
DIS=$(mktemp)
DIS_ERR=$(mktemp)
trap 'rm -f "$SYMS" "$EVEX" "$NM_OUT" "$DIS" "$DIS_ERR"' EXIT
# nm exits non-zero both for a stripped binary ("no symbols") and for one it cannot read. LC_ALL=C
# keeps that message in English so the two can be told apart.
if ! LC_ALL=C "$NM" --defined-only "$BIN" > "$NM_OUT" 2>&1; then
    grep -q 'no symbols' "$NM_OUT" && die "no function symbols in $BIN (stripped?); cannot check it"
    die "$NM failed on $BIN: $(head -1 "$NM_OUT")"
fi
mawk '$2 ~ /^[TtWw]$/ { print $1 "\t" $2 "\t" $3 }' "$NM_OUT" | sort -u > "$SYMS"
[ -s "$SYMS" ] || die "no function symbols in $BIN (stripped?); cannot check it"

# -w keeps each instruction's bytes on one line: "<addr>:\t<hex bytes>\t<mnemonic operands>". A
# failing objdump is a tool problem (2), never "leaks" (1), so it runs on its own.
"$OBJDUMP" -d -w "$BIN" > "$DIS" 2> "$DIS_ERR" || die "$OBJDUMP failed on $BIN: $(head -1 "$DIS_ERR")"
mawk -F'\t' '
    /^[0-9a-f]+ <.*>:$/ { split($0, h, " "); addr = h[1]; next }
    NF >= 3 { total++ }
    NF >= 3 && $3 !~ /\(bad\)/ {
        nb = split($2, b, " ")
        i = 1
        while (i <= nb && b[i] ~ /^(66|67|f0|f2|f3|2e|36|3e|26|64|65)$/) i++
        if ((i <= nb && b[i] == "62") || $3 ~ /^k(mov|and|andn|or|xor|xnor|not|test|ortest|shift|unpck|add)/) n[addr]++
    }
    END { for (a in n) print a "\t" n[a]; print "total\t" total + 0 }
' "$DIS" | sort -u > "$EVEX"
# A disassembly with no instruction lines means the layout was not understood: fail, do not pass.
grep -q $'^total\t[1-9]' "$EVEX" || die "no instructions parsed from $OBJDUMP's output for $BIN"

# Join on address; report shared (T/W) symbols outside the avx512bw tier.
LEAKS=$(join -t "$(printf '\t')" "$EVEX" "$SYMS" | mawk -F'\t' '$3 ~ /^[TW]$/ && $4 !~ /_avx512bw/ { print $2 "\t" $4 }')
if [ -z "$LEAKS" ]; then
    echo "check-isa-leak: no AVX-512 code outside the avx512bw tier in $BIN"
    exit 0
fi
echo "check-isa-leak: AVX-512 code in shared symbols outside the avx512bw tier in $BIN:"
echo "$LEAKS" | while IFS=$'\t' read -r count sym; do
    printf '  %5s AVX-512 instructions  %s\n' "$count" "$(echo "$sym" | c++filt 2>/dev/null || echo "$sym")"
done
exit 1
