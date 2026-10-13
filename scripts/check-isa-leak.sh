#!/usr/bin/env bash
# Fail if an x86_64 binary can run AVX-512 code from outside the avx512bw kernel tier.
#
# Usage: [ISA_LEAK_ALLOW=<ere>] scripts/check-isa-leak.sh <binary>
#
# Only the avx512bw tier's code may use AVX-512: simd_dispatch.cpp picks that tier only on CPUs
# that have it. Code shared between tiers that carries AVX-512 instructions makes the program die
# with SIGILL on CPUs without it. Why that can happen, and how bwa-mem3-sys's link order prevents
# it, is in bwa-mem3-sys/build.rs (step 4a) and CLAUDE.md gotcha #9.
#
# What is reported. A function is tier-only when its name carries the tier marker, or when every
# direct call/jump to it comes from a tier-only function and nothing takes its address, computed to
# a fixed point. The marker is `_avx512bw` ending an identifier or followed by `_`, so
# `ksw_extend2_avx512bw`, `make_kswv_kernel_avx512bw_mat` and members of class `kswv_avx512bw`
# carry it but `ksw_extend2_avx512bw2` does not; each length-prefixed identifier of a mangled C++
# name is checked, and a C name up to a compiler clone suffix (`.cold`, `.isra.0`). Identifiers in
# template arguments count too, so the members of `std::vector<kswv_avx512bw *>` are tier code;
# such instantiations exist only where the tier's own types are named. Every function that
# contains an AVX-512 instruction and is not tier-only is reported, whatever its symbol binding: a
# weak or hidden copy that the linker made local in the output, or a local helper behind a leaked
# entry point, is as reachable as a global one. Functions are keyed by address, since several
# static functions can share a name.
#
# A pointer to a function can be called from anywhere, so a function with no direct callers, or
# whose address is taken, cannot be proven tier-only and is reported unless its own name makes it a
# tier root. An address counts as taken when code refers to it other than by a direct call or jump
# (a rip-relative lea, an immediate), or when it is an aligned 8-byte word in an allocated data
# section, which covers both a pointer stored in data and a PIE's dynamic relocation addends. This
# errs toward a false alarm rather than a missed leak.
#
# Other libraries that do their own runtime CPU dispatch (crc32fast, zlib-ng, ...) carry legitimate
# AVX-512 code outside the tier. ISA_LEAK_ALLOW, an extended regex matched against mangled names,
# makes the functions it names tier roots too; name each such entry point explicitly, so that it
# cannot hide a real leak.
#
# What counts as AVX-512: an EVEX-encoded instruction (in 64-bit mode an opcode byte of 0x62 after
# any legacy prefixes is always the EVEX prefix, whatever the registers: `vpxord %ymm1,%ymm1,%ymm1`
# is AVX-512VL) or a VEX-encoded opmask (k-register) instruction.
#
# Scope. Run it on an unstripped binary built with the default features: `--features avx512`, or
# `native` on an AVX-512 host, compiles the baseline itself with AVX-512 on purpose, and a stripped
# binary (such as the PyPI wheel's extension) has no names to check and is rejected.
#
# Needs GNU binutils' objdump and nm for x86_64 (on a non-x86_64 host, set OBJDUMP and NM to the
# cross tools, e.g. x86_64-linux-gnu-objdump), mawk, and c++filt for readable names. The parser reads
# GNU objdump's layout, so another objdump (macOS's default is llvm-objdump), or output with no
# instructions in it, is rejected rather than read as "no leak". Exit 0: no leak; 1: leak(s), listed
# on stdout; 2: misuse or a tool problem.
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

NM_OUT=$(mktemp)  # nm's output, for the stripped-binary check
DIS=$(mktemp)     # objdump's disassembly
DIS_ERR=$(mktemp) # objdump's errors
WORDS=$(mktemp)   # aligned 8-byte words of the allocated data sections, one hex value per line
REPORT=$(mktemp)  # "<count>\t<name>" per leaking function, then "total\t<instructions parsed>"
trap 'rm -f "$NM_OUT" "$DIS" "$DIS_ERR" "$WORDS" "$REPORT"' EXIT

# GNU nm exits 0 on a stripped binary (printing "no symbols"); non-zero means it could not read it.
"$NM" --defined-only "$BIN" > "$NM_OUT" 2>&1 || die "$NM failed on $BIN: $(head -1 "$NM_OUT")"
mawk '$2 ~ /^[TtWw]$/ { found = 1 } END { exit !found }' "$NM_OUT" \
    || die "no function symbols in $BIN (stripped?); cannot check it"

# -w keeps each instruction's bytes on one line: "<addr>:\t<hex bytes>\t<mnemonic operands>". A
# failing objdump is a tool problem (2), never "leaks" (1), so it runs on its own.
"$OBJDUMP" -d -w "$BIN" > "$DIS" 2> "$DIS_ERR" || die "$OBJDUMP failed on $BIN: $(head -1 "$DIS_ERR")"

# Every aligned 8-byte word of the allocated, non-code sections with contents, as a hex value: a
# function address among them is a pointer to it (in a PIE, the addend of its dynamic relocation).
# -s prints 16 bytes per line as " <addr> <4 groups of 4 bytes>  <ascii>", in memory order.
SECTIONS=$("$OBJDUMP" -h -w "$BIN" 2> "$DIS_ERR" \
    | mawk '/CONTENTS/ && /ALLOC/ && !/CODE/ { printf " -j %s", $2 }') \
    || die "$OBJDUMP -h failed on $BIN: $(head -1 "$DIS_ERR")"
if [ -n "$SECTIONS" ]; then
    # shellcheck disable=SC2086 # SECTIONS is a list of "-j <section>" arguments
    "$OBJDUMP" -s $SECTIONS "$BIN" 2> "$DIS_ERR" | mawk '
        /^Contents of section / { buf = ""; next }
        /^ [0-9a-f]+ / {
            low = index("0123456789abcdef", substr($1, length($1), 1)) - 1
            hex = substr($0, length($1) + 3, 35); gsub(/ /, "", hex)
            for (j = 0; 2 * j < length(hex); j++) {
                pos = (low + j) % 8
                if (pos == 0) buf = ""
                buf = buf substr(hex, 2 * j + 1, 2)
                if (pos == 7 && length(buf) == 16) {
                    w = ""; for (k = 15; k >= 1; k -= 2) w = w substr(buf, k, 2)
                    print w
                }
            }
        }' > "$WORDS" || die "$OBJDUMP -s failed on $BIN: $(head -1 "$DIS_ERR")"
fi

mawk -F'\t' -v allow="${ISA_LEAK_ALLOW:-}" -v words="$WORDS" '
    function key(a) { sub(/^0+/, "", a); return a }
    # Does the name carry the tier marker? See the header for the rule.
    function marked(n,   s, i, j, len, id) {
        if (n !~ /^_Z/) { s = n; sub(/[.@].*$/, "", s); return s ~ /_avx512bw(_|$)/ }
        for (i = 1; i <= length(n); i++) {
            if (substr(n, i, 1) !~ /[0-9]/ || (i > 1 && substr(n, i - 1, 1) ~ /[0-9]/)) continue
            for (j = i; substr(n, j, 1) ~ /[0-9]/; j++) {}
            len = substr(n, i, j - i) + 0
            id = substr(n, j, len)
            if (len > 0 && length(id) == len && id ~ /_avx512bw(_|$)/) return 1
            i = j - 1
        }
        return 0
    }
    # Function header: "<16-hex-digit address> <name>:".
    /^[0-9a-f]+ <.*>:$/ {
        cur = key(substr($0, 1, index($0, " ") - 1))
        nm = $0; sub(/^[0-9a-f]+ </, "", nm); sub(/>:$/, "", nm)
        name[cur] = nm
        next
    }
    NF >= 3 {
        total++
        if ($3 !~ /\(bad\)/) {
            nb = split($2, b, " ")
            i = 1
            while (i <= nb && b[i] ~ /^(66|67|f0|f2|f3|2e|36|3e|26|64|65)$/) i++
            if ((i <= nb && b[i] == "62") || $3 ~ /^k(mov|and|andn|or|xor|xnor|not|test|ortest|shift|unpck|add)/) evex[cur]++
        }
        # A direct call or jump to the start of another function: "call   401136 <body>".
        if ($3 ~ /(^|[ ])(call|j[a-z]+)[a-z]*[ ]+[0-9a-f]+ <[^>+]+>$/) {
            n = split($3, o, " ")
            t = key(o[n - 1])
            if (t != cur) { ncaller[t]++; caller[t, ncaller[t]] = cur }
        } else {
            # Any other reference takes an address: "lea 0x7(%rip),%rax  # 401136 <body>", or an
            # immediate "$0x401136". Those that are not a function start are ignored in END.
            if (match($3, /# [0-9a-f]+ <[^>+]+>$/)) {
                t = substr($3, RSTART + 2)
                taken[key(substr(t, 1, index(t, " ") - 1))] = 1
            }
            ops = $3
            while (match(ops, /\$0x[0-9a-f]+/)) { taken[key(substr(ops, RSTART + 3, RLENGTH - 3))] = 1; ops = substr(ops, RSTART + RLENGTH) }
        }
    }
    END {
        while ((getline w < words) > 0) taken[key(w)] = 1
        for (f in name) tier[f] = (marked(name[f]) || (allow != "" && name[f] ~ allow))
        do {
            changed = 0
            for (f in name) {
                if (tier[f] || !(f in ncaller) || (f in taken)) continue
                ok = 1
                for (i = 1; i <= ncaller[f]; i++) if (!tier[caller[f, i]]) { ok = 0; break }
                if (ok) { tier[f] = 1; changed = 1 }
            }
        } while (changed)
        for (f in evex) if (!tier[f]) print evex[f] "\t" name[f]
        print "total\t" total + 0
    }
' "$DIS" > "$REPORT"
# A disassembly with no instruction lines means the layout was not understood: fail, do not pass.
grep -q $'^total\t[1-9]' "$REPORT" || die "no instructions parsed from $OBJDUMP's output for $BIN"

LEAKS=$(grep -v $'^total\t' "$REPORT" | sort -t $'\t' -k2,2 || true)
if [ -z "$LEAKS" ]; then
    echo "check-isa-leak: no AVX-512 code outside the avx512bw tier in $BIN"
    exit 0
fi
echo "check-isa-leak: AVX-512 code reachable from outside the avx512bw tier in $BIN:"
echo "$LEAKS" | while IFS=$'\t' read -r count sym; do
    printf '  %5s AVX-512 instructions  %s\n' "$count" "$(echo "$sym" | c++filt 2> /dev/null || echo "$sym")"
done
exit 1
