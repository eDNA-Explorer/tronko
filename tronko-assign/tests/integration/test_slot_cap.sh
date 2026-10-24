#!/usr/bin/env bash
# The candidate-slot arrays at a small cap. Builds a copy of tronko-assign with MAX_NUM_BWA_MATCHES
# 2 and AddressSanitizer, and runs it on chimeric single-end reads (three segments of three
# different leaves of the example Charadriiformes reference, out of reference order, so that BWA
# reports each segment as its own record and each read has three candidate leaves, one more than
# the slots). Passes when the run finishes with no AddressSanitizer report: the SAM parse
# (fastmap.c) must stop at the last slot instead of writing past the arrays, and the placement
# worker must not read past them.
#
# Usage: tests/integration/test_slot_cap.sh     (about 15 s: it compiles the whole program once)
# Exit 77 (skip) when the compiler cannot build with -fsanitize=address or the example reference
# (tronko-build/example_datasets/single_tree) is not beside tronko-assign/.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TA=$(cd "$HERE/../.." && pwd)
# canonical, without "..": tronko-assign refuses a path option longer than its buffer (199 characters)
ROOT=$(cd "$TA/.." && pwd)
EX=$ROOT/tronko-build/example_datasets/single_tree
CC=${CC:-gcc}
if [ ! -f "$EX/Charadriiformes.fasta" ] || [ ! -f "$EX/reference_tree.txt" ]; then
	echo "SKIP: no example reference at $EX"; exit 77
fi
T=$(mktemp -d "${TMPDIR:-/tmp}/tronko-slot-cap.XXXXXX")
if ! echo 'int main(void){return 0;}' | $CC -fsanitize=address -x c - -o "$T/probe" 2> /dev/null; then
	echo "SKIP: $CC cannot build with -fsanitize=address"; exit 77
fi
mkdir -p "$T/src"
tar -C "$TA" --exclude=./carquet --exclude=./tests --exclude=./tronko-assign -cf - . | tar -C "$T/src" -xf -
sed -i 's/^#define MAX_NUM_BWA_MATCHES [0-9]*/#define MAX_NUM_BWA_MATCHES 2/' "$T/src/global.h"
grep -q '^#define MAX_NUM_BWA_MATCHES 2$' "$T/src/global.h" || { echo "FAIL: could not set the cap"; exit 1; }
if ! (cd "$T/src" && make ARCH_FLAGS="-fsanitize=address -fno-omit-frame-pointer" > "$T/build.log" 2>&1); then
	echo "FAIL: build (log $T/build.log)"; exit 1
fi
# chimeric reads: for leaves i, i+7, i+13 take bases 200-259, 20-79 and 120-179, joined in that order
awk -v out="$T/chimeras.fasta" '
	/^>/ { if (seq != "") s[n++] = seq; seq = ""; next }
	{ seq = seq toupper($0) }
	END {
		if (seq != "") s[n++] = seq
		m = 0
		for (i = 0; i + 13 < n && m < 200; i += 3) {
			a = s[i]; b = s[i + 7]; c = s[i + 13]
			if (length(a) < 260 || length(b) < 80 || length(c) < 180) continue
			printf(">chimera_%d\n%s%s%s\n", m++, substr(a, 201, 60), substr(b, 21, 60), substr(c, 121, 60)) > out
		}
	}' "$EX/Charadriiformes.fasta"
n=$(grep -c '^>' "$T/chimeras.fasta")
[ "$n" -ge 50 ] || { echo "FAIL: only $n chimeric reads"; exit 1; }
(cd "$T" && ASAN_OPTIONS=detect_leaks=0 "$T/src/tronko-assign" -r -f "$EX/reference_tree.txt" -a "$EX/Charadriiformes.fasta" \
	-s -g "$T/chimeras.fasta" -w -6 --Cinterval 10 --number-of-cores 1 -o "$T/out.tsv" > "$T/stdout" 2> "$T/stderr")
rc=$?
if [ $rc -ne 0 ] || grep -q 'ERROR: AddressSanitizer' "$T/stderr"; then
	echo "FAIL: exit $rc; $(grep -m1 -A3 'ERROR: AddressSanitizer' "$T/stderr")"
	echo "(work directory kept: $T)"
	exit 1
fi
rows=$(($(wc -l < "$T/out.tsv") - 1))
echo "test_slot_cap: $n chimeric reads at cap 2, $rows rows, no AddressSanitizer report: PASS"
exit 0
