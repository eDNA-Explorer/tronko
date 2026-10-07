#!/usr/bin/env bash
# The path options of tronko-assign against the size of the buffer each one is read into
# (global.h, Options). For each option: a path of exactly the largest length that fits is taken
# as given, and a path one character longer stops the run, exit 1, with a message naming the
# option, before anything else happens. Also checks that the word taken is the one
# sscanf("%s") took (leading whitespace skipped, cut at the first whitespace).
#
# Usage: tests/integration/test_path_options.sh [tronko-assign binary]
# Default binary: tronko-assign/tronko-assign next to tests/ (make). Exit 77 (skip) without it.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN=${1:-${TRONKO_ASSIGN:-$HERE/../../tronko-assign}}
if [ ! -x "$BIN" ]; then
	echo "SKIP: no tronko-assign binary at $BIN (run make first)"
	exit 77
fi
fail=0
checks=0
path_of_length() { # a path of exactly $1 characters that does not exist
	local n=$1 p="/nonexistent-tronko-test/"
	while [ ${#p} -lt "$n" ]; do p="${p}x"; done
	printf '%s' "${p:0:$n}"
}
# option, buffer size (bytes, terminator included), long name
cases=(
	"-f 200 reference-file"
	"-m 200 msa-file"
	"-x 200 tax-file"
	"-d 200 partitions-directory"
	"-o 200 results"
	"-a 200 fasta-file"
	"-g 2000 single-read-file"
	"-1 2000 paired-read-file1"
	"-2 2000 paired-read-file2"
	"-t 2000 tree-file"
	"-5 1000 print-node-info"
	"-3 1000 print-alignments-dir"
	"-4 1000 tree-dir"
)
ref=/nonexistent-tronko-test/reference_tree.trkb
for c in "${cases[@]}"; do
	read -r opt size long <<< "$c"
	max=$((size - 1))
	for len in "$max" "$size"; do
		p=$(path_of_length "$len")
		if [ "$opt" = -f ]; then args=(-r -f "$p"); else args=(-r -f "$ref" "$opt" "$p"); fi
		out=$("$BIN" "${args[@]}" 2>&1)
		rc=$?
		checks=$((checks + 1))
		msg="Path given to $opt (--$long) is $len characters long; at most $max are supported."
		if [ "$len" -eq "$max" ]; then
			if grep -q "Path given to" <<< "$out"; then
				echo "FAIL $opt: a path of $len characters was refused"; fail=1
			elif ! grep -q "Cannot find reference_tree.txt file" <<< "$out"; then
				echo "FAIL $opt: a path of $len characters: the run did not reach the reference check: $out"; fail=1
			fi
		else
			if [ "$rc" -ne 1 ] || ! grep -qF "$msg" <<< "$out"; then
				echo "FAIL $opt: a path of $len characters: exit $rc, output: $out"; fail=1
			elif grep -q "Cannot find reference_tree.txt file" <<< "$out"; then
				echo "FAIL $opt: a path of $len characters: the run went on after the refusal"; fail=1
			fi
		fi
	done
done
# the word taken, as sscanf("%s") took it
out=$("$BIN" -r -f "   /nonexistent-tronko-test/a b" 2>&1)
checks=$((checks + 1))
if ! grep -qF "Cannot find reference_tree.txt file: /nonexistent-tronko-test/a. Exiting..." <<< "$out"; then
	echo "FAIL whitespace: $out"; fail=1
fi
# only the word counts against the limit: a long tail after whitespace is ignored, as before
tail=$(path_of_length 500)
out=$("$BIN" -r -f "/nonexistent-tronko-test/b $tail" 2>&1)
checks=$((checks + 1))
if ! grep -qF "Cannot find reference_tree.txt file: /nonexistent-tronko-test/b. Exiting..." <<< "$out"; then
	echo "FAIL whitespace tail: $out"; fail=1
fi
echo "test_path_options: $checks checks, $([ $fail = 0 ] && echo PASS || echo FAIL)"
exit $fail
