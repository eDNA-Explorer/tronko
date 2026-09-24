#!/usr/bin/env bash
# The two-pass Needleman-Wunsch fill against the original, on every alignment of the repository
# fixtures, with the production invocation.
#
# TRONKO_ASSIGN_BIN must be tronko-assign built with -DNW_SHADOW (make -C tronko-assign/tests
# test-shadow builds one and runs this script). In that build every alignment the two-pass fill
# computes is recomputed with the original fill, the three score matrices are compared, and any
# difference aborts the process; at exit it prints "NW_SHADOW alignments_compared=<n>
# mismatches=0". Each case must exit 0, print that line with n > 0, and reproduce the
# production golden TSV byte for byte.
#
# Environment: TRONKO_ASSIGN_BIN (required), TRONKO_ASSIGN_CORES (thread counts, default "1 4"),
# TRONKO_FIXTURE_DIR (default tests/data/assignment at the repository root, the fixtures of
# test_assignment_production_parity.sh). Exit 77 if those fixtures are absent.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
ASSIGN=${TRONKO_ASSIGN_BIN:-}
CORES=${TRONKO_ASSIGN_CORES:-"1 4"}
DATA=${TRONKO_FIXTURE_DIR:-"$REPO_ROOT/tests/data/assignment"}
FASTA="$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/tronko-nw-shadow.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

[[ -n "$ASSIGN" && -x "$ASSIGN" ]] || { echo "Set TRONKO_ASSIGN_BIN to a tronko-assign built with -DNW_SHADOW" >&2; exit 2; }
for path in "$DATA/reference_tree.trkb" "$DATA/single_4000.fasta" "$DATA/paired_2000_1.fasta" \
	"$DATA/paired_2000_2.fasta" "$DATA/expected_single_nw_trkb.tsv" \
	"$DATA/expected_unpaired_r_nw_trkb.tsv" "$DATA/expected_paired_nw_trkb.tsv"; do
	[[ -f "$path" ]] || { echo "SKIP: missing fixture $path (set TRONKO_FIXTURE_DIR)"; exit 77; }
done

COMMON=(-r -f "$DATA/reference_tree.trkb" -a "$FASTA" -w -6 --Cinterval 10)
total=0

# run_case <label> <tronko args...>
run_case() {
	local label=$1 n line count
	shift
	for n in $CORES; do
		if ! "$ASSIGN" "${COMMON[@]}" --number-of-cores "$n" "$@" -o "$TMP_DIR/$label-c$n.tsv" \
			>/dev/null 2>"$TMP_DIR/$label-c$n.log"; then
			echo "NW shadow check failed: $label, --number-of-cores $n exited non-zero" >&2
			grep -a 'NW_SHADOW' "$TMP_DIR/$label-c$n.log" >&2 || tail -5 "$TMP_DIR/$label-c$n.log" >&2
			return 1
		fi
		line=$(grep -a '^NW_SHADOW alignments_compared=' "$TMP_DIR/$label-c$n.log" | tail -1 || true)
		count=$(sed -n 's/^NW_SHADOW alignments_compared=\([0-9]*\) .*/\1/p' <<<"$line")
		if [[ -z "$count" || "$count" -eq 0 ]]; then
			echo "NW shadow check failed: $label, --number-of-cores $n: no comparison reported (is $ASSIGN a -DNW_SHADOW build?)" >&2
			return 1
		fi
		if ! cmp -s "$DATA/expected_$label.tsv" "$TMP_DIR/$label-c$n.tsv"; then
			echo "NW shadow check failed: $label, --number-of-cores $n: output differs from expected_$label.tsv" >&2
			diff -u "$DATA/expected_$label.tsv" "$TMP_DIR/$label-c$n.tsv" | head -20 >&2 || true
			return 1
		fi
		total=$((total + count))
		echo "PASS: $label, --number-of-cores $n: $count alignments, matrices identical, output identical"
	done
}

run_case single_nw_trkb -s -g "$DATA/single_4000.fasta"
run_case unpaired_r_nw_trkb -s -v -g "$DATA/paired_2000_2.fasta"
run_case paired_nw_trkb -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
echo "NW shadow check: $total alignments compared, 0 differ"
