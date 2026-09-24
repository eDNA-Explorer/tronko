#!/usr/bin/env bash
# Assignment parity with the production invocation of tronko-assign.
#
# Complements test_assignment_parity.sh (tiny fixtures, -C 1 -c 5, text reference)
# with the flags the eDNA Explorer pipeline uses in production:
#   paired:     -p -z -w -1 F -2 R -6 --number-of-cores N --Cinterval 10
#   unpaired_F: -s -w -g F ...       unpaired_R: -s -w -v -g R ...
# against a .trkb reference, on fixtures large enough to contain reads whose
# best BWA hits tie (about 0.5% of single-end and 3% of paired reads in these sets).
#
# The goldens were generated at --number-of-cores 1 from the production branch.
# TRONKO_ASSIGN_CORES lists the core counts that must reproduce them; every count
# in the list must produce a file byte-identical to the golden.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
ASSIGN=${TRONKO_ASSIGN_BIN:-"$REPO_ROOT/tronko-assign/tronko-assign"}
BASELINE_ASSIGN=${TRONKO_ASSIGN_BASELINE_BIN:-}
CORES=${TRONKO_ASSIGN_CORES:-"1"}
DATA="$REPO_ROOT/tests/data/assignment"
REFERENCE="$DATA/reference_tree.trkb"
FASTA="$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/tronko-prod-parity.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

compare_output() {
	local expected=$1
	local actual=$2
	local label=$3
	if ! cmp -s "$expected" "$actual"; then
		echo "Assignment parity failed: $label" >&2
		diff -u "$expected" "$actual" | head -40 >&2 || true
		return 1
	fi
	echo "PASS: $label"
}

if [[ ! -x "$ASSIGN" ]]; then
	echo "Missing executable: $ASSIGN" >&2
	exit 1
fi
for path in "$REFERENCE" "$FASTA" "$DATA/single_4000.fasta" \
	"$DATA/paired_2000_1.fasta" "$DATA/paired_2000_2.fasta"; do
	if [[ ! -f "$path" ]]; then
		echo "Missing parity fixture: $path" >&2
		exit 1
	fi
done

COMMON=(-r -f "$REFERENCE" -a "$FASTA" -w -6 --Cinterval 10)
echo "reference=$REFERENCE cores=$CORES"

# case <label> <expected-golden> <tronko args...>
run_case() {
	local label=$1 expected=$2
	shift 2
	if [[ -n "$BASELINE_ASSIGN" ]]; then
		expected="$TMP_DIR/baseline-$label.tsv"
		"$BASELINE_ASSIGN" "${COMMON[@]}" --number-of-cores 1 "$@" -o "$expected" >/dev/null 2>&1
	fi
	local n
	for n in $CORES; do
		"$ASSIGN" "${COMMON[@]}" --number-of-cores "$n" "$@" -o "$TMP_DIR/$label-c$n.tsv" >/dev/null 2>&1
		compare_output "$expected" "$TMP_DIR/$label-c$n.tsv" "$label, --number-of-cores $n"
	done
}

run_case single_nw_trkb "$DATA/expected_single_nw_trkb.tsv" \
	-s -g "$DATA/single_4000.fasta"
run_case unpaired_r_nw_trkb "$DATA/expected_unpaired_r_nw_trkb.tsv" \
	-s -v -g "$DATA/paired_2000_2.fasta"
run_case paired_nw_trkb "$DATA/expected_paired_nw_trkb.tsv" \
	-p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
