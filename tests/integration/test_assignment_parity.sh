#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
ASSIGN=${TRONKO_ASSIGN_BIN:-"$REPO_ROOT/tronko-assign/tronko-assign"}
REFERENCE="$REPO_ROOT/tronko-build/example_datasets/single_tree/reference_tree.txt"
FASTA="$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
DATA="$REPO_ROOT/tests/data/assignment"
BASELINE_COMMIT=71f6ec3eed6ca7c0e11de256deb7d19e10d174fc
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/tronko-parity.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

hash_file() {
	if command -v sha256sum >/dev/null 2>&1; then
		sha256sum "$1" | awk '{print $1}'
	else
		shasum -a 256 "$1" | awk '{print $1}'
	fi
}

compare_output() {
	local expected=$1
	local actual=$2
	local label=$3
	if ! cmp -s "$expected" "$actual"; then
		echo "Assignment parity failed: $label" >&2
		diff -u "$expected" "$actual" >&2 || true
		return 1
	fi
	echo "PASS: $label (all seven TSV columns)"
}

if [[ ! -x "$ASSIGN" ]]; then
	echo "Missing executable: $ASSIGN" >&2
	exit 1
fi
for path in "$REFERENCE" "$FASTA" "$DATA/single.fasta" \
	"$DATA/paired_1.fasta" "$DATA/paired_2.fasta"; do
	if [[ ! -f "$path" ]]; then
		echo "Missing parity fixture: $path" >&2
		exit 1
	fi
done

echo "baseline_commit=$BASELINE_COMMIT"
echo "reference_sha256=$(hash_file "$REFERENCE")"
echo "fasta_sha256=$(hash_file "$FASTA")"
echo "single_query_sha256=$(hash_file "$DATA/single.fasta")"
echo "paired_1_query_sha256=$(hash_file "$DATA/paired_1.fasta")"
echo "paired_2_query_sha256=$(hash_file "$DATA/paired_2.fasta")"
echo "cores=1 batch_lines=50000 match_cap=10 cinterval=5 score_constant=0.01"

"$ASSIGN" -r -f "$REFERENCE" -a "$FASTA" -6 -C 1 -c 5 \
	-s -g "$DATA/single.fasta" -o "$TMP_DIR/single-nw.tsv" -w
compare_output "$DATA/expected_single_nw.tsv" "$TMP_DIR/single-nw.tsv" \
	"single-end Needleman-Wunsch"

"$ASSIGN" -r -f "$REFERENCE" -a "$FASTA" -6 -C 1 -c 5 \
	-s -g "$DATA/single.fasta" -o "$TMP_DIR/single-wfa.tsv"
compare_output "$DATA/expected_single_wfa.tsv" "$TMP_DIR/single-wfa.tsv" \
	"single-end WFA"

"$ASSIGN" -r -f "$REFERENCE" -a "$FASTA" -6 -C 1 -c 5 \
	-p -1 "$DATA/paired_1.fasta" -2 "$DATA/paired_2.fasta" \
	-o "$TMP_DIR/paired-wfa.tsv"
compare_output "$DATA/expected_paired_wfa.tsv" "$TMP_DIR/paired-wfa.tsv" \
	"paired-end WFA"
