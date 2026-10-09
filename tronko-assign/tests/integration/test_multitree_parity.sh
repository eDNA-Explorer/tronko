#!/usr/bin/env bash
# Parity of tronko-assign on a three-tree reference, where reads have several candidate trees.
#
# The repository's example reference is one tree, so every read there has one candidate tree and
# the multi-candidate code (candidate gathering, the vote tally over several trees, the
# multi-tree LCA) never runs. The three-tree fixture (tests/data/multitree) holds three copies of
# that tree: tree 0 is the original, tree 1 has only every second leaf in the BWA FASTA, tree 2
# has its taxonomy rotated by one line. The same sequence therefore sits in two or three trees:
# 1,247 of the 2,000 pairs have two candidate trees and 1,184 have votes in two, and the chimeric
# single-end reads (two or three segments of different leaves) have up to three. The goldens were
# made by production tronko-assign (71f6ec3) at --number-of-cores 1; every thread count in
# TRONKO_ASSIGN_CORES must reproduce them byte for byte.
#
# The single-end cases depend, before the SAM SEQ fix, on the bytes the linker places after two
# string literals in BWA's SAM writer (bwamem.c formats SEQ from unconverted bases and a NUL there
# ends the SAM text, dropping the read's later candidate trees). The reads align to the reverse
# strand, which is read correctly by production's Debian build and by aarch64 gcc builds; Ubuntu
# gcc 13.3 x86-64 builds read a NUL for a reverse-strand G and differ in about 620 rows.
#
# Fixture: tests/data/multitree/ (reference_tree.trkb, multitree.fasta (its BWA-MEM3 index is built at test time),
# single_F.fasta, single_R.fasta; the pairs are tests/data/assignment/paired_2000_{1,2}.fasta
# at the repository root). Goldens: tests/data/multitree/goldens[-cap25]/expected_<case>.tsv.
# Environment: see golden_lib.sh; TRONKO_MULTITREE_DIR replaces the fixture directory.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=golden_lib.sh
source "$SCRIPT_DIR/golden_lib.sh"
MT=${TRONKO_MULTITREE_DIR:-"$TA/tests/data/multitree"}
DATA="$REPO_ROOT/tests/data/assignment"

for f in reference_tree.trkb multitree.fasta single_F.fasta single_R.fasta; do
	[[ -f "$MT/$f" ]] || { echo "SKIP: missing fixture file $MT/$f"; exit 77; }
done
for f in "$DATA/paired_2000_1.fasta" "$DATA/paired_2000_2.fasta"; do
	[[ -f $f ]] || { echo "SKIP: missing fixture file $f"; exit 77; }
done

mem3_index "$MT/multitree.fasta" || exit 2
COMMON=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 --Cinterval 10)
GOLDENS=$(golden_dir "$MT/goldens")
echo "multi-tree fixture: cap $MATCH_CAP, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# run_case <case> <tronko mode args...>
run_case() {
	local c=$1 n out
	shift
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		run_assign "$TMP_DIR/$c-c$n.log" "${COMMON[@]}" --number-of-cores "$n" "$@" -o "$out" || { FAILS=$((FAILS + 1)); continue; }
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, --number-of-cores $n"
	done
}

run_case mt_single_F -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R -s -v -g "$MT/single_R.fasta"
run_case mt_paired -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
((FAILS == 0)) || { echo "multi-tree fixture: $FAILS failure(s)" >&2; exit 1; }
