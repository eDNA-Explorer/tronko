#!/usr/bin/env bash
# Multi-batch parity of tronko-assign with the production invocation.
#
# tronko-assign reads its input in batches of -L/2 reads or pairs and runs BWA once per batch.
# The goldens here were produced by production tronko-assign (71f6ec3) at --number-of-cores 1
# with -L 400, which splits the 2,000-pair fixtures of tests/data/assignment into ten batches.
# A change to what is kept across batches (the BWA index, the leaf map, per-batch estimates)
# must reproduce them byte for byte. The output depends on the batch size (BWA's insert-size
# estimate and tie-break ids are per batch), so these goldens are valid only at -L 400; the
# default-batch goldens are those of tests/integration/test_assignment_production_parity.sh.
#
# Goldens: tests/data/multibatch[-cap25]/expected_<case>_L400.tsv. Environment: see
# golden_lib.sh; TRONKO_FIXTURE_DIR replaces tests/data/assignment at the repository root.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=golden_lib.sh
source "$SCRIPT_DIR/golden_lib.sh"
DATA=${TRONKO_FIXTURE_DIR:-"$REPO_ROOT/tests/data/assignment"}
GOLDENS=$(golden_dir "$SCRIPT_DIR/../data/multibatch")
FASTA="$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
LINES=400

for path in "$DATA/reference_tree.trkb" "$DATA/paired_2000_1.fasta" "$DATA/paired_2000_2.fasta"; do
	[[ -f "$path" ]] || { echo "SKIP: missing fixture $path (set TRONKO_FIXTURE_DIR)"; exit 77; }
done

mem3_index "$FASTA" || exit 2
COMMON=(-r -f "$DATA/reference_tree.trkb" -a "$FASTA" -w -6 --Cinterval 10 -L "$LINES")
echo "multi-batch fixture: cap $MATCH_CAP, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# run_case <label> <min-batches> <tronko args...>
run_case() {
	local label=$1 min_batches=$2 n batches out log
	shift 2
	for n in $CORES; do
		out=$TMP_DIR/$label-c$n.tsv log=$TMP_DIR/$label-c$n.log
		run_assign "$log" "${COMMON[@]}" --number-of-cores "$n" "$@" -o "$out" || { FAILS=$((FAILS + 1)); continue; }
		# In paired mode the aligner reports its insert-size estimate once per batch (the vendored
		# BWA as "[M::mem_pestat] ... orientation FF", BWA-MEM3 as "[PE] ... orientation FF").
		batches=$(grep -c 'orientation FF' "$log" || true)
		if ((batches < min_batches)); then
			echo "FAIL: $label ran $batches batches, expected at least $min_batches" >&2
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$label" "$GOLDENS/expected_${label}_L$LINES.tsv" "$out" "$label, -L $LINES, --number-of-cores $n"
	done
}

run_case paired_nw_trkb 10 -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
run_case unpaired_r_nw_trkb 0 -s -v -g "$DATA/paired_2000_2.fasta"
((FAILS == 0)) || { echo "multi-batch fixture: $FAILS failure(s)" >&2; exit 1; }
