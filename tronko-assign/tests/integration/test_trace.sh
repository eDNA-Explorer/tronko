#!/usr/bin/env bash
# --trace-read is read-only: on the multi-tree fixture (reads with several candidate trees),
# every case run with --trace-read '*' still reproduces its golden byte for byte at every
# thread count, writes exactly one well-formed record per output row, and a named read is the
# only one traced. Environment: golden_lib.sh.
set -uo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
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
TRACE_COLUMNS=16
echo "trace: multi-tree fixture, cap $MATCH_CAP, threads $CORES"

run_case() {
	local c=$1 n out trace rows records bad
	shift
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		trace=$TMP_DIR/$c-c$n.trace
		if ! run_assign "$TMP_DIR/$c-c$n.log" "${COMMON[@]}" --number-of-cores "$n" "$@" \
			--trace-read '*' --trace-file "$trace" -o "$out"; then
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c traced, --number-of-cores $n"
		rows=$(grep -vc '^Readname' "$out")
		records=$(grep -c '^trace	' "$trace")
		# 15 trace columns, then the output row (its own tab-separated fields)
		bad=$(awk -F'\t' -v min=$TRACE_COLUMNS '$1=="trace" && NF<min' "$trace" | wc -l)
		if [[ $records -ne $rows || $bad -ne 0 ]]; then
			echo "FAIL: $c, -C $n: $records trace records for $rows rows, $bad short" >&2
			FAILS=$((FAILS + 1))
		fi
	done
}

run_case mt_single_F -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R -s -v -g "$MT/single_R.fasta"
run_case mt_paired -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"

# A named read is the only one traced.
name=$(sed -n '1s/^>\([^ ]*\).*/\1/p' "$MT/single_F.fasta")
trace=$TMP_DIR/named.trace
if run_assign "$TMP_DIR/named.log" "${COMMON[@]}" --number-of-cores 4 -s -g "$MT/single_F.fasta" \
	--trace-read "$name" --trace-file "$trace" -o "$TMP_DIR/named.tsv"; then
	if [[ $(grep -c '^trace	' "$trace") -ne 1 ]] || ! grep -q "^trace	$name	" "$trace"; then
		echo "FAIL: --trace-read $name did not trace exactly that read" >&2
		FAILS=$((FAILS + 1))
	fi
else
	FAILS=$((FAILS + 1))
fi

((FAILS == 0)) || { echo "trace: $FAILS failure(s)" >&2; exit 1; }
echo "trace: PASS"
