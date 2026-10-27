#!/usr/bin/env bash
# --score-full-read on the multi-tree fixture: each case reproduces its golden
# (goldens-score-full-read, recorded at one thread) at every thread count and on the legacy path
# (TRONKO_NODESTORE=0), whose sentinel test is getscore_Arr's own; and, read by read against a
# run without the option (both from --trace-read), the option only adds leading overhang rows of
# log(0.01) each: no read's best score rises, and a read with one candidate tree loses a whole
# number of them (k x log(0.01), k >= 0, the same k on every node of its tree). Some reads must
# lose rows and some none. The NW (-w) and WFA paths both build the rows. Environment:
# golden_lib.sh (TRONKO_TESTS_RECORD=1 records the goldens).
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
# Recorded from a cap-10 build only.
[[ $MATCH_CAP == 10 ]] || { echo "SKIP: --score-full-read goldens exist for cap 10 only"; exit 77; }
mem3_index "$MT/multitree.fasta" || exit 2

BASE=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -6 --Cinterval 10)
GOLDENS="$MT/goldens-score-full-read"
mkdir -p "$GOLDENS"
echo "score-full-read: multi-tree fixture, -c 10, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# check_rows <trace without> <trace with>: prints "<bad> <reads with rows added> <reads without>"
check_rows() {
	awk -F'\t' '
		BEGIN { L = log(0.01) }
		$1 != "trace" || $3 == 0 { next }
		FNR == NR { m[$2] = $5; n[$2] = $3; next }
		!($2 in m) { bad++; next }
		{
			d = $5 - m[$2]; tol = 1e-9 * ((m[$2] < 0 ? -m[$2] : m[$2]) + 1)
			if (d > tol) { bad++; next }
			if ($3 != 1 || n[$2] != 1) next
			k = d / L; r = int(k + 0.5)
			if (k - r > 1e-6 || r - k > 1e-6) { bad++; next }
			if (r > 0) added++; else none++
		}
		END { print bad + 0, added + 0, none + 0 }' "$1" "$2"
}

run_case() {
	local c=$1 n out trace base bad added none
	shift
	base=$TMP_DIR/$c-without.trace
	if ! run_assign "$TMP_DIR/$c-without.log" "${BASE[@]}" --number-of-cores 4 "$@" \
		--trace-read '*' --trace-file "$base" -o "$TMP_DIR/$c-without.tsv"; then
		FAILS=$((FAILS + 1))
		return
	fi
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		trace=$TMP_DIR/$c-c$n.trace
		if ! run_assign "$TMP_DIR/$c-c$n.log" "${BASE[@]}" --score-full-read --number-of-cores "$n" "$@" \
			--trace-read '*' --trace-file "$trace" -o "$out"; then
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, --number-of-cores $n"
		read -r bad added none < <(check_rows "$base" "$trace")
		echo "  $c, -C $n: one-tree reads with overhang rows added $added, without $none"
		if [[ $bad -ne 0 || $added -eq 0 || $none -eq 0 ]]; then
			echo "FAIL: $c, -C $n: $bad reads changed otherwise ($added with rows added, $none without)" >&2
			FAILS=$((FAILS + 1))
		fi
	done
	[[ $RECORD == 1 ]] && return
	out=$TMP_DIR/$c-legacy.tsv
	if TRONKO_NODESTORE=0 run_assign "$TMP_DIR/$c-legacy.log" "${BASE[@]}" --score-full-read --number-of-cores 4 "$@" -o "$out"; then
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, legacy path, --number-of-cores 4"
	else
		FAILS=$((FAILS + 1))
	fi
}

run_case mt_single_F -w -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R -w -s -v -g "$MT/single_R.fasta"
run_case mt_paired -w -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
run_case mt_single_F_wfa -s -g "$MT/single_F.fasta"

((FAILS == 0)) || { echo "score-full-read: $FAILS failure(s)" >&2; exit 1; }
echo "score-full-read: PASS"
