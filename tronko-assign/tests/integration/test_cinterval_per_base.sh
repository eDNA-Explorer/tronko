#!/usr/bin/env bash
# --cinterval-per-base on the multi-tree fixture: each case reproduces its golden
# (goldens-cinterval-per-base, recorded at one thread) at every thread count; every read's
# window, as --trace-read reports it, is -c times the read's own A/C/G/T bases (counted here,
# independently); and the option is refused with pruning, whose threshold is a multiple of
# the raw -c. Environment: golden_lib.sh (TRONKO_TESTS_RECORD=1 records the goldens).
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
[[ $MATCH_CAP == 10 ]] || { echo "SKIP: --cinterval-per-base goldens exist for cap 10 only"; exit 77; }
mem3_index "$MT/multitree.fasta" || exit 2

C=0.06
COMMON=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 --Cinterval "$C" --cinterval-per-base)
GOLDENS="$MT/goldens-cinterval-per-base"
mkdir -p "$GOLDENS"
echo "cinterval-per-base: multi-tree fixture, -c $C, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# check_windows <trace> <fasta> [<fasta 2>]: every record's window is C x the read's bases
check_windows() {
	local trace=$1 bad
	shift
	bad=$(awk -F'\t' -v c="$C" '
		FNR == 1 { file++ }
		file <= nf && /^>/ { name = substr($1, 2); next }
		file <= nf { n = gsub(/[ACGTacgt]/, ""); bases[name] += n; next }
		$1 == "trace" && $3 > 0 {
			want = c * bases[$2]
			if ($7 - want > 1e-9 * (want + 1) || want - $7 > 1e-9 * (want + 1)) bad++
		}
		END { print bad + 0 }' nf=$# "$@" "$trace")
	echo "$bad"
}

run_case() {
	local c=$1 n out trace bad
	shift
	local fastas=()
	for a in "$@"; do [[ $a == *.fasta ]] && fastas+=("$a"); done
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		trace=$TMP_DIR/$c-c$n.trace
		if ! run_assign "$TMP_DIR/$c-c$n.log" "${COMMON[@]}" --number-of-cores "$n" "$@" \
			--trace-read '*' --trace-file "$trace" -o "$out"; then
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, --number-of-cores $n"
		bad=$(check_windows "$trace" "${fastas[@]}")
		if [[ $bad -ne 0 ]]; then
			echo "FAIL: $c, -C $n: $bad reads with a window other than $C x their bases" >&2
			FAILS=$((FAILS + 1))
		fi
	done
}

run_case mt_single_F -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R -s -v -g "$MT/single_R.fasta"
run_case mt_paired -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"

if run_assign "$TMP_DIR/refused.log" "${COMMON[@]}" --enable-pruning -s -g "$MT/single_F.fasta" \
	-o "$TMP_DIR/refused.tsv" 2>/dev/null; then
	echo "FAIL: --cinterval-per-base with --enable-pruning was not refused" >&2
	FAILS=$((FAILS + 1))
fi

((FAILS == 0)) || { echo "cinterval-per-base: $FAILS failure(s)" >&2; exit 1; }
echo "cinterval-per-base: PASS"
