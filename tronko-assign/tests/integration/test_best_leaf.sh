#!/usr/bin/env bash
# --best-leaf-margin on the multi-tree fixture (three copies of one tree, so every species has a
# leaf in each candidate tree): each case reproduces its golden (goldens-best-leaf, recorded at
# one thread) at every thread count, and for every read, from its --trace-read record:
#   - the rule fired exactly when L1 is voted (S(L1) >= M - W), leads the best leaf of every other
#     species by the margin (S(L2) is -inf when there is none) and, with --best-leaf-max-deficit,
#     is within it of a perfect match (the read's A/C/G/T bases, counted here, times log(1 - u));
#   - a read it fired on got one vote and is assigned to L1 (the row's tree and node);
#   - a read it did not fire on has the row of the default golden (goldens/, same -c).
# The option is refused with pruning, and --best-leaf-max-deficit without a margin. Environment:
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
[[ $MATCH_CAP == 10 ]] || { echo "SKIP: --best-leaf-margin goldens exist for cap 10 only"; exit 77; }
mem3_index "$MT/multitree.fasta" || exit 2

C=10      # the default goldens' -c
D=5.7     # about one mismatch at -u 0.01: log(3 x 0.99 / 0.01) = 5.69
T=40      # the deficit case's bound
U=0.01    # -u, tronko-assign's default
COMMON=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 --Cinterval "$C" --best-leaf-margin "$D")
GOLDENS="$MT/goldens-best-leaf"
DEFAULTS="$MT/goldens"
mkdir -p "$GOLDENS"
echo "best-leaf: multi-tree fixture, -c $C, margin $D, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# check_rule <trace> <default golden> <deficit bound or inf> <fasta>...: prints
# "<bad> <fired> <not fired> <near the boundary>"
check_rule() {
	local trace=$1 golden=$2 bound=$3
	shift 3
	awk -F'\t' -v d="$D" -v t="$bound" -v u="$U" -v nf=$# '
		FNR == 1 { file++ }
		file <= nf && /^>/ { name = substr($1, 2); next }
		file <= nf { bases[name] += gsub(/[ACGTacgt]/, ""); next }
		file == nf + 1 { row[$1] = $0; next }
		$1 != "trace" || $3 == 0 { next }
		{
			m = $5; w = $7; votes = $8; l1 = $11; l2 = $14; fired = $15
			lead = (l2 == "-") ? "inf" : l1 - l2
			deficit = bases[$2] * log(1 - u) - l1
			want = l1 >= m - w && (lead == "inf" || lead >= d) && (t == "inf" || deficit <= t)
			# a margin or deficit within 1e-9 of its bound could go either way in the last bit
			if ((lead != "inf" && lead - d < 1e-9 && d - lead < 1e-9) || (t != "inf" && deficit - t < 1e-9 && t - deficit < 1e-9)) { near++; next }
			if (fired != want) { bad++; next }
			out = ""
			for (i = 16; i <= NF; i++) out = out (i > 16 ? "\t" : "") $i
			split(out, r, "\t")
			if (fired) {
				n_fired++
				if (votes != 1 || r[6] ":" r[7] != $9) bad++
			} else {
				n_kept++
				if (out != row[$2]) bad++
			}
		}
		END { print bad + 0, n_fired + 0, n_kept + 0, near + 0 }' "$@" "$golden" "$trace"
}

run_case() {
	local c=$1 golden_case=$2 bound=$3 n out trace res
	shift 3
	local fastas=()
	for a in "$@"; do [[ $a == *.fasta ]] && fastas+=("$a"); done
	local extra=()
	[[ $bound != inf ]] && extra=(--best-leaf-max-deficit "$bound")
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		trace=$TMP_DIR/$c-c$n.trace
		if ! run_assign "$TMP_DIR/$c-c$n.log" "${COMMON[@]}" "${extra[@]}" --number-of-cores "$n" "$@" \
			--trace-read '*' --trace-file "$trace" -o "$out"; then
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, --number-of-cores $n"
		read -r bad fired kept near < <(check_rule "$trace" "$DEFAULTS/expected_$golden_case.tsv" "$bound" "${fastas[@]}")
		echo "  $c, -C $n: fired on $fired reads, kept $kept, $near at a bound"
		if [[ $bad -ne 0 || $fired -eq 0 || $kept -eq 0 ]]; then
			echo "FAIL: $c, -C $n: $bad reads break the rule ($fired fired, $kept kept)" >&2
			FAILS=$((FAILS + 1))
		fi
	done
}

run_case mt_single_F mt_single_F inf -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R mt_unpaired_R inf -s -v -g "$MT/single_R.fasta"
run_case mt_paired mt_paired inf -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"
run_case mt_single_F_deficit mt_single_F "$T" -s -g "$MT/single_F.fasta"

if run_assign "$TMP_DIR/refused.log" "${COMMON[@]}" --enable-pruning -s -g "$MT/single_F.fasta" \
	-o "$TMP_DIR/refused.tsv" 2>/dev/null; then
	echo "FAIL: --best-leaf-margin with --enable-pruning was not refused" >&2
	FAILS=$((FAILS + 1))
fi
if run_assign "$TMP_DIR/refused2.log" -r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 \
	--Cinterval "$C" --best-leaf-max-deficit "$T" -s -g "$MT/single_F.fasta" -o "$TMP_DIR/refused2.tsv" 2>/dev/null; then
	echo "FAIL: --best-leaf-max-deficit without --best-leaf-margin was not refused" >&2
	FAILS=$((FAILS + 1))
fi

((FAILS == 0)) || { echo "best-leaf: $FAILS failure(s)" >&2; exit 1; }
echo "best-leaf: PASS"
