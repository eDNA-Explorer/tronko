#!/usr/bin/env bash
# The pinned aligner: BWA-MEM3 must give the stored SAM for each fixture batch below, so that a
# different BWA-MEM3 (another release, another build, another SIMD tier or machine) cannot change
# tronko-assign's candidate search unnoticed. tronko-assign runs as in the golden tests, through a
# stand-in bwa-mem3 that passes every command to the real one and keeps each batch's SAM, so the
# command line is tronko-assign's own. For every case and thread count, two checks:
#   fields  the Tronko fields of every record (QNAME FLAG RNAME POS CIGAR RNEXT: SAM columns 1-4, 6,
#           7; what the SAM parse reads) equal tests/data/pinned-aligner/<case>.fields.tsv
#   sam     the SHA-256 of the whole SAM without its header lines equals <case>.sam.sha256
# A difference in the fields changes Tronko's input; one in the SAM alone (a tag, MAPQ, ...)
# does not, but still means the aligner is not the one pinned.
#
#   mr_paired     620 pairs, one mate of each unseedable: mate rescue
#   n1_paired     the repository's 2,000 pairs read as one batch: every pair has the same spacing,
#                 so the insert-size spread is 0 (BWA-MEM3's guard, details N1)
#   mt_paired     the same pairs on the three-tree reference: hits in several trees
#   rc_single     460 single reads with one N, random sequence, or 30 to 60 bases
#
# Record: TRONKO_TESTS_RECORD=1 TRONKO_ASSIGN_CORES=1 writes the expected files. Environment: see
# golden_lib.sh; TRONKO_PINNED_CASES limits the cases.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=golden_lib.sh
source "$SCRIPT_DIR/golden_lib.sh"
EXPECT=$TA/tests/data/pinned-aligner
DATA=$REPO_ROOT/tests/data/assignment
GAPS=$TA/tests/data/gaps
MT=$TA/tests/data/multitree
FASTA=$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta
REFERENCE=$DATA/reference_tree.trkb
CASES=${TRONKO_PINNED_CASES:-"mr_paired n1_paired mt_paired rc_single"}
for f in "$FASTA" "$REFERENCE" "$DATA/paired_2000_1.fasta" "$GAPS/mate-rescue/mr_1.fasta" "$GAPS/read-content/rc_single.fasta" \
	"$MT/multitree.fasta" "$MT/reference_tree.trkb"; do
	[[ -f $f ]] || { echo "SKIP: missing fixture file $f"; exit 77; }
done
mem3_index "$FASTA" "$MT/multitree.fasta" || exit 2

REAL=$TRONKO_BWA_MEM3
FAKE=$TMP_DIR/bin/bwa-mem3
mkdir -p "$TMP_DIR/bin"
cat >"$FAKE" <<'EOF'
#!/usr/bin/env bash
[[ $1 == mem ]] || exec "$REAL_BWA_MEM3" "$@"
out=
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do [[ ${args[i]} == -o ]] && out=${args[i+1]}; done
"$REAL_BWA_MEM3" "$@" || exit $?
n=$(find "$KEEP_DIR" -name 'batch*.sam' | wc -l)
cp "$out" "$KEEP_DIR/batch$((n + 1)).sam"
EOF
chmod +x "$FAKE"
export REAL_BWA_MEM3=$REAL

case_args() {
	case $1 in
	mr_paired) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -w -6 --Cinterval 10 -p -z -1 "$GAPS/mate-rescue/mr_1.fasta" -2 "$GAPS/mate-rescue/mr_2.fasta") ;;
	n1_paired) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -w -6 --Cinterval 10 -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta") ;;
	mt_paired) ARGS=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 --Cinterval 10 -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta") ;;
	rc_single) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -w -6 --Cinterval 10 -s -g "$GAPS/read-content/rc_single.fasta") ;;
	*) echo "unknown pinned-aligner case: $1" >&2; return 1 ;;
	esac
}

echo "pinned aligner: $("$REAL" version | head -1), $REAL (SHA-256 $(sha256sum "$REAL" | cut -c1-16)...), threads $CORES$([[ $RECORD == 1 ]] && echo ', recording')"
for c in $CASES; do
	case_args "$c" || exit 2
	for n in $CORES; do
		keep=$TMP_DIR/keep-$c-$n
		mkdir -p "$keep"
		KEEP_DIR=$keep TRONKO_BWA_MEM3=$FAKE run_assign "$keep/log" "${ARGS[@]}" --number-of-cores "$n" -o "$keep/out.tsv" ||
			{ FAILS=$((FAILS + 1)); continue; }
		if [[ $(find "$keep" -name 'batch*.sam' | wc -l) != 1 ]]; then
			echo "FAIL: $c, --number-of-cores $n: expected one batch" >&2; FAILS=$((FAILS + 1)); continue
		fi
		grep -v '^@' "$keep/batch1.sam" >"$keep/body.sam"
		cut -f1-4,6,7 "$keep/body.sam" >"$keep/fields.tsv"
		sha=$(sha256sum <"$keep/body.sam" | cut -d' ' -f1)
		if [[ $RECORD == 1 ]]; then
			mkdir -p "$EXPECT"
			cp "$keep/fields.tsv" "$EXPECT/$c.fields.tsv"
			echo "$sha" >"$EXPECT/$c.sam.sha256"
			echo "RECORDED: $c ($(wc -l <"$keep/body.sam") records)"
			continue
		fi
		label="$c, --number-of-cores $n"
		if cmp -s "$EXPECT/$c.fields.tsv" "$keep/fields.tsv"; then
			echo "PASS: $label: Tronko fields"
		else
			echo "FAIL: $label: Tronko fields differ from $EXPECT/$c.fields.tsv ($(diff "$EXPECT/$c.fields.tsv" "$keep/fields.tsv" | grep -c '^<') records)" >&2
			diff "$EXPECT/$c.fields.tsv" "$keep/fields.tsv" | head -6 >&2
			FAILS=$((FAILS + 1))
		fi
		if [[ $sha == "$(cat "$EXPECT/$c.sam.sha256" 2>/dev/null)" ]]; then
			echo "PASS: $label: whole SAM"
		else
			echo "FAIL: $label: SAM SHA-256 $sha, expected $(cat "$EXPECT/$c.sam.sha256" 2>/dev/null || echo none)" >&2
			FAILS=$((FAILS + 1))
		fi
	done
done
((FAILS == 0)) || { echo "pinned aligner: $FAILS failure(s)" >&2; exit 1; }
