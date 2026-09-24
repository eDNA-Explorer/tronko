#!/usr/bin/env bash
# Records tronko-assign/tests/data/nw_fixture_alignments.txt.gz: the leaf and read of every fifth
# needleman_wunsch_align() call (the 1st, 6th, 11th, ...) made by the three fixture cases of
# tests/integration/test_assignment_production_parity.sh at the repository root (single, unpaired
# reverse, paired; 4,000, 2,000 and 4,000 calls) at one thread, one "leaf TAB read" line per
# call, in that order.
#
# Usage, from the top of a checkout:
#   bash tronko-assign/tests/tools/record_nw_alignments.sh <commit> <output.txt.gz> [<golden-commit>]
# <commit> is a commit with the original seq-align fill: 71f6ec3 or any later commit up to this
# change's parent (the committed file is 71f6ec3's; from the commit that removes F1 on, the
# recorded leaves differ, tests/data/PROVENANCE.md). The script exports that commit's tronko-assign/ and example reference and HEAD's
# fixture reads into a temporary directory, adds a recording hook to needleman_wunsch_align()
# there, builds it, runs the three cases, checks each output against the golden of
# <golden-commit> (default HEAD) and writes <output.txt.gz>. The goldens of HEAD are those of the
# aligner HEAD runs; for a commit that still runs the vendored BWA (71f6ec3), pass a commit whose
# fixture goldens are production's, such as d019e28 (the head of pr2/01-tests). For a
# commit that runs BWA-MEM3, the example reference's BWA-MEM3 index is built with HEAD's
# integration/mem3_lib.sh (TRONKO_BWA_MEM3, default HEAD's tronko-assign/bwa-mem3). The checkout
# itself is not touched.
set -euo pipefail

[[ $# -eq 2 || $# -eq 3 ]] || { echo "usage: $0 <commit> <output.txt.gz> [<golden-commit>]" >&2; exit 2; }
COMMIT=$1
OUT=$(realpath -m "$2")
GOLDEN_COMMIT=${3:-HEAD}
TOP=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/record-nw.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# the code and the example reference (with the vendored BWA's index where <commit> has one) from
# <commit>; the fixture reads from HEAD; the goldens from <golden-commit>
git -C "$TOP" archive "$COMMIT" tronko-assign tronko-build/example_datasets/single_tree | tar -x -C "$WORK"
git -C "$TOP" archive HEAD tests/data/assignment | tar -x -C "$WORK"
mkdir -p "$WORK/goldens"
git -C "$TOP" archive "$GOLDEN_COMMIT" tests/data/assignment | tar -x -C "$WORK/goldens"

# The hook: write "a TAB b" for every fifth call to the file named by NW_RECORD. One thread only.
SRC=$WORK/tronko-assign/needleman_wunsch.c
CALL='needleman_wunsch_align2(a, b, strlen(a), strlen(b), scoring, nw, result);'
[[ $(grep -cF "$CALL" "$SRC") -eq 1 ]] || { echo "needleman_wunsch_align() not as expected in $COMMIT" >&2; exit 1; }
sed -i "s/^\(\s*\)needleman_wunsch_align2(a, b, strlen(a), strlen(b), scoring, nw, result);/\1nw_record(a, b);\n&/" "$SRC"
sed -i '0,/^void needleman_wunsch_align(/s//static void nw_record(const char *a, const char *b);\n&/' "$SRC"
cat >>"$SRC" <<'EOF'

static void nw_record(const char *a, const char *b)
{
  static FILE *fp = NULL;
  static long calls = 0;
  if(calls == 0 && getenv("NW_RECORD")) fp = fopen(getenv("NW_RECORD"), "w");
  if(fp && calls % 5 == 0) fprintf(fp, "%s\t%s\n", a, b);
  calls++;
}
EOF
make -C "$WORK/tronko-assign" >"$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log" >&2; exit 1; }
if [[ -f $WORK/tronko-assign/bwamem3.c ]]; then
	# shellcheck source=../integration/mem3_lib.sh
	source "$TOP/tronko-assign/tests/integration/mem3_lib.sh"
	mem3_index "$WORK/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
fi

D=$WORK/tests/data/assignment
COMMON=(-r -f "$D/reference_tree.trkb" -a "$WORK/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
	-w -6 --Cinterval 10 --number-of-cores 1)
run() { # run <label> <tronko args...>
	local label=$1; shift
	NW_RECORD="$WORK/$label.txt" "$WORK/tronko-assign/tronko-assign" "${COMMON[@]}" "$@" \
		-o "$WORK/$label.tsv" >/dev/null 2>"$WORK/$label.log"
	cmp -s "$WORK/$label.tsv" "$WORK/goldens/tests/data/assignment/expected_$label.tsv" || { echo "$label: output differs from its golden" >&2; exit 1; }
}
run single_nw_trkb -s -g "$D/single_4000.fasta"
run unpaired_r_nw_trkb -s -v -g "$D/paired_2000_2.fasta"
run paired_nw_trkb -p -z -1 "$D/paired_2000_1.fasta" -2 "$D/paired_2000_2.fasta"

cat "$WORK"/single_nw_trkb.txt "$WORK"/unpaired_r_nw_trkb.txt "$WORK"/paired_nw_trkb.txt | gzip -9n >"$OUT"
echo "$OUT: $(gzip -dc "$OUT" | wc -l) alignments"
