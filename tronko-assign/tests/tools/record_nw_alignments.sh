#!/usr/bin/env bash
# Records tronko-assign/tests/data/nw_fixture_alignments.txt.gz: the leaf and read of every fifth
# needleman_wunsch_align() call (the 1st, 6th, 11th, ...) made by the three fixture cases of
# tests/integration/test_assignment_production_parity.sh at the repository root (single, unpaired
# reverse, paired; 4,000, 2,000 and 4,000 calls) at one thread, one "leaf TAB read" line per
# call, in that order.
#
# Usage, from the top of a checkout:
#   bash tronko-assign/tests/tools/record_nw_alignments.sh <commit> <output.txt.gz>
# <commit> is a commit with the original seq-align fill: 71f6ec3 or any later commit up to this
# change's parent. The script exports that commit's tronko-assign/ and HEAD's fixtures into a
# temporary directory, adds a recording hook to needleman_wunsch_align() there, builds it, runs
# the three cases, checks each output against its golden and writes <output.txt.gz>. The checkout
# itself is not touched.
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 <commit> <output.txt.gz>" >&2; exit 2; }
COMMIT=$1
OUT=$(realpath -m "$2")
TOP=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/record-nw.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# the code from <commit>; the fixtures and the example reference with its BWA index from HEAD
git -C "$TOP" archive "$COMMIT" tronko-assign | tar -x -C "$WORK"
git -C "$TOP" archive HEAD tests/data/assignment tronko-build/example_datasets/single_tree | tar -x -C "$WORK"

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

D=$WORK/tests/data/assignment
COMMON=(-r -f "$D/reference_tree.trkb" -a "$WORK/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
	-w -6 --Cinterval 10 --number-of-cores 1)
run() { # run <label> <tronko args...>
	local label=$1; shift
	NW_RECORD="$WORK/$label.txt" "$WORK/tronko-assign/tronko-assign" "${COMMON[@]}" "$@" \
		-o "$WORK/$label.tsv" >/dev/null 2>"$WORK/$label.log"
	cmp -s "$WORK/$label.tsv" "$D/expected_$label.tsv" || { echo "$label: output differs from its golden" >&2; exit 1; }
}
run single_nw_trkb -s -g "$D/single_4000.fasta"
run unpaired_r_nw_trkb -s -v -g "$D/paired_2000_2.fasta"
run paired_nw_trkb -p -z -1 "$D/paired_2000_1.fasta" -2 "$D/paired_2000_2.fasta"

cat "$WORK"/single_nw_trkb.txt "$WORK"/unpaired_r_nw_trkb.txt "$WORK"/paired_nw_trkb.txt | gzip -9n >"$OUT"
echo "$OUT: $(gzip -dc "$OUT" | wc -l) alignments"
