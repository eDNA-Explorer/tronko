#!/usr/bin/env bash
# Regenerate the input files of the tronko-assign test fixtures from the repository's example
# dataset, tronko-build/example_datasets/single_tree/, and compare each with the committed file.
#
# Usage: tronko-assign/tests/tools/make_fixture_inputs.sh <out-dir>
#
# Writes under <out-dir>, with the repository's layout:
#   tests/data/assignment/
#       reference_tree.trkb            tronko-convert of the example reference_tree.txt
#       single_4000.fasta, paired_2000_1.fasta, paired_2000_2.fasta
#                                      sim_reads.py (seed 7)
#   tronko-assign/tests/data/multitree/
#       reference_tree.trkb            tronko-convert of multitree_fixture.py's three-tree reference
#       multitree.fasta, single_F.fasta, single_R.fasta
#                                      multitree_fixture.py (seed 7)
#       multitree.fasta.amb .ann .bwt .pac .sa
#                                      BWA index, built by tronko-assign in a run without -6
#   tronko-assign/tests/data/gaps/
#       mate-rescue/mr_1.fasta, mr_2.fasta              gap_fixtures.py mate-rescue (seed 11)
#       read-content/rc_single.fasta, rc_1.fasta, rc_2.fasta
#                                                       gap_fixtures.py read-content (seed 13)
#       prod-names/pn_*.fasta          the tests/data/assignment/ reads renamed as the pipeline
#                                      names them for primer 12S_MiFish_U (index from 0)
#       prod-cmdline/*.fasta.zst       the tests/data/assignment/ reads compressed as the pipeline
#                                      writes them: zstd -3 --no-check from standard input (one
#                                      frame, no checksum, no content size)
# It then compares every file with the committed one and prints IDENTICAL or DIFFERENT for each;
# it exits 1 if any differs.
#
# tronko-convert and tronko-assign are built from this checkout in a temporary directory (plain
# `make`; nothing is built or written inside the checkout). Needs bash, gcc, make, zlib and libzstd
# headers, the zstd command, python3 and awk.
#
# Goldens are not made here. They are made by the test scripts' record mode
# (TRONKO_TESTS_RECORD=1, tronko-assign/tests/README.md, "Regenerating the goldens") with
# tronko-assign built from production commit 71f6ec3, at --number-of-cores 1.
set -euo pipefail

[[ $# -eq 1 ]] || { sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 2; }
OUT=$1
TOOLS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$TOOLS/../../.." && pwd)
EXAMPLE=$ROOT/tronko-build/example_datasets/single_tree
[[ -f $EXAMPLE/reference_tree.txt && -f $EXAMPLE/Charadriiformes.fasta ]] ||
	{ echo "example dataset not found in $EXAMPLE" >&2; exit 2; }

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
A=$OUT/tests/data/assignment
MT=$OUT/tronko-assign/tests/data/multitree
GAPS=$OUT/tronko-assign/tests/data/gaps
mkdir -p "$A" "$MT" "$GAPS"/{mate-rescue,read-content,prod-names,prod-cmdline}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "== build tronko-convert and tronko-assign from $ROOT (in $WORK)"
cp -a "$ROOT/tronko-convert" "$WORK/"
make -C "$WORK/tronko-convert" clean >/dev/null
make -C "$WORK/tronko-convert" >"$WORK/build-convert.log" 2>&1 ||
	{ cat "$WORK/build-convert.log" >&2; exit 1; }
CONVERT=$WORK/tronko-convert/tronko-convert
tar -C "$ROOT" --exclude=tronko-assign/tests --exclude=tronko-assign/carquet -cf - tronko-assign | tar -C "$WORK" -xf -
make -C "$WORK/tronko-assign" clean >/dev/null
make -C "$WORK/tronko-assign" >"$WORK/build-assign.log" 2>&1 ||
	{ cat "$WORK/build-assign.log" >&2; exit 1; }
ASSIGN=$WORK/tronko-assign/tronko-assign

echo "== tests/data/assignment"
"$CONVERT" -i "$EXAMPLE/reference_tree.txt" -o "$A/reference_tree.trkb" >"$WORK/convert-single.log" 2>&1
python3 "$TOOLS/sim_reads.py" "$EXAMPLE/Charadriiformes.fasta" "$A"

echo "== tronko-assign/tests/data/multitree"
python3 "$TOOLS/multitree_fixture.py" "$EXAMPLE" "$WORK/multitree"
"$CONVERT" -i "$WORK/multitree/reference_tree.txt" -o "$MT/reference_tree.trkb" >"$WORK/convert-multitree.log" 2>&1
cp "$WORK/multitree"/{multitree.fasta,single_F.fasta,single_R.fasta} "$MT/"
# The BWA index: tronko-assign builds it next to the FASTA when run without -6. Build it in a run
# of its own on ten reads (building it in-process advances the process's random state; the tests
# pass -6 and never build it).
rm -f "$MT"/multitree.fasta.{amb,ann,bwt,pac,sa}
mkdir -p "$WORK/idx"
head -20 "$MT/single_F.fasta" >"$WORK/idx/tiny.fasta"
(cd "$WORK/idx" && "$ASSIGN" -r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w \
	--Cinterval 10 --number-of-cores 1 -s -g tiny.fasta -o tiny.tsv >index-build.log 2>&1) ||
	{ cat "$WORK/idx/index-build.log" >&2; exit 1; }

echo "== tronko-assign/tests/data/gaps"
python3 "$TOOLS/gap_fixtures.py" mate-rescue "$EXAMPLE/Charadriiformes.fasta" "$GAPS/mate-rescue"
python3 "$TOOLS/gap_fixtures.py" read-content "$EXAMPLE/Charadriiformes.fasta" "$GAPS/read-content"
rename_reads() { awk -v pfx="$2" '/^>/ {print ">" pfx n++; next} {print}' "$1"; }
rename_reads "$A/paired_2000_1.fasta" 12S_MiFish_U_paired_F_ >"$GAPS/prod-names/pn_paired_F.fasta"
rename_reads "$A/paired_2000_2.fasta" 12S_MiFish_U_paired_R_ >"$GAPS/prod-names/pn_paired_R.fasta"
rename_reads "$A/single_4000.fasta" 12S_MiFish_U_unpaired_F_ >"$GAPS/prod-names/pn_unpaired_F.fasta"
rename_reads "$A/paired_2000_2.fasta" 12S_MiFish_U_unpaired_R_ >"$GAPS/prod-names/pn_unpaired_R.fasta"
for f in single_4000 paired_2000_1 paired_2000_2; do
	zstd -q -3 --no-check -c <"$A/$f.fasta" >"$GAPS/prod-cmdline/$f.fasta.zst"
done

echo "== compare with the committed files ($(zstd --version | head -1))"
different=0
while IFS= read -r p; do
	if cmp -s "$OUT/$p" "$ROOT/$p"; then r=IDENTICAL; else r=DIFFERENT; different=$((different + 1)); fi
	printf '%-10s %s\n' "$r" "$p"
done < <(cd "$OUT" && find tests tronko-assign -type f | LC_ALL=C sort)
((different == 0)) || { echo "$different file(s) differ from the committed ones" >&2; exit 1; }
echo "every file is identical to the committed one"
