#!/bin/bash
# Standalone equivalent of build-tronko-db.sh Step 5/5 (Finalize).
#
# Needed when tronko-build was run directly rather than through build-tronko-db.sh --
# e.g. after resuming an interrupted partition build with -b/--sequential-clusters,
# which the wrapper cannot express (its getopts has no 'b' and it never passes -b
# through). Killing the wrapper to relaunch tronko-build directly means Step 5 never
# runs; this script performs it.
#
# Usage: finalize-db.sh <output-dir> <primer-name>
set -euo pipefail

OUTPUT_DIR="${1:?usage: finalize-db.sh <output-dir> <primer-name>}"
PRIMER="${2:?usage: finalize-db.sh <output-dir> <primer-name>}"

[ -d "$OUTPUT_DIR" ] || { echo "not a directory: $OUTPUT_DIR" >&2; exit 1; }
[ -f "$OUTPUT_DIR/final_partitions.txt" ] || {
  echo "ERROR: $OUTPUT_DIR/final_partitions.txt not found." >&2
  echo "tronko-build did not finish -- do not finalize a partial build." >&2
  exit 1
}
[ -f "$OUTPUT_DIR/reference_tree.txt" ] || {
  echo "ERROR: $OUTPUT_DIR/reference_tree.txt not found (tronko-build incomplete)." >&2
  exit 1
}

echo "=== Finalize: FASTA concat + BWA index ==="
echo "Concatenating raw partition FASTAs per final_partitions.txt..."

: > "$OUTPUT_DIR/${PRIMER}.fasta"
: > "$OUTPUT_DIR/${PRIMER}_taxonomy.txt"
missing=0
count=0
while read -r partnum; do
  [ -n "$partnum" ] || continue
  fa="$OUTPUT_DIR/partition${partnum}.fasta"
  tx="$OUTPUT_DIR/partition${partnum}_taxonomy.txt"
  if [ ! -f "$fa" ] || [ ! -f "$tx" ]; then
    echo "  MISSING: partition${partnum}" >&2
    missing=$((missing+1))
    continue
  fi
  cat "$fa" >> "$OUTPUT_DIR/${PRIMER}.fasta"
  cat "$tx" >> "$OUTPUT_DIR/${PRIMER}_taxonomy.txt"
  count=$((count+1))
done < "$OUTPUT_DIR/final_partitions.txt"

echo "  concatenated $count partitions ($missing missing)"
if [ "$missing" -gt 0 ]; then
  echo "ERROR: $missing partition(s) listed in final_partitions.txt are absent." >&2
  echo "Refusing to build an index over an incomplete database." >&2
  exit 1
fi

nseq=$(grep -c '^>' "$OUTPUT_DIR/${PRIMER}.fasta")
ntax=$(grep -vc '^$' "$OUTPUT_DIR/${PRIMER}_taxonomy.txt")
echo "  ${PRIMER}.fasta: $nseq sequences; taxonomy: $ntax lines"
[ "$nseq" = "$ntax" ] || { echo "ERROR: FASTA/taxonomy count mismatch." >&2; exit 1; }

echo "Building BWA index..."
bwa index "$OUTPUT_DIR/${PRIMER}.fasta"

# tronko-convert is built in-tree but is not necessarily on PATH; the wrapper's
# `command -v` check silently skips the conversion in that case, leaving a database
# tronko-assign cannot load. Look in the repo before giving up.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERT=""
if command -v tronko-convert &> /dev/null; then
  CONVERT="tronko-convert"
elif [ -x "$HERE/tronko-convert/tronko-convert" ]; then
  CONVERT="$HERE/tronko-convert/tronko-convert"
fi

if [ -n "$CONVERT" ]; then
  echo "Converting reference tree to .trkb format ($CONVERT)..."
  "$CONVERT" -i "$OUTPUT_DIR/reference_tree.txt" -o "$OUTPUT_DIR/reference_tree.trkb"
  [ -s "$OUTPUT_DIR/reference_tree.trkb" ] || { echo "ERROR: .trkb conversion produced no output" >&2; exit 1; }
  gzip -f "$OUTPUT_DIR/reference_tree.txt"
else
  echo "WARNING: tronko-convert not found (not on PATH, not at $HERE/tronko-convert/)." >&2
  echo "         Keeping .txt.gz only -- tronko-assign may require .trkb." >&2
  gzip -f "$OUTPUT_DIR/reference_tree.txt"
fi

echo "Finalize complete."
