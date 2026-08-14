#!/bin/bash
# Run one {pre,post} x cell tronko assignment at production parameters.
#   run_cell.sh <cell_dir> <db_dir> <marker_fasta_name> <binary: pre|post>
# Production parameters are also the compiled-in defaults at both commits;
# cores=1 because tronko is ~3.2% non-deterministic above -C 1.
set -uo pipefail

CELL=$1; DB=$2; FASTA=$3; BIN=$4
OUT="$CELL/${BIN}.tsv"

# A non-empty output is not proof of completion: a killed run leaves a partial
# TSV behind. Only skip when the row count matches the query count exactly
# (+1 for the header), otherwise discard the partial and redo the cell.
WANT=$(( $(grep -c '^>' "$CELL/pF.fasta") + 1 ))
if [ -s "$OUT" ]; then
  HAVE=$(wc -l < "$OUT" | tr -d ' ')
  if [ "$HAVE" -eq "$WANT" ]; then
    echo "### $(basename "$CELL")/$BIN already done ($HAVE lines)"
    exit 0
  fi
  echo "### $(basename "$CELL")/$BIN partial ($HAVE/$WANT lines) -- redoing"
  mv "$OUT" "$OUT.partial.$(date '+%H%M%S')"
fi

echo "### $(basename "$CELL")/$BIN start $(date '+%H:%M:%S')"
/Users/ryanmartin/tronko-ab/bin/tronko-assign."$BIN" -r \
  -f "$DB/reference_tree.trkb" \
  -a "$DB/$FASTA" \
  -p -z -w \
  -1 "$CELL/pF.fasta" \
  -2 "$CELL/pR.fasta" \
  -6 --number-of-cores 1 --Cinterval 0.02 --aligner minimap2 \
  -u 0.0001 --max-leaf-matches 10 \
  --best-leaf-threshold -0.1 --best-leaf-max-votes 10 \
  -o "$OUT" > "$CELL/${BIN}.log" 2>&1
rc=$?
HAVE=$(wc -l < "$OUT" 2>/dev/null | tr -d ' ' || echo 0)
echo "### $(basename "$CELL")/$BIN done $(date '+%H:%M:%S') rc=$rc lines=$HAVE/$WANT"
if [ "$rc" -eq 0 ] && [ "$HAVE" -ne "$WANT" ]; then
  echo "### $(basename "$CELL")/$BIN INCOMPLETE despite rc=0 -- treating as failure"
  exit 1
fi
exit $rc
