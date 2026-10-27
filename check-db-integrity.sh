#!/bin/bash
# Integrity spot-check for a tronko-build -y output directory.
#
# Intended for periodic sampling during a long concurrent (-J N) build, and as the
# accept/reject gate at the end. These are invariants that must hold regardless of
# thread scheduling -- unlike exact output equality, which does NOT hold once
# VeryFastTree runs with more than ~2 threads (its multithreaded floating-point
# reduction order varies run to run, changing tree topology and therefore where the
# min-variance split lands). That nondeterminism is pre-existing and independent of
# cluster-level concurrency; do not use "output differs from a -J 1 rerun" as a bug
# signal. Use these invariants instead.
#
# Usage: check-db-integrity.sh <partitions-directory> [expected-total-seqs] [sample-n]
#   sample-n: check only the N most recent partitions (default: all)
set -u

DB="${1:?usage: check-db-integrity.sh <partitions-directory> [expected-total-seqs] [sample-n]}"
EXPECTED="${2:-}"
SAMPLE="${3:-0}"

[ -f "$DB/final_partitions.txt" ] || { echo "no final_partitions.txt in $DB (build not finalized)" >&2; exit 1; }

mapfile -t PARTS < "$DB/final_partitions.txt"
if [ "$SAMPLE" -gt 0 ] 2>/dev/null; then
  mapfile -t PARTS < <(printf '%s\n' "${PARTS[@]}" | sort -n | tail -"$SAMPLE")
fi

mismatch=0; small=0; missing=0; ok=0
acc_file=$(mktemp)
trap 'rm -f "$acc_file"' EXIT

for n in "${PARTS[@]}"; do
  [ -n "$n" ] || continue
  msa="$DB/partition${n}_MSA.fasta"
  tax="$DB/partition${n}_taxonomy.txt"
  tre="$DB/RAxML_bestTree.partition${n}.reroot"
  fa="$DB/partition${n}.fasta"
  if [ ! -s "$msa" ] || [ ! -s "$tax" ] || [ ! -s "$tre" ]; then
    echo "MISSING/EMPTY triplet: partition$n"; missing=$((missing+1)); continue
  fi
  a=$(grep -c '^>' "$msa")
  b=$(grep -vc '^$' "$tax")
  # leaf count == commas + 1 for a Newick string
  c=$(( $(tr -cd ',' < "$tre" | wc -c) + 1 ))
  if [ "$a" != "$b" ] || [ "$a" != "$c" ]; then
    echo "MISMATCH partition$n: msa=$a taxonomy=$b tree_leaves=$c"
    mismatch=$((mismatch+1))
  else
    ok=$((ok+1))
  fi
  # createNewRoots refuses to split below 4 sequences (tronko-build.c:1097)
  [ "$a" -lt 4 ] && { echo "UNDERSIZED partition$n: $a seqs"; small=$((small+1)); }
  [ -s "$fa" ] && grep '^>' "$fa" >> "$acc_file"
done

total=$(wc -l < "$acc_file")
uniq=$(sort -u "$acc_file" | wc -l)

echo
echo "partitions checked : ${#PARTS[@]}"
echo "consistent         : $ok"
echo "mismatched         : $mismatch"
echo "missing/empty      : $missing"
echo "undersized (<4)    : $small"
echo "sequences          : $total total, $uniq unique"

rc=0
if [ "$mismatch" -gt 0 ] || [ "$missing" -gt 0 ] || [ "$small" -gt 0 ]; then rc=1; fi
if [ "$total" != "$uniq" ]; then
  echo "DUPLICATE SEQUENCES: $(( total - uniq )) duplicated across partitions"
  rc=1
fi
if [ -n "$EXPECTED" ] && [ "$SAMPLE" = 0 ]; then
  if [ "$uniq" != "$EXPECTED" ]; then
    echo "SEQUENCE COUNT MISMATCH: found $uniq, expected $EXPECTED"
    rc=1
  else
    echo "sequence count matches input ($EXPECTED)"
  fi
fi

[ "$rc" = 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"
exit "$rc"
