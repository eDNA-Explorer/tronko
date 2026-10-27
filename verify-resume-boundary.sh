#!/bin/bash
# Compute a safe -b value for resuming an interrupted tronko-build -y run.
#
# tronko-build's -b skip logic is a pure integer comparison (tronko-build.c:1124 gates
# the disk write, :1160 gates the fork/compute). It never stats or validates what is
# already on disk. If -b is set at or above a partition whose triplet is incomplete --
# e.g. VeryFastTree was mid-write when the process was killed -- that partition is
# silently dropped from the database: exit code 0, one stderr line, N-1 trees.
#
# This script finds the highest partition number whose triplet is complete, and every
# triplet below it in the scan window is complete too, then applies a safety margin.
#
# Usage: verify-resume-boundary.sh <partitions-directory> [window] [margin]
set -u

DB="${1:?usage: verify-resume-boundary.sh <partitions-directory> [window] [margin]}"
WINDOW="${2:-50}"   # how many of the top partitions to inspect
MARGIN="${3:-10}"   # extra slots to back off past the last good one

[ -d "$DB" ] || { echo "not a directory: $DB" >&2; exit 1; }

# Note: a plain `ls partition*_MSA.fasta` glob dies with "Argument list too long"
# once the directory holds ~90k+ matches. find|grep is the pattern that works
# (same as repair-orphaned-partitions.sh:59-60).
mapfile -t TOP < <(find "$DB" -maxdepth 1 -name 'partition*_MSA.fasta' \
  | grep -oE 'partition[0-9]+' | grep -oE '[0-9]+' | sort -n | tail -"$WINDOW")

[ "${#TOP[@]}" -gt 0 ] || { echo "no partition*_MSA.fasta found in $DB" >&2; exit 1; }

complete(){ # $1 = partition number
  local n="$1" msa tax tre
  msa="$DB/partition${n}_MSA.fasta"
  tax="$DB/partition${n}_taxonomy.txt"
  tre="$DB/RAxML_bestTree.partition${n}.reroot"
  [ -s "$msa" ] && [ -s "$tax" ] && [ -s "$tre" ] || return 1
  # Newick must be non-trivial and terminated
  [ "$(wc -c < "$tre")" -ge 5 ] || return 1
  tr -d ' \n\r' < "$tre" | grep -q ';$' || return 1
  # MSA sequence count must match taxonomy line count
  local a b
  a=$(grep -c '^>' "$msa")
  b=$(grep -vc '^$' "$tax")
  [ "$a" = "$b" ] || return 1
  return 0
}

echo "Inspecting top ${#TOP[@]} partitions in $DB"
last_good=""
first_bad=""
for n in "${TOP[@]}"; do
  if complete "$n"; then
    last_good="$n"
  else
    [ -n "$first_bad" ] || first_bad="$n"
    echo "  INCOMPLETE: partition$n"
    last_good=""   # anything at or above an incomplete triplet is untrustworthy
  fi
done

if [ -n "$first_bad" ]; then
  # Back off below the lowest incomplete partition seen.
  base=$(( first_bad - 1 ))
else
  base="${TOP[-1]}"
fi

safe=$(( base - MARGIN ))
[ "$safe" -lt 0 ] && safe=0

echo
echo "highest partition on disk : ${TOP[-1]}"
[ -n "$first_bad" ] && echo "lowest incomplete         : $first_bad"
echo "safety margin             : $MARGIN"
echo
echo "  -b $safe"
echo
echo "Re-run of the ~$MARGIN partitions above this boundary is expected and cheap."
echo "Pair with --sequential-clusters = (clusters started), which you should ROUND UP:"
echo "  grep -c '^m->numspec:' build.log"
echo "Undercounting is the dangerous direction (an already-started cluster would be"
echo "dispatched to the worker pool and could collide with slots already on disk);"
echo "overcounting only costs a few clusters replayed single-threaded."
