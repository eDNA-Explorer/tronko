#!/bin/bash
# Repair orphaned top-level clusters in a FINISHED tronko database.
#
# See check-orphaned-partitions.sh for the root cause. This script performs exactly the
# materialisation that tronko-build's two-step branch already does at
# tronko-build.c:2125-2141 (`if (strcmp(directory, opt.readdir)==0)` -> copy the
# cluster's files into the partitions dir), and then additionally gives them real
# partition numbers so they land in final_partitions.txt and the DB FASTA.
#
# It does NOT touch reference_tree.txt/.trkb: those already contain every leaf
# tree at the correct index (tree index == tree_list.txt line number), so no
# rebuild of the expensive partitioning stage is ever required.
#
# Safe to run only after the build has exited. Idempotent: a repaired DB reports
# 0 orphans and the script becomes a no-op.
#
# Usage: repair-orphaned-partitions.sh <db_dir> [--dry-run]

set -euo pipefail

DB="${1:?usage: $0 <db_dir> [--dry-run]}"
DRY=0
[[ "${2:-}" == "--dry-run" ]] && DRY=1

# Absolute from here on. The rewritten tree_list.txt lines below embed $DB, and
# tronko-build writes every other line as an absolute path -- so invoking this
# with a relative path (e.g. from the repo root) leaves a handful of
# repo-root-relative entries among thousands of absolute ones, which resolve
# only from the directory the repair happened to be run from.
DB=$(cd "$DB" 2>/dev/null && pwd) || { echo "ERROR: no such directory: $1" >&2; exit 1; }

TL="$DB/tree_list.txt"
FP="$DB/final_partitions.txt"

[[ -f "$TL" ]] || { echo "ERROR: no tree_list.txt in $DB" >&2; exit 1; }
[[ -f "$FP" ]] || { echo "ERROR: no final_partitions.txt in $DB" >&2; exit 1; }

if pgrep -f "build-tronko-db.sh.*$DB" >/dev/null 2>&1 || pgrep -f "tronko-build .*$DB" >/dev/null 2>&1; then
    echo "ERROR: a build is still running against $DB — refusing to modify it." >&2
    exit 1
fi

# The DB FASTA is the one that is neither a partition file nor the marker alias.
PRIMER_FA=$(find "$DB" -maxdepth 1 -name '*.fasta' ! -name 'partition*' ! -name 'marker.fasta' ! -name 'input.fasta' | head -1)
[[ -n "$PRIMER_FA" ]] || { echo "ERROR: cannot find the DB FASTA in $DB" >&2; exit 1; }
PRIMER=$(basename "$PRIMER_FA" .fasta)
PRIMER_TAX="$DB/${PRIMER}_taxonomy.txt"
echo "DB     : $DB"
echo "primer : $PRIMER"

mapfile -t ORPHANS < <(grep -nvE 'RAxML_bestTree\.partition[0-9]+\.reroot' "$TL" || true)
if (( ${#ORPHANS[@]} == 0 )); then
    echo "Nothing to do: 0 orphaned clusters."
    exit 0
fi
echo "orphans: ${#ORPHANS[@]}"

# Fresh partition numbers, above every index already in use, so nothing collides.
NEXT=$(find "$DB" -maxdepth 1 -name 'partition*_MSA.fasta' \
        | grep -oE 'partition[0-9]+' | grep -oE '[0-9]+' | sort -n | tail -1)
NEXT=$((NEXT + 1))
echo "new partition numbers start at: $NEXT"

MAP=$(mktemp)
PLAN=$(mktemp)
trap 'rm -f "$MAP" "$PLAN"' EXIT

idx=$NEXT
for entry in "${ORPHANS[@]}"; do
    lineno=${entry%%:*}
    path=${entry#*:}
    [[ -z "$path" ]] && continue
    id=$(basename "$path"); id=${id#RAxML_bestTree.}; id=${id%.reroot}
    dir=$(dirname "$path")
    for suffix in ".fasta" "_MSA.fasta" "_taxonomy.txt"; do
        [[ -f "$dir/$id$suffix" ]] || { echo "ERROR: missing source $dir/$id$suffix" >&2; exit 1; }
    done
    [[ -f "$dir/RAxML_bestTree.$id.reroot" ]] || { echo "ERROR: missing source tree for $id" >&2; exit 1; }
    printf "%s\t%s\t%s\t%s\n" "$lineno" "$id" "$dir" "$idx" >> "$PLAN"
    printf "%s\t%s\n" "$lineno" "$DB/RAxML_bestTree.partition$idx.reroot" >> "$MAP"
    printf "  line %-7s cluster %-8s -> partition%s (%s seqs)\n" \
           "$lineno" "$id" "$idx" "$(grep -c '^>' "$dir/$id.fasta")"
    idx=$((idx + 1))
done

if (( DRY )); then echo "(dry run — nothing written)"; exit 0; fi

echo "--- backing up list files ---"
cp -n "$TL" "$TL.orig" 2>/dev/null || true
cp -n "$FP" "$FP.orig" 2>/dev/null || true

echo "--- 1/5 copying cluster files into the DB as partition{N} ---"
while IFS=$'\t' read -r lineno id dir n; do
    cp "$dir/$id.fasta"                  "$DB/partition$n.fasta"
    cp "$dir/${id}_MSA.fasta"            "$DB/partition${n}_MSA.fasta"
    cp "$dir/${id}_taxonomy.txt"         "$DB/partition${n}_taxonomy.txt"
    cp "$dir/RAxML_bestTree.$id.reroot"  "$DB/RAxML_bestTree.partition$n.reroot"
done < "$PLAN"

echo "--- 2/5 rewriting tree_list.txt in place (line order is the tree index) ---"
awk -v mapfile="$MAP" '
  BEGIN{ while((getline l < mapfile) > 0){ split(l, a, "\t"); m[a[1]] = a[2] } }
  { if (FNR in m) print m[FNR]; else print }
' "$TL" > "$TL.new"
[[ "$(wc -l < "$TL.new")" == "$(wc -l < "$TL")" ]] || { echo "ERROR: tree_list line count changed" >&2; exit 1; }
mv "$TL.new" "$TL"

echo "--- 3/5 regenerating final_partitions.txt from tree_list.txt ---"
grep -oE 'RAxML_bestTree\.partition[0-9]+\.reroot' "$TL" | grep -oE '[0-9]+' > "$FP.new"
mv "$FP.new" "$FP"

echo "--- 4/5 rebuilding ${PRIMER}.fasta + taxonomy in final_partitions order ---"
awk -v d="$DB" '{print d"/partition"$1".fasta"}'         "$FP" | xargs cat > "$PRIMER_FA.new"
awk -v d="$DB" '{print d"/partition"$1"_taxonomy.txt"}'  "$FP" | xargs cat > "$PRIMER_TAX.new"
mv "$PRIMER_FA.new"  "$PRIMER_FA"
mv "$PRIMER_TAX.new" "$PRIMER_TAX"

echo "--- 5/5 rebuilding BWA index ---"
bwa index "$PRIMER_FA"

echo "=== verification ==="
printf "  tree_list.txt        : %s\n" "$(wc -l < "$TL")"
printf "  final_partitions.txt : %s\n" "$(wc -l < "$FP")"
printf "  %s.fasta seqs : %s\n" "$PRIMER" "$(grep -c '^>' "$PRIMER_FA")"
printf "  %s_taxonomy lines: %s\n" "$PRIMER" "$(wc -l < "$PRIMER_TAX")"
printf "  remaining orphans    : %s\n" \
       "$(grep -cvE 'RAxML_bestTree\.partition[0-9]+\.reroot' "$TL" || true)"
