#!/usr/bin/env bash
# Run ONLY Step 4 (tronko-build partitioning) against an already-cached merged/ dir.
#
# build-tronko-db.sh cannot express this: its getopts has no 'b', it never passes
# -b or --sequential-clusters through, and it hardcodes -c $((THREADS/3)). Steps
# 1-3 (AncestralClust, per-cluster FAMSA/trees, merge+renumber) must already be
# cached in the merged dir. Run finalize-db.sh afterwards for Step 5.
#
# Usage: run-step4.sh <merged-dir> <output-dir> [sp-threshold] [-J N] [-c N] [extra tronko-build args...]
set -uo pipefail

MERGED="${1:?usage: run-step4.sh <merged-dir> <output-dir> [sp] [-J N] [-c N] ...}"
OUT="${2:?usage: run-step4.sh <merged-dir> <output-dir> [sp] [-J N] [-c N] ...}"
SP="${3:-0.10}"
shift 3 2>/dev/null || shift $#

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$REPO/tronko-build/tronko-build"
export PATH="$REPO/bin:$REPO/tronko-build:$REPO/tronko-convert:/usr/local/bin:$PATH"

[ -x "$BIN" ] || { echo "ERROR: $BIN not built (run: cd tronko-build && make)" >&2; exit 1; }
[ -d "$MERGED" ] || { echo "ERROR: merged dir not found: $MERGED" >&2; exit 1; }

# -n must equal the MSA count EXACTLY. readFilesInDir writes msa_files[file_count]
# with file_count unbounded by number_of_partitions (readreference.c:58,95-107), so a
# count higher than -n is a heap buffer overflow, and a lower one leaves uninitialised
# filenames. This is a correctness guard, not a convenience check.
NCLUST=$(find "$MERGED" -maxdepth 1 -name '*_MSA.fasta' | wc -l)
[ "$NCLUST" -gt 0 ] || { echo "ERROR: no *_MSA.fasta in $MERGED" >&2; exit 1; }

# Every MSA needs its taxonomy + tree companion, matched by stem.
missing=0
while IFS= read -r m; do
  n=$(basename "$m" _MSA.fasta)
  [ -s "$MERGED/${n}_taxonomy.txt" ] || { echo "  missing taxonomy: $n" >&2; missing=1; }
  [ -s "$MERGED/RAxML_bestTree.${n}.reroot" ] || { echo "  missing tree: $n" >&2; missing=1; }
done < <(find "$MERGED" -maxdepth 1 -name '*_MSA.fasta')
[ "$missing" = 0 ] || { echo "ERROR: incomplete input triplets in $MERGED" >&2; exit 1; }

mkdir -p "$OUT"
LOG="$OUT/build.log"

echo "=== Step 4/5: tronko-build partitioning ==="
echo "  input   : $MERGED ($NCLUST clusters)"
echo "  output  : $OUT"
echo "  binary  : $BIN"
echo "  extra   : $*"
printf 'STEP4_START %s\nclusters %s\nbinary %s\nargs %s\n' \
  "$(date -u +%FT%TZ)" "$NCLUST" "$BIN" "$*" >> "$OUT/status.txt"

# stdbuf -oL: redirected stdout is block-buffered, so a killed build loses its last
# few KB of log -- exactly the lines needed to derive --sequential-clusters for a
# resume. Line-buffering costs nothing here and makes the log trustworthy.
stdbuf -oL -eL "$BIN" -y \
  -e "$MERGED" \
  -n "$NCLUST" \
  -d "$OUT" \
  -s -u "$SP" \
  --tree-tool veryfasttree \
  -E --legacy-sp \
  "$@" \
  >> "$LOG" 2>&1
rc=$?

printf 'STEP4_EXIT %s rc=%s\n' "$(date -u +%FT%TZ)" "$rc" >> "$OUT/status.txt"
echo "Step 4 exited rc=$rc"
[ "$rc" = 0 ] && echo "Next: check-db-integrity.sh '$OUT' <expected-seqs> && finalize-db.sh '$OUT' <primer>"
exit $rc
