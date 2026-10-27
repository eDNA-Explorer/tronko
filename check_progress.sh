#!/usr/bin/env bash
set -euo pipefail

MARKER=CO1_mlCOIintF_Fol-degen-rev
BASE="$HOME/tronko-build-branch/databases/$MARKER/lca/ac/default"
OUT="$BASE/sp0.10"
LOG="$OUT/build.log"
NCLUSTERS=232
TOTAL_SEQS=3537341   # sum of seqs across all 232 merged clusters, fixed for this build

# Expected final partition count, calibrated from the archived sequential run:
# 95,889 partition MSAs for 820,964 sequences = 0.1168 partitions/seq.
EXPECTED_PARTITIONS=413000

bar(){ # $1 = percent
  local pct=$1 filled empty
  [ "$pct" -gt 100 ] && pct=100
  filled=$((pct / 2)); empty=$((50 - filled))
  printf '%s%s' "$(printf '#%.0s' $(seq 1 $filled 2>/dev/null))" \
                "$(printf '.%.0s' $(seq 1 $empty 2>/dev/null))"
}

started=$(grep -c '^m->numspec:' "$LOG" 2>/dev/null || true)
jobs=$(grep -oE '\-J [0-9]+' "$OUT/status.txt" 2>/dev/null | tail -1 | awk '{print $2}')
jobs=${jobs:-1}

# Partition count on disk is the honest progress signal under -J: clusters complete
# out of order, so "clusters started" tells you nothing about how much work is done.
parts=$(find "$OUT" -maxdepth 1 -name 'partition*_MSA.fasta' 2>/dev/null | wc -l)
ppct=$((parts * 100 / EXPECTED_PARTITIONS))

cpct=$((started * 100 / NCLUSTERS))

echo "Step 4 (tronko-build partitioning, -J $jobs concurrent):"
echo "[$(bar $ppct)] $parts/~$EXPECTED_PARTITIONS partitions written (~$ppct%)  <- primary progress"
echo "[$(bar $cpct)] $started/$NCLUSTERS clusters started ($cpct%)"

# Up to $jobs clusters are in flight at once, so "started" overstates completion.
# Report the in-flight set rather than pretending there is a single current cluster.
if [ "$started" -gt 0 ]; then
  inflight=$(grep '^m->numspec:' "$LOG" | tail -"$jobs" | awk '{print $2}' | paste -sd, - )
  echo "Up to $jobs clusters in flight; most recently started sizes: $inflight"
fi

# Live utilisation — the whole point of the -J change; if this drops, something stalled.
main_pid=$(pgrep -f "tronko-build -y -e .*merged" | head -1 || true)
if [ -n "$main_pid" ] && [ -r "/proc/$main_pid/status" ]; then
  threads=$(awk '/^Threads/{print $2}' "/proc/$main_pid/status")
  rss=$(awk '/^VmRSS/{printf "%.1f", $2/1048576}' "/proc/$main_pid/status")
  idle=$(vmstat 1 2 | tail -1 | awk '{print $15}')
  echo "Worker threads: $threads   main RSS: ${rss} GB   CPU idle: ${idle}%"
  echo "Pipeline children: $(pgrep -c famsa 2>/dev/null || echo 0) famsa, $(pgrep -c VeryFastTree 2>/dev/null || echo 0) VeryFastTree"
else
  echo "tronko-build not running (finished, or died -- check $OUT/status.txt)"
fi

# grep -c prints "0" and exits 1 when there are no matches, so use || true, not || echo 0
errs=$(grep -cE '^ERROR|Cannot open|Out of memory|slot limit' "$LOG" 2>/dev/null || true)
echo "Errors so far: ${errs:-0}"
