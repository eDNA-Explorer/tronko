#!/usr/bin/env bash
# benchmark.sh — time ancestralclust on a dataset, report per-phase timings
# Usage: bash ancestralclust/benchmark.sh <input.fasta> [threads] [label]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BINARY="$SCRIPT_DIR/ancestralclust"
INPUT_FASTA="${1:?Usage: benchmark.sh <input.fasta> [threads] [label]}"
THREADS="${2:-1}"
LABEL="${3:-benchmark}"
OUTPUT_DIR="/tmp/ancestralclust_bench_$$"

if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: Build ancestralclust first (make)" >&2
    exit 1
fi

NUM_SEQS=$(grep -c "^>" "$INPUT_FASTA")
echo "=== AncestralClust Benchmark: $LABEL ==="
echo "Input:   $INPUT_FASTA ($NUM_SEQS sequences)"
echo "Threads: $THREADS"
echo "Output:  $OUTPUT_DIR"
echo ""

mkdir -p "$OUTPUT_DIR"

# Time the full run
START=$(date +%s%N)
# Default -l now reads all sequences in one pass (no need to specify)
"$BINARY" -f -u \
    -i "$INPUT_FASTA" \
    -b 50 -r 100 -p 75 \
    -c "$THREADS" \
    -d "$OUTPUT_DIR" 2>&1
END=$(date +%s%N)

ELAPSED_MS=$(( (END - START) / 1000000 ))
echo ""
echo "=== $LABEL Results ==="
echo "Total wall time: ${ELAPSED_MS}ms ($(echo "scale=1; $ELAPSED_MS / 1000" | bc)s)"
echo "Sequences: $NUM_SEQS"
echo "Clusters produced: $(ls "$OUTPUT_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')"
echo ""

# Cleanup
rm -rf "$OUTPUT_DIR"
