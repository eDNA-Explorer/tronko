#!/bin/bash
# benchmark_optimization.sh — Time tronko-assign with 10K reads to measure optimization impact
#
# Usage:
#   ./tests/benchmark_optimization.sh [path-to-binary] [iterations]
#
# Defaults: ./tronko-assign/tronko-assign, 3 iterations
# Runs NW and WFA aligners, reports wall time and RSS.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
GOLDEN_DIR="$SCRIPT_DIR/golden_outputs"

BINARY="${1:-$REPO_DIR/tronko-assign/tronko-assign}"
ITERATIONS="${2:-3}"

SINGLE_REF="$REPO_DIR/tronko-build/example_datasets/single_tree/reference_tree.txt"
SINGLE_FASTA="$REPO_DIR/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
BENCH_10K="$GOLDEN_DIR/benchmark_reads_10000.fasta"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary not found: $BINARY"
    exit 1
fi

if [ ! -f "$BENCH_10K" ]; then
    echo "ERROR: Benchmark reads not found: $BENCH_10K"
    echo "Generate them first (see tests/golden_outputs/)"
    exit 1
fi

echo "=============================================="
echo "Optimization Benchmark (10K reads)"
echo "=============================================="
echo "Binary:     $BINARY"
echo "Iterations: $ITERATIONS"
echo ""

run_benchmark() {
    local label="$1"
    shift
    local output_file="$TMPDIR/bench_output.txt"

    echo "--- $label ---"
    local times=()
    local rss_values=()

    for i in $(seq 1 "$ITERATIONS"); do
        # /usr/bin/time on macOS outputs to stderr with -l flag
        local timing_file="$TMPDIR/timing_${i}.txt"
        /usr/bin/time -l "$BINARY" "$@" -o "$output_file" -6 -C 1 2>"$timing_file" || true

        local wall=$(grep "real" "$timing_file" | awk '{print $1}')
        local rss=$(grep "maximum resident" "$timing_file" | awk '{print $1}')
        local rss_mb=$(echo "scale=1; $rss / 1048576" | bc)

        times+=("$wall")
        rss_values+=("$rss_mb")
        echo "  Run $i: ${wall}s  RSS: ${rss_mb}MB"
    done

    # Calculate average
    local sum=0
    for t in "${times[@]}"; do
        sum=$(echo "$sum + $t" | bc)
    done
    local avg=$(echo "scale=2; $sum / ${#times[@]}" | bc)

    local rss_sum=0
    for r in "${rss_values[@]}"; do
        rss_sum=$(echo "$rss_sum + $r" | bc)
    done
    local rss_avg=$(echo "scale=1; $rss_sum / ${#rss_values[@]}" | bc)

    echo "  Average: ${avg}s  RSS: ${rss_avg}MB"

    # Verify correctness against golden output
    local golden_file="$GOLDEN_DIR/benchmark_10k_singleend_${label}.txt"
    if [ -f "$golden_file" ]; then
        if diff -q "$golden_file" "$output_file" >/dev/null 2>&1; then
            echo "  Correctness: PASS (matches golden output)"
        else
            echo "  Correctness: FAIL (output differs from golden!)"
            diff "$golden_file" "$output_file" | head -5 | sed 's/^/    /'
        fi
    else
        echo "  Correctness: (no golden file for $label)"
    fi
    echo ""
}

run_benchmark "nw" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -s -g "$BENCH_10K" -w

run_benchmark "wfa" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -s -g "$BENCH_10K"

echo "=============================================="
echo "Done. Compare these numbers after applying optimizations."
echo "=============================================="
