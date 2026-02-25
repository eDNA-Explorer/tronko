#!/bin/bash
# test_golden_outputs.sh — Validate optimized tronko-assign produces identical output to baseline
#
# Usage:
#   ./tests/test_golden_outputs.sh [path-to-tronko-assign-binary]
#
# If no binary is specified, uses ./tronko-assign/tronko-assign
# Exit code: 0 = all tests pass, 1 = any test fails

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
GOLDEN_DIR="$SCRIPT_DIR/golden_outputs"

BINARY="${1:-$REPO_DIR/tronko-assign/tronko-assign}"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# Paths to test data
SINGLE_REF="$REPO_DIR/tronko-build/example_datasets/single_tree/reference_tree.txt"
SINGLE_FASTA="$REPO_DIR/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
SINGLE_SE_READS="$REPO_DIR/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta"
SINGLE_PE_R1="$REPO_DIR/example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta"
SINGLE_PE_R2="$REPO_DIR/example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta"
BENCH_10K="$GOLDEN_DIR/benchmark_reads_10000.fasta"

PASSED=0
FAILED=0
TOTAL=0

run_test() {
    local test_name="$1"
    local golden_file="$2"
    shift 2
    local output_file="$TMPDIR/${test_name}.txt"

    TOTAL=$((TOTAL + 1))
    echo -n "  [$TOTAL] $test_name ... "

    # Run tronko-assign with the given arguments
    if ! "$BINARY" "$@" -o "$output_file" -6 -C 1 2>/dev/null; then
        echo "FAIL (tronko-assign exited non-zero)"
        FAILED=$((FAILED + 1))
        return
    fi

    # Compare output
    if diff -q "$golden_file" "$output_file" >/dev/null 2>&1; then
        echo "PASS"
        PASSED=$((PASSED + 1))
    else
        echo "FAIL"
        echo "    Differences found:"
        diff "$golden_file" "$output_file" | head -20 | sed 's/^/    /'
        FAILED=$((FAILED + 1))
    fi
}

echo "=============================================="
echo "Golden Output Regression Tests"
echo "=============================================="
echo "Binary: $BINARY"
echo "Golden: $GOLDEN_DIR"
echo ""

# Verify binary exists
if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary not found or not executable: $BINARY"
    exit 1
fi

# Verify golden outputs exist
if [ ! -d "$GOLDEN_DIR" ]; then
    echo "ERROR: Golden outputs directory not found: $GOLDEN_DIR"
    echo "Run 'make golden' or generate golden outputs first."
    exit 1
fi

echo "--- Quick Tests (164 reads) ---"

run_test "single_tree_singleend_nw" \
    "$GOLDEN_DIR/single_tree_singleend_nw.txt" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -s -g "$SINGLE_SE_READS" -w

run_test "single_tree_singleend_wfa" \
    "$GOLDEN_DIR/single_tree_singleend_wfa.txt" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -s -g "$SINGLE_SE_READS"

run_test "single_tree_pairedend_nw" \
    "$GOLDEN_DIR/single_tree_pairedend_nw.txt" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -p -1 "$SINGLE_PE_R1" -2 "$SINGLE_PE_R2" -w

run_test "single_tree_pairedend_wfa" \
    "$GOLDEN_DIR/single_tree_pairedend_wfa.txt" \
    -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
    -p -1 "$SINGLE_PE_R1" -2 "$SINGLE_PE_R2"

echo ""
echo "--- Benchmark Tests (10K reads) ---"

if [ -f "$BENCH_10K" ]; then
    run_test "benchmark_10k_singleend_nw" \
        "$GOLDEN_DIR/benchmark_10k_singleend_nw.txt" \
        -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
        -s -g "$BENCH_10K" -w

    run_test "benchmark_10k_singleend_wfa" \
        "$GOLDEN_DIR/benchmark_10k_singleend_wfa.txt" \
        -r -f "$SINGLE_REF" -a "$SINGLE_FASTA" \
        -s -g "$BENCH_10K"
else
    echo "  (skipped — benchmark_reads_10000.fasta not found)"
fi

echo ""
echo "=============================================="
echo "Results: $PASSED passed, $FAILED failed, $TOTAL total"
echo "=============================================="

if [ "$FAILED" -gt 0 ]; then
    exit 1
else
    echo "All golden output tests passed!"
    exit 0
fi
