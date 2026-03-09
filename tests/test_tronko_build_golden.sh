#!/bin/bash
# tests/test_tronko_build_golden.sh
# Golden output regression tests for tronko-build optimizations
#
# Usage:
#   Generate golden: ./tests/test_tronko_build_golden.sh --generate
#   Run tests:       ./tests/test_tronko_build_golden.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GOLDEN_DIR="$SCRIPT_DIR/golden_outputs_build"
TRONKO_BUILD="$PROJECT_DIR/tronko-build/tronko-build"
EXAMPLE_DIR="$PROJECT_DIR/tronko-build/example_datasets"

GENERATE=0
if [ "$1" = "--generate" ]; then
    GENERATE=1
    mkdir -p "$GOLDEN_DIR"
fi

PASS=0
FAIL=0

# Exact match test
run_test() {
    local test_name="$1"
    local golden_file="$GOLDEN_DIR/${test_name}_reference_tree.txt"
    local tmpdir
    tmpdir=$(mktemp -d)

    shift
    echo "--- Test: $test_name ---"
    local start_time=$(date +%s)

    "$TRONKO_BUILD" "$@" -d "$tmpdir" 2>&1 | tail -5
    local end_time=$(date +%s)
    echo "  Time: $((end_time - start_time))s"

    local output="$tmpdir/reference_tree.txt"
    if [ ! -f "$output" ]; then
        echo "  ERROR: reference_tree.txt not produced!"
        FAIL=$((FAIL + 1))
        rm -rf "$tmpdir"
        return
    fi

    if [ "$GENERATE" -eq 1 ]; then
        cp "$output" "$golden_file"
        echo "  Golden saved: $golden_file ($(wc -l < "$golden_file") lines)"
    else
        if [ ! -f "$golden_file" ]; then
            echo "  SKIP: no golden file"
            rm -rf "$tmpdir"
            return
        fi
        if diff -q "$golden_file" "$output" > /dev/null 2>&1; then
            echo "  PASS"
            PASS=$((PASS + 1))
        else
            echo "  FAIL: output differs from golden"
            diff "$golden_file" "$output" | head -20
            FAIL=$((FAIL + 1))
        fi
    fi

    rm -rf "$tmpdir"
}

# Smoke test - just checks output produced and has expected structure
run_smoke_test() {
    local test_name="$1"
    local tmpdir
    tmpdir=$(mktemp -d)

    shift
    echo "--- Smoke Test: $test_name ---"
    local start_time=$(date +%s)

    "$TRONKO_BUILD" "$@" -d "$tmpdir" 2>&1 | tail -5
    local end_time=$(date +%s)
    echo "  Time: $((end_time - start_time))s"

    local output="$tmpdir/reference_tree.txt"
    if [ ! -f "$output" ]; then
        echo "  ERROR: reference_tree.txt not produced!"
        FAIL=$((FAIL + 1))
        rm -rf "$tmpdir"
        return
    fi

    local lines=$(wc -l < "$output")
    if [ "$lines" -gt 100 ]; then
        echo "  PASS ($lines lines)"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: output too small ($lines lines)"
        FAIL=$((FAIL + 1))
    fi

    rm -rf "$tmpdir"
}

# Test 1: Single tree mode (deterministic - exact match)
run_test "single_tree" \
    -l \
    -m "$EXAMPLE_DIR/single_tree/Charadriiformes_MSA.fasta" \
    -x "$EXAMPLE_DIR/single_tree/Charadriiformes_taxonomy.txt" \
    -t "$EXAMPLE_DIR/single_tree/RAxML_bestTree.Charadriiformes.reroot"

# Test 2: Partition mode with VeryFastTree (non-deterministic due to srand(time) - smoke test)
run_smoke_test "partition_vft" \
    -y \
    -e "$EXAMPLE_DIR/multiple_trees/multiple_MSA" \
    -n 5 \
    -v -f 50 -a

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
