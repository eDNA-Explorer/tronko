#!/usr/bin/env bash
# test_consistency.sh — verify ancestralclust produces identical output to golden reference
# Run from the tronko-fork root directory: bash ancestralclust/test_consistency.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GOLDEN_DIR="$SCRIPT_DIR/test_data/golden_output"
TEST_OUTPUT="/tmp/ancestralclust_test_$$"
INPUT_FASTA="$ROOT_DIR/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
BINARY="$SCRIPT_DIR/ancestralclust"
PASS=0
FAIL=0

cleanup() { rm -rf "$TEST_OUTPUT"; }
trap cleanup EXIT

# --- Checks ---
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: ancestralclust binary not found. Run 'make' in ancestralclust/ first." >&2
    exit 1
fi
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA not found: $INPUT_FASTA" >&2
    exit 1
fi
if [[ ! -d "$GOLDEN_DIR" ]]; then
    echo "ERROR: Golden output not found: $GOLDEN_DIR" >&2
    exit 1
fi

mkdir -p "$TEST_OUTPUT"

echo "=== AncestralClust Consistency Test ==="
echo "Binary:  $BINARY"
echo "Input:   $INPUT_FASTA ($(grep -c '^>' "$INPUT_FASTA") sequences)"
echo "Golden:  $GOLDEN_DIR"
echo "Output:  $TEST_OUTPUT"
echo ""

# --- Run ancestralclust with identical parameters ---
echo "Running ancestralclust (single-threaded for determinism)..."
"$BINARY" -f -u \
    -i "$INPUT_FASTA" \
    -b 50 -r 100 -p 10 \
    -c 1 \
    -d "$TEST_OUTPUT" 2>&1

echo ""
echo "=== Comparing outputs ==="

# --- Compare .clstr file ---
if diff -q "$GOLDEN_DIR/output.clstr" "$TEST_OUTPUT/output.clstr" > /dev/null 2>&1; then
    echo "PASS: output.clstr matches golden reference"
    ((PASS++))
else
    echo "FAIL: output.clstr differs from golden reference"
    diff "$GOLDEN_DIR/output.clstr" "$TEST_OUTPUT/output.clstr" | head -20
    ((FAIL++))
fi

# --- Compare cluster FASTA files ---
golden_fastas=$(ls "$GOLDEN_DIR"/*.fasta 2>/dev/null | sort)
test_fastas=$(ls "$TEST_OUTPUT"/*.fasta 2>/dev/null | sort)

golden_count=$(echo "$golden_fastas" | wc -l | tr -d ' ')
test_count=$(echo "$test_fastas" | wc -l | tr -d ' ')

if [[ "$golden_count" != "$test_count" ]]; then
    echo "FAIL: different number of cluster files (golden=$golden_count, test=$test_count)"
    ((FAIL++))
else
    echo "PASS: same number of cluster files ($golden_count)"
    ((PASS++))
fi

for golden_fasta in $golden_fastas; do
    base=$(basename "$golden_fasta")
    test_fasta="$TEST_OUTPUT/$base"

    if [[ ! -f "$test_fasta" ]]; then
        echo "FAIL: missing cluster file $base"
        ((FAIL++))
        continue
    fi

    # Compare sequence headers (sorted, order-independent within cluster)
    golden_headers=$(grep "^>" "$golden_fasta" | sort)
    test_headers=$(grep "^>" "$test_fasta" | sort)

    if [[ "$golden_headers" == "$test_headers" ]]; then
        n_seqs=$(grep -c "^>" "$golden_fasta")
        echo "PASS: $base — $n_seqs sequences match"
        ((PASS++))
    else
        echo "FAIL: $base — sequence headers differ"
        diff <(echo "$golden_headers") <(echo "$test_headers") | head -10
        ((FAIL++))
    fi

    # Compare actual sequence content (byte-identical)
    if diff -q "$golden_fasta" "$test_fasta" > /dev/null 2>&1; then
        echo "PASS: $base — byte-identical"
        ((PASS++))
    else
        echo "FAIL: $base — content differs (sequences may be reordered)"
        ((FAIL++))
    fi
done

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    echo "CONSISTENCY CHECK FAILED"
    exit 1
else
    echo "ALL CHECKS PASSED"
    exit 0
fi
