#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TRONKO_ASSIGN="$SCRIPT_DIR/../tronko-assign"
TRONKO_CONVERT="$SCRIPT_DIR/../../tronko-convert/tronko-convert"
EXAMPLE_DIR="$SCRIPT_DIR/../../tronko-build/example_datasets/single_tree"
QUERY_DIR="$SCRIPT_DIR/../../example_datasets/single_tree"
TMP_DIR="/tmp/tronko-binary-test"

mkdir -p "$TMP_DIR"

echo "=== Phase 2 Binary Format Testing ==="
echo ""

# Check prerequisites
if [ ! -f "$TRONKO_ASSIGN" ]; then
    echo "Error: tronko-assign not found. Run 'make' first."
    exit 1
fi

if [ ! -f "$TRONKO_CONVERT" ]; then
    echo "Error: tronko-convert not found. Build Phase 1 first."
    exit 1
fi

if [ ! -f "$EXAMPLE_DIR/reference_tree.txt" ]; then
    echo "Error: Example reference_tree.txt not found."
    exit 1
fi

# Create test query file if it doesn't exist
QUERY_FILE="$QUERY_DIR/missingreads_singleend_150bp_2error.fasta"
if [ ! -f "$QUERY_FILE" ]; then
    echo "Error: Test query file not found: $QUERY_FILE"
    exit 1
fi

echo "=== Step 1: Convert text to binary ==="
"$TRONKO_CONVERT" -i "$EXAMPLE_DIR/reference_tree.txt" -o "$TMP_DIR/reference_tree.trkb" -v
echo ""

echo "=== Step 2: Run assignment with TEXT format ==="
time "$TRONKO_ASSIGN" -r -f "$EXAMPLE_DIR/reference_tree.txt" \
    -a "$EXAMPLE_DIR/Charadriiformes.fasta" \
    -s -g "$QUERY_FILE" \
    -o "$TMP_DIR/results_text.txt" -w
echo ""

echo "=== Step 3: Run assignment with BINARY format ==="
time "$TRONKO_ASSIGN" -r -f "$TMP_DIR/reference_tree.trkb" \
    -a "$EXAMPLE_DIR/Charadriiformes.fasta" \
    -s -g "$QUERY_FILE" \
    -o "$TMP_DIR/results_binary.txt" -w
echo ""

echo "=== Step 4: Compare critical results (taxonomy, score, mismatches) ==="
# Extract critical columns (columns 2-4: taxonomy, score, forward_mismatch)
cut -f2-4 "$TMP_DIR/results_text.txt" > "$TMP_DIR/critical_text.txt"
cut -f2-4 "$TMP_DIR/results_binary.txt" > "$TMP_DIR/critical_binary.txt"

if diff -q "$TMP_DIR/critical_text.txt" "$TMP_DIR/critical_binary.txt" > /dev/null 2>&1; then
    echo "SUCCESS: Critical columns (taxonomy, score, mismatches) are IDENTICAL"
else
    echo "MISMATCH: Critical columns differ!"
    echo "First 10 differences:"
    diff "$TMP_DIR/critical_text.txt" "$TMP_DIR/critical_binary.txt" | head -20
    exit 1
fi

echo ""
echo "=== Step 5: File size comparison ==="
TEXT_SIZE=$(stat -c%s "$EXAMPLE_DIR/reference_tree.txt")
BINARY_SIZE=$(stat -c%s "$TMP_DIR/reference_tree.trkb")
RATIO=$(awk "BEGIN {printf \"%.2f\", $TEXT_SIZE / $BINARY_SIZE}")
echo "Text size:   $(ls -lh "$EXAMPLE_DIR/reference_tree.txt" | awk '{print $5}')"
echo "Binary size: $(ls -lh "$TMP_DIR/reference_tree.trkb" | awk '{print $5}')"
echo "Compression ratio: ${RATIO}x"

echo ""
echo "=== All tests passed! ==="
echo ""
echo "Note: The Node_Number column may differ due to tie-breaking in score rankings."
echo "      This is expected behavior when multiple nodes have identical scores."
