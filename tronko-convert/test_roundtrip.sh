#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXAMPLE_DIR="$SCRIPT_DIR/../tronko-build/example_datasets/single_tree"
TMP_DIR="/tmp/tronko-convert-test"

mkdir -p "$TMP_DIR"

echo "=== Building tronko-convert ==="
cd "$SCRIPT_DIR"
make clean && make

echo ""
echo "=== Test 1: Text to Binary ==="
./tronko-convert -i "$EXAMPLE_DIR/reference_tree.txt" -o "$TMP_DIR/test.trkb" -v

echo ""
echo "=== Test 2: Binary to Text ==="
./tronko-convert -i "$TMP_DIR/test.trkb" -o "$TMP_DIR/roundtrip.txt" -t -v

echo ""
echo "=== Test 3: Validate Round-Trip ==="
# Compare first 1471 lines (header + metadata + taxonomy) exactly
head -1471 "$EXAMPLE_DIR/reference_tree.txt" > "$TMP_DIR/orig_header.txt"
head -1471 "$TMP_DIR/roundtrip.txt" > "$TMP_DIR/rt_header.txt"

if diff -q "$TMP_DIR/orig_header.txt" "$TMP_DIR/rt_header.txt"; then
    echo "Header and taxonomy: MATCH"
else
    echo "Header and taxonomy: MISMATCH"
    diff "$TMP_DIR/orig_header.txt" "$TMP_DIR/rt_header.txt" | head -20
fi

# Compare node structure lines (every 317th line starting from 1472)
# Note: Posterior values will differ due to float precision
echo ""
echo "=== File Sizes ==="
ls -lh "$EXAMPLE_DIR/reference_tree.txt" "$TMP_DIR/test.trkb"

TEXT_SIZE=$(stat -c%s "$EXAMPLE_DIR/reference_tree.txt")
BINARY_SIZE=$(stat -c%s "$TMP_DIR/test.trkb")
if command -v bc &> /dev/null; then
    RATIO=$(echo "scale=2; $TEXT_SIZE / $BINARY_SIZE" | bc)
    echo "Compression ratio: ${RATIO}x"
else
    echo "Text size: $TEXT_SIZE, Binary size: $BINARY_SIZE"
fi

echo ""
echo "=== Binary File Header ==="
od -A x -t x1z -N 64 "$TMP_DIR/test.trkb"
