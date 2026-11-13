#!/bin/bash
# tronko-build regression tests.
#   1. Single-tree mode must reproduce the committed example reference_tree.txt byte for byte.
#   2. Partition mode (-y, default flags) must produce a database; it needs raxmlHPC-PTHREADS,
#      famsa and nw_reroot on PATH and is skipped without them.
# Usage: tests/test_tronko_build_golden.sh

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TRONKO_BUILD="$PROJECT_DIR/tronko-build/tronko-build"
EXAMPLE_DIR="$PROJECT_DIR/tronko-build/example_datasets"
PASS=0
FAIL=0

echo "--- single_tree (exact match) ---"
out=$(mktemp -d)
start=$(date +%s)
"$TRONKO_BUILD" -l \
    -m "$EXAMPLE_DIR/single_tree/Charadriiformes_MSA.fasta" \
    -x "$EXAMPLE_DIR/single_tree/Charadriiformes_taxonomy.txt" \
    -t "$EXAMPLE_DIR/single_tree/RAxML_bestTree.Charadriiformes.reroot" \
    -d "$out" > "$out/log" 2>&1
echo "  $(( $(date +%s) - start ))s"
if cmp -s "$EXAMPLE_DIR/single_tree/reference_tree.txt" "$out/reference_tree.txt"; then
    echo "  PASS"; PASS=$((PASS + 1))
else
    echo "  FAIL: differs from example_datasets/single_tree/reference_tree.txt"
    tail -5 "$out/log"; FAIL=$((FAIL + 1))
fi
rm -rf "$out"

echo "--- partition (smoke) ---"
if command -v raxmlHPC-PTHREADS >/dev/null && command -v famsa >/dev/null && command -v nw_reroot >/dev/null; then
    out=$(mktemp -d)
    start=$(date +%s)
    "$TRONKO_BUILD" -y -e "$EXAMPLE_DIR/multiple_trees/multiple_MSA" -n 5 -d "$out" -s -u 0.1 > "$out/log" 2>&1
    echo "  $(( $(date +%s) - start ))s"
    if [ "$(wc -l < "$out/reference_tree.txt" 2>/dev/null || echo 0)" -gt 100 ]; then
        echo "  PASS"; PASS=$((PASS + 1))
    else
        echo "  FAIL"; tail -5 "$out/log"; FAIL=$((FAIL + 1))
    fi
    rm -rf "$out"
else
    echo "  SKIP: raxmlHPC-PTHREADS, famsa or nw_reroot not on PATH"
fi

echo "=== $PASS passed, $FAIL failed ==="
[ "$FAIL" -eq 0 ]
