#!/bin/bash
# tronko-build regression tests.
#   1. Single-tree mode must reproduce the committed example reference_tree.txt byte for byte.
#   2. Partition mode (-y) must give main's 9 trees for the example, the same database and
#      partition numbers at -J 1 and -J 3, and leave its input untouched. It needs
#      raxmlHPC-PTHREADS, famsa and nw_reroot on PATH and is skipped without them.
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

echo "--- partition (9 trees as on main; same output and partition numbers at -J 1 and -J 3) ---"
if command -v raxmlHPC-PTHREADS >/dev/null && command -v famsa >/dev/null && command -v nw_reroot >/dev/null; then
    out=$(mktemp -d)
    in="$EXAMPLE_DIR/multiple_trees/multiple_MSA"
    before=$(cat "$in"/* | cksum)
    start=$(date +%s)
    "$TRONKO_BUILD" -y -e "$in" -n 5 -d "$out/j1" -s -u 0.1 > "$out/j1.log" 2>&1; rc1=$?
    "$TRONKO_BUILD" -y -e "$in" -n 5 -d "$out/j3" -s -u 0.1 -J 3 > "$out/j3.log" 2>&1; rc3=$?
    echo "  $(( $(date +%s) - start ))s"
    if [ "$rc1" -eq 0 ] && [ "$rc3" -eq 0 ] \
        && [ "$(head -1 "$out/j1/reference_tree.txt")" = 9 ] \
        && cmp -s "$out/j1/reference_tree.txt" "$out/j3/reference_tree.txt" \
        && cmp -s "$out/j1/final_partitions.txt" "$out/j3/final_partitions.txt" \
        && [ "$(ls "$out/j1")" = "$(ls "$out/j3")" ] \
        && [ "$(cat "$in"/* | cksum)" = "$before" ]; then
        echo "  PASS"; PASS=$((PASS + 1))
    else
        echo "  FAIL (exit $rc1 at -J 1, $rc3 at -J 3)"; tail -5 "$out/j1.log" "$out/j3.log"; FAIL=$((FAIL + 1))
    fi
    rm -rf "$out"
else
    echo "  SKIP: raxmlHPC-PTHREADS, famsa or nw_reroot not on PATH"
fi

echo "=== $PASS passed, $FAIL failed ==="
[ "$FAIL" -eq 0 ]
