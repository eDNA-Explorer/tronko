#!/bin/bash
# tests/verify_no_logic_change.sh
#
# Builds tronko-build from both main and the current optimized branch,
# runs the single-tree example dataset with each, and compares the
# reference_tree.txt output byte-for-byte.
#
# Usage:
#   ./tests/verify_no_logic_change.sh
#
# Requirements:
#   - Both branches must compile (gcc with appropriate version)
#   - The example dataset must exist at tronko-build/example_datasets/single_tree/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
EXAMPLE_DIR="$PROJECT_DIR/tronko-build/example_datasets/single_tree"

# Temp dirs
MAIN_BUILD_DIR=$(mktemp -d)
OPT_BUILD_DIR=$(mktemp -d)
MAIN_OUTPUT_DIR=$(mktemp -d)
OPT_OUTPUT_DIR=$(mktemp -d)
cleanup() {
    rm -rf "$MAIN_BUILD_DIR" "$OPT_BUILD_DIR" "$MAIN_OUTPUT_DIR" "$OPT_OUTPUT_DIR"
}
trap cleanup EXIT

CURRENT_BRANCH=$(git -C "$PROJECT_DIR" rev-parse --abbrev-ref HEAD)
echo "Current branch: $CURRENT_BRANCH"
echo ""

##############################################################################
# Step 1: Build tronko-build from main branch
##############################################################################
echo "=== Building tronko-build from main branch ==="
git -C "$PROJECT_DIR" worktree add "$MAIN_BUILD_DIR" main --quiet 2>/dev/null || {
    echo "ERROR: Could not create worktree for main branch"
    exit 1
}
(
    cd "$MAIN_BUILD_DIR/tronko-build"
    # Main branch uses original Makefile (no -fopenmp, plain gcc)
    make clean 2>/dev/null || true
    make 2>&1 | tail -3
)
MAIN_BINARY="$MAIN_BUILD_DIR/tronko-build/tronko-build"
if [ ! -x "$MAIN_BINARY" ]; then
    echo "ERROR: main branch tronko-build binary not found"
    exit 1
fi
echo "  Built: $MAIN_BINARY"
echo ""

##############################################################################
# Step 2: Build tronko-build from optimized branch
##############################################################################
echo "=== Building tronko-build from $CURRENT_BRANCH ==="
(
    cd "$PROJECT_DIR/tronko-build"
    make clean 2>/dev/null || true
    make 2>&1 | tail -3
)
OPT_BINARY="$PROJECT_DIR/tronko-build/tronko-build"
if [ ! -x "$OPT_BINARY" ]; then
    echo "ERROR: optimized branch tronko-build binary not found"
    exit 1
fi
echo "  Built: $OPT_BINARY"
echo ""

##############################################################################
# Step 3: Run single-tree test with MAIN branch binary
##############################################################################
echo "=== Running single-tree test with MAIN branch binary ==="
echo "  Output dir: $MAIN_OUTPUT_DIR"
time "$MAIN_BINARY" -l \
    -m "$EXAMPLE_DIR/Charadriiformes_MSA.fasta" \
    -x "$EXAMPLE_DIR/Charadriiformes_taxonomy.txt" \
    -t "$EXAMPLE_DIR/RAxML_bestTree.Charadriiformes.reroot" \
    -d "$MAIN_OUTPUT_DIR" 2>&1 | tail -5
echo ""

##############################################################################
# Step 4: Run single-tree test with OPTIMIZED branch binary
##############################################################################
echo "=== Running single-tree test with OPTIMIZED branch binary ==="
echo "  Output dir: $OPT_OUTPUT_DIR"
# Force single-threaded to eliminate any floating-point ordering differences
OMP_NUM_THREADS=1 time "$OPT_BINARY" -l \
    -m "$EXAMPLE_DIR/Charadriiformes_MSA.fasta" \
    -x "$EXAMPLE_DIR/Charadriiformes_taxonomy.txt" \
    -t "$EXAMPLE_DIR/RAxML_bestTree.Charadriiformes.reroot" \
    -d "$OPT_OUTPUT_DIR" 2>&1 | tail -5
echo ""

##############################################################################
# Step 5: Compare reference_tree.txt outputs
##############################################################################
echo "=== Comparing reference_tree.txt ==="
MAIN_REF="$MAIN_OUTPUT_DIR/reference_tree.txt"
OPT_REF="$OPT_OUTPUT_DIR/reference_tree.txt"

if [ ! -f "$MAIN_REF" ]; then
    echo "ERROR: main branch did not produce reference_tree.txt"
    exit 1
fi
if [ ! -f "$OPT_REF" ]; then
    echo "ERROR: optimized branch did not produce reference_tree.txt"
    exit 1
fi

MAIN_LINES=$(wc -l < "$MAIN_REF")
OPT_LINES=$(wc -l < "$OPT_REF")
MAIN_SIZE=$(wc -c < "$MAIN_REF")
OPT_SIZE=$(wc -c < "$OPT_REF")

echo "  Main:      $MAIN_LINES lines, $MAIN_SIZE bytes"
echo "  Optimized: $OPT_LINES lines, $OPT_SIZE bytes"

if diff -q "$MAIN_REF" "$OPT_REF" > /dev/null 2>&1; then
    echo ""
    echo "  *** PASS: reference_tree.txt files are BYTE-IDENTICAL ***"
    echo ""
else
    echo ""
    echo "  *** FAIL: reference_tree.txt files DIFFER ***"
    echo ""
    echo "  First differences:"
    diff --unified=3 "$MAIN_REF" "$OPT_REF" | head -40
    echo ""
    echo "  Checking if differences are only floating-point precision..."
    # Compare structure (non-posterior lines: lines with text names, header)
    MAIN_HEADER=$(head -4 "$MAIN_REF")
    OPT_HEADER=$(head -4 "$OPT_REF")
    if [ "$MAIN_HEADER" = "$OPT_HEADER" ]; then
        echo "  Header (tree count, dimensions): MATCH"
    else
        echo "  Header: MISMATCH"
        echo "    Main: $(head -4 "$MAIN_REF")"
        echo "    Opt:  $(head -4 "$OPT_REF")"
    fi

    # Count lines that differ
    DIFF_COUNT=$(diff "$MAIN_REF" "$OPT_REF" | grep -c '^[<>]' || true)
    echo "  Total differing lines: $DIFF_COUNT (out of $MAIN_LINES)"

    # Check if structural lines (node headers with names) match
    echo ""
    echo "  Comparing node structure (non-numeric lines)..."
    grep -E '[a-zA-Z]' "$MAIN_REF" > /tmp/main_structure.txt || true
    grep -E '[a-zA-Z]' "$OPT_REF" > /tmp/opt_structure.txt || true
    if diff -q /tmp/main_structure.txt /tmp/opt_structure.txt > /dev/null 2>&1; then
        echo "  Node structure: MATCH (differences are only in posterior values)"
    else
        echo "  Node structure: MISMATCH"
        diff /tmp/main_structure.txt /tmp/opt_structure.txt | head -20
    fi

    # Compute max relative difference in posterior values
    echo ""
    echo "  Computing max relative difference in posterior values..."
    python3 -c "
import sys
max_rel = 0.0
max_abs = 0.0
diff_count = 0
total_values = 0
with open('$MAIN_REF') as f1, open('$OPT_REF') as f2:
    for line1, line2 in zip(f1, f2):
        parts1 = line1.strip().split('\t')
        parts2 = line2.strip().split('\t')
        for v1s, v2s in zip(parts1, parts2):
            try:
                v1 = float(v1s)
                v2 = float(v2s)
                total_values += 1
                abs_diff = abs(v1 - v2)
                if abs_diff > 0:
                    diff_count += 1
                    if abs_diff > max_abs:
                        max_abs = abs_diff
                    denom = max(abs(v1), abs(v2), 1e-300)
                    rel = abs_diff / denom
                    if rel > max_rel:
                        max_rel = rel
            except ValueError:
                if v1s != v2s:
                    diff_count += 1
print(f'  Total numeric values compared: {total_values}')
print(f'  Values that differ: {diff_count}')
print(f'  Max absolute difference: {max_abs:.2e}')
print(f'  Max relative difference: {max_rel:.2e}')
if max_rel < 1e-10:
    print('  => Differences are within floating-point precision (ACCEPTABLE)')
else:
    print('  => Differences exceed floating-point precision (INVESTIGATE)')
" 2>&1
fi

##############################################################################
# Step 6: Also test with multiple threads to check for threading issues
##############################################################################
echo ""
echo "=== Running optimized binary with multiple threads ==="
MULTI_OUTPUT_DIR=$(mktemp -d)
OMP_NUM_THREADS=4 time "$OPT_BINARY" -l \
    -m "$EXAMPLE_DIR/Charadriiformes_MSA.fasta" \
    -x "$EXAMPLE_DIR/Charadriiformes_taxonomy.txt" \
    -t "$EXAMPLE_DIR/RAxML_bestTree.Charadriiformes.reroot" \
    -d "$MULTI_OUTPUT_DIR" 2>&1 | tail -5

MULTI_REF="$MULTI_OUTPUT_DIR/reference_tree.txt"
echo ""
echo "  Comparing single-thread vs multi-thread optimized output..."
if diff -q "$OPT_REF" "$MULTI_REF" > /dev/null 2>&1; then
    echo "  *** PASS: Single-thread and multi-thread outputs are BYTE-IDENTICAL ***"
else
    echo "  *** WARN: Single-thread vs multi-thread outputs DIFFER ***"
    DIFF_COUNT=$(diff "$OPT_REF" "$MULTI_REF" | grep -c '^[<>]' || true)
    echo "  Differing lines: $DIFF_COUNT"
    echo "  (This may indicate thread-safety issues in posterior computation)"
fi
rm -rf "$MULTI_OUTPUT_DIR"

##############################################################################
# Cleanup worktree
##############################################################################
git -C "$PROJECT_DIR" worktree remove "$MAIN_BUILD_DIR" --force 2>/dev/null || true

echo ""
echo "=== Verification complete ==="
