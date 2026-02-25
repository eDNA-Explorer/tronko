#!/bin/bash
#
# reproducibility_examples.sh - Examples for achieving reproducibility with tronko-assign
#
# This script demonstrates three approaches to reproducibility:
# 1. Single-threaded mode (100% deterministic)
# 2. Standard multi-threaded mode (~3% variance)
# 3. Variance testing

set -e

# Configuration (UPDATE THESE PATHS)
REF_TRKB="path/to/reference_tree.trkb"
REF_FASTA="path/to/reference.fasta"
FORWARD_READS="path/to/forward.fasta"
REVERSE_READS="path/to/reverse.fasta"

# Check if paths are configured
if [[ "$REF_TRKB" == "path/to/"* ]]; then
    echo "ERROR: Please edit this script and set the file paths at the top"
    echo ""
    echo "You need to set:"
    echo "  REF_TRKB      - Path to tronko-build reference database (.trkb file)"
    echo "  REF_FASTA     - Path to reference sequences (.fasta file)"
    echo "  FORWARD_READS - Path to forward reads"
    echo "  REVERSE_READS - Path to reverse reads"
    exit 1
fi

echo "=========================================="
echo "  tronko-assign Reproducibility Examples"
echo "=========================================="
echo ""

# Example 1: Single-threaded deterministic mode
echo "Example 1: Single-threaded mode (100% reproducible)"
echo "----------------------------------------------------"
echo "This mode is SLOW but produces identical results every time."
echo "Use for: validation, testing, regression checks, publications requiring exact reproducibility"
echo ""
echo "Command:"
echo "  tronko-assign -C 1 -r -f $REF_TRKB -a $REF_FASTA -p \\"
echo "    -1 $FORWARD_READS -2 $REVERSE_READS -o deterministic_results.txt"
echo ""
read -p "Press Enter to run, or Ctrl+C to skip..."

tronko-assign -C 1 -r -f "$REF_TRKB" -a "$REF_FASTA" -p \
  -1 "$FORWARD_READS" -2 "$REVERSE_READS" -o deterministic_results.txt

echo "✓ Results saved to: deterministic_results.txt"
echo ""

# Example 2: Multi-threaded mode
echo "Example 2: Multi-threaded mode (fast, ~3% variance)"
echo "----------------------------------------------------"
echo "This mode is FAST but has ~3% variance between runs."
echo "Use for: standard workflows, exploratory analysis, production pipelines"
echo ""
echo "Command:"
echo "  tronko-assign -C 16 -r -f $REF_TRKB -a $REF_FASTA -p \\"
echo "    -1 $FORWARD_READS -2 $REVERSE_READS -o fast_results.txt"
echo ""
read -p "Press Enter to run, or Ctrl+C to skip..."

tronko-assign -C 16 -r -f "$REF_TRKB" -a "$REF_FASTA" -p \
  -1 "$FORWARD_READS" -2 "$REVERSE_READS" -o fast_results.txt

echo "✓ Results saved to: fast_results.txt"
echo ""

# Example 3: Test variance
echo "Example 3: Test variance on your dataset"
echo "-----------------------------------------"
echo "Run the determinism test script to measure variance."
echo ""
echo "Command:"
echo "  cd tronko-assign"
echo "  ./scripts/test_determinism.sh \\"
echo "    -r $REF_TRKB \\"
echo "    -a $REF_FASTA \\"
echo "    -1 $FORWARD_READS \\"
echo "    -2 $REVERSE_READS \\"
echo "    -n 3 -c 16 -k"
echo ""
read -p "Press Enter to run, or Ctrl+C to skip..."

cd "$(dirname "$0")/.."  # Go to tronko-assign directory
./scripts/test_determinism.sh \
  -r "$REF_TRKB" \
  -a "$REF_FASTA" \
  -1 "$FORWARD_READS" \
  -2 "$REVERSE_READS" \
  -n 3 -c 16 -k

echo ""
echo "=========================================="
echo "  Examples Complete"
echo "=========================================="
echo ""
echo "Files created:"
echo "  deterministic_results.txt - Single-threaded (reproducible)"
echo "  fast_results.txt          - Multi-threaded (fast)"
echo "  /tmp/determinism_test/    - Variance test results"
echo ""
echo "Compare the results:"
echo "  diff <(cut -f1,2 deterministic_results.txt | sort) \\"
echo "       <(cut -f1,2 fast_results.txt | sort) | wc -l"
echo ""
echo "This will show how many reads differ between modes."
