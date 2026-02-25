#!/bin/bash
#
# test_determinism.sh - Test tronko-assign output for non-determinism
#
# This script runs tronko-assign multiple times with identical input and
# compares the taxonomic path assignments to measure variance.
#
# Usage:
#   ./scripts/test_determinism.sh [OPTIONS]
#
# Options:
#   -n NUM_RUNS     Number of test runs (default: 3)
#   -c CORES        Number of CPU cores (default: 4)
#   -r REF_TRKB     Path to reference .trkb file (required)
#   -a REF_FASTA    Path to reference .fasta file (required)
#   -1 FORWARD      Path to forward reads file (required)
#   -2 REVERSE      Path to reverse reads file (required)
#   -o OUTPUT_DIR   Output directory for results (default: /tmp/determinism_test)
#   -k              Keep intermediate files after completion
#   -h              Show this help message
#
# Example:
#   ./scripts/test_determinism.sh \
#     -r /path/to/reference_tree.trkb \
#     -a /path/to/reference.fasta \
#     -1 /path/to/reads_F.fasta \
#     -2 /path/to/reads_R.fasta \
#     -n 5 -c 4

set -e

# Default values
NUM_RUNS=3
CORES=4
OUTPUT_DIR="/tmp/determinism_test"
KEEP_FILES=false
TRONKO_ASSIGN="$(dirname "$0")/../tronko-assign"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

usage() {
    head -28 "$0" | tail -25 | sed 's/^# //' | sed 's/^#//'
    exit 0
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[OK]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse arguments
while getopts "n:c:r:a:1:2:o:kh" opt; do
    case $opt in
        n) NUM_RUNS="$OPTARG" ;;
        c) CORES="$OPTARG" ;;
        r) REF_TRKB="$OPTARG" ;;
        a) REF_FASTA="$OPTARG" ;;
        1) FORWARD_READS="$OPTARG" ;;
        2) REVERSE_READS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        k) KEEP_FILES=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$REF_TRKB" || -z "$REF_FASTA" || -z "$FORWARD_READS" || -z "$REVERSE_READS" ]]; then
    log_error "Missing required arguments"
    echo ""
    usage
fi

# Validate files exist
for f in "$REF_TRKB" "$REF_FASTA" "$FORWARD_READS" "$REVERSE_READS"; do
    if [[ ! -f "$f" ]]; then
        log_error "File not found: $f"
        exit 1
    fi
done

# Validate tronko-assign exists
if [[ ! -x "$TRONKO_ASSIGN" ]]; then
    log_error "tronko-assign not found or not executable: $TRONKO_ASSIGN"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count reads in input
READ_COUNT=$(grep -c "^>" "$FORWARD_READS" 2>/dev/null || zstdcat "$FORWARD_READS" 2>/dev/null | grep -c "^>" || echo "unknown")

log_info "Starting determinism test"
echo "  Reference TRKB: $REF_TRKB"
echo "  Reference FASTA: $REF_FASTA"
echo "  Forward reads: $FORWARD_READS"
echo "  Reverse reads: $REVERSE_READS"
echo "  Read count: $READ_COUNT"
echo "  Number of runs: $NUM_RUNS"
echo "  CPU cores: $CORES"
echo "  Output directory: $OUTPUT_DIR"
echo ""

# Run tronko-assign multiple times
declare -a RUN_FILES
declare -a RUN_TIMES

for i in $(seq 1 $NUM_RUNS); do
    OUTPUT_FILE="$OUTPUT_DIR/run_${i}.txt"
    RUN_FILES+=("$OUTPUT_FILE")

    log_info "Running test $i of $NUM_RUNS..."

    START_TIME=$(date +%s)

    "$TRONKO_ASSIGN" -r \
        -f "$REF_TRKB" \
        -a "$REF_FASTA" \
        -6 -C "$CORES" -c 10 \
        -p -z -w \
        -1 "$FORWARD_READS" \
        -2 "$REVERSE_READS" \
        -o "$OUTPUT_FILE" 2>/dev/null

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    RUN_TIMES+=("$ELAPSED")

    LINE_COUNT=$(wc -l < "$OUTPUT_FILE")
    log_success "Run $i completed: $LINE_COUNT lines, ${ELAPSED}s"
done

echo ""
log_info "Extracting taxonomic paths..."

# Extract readname + taxonomic_path from each run
declare -a TAXA_FILES
for i in $(seq 1 $NUM_RUNS); do
    TAXA_FILE="$OUTPUT_DIR/taxa_${i}.txt"
    TAXA_FILES+=("$TAXA_FILE")
    cut -f1,2 "${RUN_FILES[$((i-1))]}" | sort > "$TAXA_FILE"
done

echo ""
log_info "Comparing runs..."
echo ""

# Create comparison AWK script
cat > "$OUTPUT_DIR/compare.awk" << 'EOF'
BEGIN { FS="\t" }
NR==FNR { taxa[$1]=$2; next }
{ if ($1 in taxa && taxa[$1] != $2) print $1 "\t" taxa[$1] "\t" $2 }
EOF

# Compare all pairs
echo "=============================================="
echo "         TAXONOMIC PATH VARIANCE REPORT       "
echo "=============================================="
echo ""

TOTAL_READS=$((LINE_COUNT - 1))  # Subtract header
MAX_DIFF=0
TOTAL_DIFF=0
PAIR_COUNT=0

# Helper function to calculate percentage without bc
calc_percent() {
    local num=$1
    local denom=$2
    local whole=$((num * 100 / denom))
    local frac=$(( (num * 10000 / denom) % 100 ))
    printf "%d.%02d" "$whole" "$frac"
}

for i in $(seq 1 $NUM_RUNS); do
    for j in $(seq $((i+1)) $NUM_RUNS); do
        DIFF_FILE="$OUTPUT_DIR/diff_${i}_${j}.txt"
        awk -f "$OUTPUT_DIR/compare.awk" "${TAXA_FILES[$((i-1))]}" "${TAXA_FILES[$((j-1))]}" > "$DIFF_FILE"

        DIFF_COUNT=$(wc -l < "$DIFF_FILE")
        PERCENT=$(calc_percent $DIFF_COUNT $TOTAL_READS)

        if [[ $DIFF_COUNT -gt $MAX_DIFF ]]; then
            MAX_DIFF=$DIFF_COUNT
        fi
        TOTAL_DIFF=$((TOTAL_DIFF + DIFF_COUNT))
        PAIR_COUNT=$((PAIR_COUNT + 1))

        if [[ $DIFF_COUNT -eq 0 ]]; then
            echo -e "Run $i vs Run $j: ${GREEN}IDENTICAL${NC} (0 differences)"
        else
            echo -e "Run $i vs Run $j: ${YELLOW}$DIFF_COUNT differences${NC} (${PERCENT}%)"
        fi
    done
done

AVG_DIFF=$((TOTAL_DIFF / PAIR_COUNT))
AVG_PERCENT=$(calc_percent $AVG_DIFF $TOTAL_READS)
MAX_PERCENT=$(calc_percent $MAX_DIFF $TOTAL_READS)

echo ""
echo "=============================================="
echo "                   SUMMARY                    "
echo "=============================================="
echo ""
echo "Total reads:        $TOTAL_READS"
echo "Number of runs:     $NUM_RUNS"
echo "Comparisons made:   $PAIR_COUNT"
echo ""
echo "Average differences: $AVG_DIFF (${AVG_PERCENT}%)"
echo "Maximum differences: $MAX_DIFF (${MAX_PERCENT}%)"
echo ""

# Statistical analysis
echo ""
echo "=============================================="
echo "           STATISTICAL ANALYSIS              "
echo "=============================================="
echo ""

# Collect all pairwise difference counts
declare -a ALL_DIFFS
MIN_DIFF=$MAX_DIFF
for i in $(seq 1 $NUM_RUNS); do
    for j in $(seq $((i+1)) $NUM_RUNS); do
        DIFF_FILE="$OUTPUT_DIR/diff_${i}_${j}.txt"
        DIFF_COUNT=$(wc -l < "$DIFF_FILE")
        ALL_DIFFS+=("$DIFF_COUNT")
        if [[ $DIFF_COUNT -lt $MIN_DIFF ]]; then
            MIN_DIFF=$DIFF_COUNT
        fi
    done
done

# Calculate median
SORTED_DIFFS=($(printf '%s\n' "${ALL_DIFFS[@]}" | sort -n))
MID_INDEX=$((PAIR_COUNT / 2))
if [[ $((PAIR_COUNT % 2)) -eq 0 ]]; then
    # Even number of elements
    MEDIAN_DIFF=$(( (${SORTED_DIFFS[$((MID_INDEX-1))]} + ${SORTED_DIFFS[$MID_INDEX]}) / 2 ))
else
    # Odd number of elements
    MEDIAN_DIFF=${SORTED_DIFFS[$MID_INDEX]}
fi

# Calculate standard deviation
SUM_SQUARED_DIFF=0
for diff in "${ALL_DIFFS[@]}"; do
    DEVIATION=$((diff - AVG_DIFF))
    SQ_DEV=$((DEVIATION * DEVIATION))
    SUM_SQUARED_DIFF=$((SUM_SQUARED_DIFF + SQ_DEV))
done
VARIANCE=$((SUM_SQUARED_DIFF / PAIR_COUNT))
# Integer square root approximation
STDDEV=$(awk "BEGIN {printf \"%.0f\", sqrt($VARIANCE)}")

MEDIAN_PERCENT=$(calc_percent $MEDIAN_DIFF $TOTAL_READS)
MIN_PERCENT=$(calc_percent $MIN_DIFF $TOTAL_READS)
STDDEV_PERCENT=$(calc_percent $STDDEV $TOTAL_READS)

echo "Pairwise Variance Statistics:"
echo "  Mean differences:    $AVG_DIFF (${AVG_PERCENT}%)"
echo "  Median differences:  $MEDIAN_DIFF (${MEDIAN_PERCENT}%)"
echo "  Std deviation:       $STDDEV (${STDDEV_PERCENT}%)"
echo "  Min differences:     $MIN_DIFF (${MIN_PERCENT}%)"
echo "  Max differences:     $MAX_DIFF (${MAX_PERCENT}%)"
echo ""

# Per-level variance analysis
echo "Analyzing differences by taxonomic level..."

# Create AWK script to analyze taxonomic level differences
cat > "$OUTPUT_DIR/analyze_levels.awk" << 'AWK_SCRIPT'
BEGIN {
    FS="\t"
    levels[0]="Phylum"
    levels[1]="Class"
    levels[2]="Order"
    levels[3]="Family"
    levels[4]="Genus"
    levels[5]="Species"
    levels[6]="Resolution"
    for (i=0; i<=6; i++) count[i]=0
}
{
    # Split both taxonomies
    split($2, tax1, ";")
    split($3, tax2, ";")

    # Find first level of difference
    for (i=1; i<=7; i++) {
        if (tax1[i] != tax2[i]) {
            count[i-1]++
            break
        }
    }

    # Check if only resolution differs (all levels same but length differs)
    if (tax1[1]==tax2[1] && tax1[2]==tax2[2] && tax1[3]==tax2[3] &&
        tax1[4]==tax2[4] && tax1[5]==tax2[5] && tax1[6]==tax2[6] &&
        (length($2) != length($3))) {
        count[6]++
    }
}
END {
    for (i=0; i<=6; i++) {
        print levels[i] "\t" count[i]
    }
}
AWK_SCRIPT

# Run analysis on first diff file (representative)
if [[ -f "$OUTPUT_DIR/diff_1_2.txt" && -s "$OUTPUT_DIR/diff_1_2.txt" ]]; then
    echo ""
    echo "Taxonomic Level Variance (Run 1 vs Run 2):"
    echo "  Level         Count   Percentage"
    echo "  ----------  -------  -----------"

    awk -f "$OUTPUT_DIR/analyze_levels.awk" "$OUTPUT_DIR/diff_1_2.txt" | \
    while IFS=$'\t' read -r level count; do
        if [[ $count -gt 0 ]]; then
            percent=$(calc_percent $count $TOTAL_READS)
            printf "  %-10s  %7d  %11s%%\n" "$level" "$count" "$percent"
        fi
    done
else
    echo "  (All taxonomic paths identical)"
fi

# Run timing statistics
echo ""
echo "Performance Statistics:"
TOTAL_TIME=0
MIN_TIME=999999
MAX_TIME=0
for time in "${RUN_TIMES[@]}"; do
    TOTAL_TIME=$((TOTAL_TIME + time))
    if [[ $time -lt $MIN_TIME ]]; then MIN_TIME=$time; fi
    if [[ $time -gt $MAX_TIME ]]; then MAX_TIME=$time; fi
done
AVG_TIME=$((TOTAL_TIME / NUM_RUNS))
echo "  Total time:   ${TOTAL_TIME}s"
echo "  Average time: ${AVG_TIME}s"
echo "  Min time:     ${MIN_TIME}s"
echo "  Max time:     ${MAX_TIME}s"
echo "  Time stddev:  $((MAX_TIME - MIN_TIME))s range"

# Text-based histogram of differences
echo ""
echo "Variance Distribution (Histogram):"
if [[ ${#ALL_DIFFS[@]} -eq 0 || $MAX_DIFF -eq 0 ]]; then
    echo "  (No variance detected)"
else
    echo "  Showing distribution of difference counts across $PAIR_COUNT comparisons"
    echo ""

    # Find bins for histogram (10 bins)
    NUM_BINS=10
    BIN_WIDTH=$(( (MAX_DIFF - MIN_DIFF + NUM_BINS - 1) / NUM_BINS ))
    if [[ $BIN_WIDTH -eq 0 ]]; then BIN_WIDTH=1; fi

    declare -A BINS
    for diff in "${ALL_DIFFS[@]}"; do
        BIN_INDEX=$(( (diff - MIN_DIFF) / BIN_WIDTH ))
        BINS[$BIN_INDEX]=$((${BINS[$BIN_INDEX]:-0} + 1))
    done

    MAX_BAR_COUNT=0
    for bin_count in "${BINS[@]}"; do
        if [[ $bin_count -gt $MAX_BAR_COUNT ]]; then
            MAX_BAR_COUNT=$bin_count
        fi
    done

    # Scale bars to 40 characters max
    for bin in $(seq 0 $((NUM_BINS-1))); do
        BIN_START=$((MIN_DIFF + bin * BIN_WIDTH))
        BIN_END=$((BIN_START + BIN_WIDTH - 1))
        BIN_COUNT=${BINS[$bin]:-0}

        if [[ $MAX_BAR_COUNT -gt 0 ]]; then
            BAR_LENGTH=$((BIN_COUNT * 40 / MAX_BAR_COUNT))
        else
            BAR_LENGTH=0
        fi

        BAR=$(printf '%.0s#' $(seq 1 $BAR_LENGTH))
        printf "  %5d-%5d: %s (%d)\n" "$BIN_START" "$BIN_END" "$BAR" "$BIN_COUNT"
    done
fi

echo ""

if [[ $MAX_DIFF -eq 0 ]]; then
    echo -e "${GREEN}RESULT: All runs produced identical taxonomic assignments!${NC}"
    EXIT_CODE=0
else
    echo -e "${YELLOW}RESULT: Non-deterministic behavior detected${NC}"
    echo ""
    echo "Sample differences (first 10 from Run 1 vs Run 2):"
    echo "------------------------------------------------"
    head -10 "$OUTPUT_DIR/diff_1_2.txt" | while IFS=$'\t' read -r readname taxa1 taxa2; do
        echo "  $readname"
        echo "    Run 1: $taxa1"
        echo "    Run 2: $taxa2"
        echo ""
    done
    EXIT_CODE=1
fi

echo ""
echo "Output files saved to: $OUTPUT_DIR"

# Cleanup if not keeping files
if [[ "$KEEP_FILES" = false ]]; then
    log_info "Cleaning up intermediate files (use -k to keep)"
    rm -f "$OUTPUT_DIR"/run_*.txt "$OUTPUT_DIR"/taxa_*.txt "$OUTPUT_DIR/compare.awk" "$OUTPUT_DIR/analyze_levels.awk"
else
    log_info "Keeping all intermediate files"
    echo "  Full results: $OUTPUT_DIR/run_*.txt"
    echo "  Taxa files:   $OUTPUT_DIR/taxa_*.txt"
    echo "  Diff files:   $OUTPUT_DIR/diff_*.txt"
fi

exit $EXIT_CODE
