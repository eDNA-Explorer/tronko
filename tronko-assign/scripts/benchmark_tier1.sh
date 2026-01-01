#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRONKO_DIR="$(dirname "$SCRIPT_DIR")"
cd "$TRONKO_DIR"

DATASET_DIR="${1:-../tronko-build/example_datasets/single_tree}"
RESULTS_FILE="${2:-benchmark_results_$(date +%Y%m%d_%H%M%S).tsv}"
ITERATIONS="${3:-3}"

echo "=============================================="
echo "Tier 1 Optimization Benchmark Suite"
echo "=============================================="
echo "Dataset: $DATASET_DIR"
echo "Results: $RESULTS_FILE"
echo "Iterations: $ITERATIONS"
echo ""

# Detect test reads file
if [ -f "$DATASET_DIR/missingreads_singleend_150bp_2error.fasta" ]; then
    READS_FILE="$DATASET_DIR/missingreads_singleend_150bp_2error.fasta"
elif [ -f "../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta" ]; then
    READS_FILE="../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta"
else
    echo "ERROR: Could not find test reads file"
    exit 1
fi

echo "Reads file: $READS_FILE"

# Detect reference files
if [ -f "$DATASET_DIR/reference_tree.txt" ]; then
    REF_FILE="$DATASET_DIR/reference_tree.txt"
elif [ -f "$DATASET_DIR/reference_tree.txt.gz" ]; then
    REF_FILE="$DATASET_DIR/reference_tree.txt.gz"
elif [ -f "$DATASET_DIR/reference_tree.trkb" ]; then
    REF_FILE="$DATASET_DIR/reference_tree.trkb"
elif [ -f "$DATASET_DIR/reference_tree.trkb.gz" ]; then
    REF_FILE="$DATASET_DIR/reference_tree.trkb.gz"
else
    echo "ERROR: Could not find reference tree file"
    exit 1
fi

# Detect fasta file
if [ -f "$DATASET_DIR/Charadriiformes.fasta" ]; then
    FASTA_FILE="$DATASET_DIR/Charadriiformes.fasta"
elif [ -f "$DATASET_DIR/16S_Bacteria.fasta" ]; then
    FASTA_FILE="$DATASET_DIR/16S_Bacteria.fasta"
else
    echo "ERROR: Could not find reference fasta file"
    exit 1
fi

echo "Reference: $REF_FILE"
echo "Fasta: $FASTA_FILE"
echo ""

# Header
echo -e "config\twall_time_avg\twall_time_min\twall_time_max\tuser_time\tsys_time\taccuracy_pct" > "$RESULTS_FILE"

# Function to run benchmark and collect timing
run_benchmark() {
    local binary="$1"
    local config_name="$2"
    local extra_flags="$3"
    local output_file="/tmp/output_${config_name}.txt"

    echo "Benchmarking: $config_name"
    echo "  Binary: $binary"
    echo "  Flags: $extra_flags"

    local times=()
    for i in $(seq 1 $ITERATIONS); do
        TIMEFORMAT='%R %U %S'
        timing=$( { time "$binary" -r -f "$REF_FILE" -a "$FASTA_FILE" -s -g "$READS_FILE" -o "$output_file" -w -6 $extra_flags >/dev/null 2>&1; } 2>&1 )
        wall=$(echo "$timing" | awk '{print $1}')
        user=$(echo "$timing" | awk '{print $2}')
        sys=$(echo "$timing" | awk '{print $3}')
        times+=("$wall")
        echo "    Run $i: ${wall}s (user: ${user}s, sys: ${sys}s)"
    done

    # Calculate min, max, avg
    min=$(printf '%s\n' "${times[@]}" | sort -n | head -1)
    max=$(printf '%s\n' "${times[@]}" | sort -n | tail -1)
    sum=0
    for t in "${times[@]}"; do
        sum=$(awk "BEGIN {print $sum + $t}")
    done
    avg=$(awk "BEGIN {printf \"%.3f\", $sum / ${#times[@]}}")

    # Calculate accuracy against baseline
    local accuracy="100.00"
    if [ -f "/tmp/output_baseline.txt" ] && [ "$config_name" != "baseline" ]; then
        total=$(tail -n +2 /tmp/output_baseline.txt | wc -l)
        if [ "$total" -gt 0 ]; then
            # Count matching lines (taxonomy assignments)
            diff_count=$(diff /tmp/output_baseline.txt "$output_file" 2>/dev/null | grep -c "^<" || true)
            accuracy=$(awk "BEGIN {printf \"%.2f\", ($total - $diff_count) * 100 / $total}")
        fi
    fi

    echo "  Average: ${avg}s (min: ${min}s, max: ${max}s), Accuracy: ${accuracy}%"
    echo -e "${config_name}\t${avg}\t${min}\t${max}\t${user}\t${sys}\t${accuracy}" >> "$RESULTS_FILE"
    echo ""
}

# ============================================
# Phase 1: Build Binaries
# ============================================
echo "=============================================="
echo "Phase 1: Building Binaries"
echo "=============================================="

echo "Building baseline..."
make clean >/dev/null 2>&1
make >/dev/null 2>&1
cp tronko-assign tronko-assign-baseline
echo "  Done"

echo "Building with NATIVE_ARCH..."
make clean >/dev/null 2>&1
make NATIVE_ARCH=1 >/dev/null 2>&1
cp tronko-assign tronko-assign-native
echo "  Done"

echo "Building with FAST_MATH..."
make clean >/dev/null 2>&1
make FAST_MATH=1 >/dev/null 2>&1
cp tronko-assign tronko-assign-fastmath
echo "  Done"

echo "Building with NATIVE_ARCH + FAST_MATH..."
make clean >/dev/null 2>&1
make NATIVE_ARCH=1 FAST_MATH=1 >/dev/null 2>&1
cp tronko-assign tronko-assign-optimized
echo "  Done"

echo ""

# ============================================
# Phase 2: Compile-Time Benchmarks
# ============================================
echo "=============================================="
echo "Phase 2: Compile-Time Optimization Benchmarks"
echo "=============================================="

run_benchmark "./tronko-assign-baseline" "baseline" ""
run_benchmark "./tronko-assign-native" "native_arch" ""
run_benchmark "./tronko-assign-fastmath" "fast_math" ""
run_benchmark "./tronko-assign-optimized" "native+fastmath" ""

# Verify compile-time outputs match baseline
echo "Verifying output correctness..."
for config in native_arch fast_math native+fastmath; do
    if diff -q /tmp/output_baseline.txt "/tmp/output_${config}.txt" >/dev/null 2>&1; then
        echo "  $config: IDENTICAL to baseline"
    else
        echo "  $config: DIFFERENT from baseline (potential issue!)"
    fi
done
echo ""

# ============================================
# Phase 3: Runtime Optimization Benchmarks
# ============================================
echo "=============================================="
echo "Phase 3: Runtime Optimization Benchmarks"
echo "=============================================="

# Use optimized binary for runtime tests
OPTIMIZED_BINARY="./tronko-assign-optimized"

run_benchmark "$OPTIMIZED_BINARY" "rt_no_opts" ""
run_benchmark "$OPTIMIZED_BINARY" "rt_early_term" "--early-termination"
run_benchmark "$OPTIMIZED_BINARY" "rt_pruning" "--enable-pruning"
run_benchmark "$OPTIMIZED_BINARY" "rt_both" "--early-termination --enable-pruning"
run_benchmark "$OPTIMIZED_BINARY" "rt_aggressive" "--early-termination --enable-pruning --strike-box 0.5 --max-strikes 3"
run_benchmark "$OPTIMIZED_BINARY" "rt_conservative" "--early-termination --enable-pruning --strike-box 2.0 --max-strikes 10 --pruning-factor 3.0"

# ============================================
# Phase 4: Parameter Sensitivity
# ============================================
echo "=============================================="
echo "Phase 4: Parameter Sensitivity Analysis"
echo "=============================================="

echo "Testing strike-box values..."
for sb in 0.5 1.0 1.5 2.0 3.0; do
    run_benchmark "$OPTIMIZED_BINARY" "sb_${sb}" "--early-termination --strike-box $sb"
done

echo "Testing max-strikes values..."
for ms in 3 6 10 15; do
    run_benchmark "$OPTIMIZED_BINARY" "ms_${ms}" "--early-termination --max-strikes $ms"
done

echo "Testing pruning-factor values..."
for pf in 1.5 2.0 2.5 3.0; do
    run_benchmark "$OPTIMIZED_BINARY" "pf_${pf}" "--enable-pruning --pruning-factor $pf"
done

# ============================================
# Summary
# ============================================
echo "=============================================="
echo "Benchmark Complete!"
echo "=============================================="
echo ""
echo "Results saved to: $RESULTS_FILE"
echo ""
echo "Summary:"
cat "$RESULTS_FILE" | column -t -s $'\t'
