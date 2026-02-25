# Tier 1 Optimizations Benchmarking Plan

## Overview

This plan outlines how to benchmark the Tier 1 optimizations implemented in tronko-assign. The goal is to measure both performance improvements and accuracy impacts for each optimization individually and in combination.

## Prerequisites

- Tier 1 optimizations implemented (see `2026-01-01-tier1-optimizations-toggleable.md`)
- Access to test datasets of varying sizes
- `hyperfine` or similar benchmarking tool (optional but recommended)

## Optimization Flags Reference

### Compile-Time Flags
```bash
make NATIVE_ARCH=1    # -march=native -mtune=native
make FAST_MATH=1      # -ffast-math
make ENABLE_OPENMP=1  # -fopenmp
```

### Runtime Flags
```bash
--early-termination      # Enable strike-based early termination
--no-early-termination   # Disable (default)
--strike-box [FLOAT]     # Strike threshold multiplier of Cinterval (default: 1.0)
--max-strikes [INT]      # Strikes before termination (default: 6)
--enable-pruning         # Enable subtree pruning
--disable-pruning        # Disable (default)
--pruning-factor [FLOAT] # Pruning threshold multiplier of Cinterval (default: 2.0)
```

---

## Phase 1: Establish Baseline Measurements

### 1.1 Build Baseline Binary
```bash
cd tronko-assign
make clean && make
cp tronko-assign tronko-assign-baseline
```

### 1.2 Select Test Datasets

Recommended dataset progression:

| Dataset | Size | Purpose |
|---------|------|---------|
| single_tree example | ~164 reads | Quick sanity check |
| 16S_Bacteria subset | ~1000 reads | Medium benchmark |
| 16S_Bacteria full | ~10000+ reads | Performance benchmark |
| Production dataset | Varies | Real-world validation |

### 1.3 Run Baseline Benchmarks

```bash
# Quick test (single_tree)
time ./tronko-assign-baseline -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/baseline_single.txt -w

# Medium test (16S_Bacteria - adjust paths as needed)
time ./tronko-assign-baseline -r \
    -f /path/to/16S_Bacteria/reference_tree.txt \
    -a /path/to/16S_Bacteria/reference.fasta \
    -s -g /path/to/test_reads.fasta \
    -o /tmp/baseline_16s.txt -w
```

### 1.4 Record Baseline Metrics

Create a results file `benchmark_results.tsv`:
```
configuration	dataset	reads	wall_time_s	user_time_s	sys_time_s	max_rss_kb	accuracy_pct
baseline	single_tree	164	X.XX	X.XX	X.XX	XXXXX	100.00
baseline	16s_subset	1000	X.XX	X.XX	X.XX	XXXXX	100.00
```

---

## Phase 2: Compile-Time Optimization Benchmarks

### 2.1 Build Optimized Binaries

```bash
# Native architecture only
make clean && make NATIVE_ARCH=1
cp tronko-assign tronko-assign-native

# Fast math only
make clean && make FAST_MATH=1
cp tronko-assign tronko-assign-fastmath

# OpenMP only
make clean && make ENABLE_OPENMP=1
cp tronko-assign tronko-assign-openmp

# All compile-time optimizations
make clean && make NATIVE_ARCH=1 FAST_MATH=1 ENABLE_OPENMP=1
cp tronko-assign tronko-assign-all-compile
```

### 2.2 Run Compile-Time Benchmarks

```bash
# Test each binary with same dataset
for binary in baseline native fastmath openmp all-compile; do
    echo "Testing tronko-assign-${binary}..."
    /usr/bin/time -v ./tronko-assign-${binary} -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_${binary}.txt -w 2>&1 | tee /tmp/timing_${binary}.txt
done
```

### 2.3 Verify Output Correctness

```bash
# All compile-time optimizations should produce identical output
for binary in native fastmath openmp all-compile; do
    diff /tmp/output_baseline.txt /tmp/output_${binary}.txt && \
        echo "${binary}: IDENTICAL" || echo "${binary}: DIFFERENT"
done
```

---

## Phase 3: Runtime Optimization Benchmarks

### 3.1 Build Single Binary for Runtime Tests

```bash
# Use all compile-time optimizations for runtime tests
make clean && make NATIVE_ARCH=1 FAST_MATH=1
cp tronko-assign tronko-assign-optimized
```

### 3.2 Test Runtime Flag Combinations

```bash
# Define test configurations
declare -A configs
configs["rt_baseline"]=""
configs["rt_early_term"]="--early-termination"
configs["rt_pruning"]="--enable-pruning"
configs["rt_both"]="--early-termination --enable-pruning"
configs["rt_aggressive"]="--early-termination --enable-pruning --strike-box 0.5 --max-strikes 3"
configs["rt_conservative"]="--early-termination --enable-pruning --strike-box 2.0 --max-strikes 10 --pruning-factor 3.0"

# Run each configuration
for config_name in "${!configs[@]}"; do
    flags="${configs[$config_name]}"
    echo "Testing ${config_name}: ${flags}"
    /usr/bin/time -v ./tronko-assign-optimized -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_${config_name}.txt -w ${flags} 2>&1 | tee /tmp/timing_${config_name}.txt
done
```

### 3.3 Calculate Accuracy Metrics

```bash
# Compare each output to baseline and calculate accuracy
baseline="/tmp/output_rt_baseline.txt"

for config_name in "${!configs[@]}"; do
    output="/tmp/output_${config_name}.txt"

    # Count total lines (excluding header)
    total=$(tail -n +2 "$baseline" | wc -l)

    # Count matching lines
    matching=$(diff "$baseline" "$output" | grep -c "^<" || true)
    different=$((total - matching))

    # Calculate accuracy
    accuracy=$(echo "scale=2; ($total - $different) * 100 / $total" | bc)

    echo "${config_name}: ${accuracy}% accuracy (${different}/${total} different)"
done
```

---

## Phase 4: Parameter Sensitivity Analysis

### 4.1 Strike Box Sensitivity

Test how `--strike-box` affects speed vs accuracy:

```bash
for strike_box in 0.25 0.5 1.0 1.5 2.0 3.0 5.0; do
    echo "Testing strike-box=${strike_box}"
    /usr/bin/time -v ./tronko-assign-optimized -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_sb_${strike_box}.txt -w \
        --early-termination --strike-box ${strike_box} 2>&1 | tee /tmp/timing_sb_${strike_box}.txt
done
```

### 4.2 Max Strikes Sensitivity

```bash
for max_strikes in 1 2 3 4 6 8 10 15 20; do
    echo "Testing max-strikes=${max_strikes}"
    /usr/bin/time -v ./tronko-assign-optimized -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_ms_${max_strikes}.txt -w \
        --early-termination --max-strikes ${max_strikes} 2>&1 | tee /tmp/timing_ms_${max_strikes}.txt
done
```

### 4.3 Pruning Factor Sensitivity

```bash
for pruning_factor in 1.0 1.5 2.0 2.5 3.0 4.0 5.0; do
    echo "Testing pruning-factor=${pruning_factor}"
    /usr/bin/time -v ./tronko-assign-optimized -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_pf_${pruning_factor}.txt -w \
        --enable-pruning --pruning-factor ${pruning_factor} 2>&1 | tee /tmp/timing_pf_${pruning_factor}.txt
done
```

---

## Phase 5: Large-Scale Benchmarking

### 5.1 Using hyperfine (Recommended)

```bash
# Install hyperfine if not present
# cargo install hyperfine  OR  apt install hyperfine

# Compare baseline vs optimized with statistical rigor
hyperfine --warmup 2 --runs 10 \
    './tronko-assign-baseline -r -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w' \
    './tronko-assign-optimized -r -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w' \
    './tronko-assign-optimized --early-termination -r -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w' \
    './tronko-assign-optimized --enable-pruning -r -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w' \
    './tronko-assign-optimized --early-termination --enable-pruning -r -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w'
```

### 5.2 Memory Profiling

```bash
# Use valgrind for detailed memory analysis
valgrind --tool=massif ./tronko-assign-baseline -r \
    -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w

ms_print massif.out.* > memory_baseline.txt

valgrind --tool=massif ./tronko-assign-optimized --early-termination --enable-pruning -r \
    -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w

ms_print massif.out.* > memory_optimized.txt
```

### 5.3 CPU Profiling

```bash
# Use perf for CPU analysis
perf record ./tronko-assign-baseline -r \
    -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w
perf report > perf_baseline.txt

perf record ./tronko-assign-optimized --early-termination --enable-pruning -r \
    -f db.txt -a ref.fasta -s -g reads.fasta -o /tmp/out.txt -w
perf report > perf_optimized.txt
```

---

## Phase 6: Results Analysis

### 6.1 Expected Results Template

```markdown
## Benchmark Results: [Dataset Name]

### Environment
- CPU: [model]
- RAM: [size]
- OS: [version]
- GCC: [version]

### Compile-Time Optimizations

| Configuration | Time (s) | Speedup | Memory (MB) | Output Match |
|--------------|----------|---------|-------------|--------------|
| baseline | X.XX | 1.00x | XXX | - |
| NATIVE_ARCH | X.XX | X.XXx | XXX | Yes |
| FAST_MATH | X.XX | X.XXx | XXX | Yes |
| ENABLE_OPENMP | X.XX | X.XXx | XXX | Yes |
| All compile | X.XX | X.XXx | XXX | Yes |

### Runtime Optimizations (with all compile-time)

| Configuration | Time (s) | Speedup | Accuracy (%) |
|--------------|----------|---------|--------------|
| No runtime opts | X.XX | 1.00x | 100.00 |
| --early-termination | X.XX | X.XXx | XX.XX |
| --enable-pruning | X.XX | X.XXx | 100.00 |
| Both | X.XX | X.XXx | XX.XX |

### Parameter Sensitivity

#### Strike Box
| Value | Time (s) | Accuracy (%) |
|-------|----------|--------------|
| 0.5 | X.XX | XX.XX |
| 1.0 | X.XX | XX.XX |
| 2.0 | X.XX | XX.XX |

[Additional tables...]
```

### 6.2 Key Metrics to Report

1. **Wall-clock time**: Primary performance metric
2. **Speedup factor**: Optimized time / baseline time
3. **Accuracy percentage**: Matching assignments / total assignments
4. **Memory usage**: Peak RSS
5. **Nodes visited**: If instrumented (requires code changes)

---

## Phase 7: Recommendations Generation

Based on benchmark results, generate recommendations:

### For Maximum Speed (Accuracy Tradeoff Acceptable)
```bash
make NATIVE_ARCH=1 FAST_MATH=1 ENABLE_OPENMP=1
./tronko-assign --early-termination --enable-pruning \
    --strike-box 0.5 --max-strikes 3 --pruning-factor 1.5 ...
```

### For Balanced Speed/Accuracy
```bash
make NATIVE_ARCH=1 FAST_MATH=1
./tronko-assign --enable-pruning --pruning-factor 2.0 ...
```

### For Maximum Accuracy (Some Speedup)
```bash
make NATIVE_ARCH=1
./tronko-assign ...  # No runtime optimizations
```

---

## Automated Benchmark Script

Save as `scripts/benchmark_tier1.sh`:

```bash
#!/bin/bash
set -e

DATASET_DIR="${1:-../tronko-build/example_datasets/single_tree}"
RESULTS_FILE="${2:-benchmark_results_$(date +%Y%m%d_%H%M%S).tsv}"

echo "Running Tier 1 Optimization Benchmarks"
echo "Dataset: $DATASET_DIR"
echo "Results: $RESULTS_FILE"

# Header
echo -e "config\ttime_real\ttime_user\ttime_sys\tmax_rss_kb\taccuracy_pct" > "$RESULTS_FILE"

# Build binaries
echo "Building binaries..."
make clean && make
cp tronko-assign tronko-assign-baseline

make clean && make NATIVE_ARCH=1 FAST_MATH=1
cp tronko-assign tronko-assign-optimized

# Run baseline
echo "Running baseline..."
/usr/bin/time -f "%e\t%U\t%S\t%M" -o /tmp/time_baseline.txt \
    ./tronko-assign-baseline -r \
    -f "$DATASET_DIR/reference_tree.txt" \
    -a "$DATASET_DIR/Charadriiformes.fasta" \
    -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/output_baseline.txt -w 2>/dev/null

echo -e "baseline\t$(cat /tmp/time_baseline.txt)\t100.00" >> "$RESULTS_FILE"

# Run optimized configurations
for config in "" "--early-termination" "--enable-pruning" "--early-termination --enable-pruning"; do
    config_name=$(echo "$config" | tr -d '-' | tr ' ' '_')
    [ -z "$config_name" ] && config_name="no_runtime_opts"

    echo "Running $config_name..."
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o /tmp/time_${config_name}.txt \
        ./tronko-assign-optimized $config -r \
        -f "$DATASET_DIR/reference_tree.txt" \
        -a "$DATASET_DIR/Charadriiformes.fasta" \
        -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
        -o /tmp/output_${config_name}.txt -w 2>/dev/null

    # Calculate accuracy
    total=$(tail -n +2 /tmp/output_baseline.txt | wc -l)
    diff_count=$(diff /tmp/output_baseline.txt /tmp/output_${config_name}.txt | grep -c "^<" || true)
    accuracy=$(echo "scale=2; ($total - $diff_count) * 100 / $total" | bc)

    echo -e "${config_name}\t$(cat /tmp/time_${config_name}.txt)\t${accuracy}" >> "$RESULTS_FILE"
done

echo "Benchmark complete. Results in $RESULTS_FILE"
cat "$RESULTS_FILE"
```

---

## Success Criteria

- [x] Baseline measurements recorded for all test datasets
- [x] All compile-time flag combinations tested and produce identical output
- [x] Runtime optimization speedups measured
- [x] Accuracy impact quantified for early termination (BROKEN - see notes below)
- [x] Parameter sensitivity documented
- [x] Recommendations generated for different use cases
- [x] Results documented in `thoughts/shared/research/` or experiment log

### Benchmark Results (2026-01-01)

**CRITICAL FINDING**: Early termination is broken and produces 0% accuracy. See `thoughts/shared/research/2026-01-01-tier1-benchmark-results.md` for full details.

Working optimizations:
- Compile-time flags (NATIVE_ARCH, FAST_MATH): ~2% speedup, 100% accuracy
- Subtree pruning: 100% accuracy, minimal speedup on small trees

Broken:
- Early termination: All parameter combinations produce incorrect results

---

## Known Behaviors

### Early Termination Accuracy Impact
The early termination optimization trades accuracy for speed. With default parameters (`--strike-box 1.0 --max-strikes 6`), expect:
- Significant speedup on large trees
- Some accuracy reduction (assignments may be less specific)
- Impact varies by dataset characteristics

### Pruning Accuracy
Subtree pruning with `--pruning-factor 2.0` should maintain 100% accuracy because it only skips subtrees that mathematically cannot contain optimal nodes.

### Compile-Time Optimizations
All compile-time flags should produce byte-identical output to baseline. Any differences indicate a bug.
