---
date: 2026-01-03T12:54:40-08:00
author: Claude (Sonnet 4.5)
git_commit: e04ddc1b694659c031a78c27f1c4d01074369950
branch: experimental
repository: tronko
status: ready_for_review
tags: [plan, documentation, reproducibility, testing]
---

# Implementation Plan: Reproducibility Documentation and Test Enhancements

**Date**: 2026-01-03T12:54:40-08:00
**Author**: Claude (Sonnet 4.5)
**Status**: Ready for Review
**Related Research**: `thoughts/shared/research/2026-01-03-non-determinism-mitigation-strategies.md`

## Overview

This plan implements the "immediate actions" from the non-determinism mitigation strategy research:
1. Update README.md with reproducibility section
2. Enhance test_determinism.sh with statistical analysis
3. Add usage examples for deterministic mode

## Goals

- **User Education**: Help users understand when determinism matters and how to achieve it
- **Test Improvements**: Provide quantitative metrics beyond simple difference counts
- **Usage Clarity**: Make `-C 1` flag discoverable and document variance expectations

## Success Criteria

- [ ] README includes "Reproducibility and Variance" section before or after "Testing" section
- [ ] Test script outputs mean/median/stddev statistics
- [ ] Usage examples show deterministic mode, standard mode, and variance testing
- [ ] Documentation is clear and accessible to non-experts
- [ ] All changes maintain existing functionality (backward compatible)

---

## Task 1: Add Reproducibility Section to README.md

### Location
File: `/home/jimjeffers/Work/tronko/README.md`
Insert after line 325 (after "# Testing" section)

### Content to Add

```markdown
# Reproducibility and Variance

`tronko-assign` uses BWA's multi-threaded alignment engine, which introduces non-deterministic ordering of candidate matches. This results in approximately **3% variance** in taxonomic assignments between runs with identical inputs when using multiple cores.

## Understanding the Variance

- **Magnitude**: ~3% of reads may receive different taxonomic assignments across runs
- **Nature**: Differences are unbiased (equal probability of more/less specific results)
- **Cause**: BWA's work-stealing thread scheduler finds matches in variable order
- **Scientific Impact**: Variance is typically smaller than biological variation in most studies

### What This Means

The variance affects individual read assignments, not overall community composition. For most metabarcoding studies:
- ✅ Community-level statistics (abundance, diversity) remain stable
- ✅ Biological replicates average out technical variance
- ✅ Differential abundance analysis is unaffected when using proper statistical methods
- ⚠️ Exact read-by-read reproducibility requires special modes (see below)

## Achieving Reproducible Results

### Option 1: Single-Threaded Mode (100% Deterministic)

Use the `-C 1` flag to run in deterministic single-threaded mode:

```bash
tronko-assign -C 1 -r -f reference.trkb -a reference.fasta -p \
  -1 forward.fasta -2 reverse.fasta -o results.txt
```

**Trade-offs**:
- ✅ Perfect reproducibility (identical results every run)
- ❌ 10-16x slower than multi-threaded mode
- 👍 **Use for**: Validation, testing, regression checks, publications requiring exact reproducibility

### Option 2: Standard Multi-Threaded Mode (Fast, ~3% Variance)

Use multiple cores for performance (default or explicit `-C N`):

```bash
tronko-assign -C 16 -r -f reference.trkb -a reference.fasta -p \
  -1 forward.fasta -2 reverse.fasta -o results.txt
```

**Trade-offs**:
- ✅ Fast execution (scales with CPU cores)
- ⚠️ ~3% variance between runs
- 👍 **Use for**: Standard workflows, exploratory analysis, production pipelines

### Option 3: Consensus Voting (Best of Both Worlds)

Run multiple times and aggregate results for both speed and reproducibility:

```bash
# Run 3 times in parallel
tronko-assign -C 16 [...] -o run1.txt &
tronko-assign -C 16 [...] -o run2.txt &
tronko-assign -C 16 [...] -o run3.txt &
wait

# Aggregate results (see consensus script in scripts/)
# (Consensus voting script to be implemented)
```

**Trade-offs**:
- ✅ Reduces variance from ~3% to ~0.3-0.6%
- ✅ Can parallelize across runs (3 runs ≈ 3x cost but parallelizable)
- 👍 **Use for**: Critical datasets, regulatory submissions, important publications

## Testing for Variance

Measure variance on your specific dataset using the test script:

```bash
cd tronko-assign

./scripts/test_determinism.sh \
  -r /path/to/reference.trkb \
  -a /path/to/reference.fasta \
  -1 /path/to/forward.fasta \
  -2 /path/to/reverse.fasta \
  -n 3 \
  -c 16 \
  -k
```

This will run 3 replicates and report variance statistics.

### Expected Results

- **Single-threaded** (`-c 1`): 0% variance (identical results)
- **Multi-threaded** (e.g., `-c 16`): 2-4% variance depending on dataset and core count
- **Higher core counts**: Slightly higher variance due to increased thread contention

## Recommendations by Use Case

| Use Case | Recommended Mode | Cores | Expected Variance |
|----------|------------------|-------|-------------------|
| Validation / Testing | Single-threaded | `-C 1` | 0% |
| Exploratory Analysis | Multi-threaded | `-C 4-16` | ~3% |
| Production Pipeline | Multi-threaded | `-C 8-16` | ~3% |
| Publication (standard) | Multi-threaded | `-C 4-16` | ~3% |
| Publication (exact reproducibility) | Single-threaded | `-C 1` | 0% |
| Regulatory Submission | Consensus voting | `-C 16`, 3 runs | <0.5% |

## Technical Details

For developers and advanced users:

- **Root cause**: BWA's `kt_for()` work-stealing scheduler in `bwa_source_files/kthread.c`
- **Propagation**: Variable match order → different tie-breaking in `placement.c:893`
- **Floating-point sensitivity**: Confidence interval comparisons (`placement.c:930`) sensitive to small score differences
- **Not affected by**: WFA2 alignment (deterministic), scoring algorithm (deterministic), result collection order (deterministic)

For detailed technical analysis, see:
- Research document: `thoughts/shared/research/2026-01-03-non-determinism-mitigation-strategies.md`
- Baseline measurement: `thoughts/shared/research/2026-01-03-non-determinism-baseline-measurement.md`

## FAQ

**Q: Is this a bug?**
A: No, this is inherent to BWA's multi-threaded design. It's a performance/reproducibility trade-off present in many bioinformatics tools.

**Q: Will this affect my biological conclusions?**
A: For most metabarcoding studies, no. The variance is small (~3%) and unbiased, similar to technical replicates.

**Q: Can I use tronko-assign in CI/CD pipelines?**
A: Yes, but use statistical comparison (e.g., expect <5% difference) rather than exact output matching. Or use `-C 1` for exact regression tests.

**Q: Does this affect accuracy?**
A: No, accuracy is identical. The algorithm is sound; only the specific reads assigned to borderline taxa vary slightly.

**Q: Will this be fixed in future versions?**
A: We are exploring options including BWA seed control and alternative aligners. For now, single-threaded mode provides perfect reproducibility when needed.
```

### Implementation Steps

1. Open `/home/jimjeffers/Work/tronko/README.md`
2. Locate line 325 (end of "# Testing" section)
3. Insert the new "# Reproducibility and Variance" section
4. Verify formatting and links

### Testing

- [ ] Verify Markdown renders correctly on GitHub
- [ ] Check all internal links work
- [ ] Ensure examples are executable
- [ ] Confirm table formatting is preserved

---

## Task 2: Enhance test_determinism.sh with Statistical Analysis

### Current Functionality

The existing script (`tronko-assign/scripts/test_determinism.sh`) provides:
- Multiple run execution
- Pairwise comparison of taxonomic paths
- Difference counting and percentage calculation
- Sample difference display

### Enhancements Needed

Add the following statistical metrics:

1. **Mean/Median/StdDev of differences** across all pairwise comparisons
2. **Per-taxonomic-level variance** (phylum, class, order, family, genus, species)
3. **Variance distribution histogram** (text-based)
4. **Confidence intervals** for variance estimate
5. **Run timing statistics** (already collected, need summary)

### Implementation

#### File: `/home/jimjeffers/Work/tronko/tronko-assign/scripts/test_determinism.sh`

**Add after line 216** (after MAX_PERCENT calculation):

```bash
# Calculate statistical metrics
echo ""
echo "=============================================="
echo "           STATISTICAL ANALYSIS              "
echo "=============================================="
echo ""

# Collect all pairwise difference counts
declare -a ALL_DIFFS
for i in $(seq 1 $NUM_RUNS); do
    for j in $(seq $((i+1)) $NUM_RUNS); do
        DIFF_FILE="$OUTPUT_DIR/diff_${i}_${j}.txt"
        DIFF_COUNT=$(wc -l < "$DIFF_FILE")
        ALL_DIFFS+=("$DIFF_COUNT")
    done
done

# Calculate mean (already have AVG_DIFF)
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
STDDEV_PERCENT=$(calc_percent $STDDEV $TOTAL_READS)

echo "Pairwise Variance Statistics:"
echo "  Mean differences:    $AVG_DIFF (${AVG_PERCENT}%)"
echo "  Median differences:  $MEDIAN_DIFF (${MEDIAN_PERCENT}%)"
echo "  Std deviation:       $STDDEV (${STDDEV_PERCENT}%)"
echo "  Min differences:     ${SORTED_DIFFS[0]}"
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
echo "  Showing distribution of difference counts across $PAIR_COUNT comparisons"
echo ""

# Find bins for histogram (10 bins)
NUM_BINS=10
BIN_WIDTH=$(( (MAX_DIFF - ${SORTED_DIFFS[0]} + NUM_BINS - 1) / NUM_BINS ))
if [[ $BIN_WIDTH -eq 0 ]]; then BIN_WIDTH=1; fi

declare -A BINS
for diff in "${ALL_DIFFS[@]}"; do
    BIN_INDEX=$(( (diff - ${SORTED_DIFFS[0]}) / BIN_WIDTH ))
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
    BIN_START=$((${SORTED_DIFFS[0]} + bin * BIN_WIDTH))
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
```

**Add cleanup for new files** (after line 254):

```bash
rm -f "$OUTPUT_DIR"/run_*.txt "$OUTPUT_DIR"/taxa_*.txt \
      "$OUTPUT_DIR/compare.awk" "$OUTPUT_DIR/analyze_levels.awk"
```

### New Output Example

```
==============================================
           STATISTICAL ANALYSIS
==============================================

Pairwise Variance Statistics:
  Mean differences:    42498 (3.21%)
  Median differences:  42450 (3.21%)
  Std deviation:       1200 (0.09%)
  Min differences:     41200
  Max differences:     43800 (3.31%)

Taxonomic Level Variance (Run 1 vs Run 2):
  Level         Count   Percentage
  ----------  -------  -----------
  Phylum           120        0.01%
  Class            450        0.03%
  Order           2100        0.16%
  Family          5800        0.44%
  Genus          12000        0.91%
  Species        18000        1.36%
  Resolution      4028        0.30%

Performance Statistics:
  Total time:   450s
  Average time: 150s
  Min time:     148s
  Max time:     152s
  Time stddev:  4s range

Variance Distribution (Histogram):
  Showing distribution of difference counts across 3 comparisons

  41200-41460: #################################### (1)
  41461-41720: (0)
  41721-41980: (0)
  41981-42240: (0)
  42241-42500: #################### (1)
  42501-42760: (0)
  42761-43020: (0)
  43021-43280: (0)
  43281-43540: (0)
  43541-43800: #################### (1)
```

### Testing

- [ ] Test with 3 runs (existing dataset)
- [ ] Test with 5 runs
- [ ] Test with identical runs (deterministic mode, expect 0%)
- [ ] Verify histogram renders correctly
- [ ] Check AWK script compatibility (bash 3.x and 4.x)

---

## Task 3: Create Usage Examples

### File: New file `tronko-assign/examples/reproducibility_examples.sh`

Create executable example script demonstrating all reproducibility modes:

```bash
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
```

Make executable:
```bash
chmod +x tronko-assign/examples/reproducibility_examples.sh
```

### Create examples directory README

**File**: `tronko-assign/examples/README.md`

```markdown
# tronko-assign Examples

This directory contains example scripts demonstrating various tronko-assign features.

## reproducibility_examples.sh

Demonstrates three approaches to reproducibility:

1. **Single-threaded mode** - 100% deterministic but slow
2. **Multi-threaded mode** - Fast but ~3% variance
3. **Variance testing** - Measure variance on your dataset

### Usage

1. Edit the script to set file paths:
   ```bash
   vim reproducibility_examples.sh
   # Update paths at the top:
   REF_TRKB="path/to/reference_tree.trkb"
   REF_FASTA="path/to/reference.fasta"
   FORWARD_READS="path/to/forward.fasta"
   REVERSE_READS="path/to/reverse.fasta"
   ```

2. Run the script:
   ```bash
   ./reproducibility_examples.sh
   ```

3. The script will run each example interactively, showing commands and waiting for confirmation.

### Requirements

- tronko-assign binary in parent directory
- Reference database (.trkb file from tronko-build)
- Reference sequences (.fasta file)
- Query reads (paired-end FASTA files)

### Expected Output

- `deterministic_results.txt` - Results from single-threaded mode
- `fast_results.txt` - Results from multi-threaded mode
- `/tmp/determinism_test/` - Variance test outputs

See the main README.md "Reproducibility and Variance" section for more details.
```

### Testing

- [ ] Create examples directory: `mkdir -p tronko-assign/examples`
- [ ] Test script with example dataset
- [ ] Verify error messages for unconfigured paths
- [ ] Check interactive prompts work correctly
- [ ] Ensure script is portable (bash 3.x compatible)

---

## Task 4: Update tronko-assign Usage Output (Optional Enhancement)

### File: `tronko-assign/options.c`

Add a note about reproducibility to the `-C` flag help text:

**Line 68** (current):
```c
"-C [INT], number of cores [default:1]\n"
```

**Change to**:
```c
"-C [INT], number of cores [default:1] (use -C 1 for reproducible results)\n"
```

This is a minor change but makes the reproducibility option more discoverable.

---

## Implementation Order

Recommended sequence:

1. **Task 1**: Add reproducibility section to README.md ⏱️ 30 min
   - High visibility, immediate benefit
   - Sets context for other changes

2. **Task 3**: Create usage examples ⏱️ 45 min
   - Provides practical demonstrations
   - Can reference in README

3. **Task 2**: Enhance test script ⏱️ 2-3 hours
   - Most complex task
   - Requires testing and validation

4. **Task 4**: Update options help text ⏱️ 5 min
   - Quick improvement
   - Optional but helpful

**Total Estimated Time**: 4-5 hours

---

## Testing Plan

### 1. Documentation Testing

- [ ] Preview README.md rendering on GitHub
- [ ] Check all code blocks for syntax highlighting
- [ ] Verify all links work
- [ ] Test example commands (if possible)

### 2. Script Testing

Run enhanced test_determinism.sh with:

**Test Case 1: Deterministic mode (expect 0% variance)**
```bash
./scripts/test_determinism.sh \
  -r reference.trkb -a reference.fasta \
  -1 forward.fasta -2 reverse.fasta \
  -n 3 -c 1 -k
```

Expected: All statistics should be 0

**Test Case 2: Multi-threaded mode (expect ~3% variance)**
```bash
./scripts/test_determinism.sh \
  -r reference.trkb -a reference.fasta \
  -1 forward.fasta -2 reverse.fasta \
  -n 3 -c 16 -k
```

Expected: Statistics should show ~2-4% variance

**Test Case 3: Many runs (test statistical stability)**
```bash
./scripts/test_determinism.sh \
  -r reference.trkb -a reference.fasta \
  -1 forward.fasta -2 reverse.fasta \
  -n 5 -c 16 -k
```

Expected: More data points, histogram should have better distribution

### 3. Example Script Testing

```bash
cd tronko-assign/examples
# Edit reproducibility_examples.sh with real paths
vim reproducibility_examples.sh

# Run examples
./reproducibility_examples.sh
```

Expected: All three examples run successfully, files created

### 4. Integration Testing

- [ ] Run full tronko-assign workflow with `-C 1` flag
- [ ] Run full workflow with `-C 16` flag
- [ ] Compare results and verify ~3% difference
- [ ] Check that both workflows complete successfully

---

## Rollback Plan

If any issues arise:

1. **README.md changes**: Easily reverted via git
2. **Test script changes**: Keep backup of original, can revert
3. **Example scripts**: New files, can be deleted
4. **Options.c changes**: Single line, easy to revert

All changes are backward compatible and non-breaking.

---

## Success Metrics

After implementation:

- [ ] Users can easily find information about reproducibility
- [ ] `-C 1` flag usage is clear and documented
- [ ] Variance testing provides actionable statistics
- [ ] Examples demonstrate all three approaches
- [ ] Documentation answers common questions (FAQ)

---

## Future Enhancements (Not in This Plan)

1. **Consensus voting script** - Aggregate multiple runs
2. **BWA seed control** - Modify BWA internals for determinism
3. **CI/CD integration examples** - GitHub Actions workflow
4. **Performance benchmarks** - Document speed vs cores

These can be addressed in future work based on user feedback.

---

## Review Checklist

Before implementation:

- [ ] Plan reviewed by maintainer
- [ ] All file paths verified
- [ ] AWK scripts tested for portability
- [ ] Bash scripts tested on bash 3.x and 4.x
- [ ] Markdown syntax validated
- [ ] Examples tested with real data

---

## Notes

- All changes preserve existing functionality
- Documentation uses accessible language (avoids jargon where possible)
- Examples are copy-paste friendly
- Test enhancements provide quantitative metrics for decision-making
- FAQ addresses common questions proactively
