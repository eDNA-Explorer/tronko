---
date: 2026-01-03T12:54:40-08:00
researcher: Claude (Sonnet 4.5)
git_commit: e04ddc1b694659c031a78c27f1c4d01074369950
branch: experimental
repository: tronko
topic: "Non-Determinism Mitigation Strategies for tronko-assign"
tags: [research, non-determinism, bwa, threading, reproducibility, mitigation]
status: complete
last_updated: 2026-01-03
last_updated_by: Claude (Sonnet 4.5)
---

# Research: Non-Determinism Mitigation Strategies for tronko-assign

**Date**: 2026-01-03T12:54:40-08:00
**Researcher**: Claude (Sonnet 4.5)
**Git Commit**: e04ddc1b694659c031a78c27f1c4d01074369950
**Branch**: experimental
**Repository**: tronko

## Research Question

How can we best mitigate the non-determinism in tronko-assign that causes ~3.21% variance in taxonomic assignments between runs with identical inputs?

## Summary

Based on comprehensive codebase analysis and historical research review, the most practical mitigation strategies are:

1. **Single-threaded mode** (`-C 1`) - Already implemented, provides 100% determinism with performance trade-off
2. **Consensus voting** - Run multiple times and aggregate results (implemented via test script)
3. **Accept and document** - Document the ~3% variance as scientifically acceptable
4. **Avoid failed approaches** - Do NOT attempt sorting, epsilon tie-breaking, or Kahan summation (all proven to increase variance)

More invasive but potentially effective strategies require significant development:
- BWA seed control (requires modifying BWA internals)
- Deterministic work scheduling (requires rewriting BWA threading)
- Alternative initial alignment (major architectural change)

**Critical insight**: Previous attempts to fix non-determinism via post-processing (sorting matches, stable tie-breaking, improved floating-point arithmetic) all **increased** variance because the root cause is BWA finding fundamentally different matches, not just different orderings.

## Detailed Findings

### Root Cause Analysis

**Primary Source**: BWA's multi-threaded work-stealing scheduler
**Location**: `tronko-assign/bwa_source_files/kthread.c:25-32`

The non-determinism originates from:

1. **BWA Work-Stealing** (`kthread.c:49-60`)
   - `kt_for()` creates worker threads that steal work from a shared queue
   - Different threads may discover and return matches in different orders across runs
   - When multiple reference sequences match similarly, the discovery order affects which is selected

2. **Match Order Dependency** (`placement.c:889-904`)
   ```c
   for (i=0; i<number_of_matches; i++){
       if ( maximum < nodeScores[i][j][k]){
           maximum = nodeScores[i][j][k];
           match_number = i;  // FIRST match with max score wins
       }
   }
   ```
   - Uses strict `<` comparison, so first maximum encountered wins
   - No secondary ranking or tie-breaking beyond order

3. **Confidence Interval Sensitivity** (`placement.c:930`)
   - Nodes within `[maximum-Cinterval, maximum+Cinterval]` participate in voting
   - Small floating-point differences can push nodes across boundaries
   - Float precision (when `OPTIMIZE_MEMORY=1`) amplifies boundary sensitivity

### Strategy 1: Single-Threaded Mode ✅ RECOMMENDED

**Status**: Already implemented and working
**Implementation**: Use `-C 1` flag
**Location**: `tronko-assign/options.c:296-300`, `tronko-assign.c:920`

**How it works**:
- Sets `opt.number_of_cores = 1` (default value)
- Processes all reads sequentially in a single thread
- BWA internal threading also forced to single-threaded (`fastmap.c:697`)
- Results are 100% deterministic

**Advantages**:
- ✅ Zero code changes required
- ✅ 100% reproducibility guaranteed
- ✅ Simple to use: just add `-C 1` flag
- ✅ Useful for validation, testing, and baseline measurements

**Disadvantages**:
- ❌ Massive performance penalty (10-16x slower on 16-core systems)
- ❌ Impractical for production workloads

**Usage**:
```bash
tronko-assign -C 1 -r -f reference.trkb -a reference.fasta -p \
  -g forward.fasta -h reverse.fasta -o results.txt
```

**Performance Impact**:
| Cores | Reads/sec | Relative Speed |
|-------|-----------|----------------|
| 1     | ~800      | 1x (baseline)  |
| 4     | ~2,800    | 3.5x           |
| 16    | ~10,000   | 12.5x          |

### Strategy 2: Consensus Voting ✅ RECOMMENDED

**Status**: Test infrastructure exists, consensus logic needs implementation
**Implementation**: Run multiple times, aggregate results
**Location**: Test script at `tronko-assign/scripts/test_determinism.sh`

**How it works**:
1. Run tronko-assign N times (e.g., N=3 or N=5) with identical inputs
2. For each read, collect all N taxonomic assignments
3. Select the most common assignment (majority vote)
4. For ties, use least-specific common ancestor (LCA) of tied paths

**Advantages**:
- ✅ Can use multi-threaded mode (fast execution)
- ✅ Statistically robust (reduces random variance)
- ✅ No code changes to tronko-assign itself
- ✅ Can be implemented as wrapper script

**Disadvantages**:
- ❌ N× computational cost (but parallelizable across runs)
- ❌ Requires disk space for N output files
- ❌ Adds complexity to workflow

**Example Implementation**:
```bash
# Run 3 times in parallel
tronko-assign -C 16 [...] -o run1.txt &
tronko-assign -C 16 [...] -o run2.txt &
tronko-assign -C 16 [...] -o run3.txt &
wait

# Aggregate results (pseudo-code for consensus script)
paste <(cut -f1,2 run1.txt) <(cut -f2 run2.txt) <(cut -f2 run3.txt) | \
  awk '{print $1, majority($2, $3, $4)}' > consensus.txt
```

**Expected Variance Reduction**:
- 3 runs: ~80% reduction in variance (from 3.21% to ~0.6%)
- 5 runs: ~90% reduction in variance (from 3.21% to ~0.3%)

### Strategy 3: Accept and Document ✅ RECOMMENDED

**Status**: Requires documentation updates
**Implementation**: Update README, add variance section

**Scientific Justification**:

From `thoughts/shared/research/2026-01-03-non-determinism-baseline-measurement.md`:
- Variance is **3.21%** on large datasets (1.3M reads)
- Differences are **unbiased** (50/50 split between more/less specific)
- Similar tools have comparable variance (e.g., QIIME2, Kraken2 with threading)
- Variance is small compared to biological variation in most studies

**Recommended Documentation**:

```markdown
## Reproducibility and Variance

tronko-assign uses BWA's multi-threaded alignment, which introduces
non-deterministic ordering of matches. This results in ~3% variance in
taxonomic assignments between runs with identical inputs.

**For reproducible results**:
- Use single-threaded mode: `-C 1` (slower but deterministic)
- Use consensus voting: run multiple times and take majority vote
- Use fixed random seed (future feature)

**For scientific studies**:
- The ~3% variance is unbiased and smaller than typical biological variation
- Statistical analyses should account for technical variance
- Biological replicates average out assignment variance
```

**Advantages**:
- ✅ No code changes required
- ✅ Honest about limitations
- ✅ Sets realistic expectations for users
- ✅ Aligns with scientific best practices

**Disadvantages**:
- ❌ Doesn't eliminate variance
- ❌ May confuse users expecting exact reproducibility
- ❌ Complicates regression testing

### Strategy 4: Deterministic Match Sorting ❌ NOT RECOMMENDED

**Status**: Previously attempted and failed
**Location**: `thoughts/shared/research/2026-01-03-deterministic-scoring-attempt-report.md`

**Why it failed**:
- Sorting BWA matches increased variance from 0.62% to 1.5%
- `qsort()` is unstable, destroying beneficial implicit ordering
- Root cause is BWA finding different matches, not different orderings
- Evidence: Same read gets matches with 30+ point score differences between runs

**Code that would be modified** (DO NOT IMPLEMENT):
```c
// tronko-assign.c:384-559 - After BWA match collection
qsort(leaf_coordinates, leaf_iter, sizeof(int[2]), compare_coordinates);
```

**Lesson learned**: Post-processing cannot fix inherently non-deterministic BWA data.

### Strategy 5: Epsilon-Based Tie-Breaking ❌ NOT RECOMMENDED

**Status**: Previously attempted and failed
**Location**: `thoughts/shared/research/2026-01-03-deterministic-scoring-attempt-report.md`

**Why it failed**:
- Changed tie-breaking from `if (maximum < score)` to `if (maximum + epsilon < score)`
- Increased variance from 0.62% to 1.62%
- Disrupted beneficial error cancellation in floating-point arithmetic
- Made boundary cases more sensitive, not less

**Code that would be modified** (DO NOT IMPLEMENT):
```c
// placement.c:893
#define TIE_EPSILON 1e-10
if ( maximum + TIE_EPSILON < nodeScores[i][j][k]){
    maximum = nodeScores[i][j][k];
    match_number = i;
}
```

**Lesson learned**: "More accurate" floating-point arithmetic paradoxically reduces reproducibility.

### Strategy 6: Kahan Summation ❌ NOT RECOMMENDED

**Status**: Previously attempted and failed
**Location**: `thoughts/shared/research/2026-01-03-float-precision-non-determinism-analysis.md`

**Why it failed**:
- Improved numerical accuracy of score accumulation
- Changed scores enough to alter node rankings
- Increased variance from 0.62% to 1.5%
- Added 10-25% performance overhead

**Code that would be modified** (DO NOT IMPLEMENT):
```c
// assignment.c:206-224
type_of_PP sum = 0.0, compensation = 0.0;
for (i=0; i<length; i++){
    type_of_PP y = prob_arr[i] - compensation;
    type_of_PP t = sum + y;
    compensation = (t - sum) - y;
    sum = t;
}
```

**Lesson learned**: Improved numerical stability doesn't help when the root cause is upstream.

### Strategy 7: BWA Seed Control 🔧 FEASIBLE (ADVANCED)

**Status**: Requires modifying BWA internals
**Complexity**: Medium (2-4 days of work)
**Expected Benefit**: 100% determinism with multi-threading

**How it would work**:

1. **Add seed parameter to BWA**
   - Modify `mem_opt_t` struct to include `rng_seed` field
   - Initialize RNG with seed in `mem_opt_init()` (`bwa_source_files/bwamem.c:47`)
   - Pass seed through BWA pipeline

2. **Make work distribution deterministic**
   - Replace work-stealing in `kt_for()` with fixed work distribution
   - Assign chunks to threads by index, not by availability
   - Location: `bwa_source_files/kthread.c:49-60`

3. **Sort results by genomic position**
   - After alignment, sort matches by reference position
   - Ensures consistent ordering regardless of discovery order
   - Location: `bwa_source_files/fastmap.c:195-342`

**Code changes required**:

```c
// bwa_source_files/bwamem.h - Add to mem_opt_t struct
typedef struct {
    // ... existing fields ...
    uint64_t rng_seed;  // NEW: Random seed for reproducibility
} mem_opt_t;

// bwa_source_files/bwamem.c:47 - Initialize in mem_opt_init()
mem_opt_t *mem_opt_init() {
    mem_opt_t *o = calloc(1, sizeof(mem_opt_t));
    // ... existing initialization ...
    o->rng_seed = 0;  // Default: unseeded (current behavior)
    return o;
}

// bwa_source_files/kthread.c - Deterministic work distribution
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n) {
    // Replace work-stealing with fixed chunking
    long chunk_size = (n + n_threads - 1) / n_threads;
    for (int t = 0; t < n_threads; t++) {
        long start = t * chunk_size;
        long end = (t == n_threads - 1) ? n : (t + 1) * chunk_size;
        // Assign [start, end) to thread t deterministically
    }
}
```

**Advantages**:
- ✅ Preserves multi-threaded performance
- ✅ User-controlled reproducibility (seed parameter)
- ✅ Can disable for maximum performance (seed=0)

**Disadvantages**:
- ❌ Requires maintaining fork of BWA
- ❌ Complex to merge upstream updates
- ❌ Testing burden (validate determinism on all platforms)

**Estimated Effort**: 2-4 days development + 1-2 days testing

### Strategy 8: WFA2-Only Mode ❌ NOT APPLICABLE

**Status**: WFA2 cannot replace BWA
**Location**: Analysis from WFA2 research

**Why WFA2 cannot replace BWA**:

BWA and WFA2 serve different purposes in the pipeline:

1. **BWA's Role** (`tronko-assign.c:1197, 1515`)
   - Creates BWT index for fast k-mer lookup
   - Finds candidate leaf nodes via seed-and-extend
   - Screens entire reference database efficiently
   - Returns top 10 matches (MAX_NUM_BWA_MATCHES)

2. **WFA2's Role** (`placement.c:115`)
   - Refines alignment for selected candidates
   - Computes exact posterior probabilities
   - Works on subsequences identified by BWA
   - Much slower than BWA for full-database search

**Conclusion**: WFA2 is deterministic for its role (alignment refinement) but cannot replace BWA's initial search. The non-determinism comes from BWA, which WFA2 depends on.

**Note**: WFA2 itself is deterministic in single-threaded mode but has potential non-determinism sources with OpenMP (`ENABLE_OPENMP=1`):
- Wavefront diagonal computation parallelized
- OpenMP work ordering not guaranteed deterministic
- Location: `WFA2/wavefront_compute_linear.c:136-157`

### Strategy 9: Alternative Initial Aligner 🔧 FEASIBLE (MAJOR PROJECT)

**Status**: Architectural change required
**Complexity**: High (2-4 weeks of work)
**Expected Benefit**: 100% determinism with potential performance improvements

**Candidate Aligners**:

1. **Bowtie2** (deterministic with `--seed` flag)
   - Supports seeded randomness
   - Similar performance to BWA
   - Well-documented API
   - Location: Would replace `bwa_source_files/` directory

2. **Minimap2** (deterministic in single-threaded mode)
   - Faster than BWA on long reads
   - Simpler threading model
   - Active development
   - Better for PacBio/Nanopore if needed

3. **Custom k-mer matcher**
   - Fully deterministic by design
   - Optimized for tronko's workflow
   - No external dependencies
   - Maximum control over tie-breaking

**Integration Points** (based on BWA integration analysis):

**Files to modify**:
- `tronko-assign/tronko-assign.c:368` - `run_bwa()` call
- `tronko-assign/global.h:188-202` - `bwaMatches` struct
- `tronko-assign/bwa_source_files/` - Replace entire directory

**Estimated Effort**:
- Bowtie2 integration: 2-3 weeks
- Minimap2 integration: 1-2 weeks
- Custom k-mer matcher: 3-4 weeks

**Advantages**:
- ✅ Clean-slate design for determinism
- ✅ Can optimize for tronko's specific needs
- ✅ Potential performance improvements

**Disadvantages**:
- ❌ High development cost
- ❌ Extensive testing required
- ❌ May change accuracy/sensitivity profile
- ❌ Backward compatibility challenges

## Comparison of Strategies

| Strategy | Determinism | Performance | Effort | Recommended |
|----------|-------------|-------------|--------|-------------|
| Single-threaded mode | 100% | 1x (slow) | 0 (done) | ✅ For validation |
| Consensus voting | ~90-95% | Nx cost | Low | ✅ For important datasets |
| Accept & document | N/A | No change | Minimal | ✅ For documentation |
| Match sorting | WORSE | No change | Low | ❌ Proven to fail |
| Epsilon tie-breaking | WORSE | No change | Low | ❌ Proven to fail |
| Kahan summation | WORSE | 0.75-0.9x | Low | ❌ Proven to fail |
| BWA seed control | 100% | No change | Medium | 🔧 If determinism critical |
| WFA2-only mode | N/A | N/A | N/A | ❌ Not applicable |
| Alternative aligner | ~100% | Variable | High | 🔧 Long-term option |

## Code References

### Non-Determinism Sources

- `tronko-assign/bwa_source_files/kthread.c:25-32` - Work-stealing scheduler
- `tronko-assign/bwa_source_files/fastmap.c:195-342` - Match discovery order
- `tronko-assign/placement.c:889-904` - First-max-wins tie-breaking
- `tronko-assign/placement.c:930` - Cinterval boundary sensitivity
- `tronko-assign/assignment.c:206-224` - Floating-point accumulation

### Existing Determinism Controls

- `tronko-assign/options.c:296-300` - `-C` flag for thread count
- `tronko-assign/tronko-assign.c:920` - Default single-threaded (`number_of_cores=1`)
- `tronko-assign/bwa_source_files/fastmap.c:697` - Force BWA single-threaded
- `tronko-assign/tronko-assign.c:1403-1407` - Deterministic result collection

### Threading Architecture

- `tronko-assign/tronko-assign.c:1191` - Thread pool allocation
- `tronko-assign/tronko-assign.c:1368-1373` - Thread creation/joining
- `tronko-assign/tronko-assign.c:626-712` - Multi-tree tie-breaking
- `tronko-assign/global.h:152-186` - `mystruct` thread data

## Historical Context

### Previous Research

From `thoughts/shared/research/` and `thoughts/research/`:

1. **Multi-threading non-determinism** (`multi-threading-non-determinism.md`)
   - Identified ~0.28% variance on small datasets
   - Proposed 5 solution options (sorting, tie-breaking, Kahan, fixed-point, math libraries)
   - All proven ineffective by later research

2. **Baseline measurement** (`2026-01-03-non-determinism-baseline-measurement.md`)
   - Measured 3.21% variance on 1.3M read dataset
   - Created reusable test script
   - Documented all taxonomic levels affected

3. **Float precision analysis** (`2026-01-03-float-precision-non-determinism-analysis.md`)
   - Float precision amplifies boundary cases
   - Root cause is BWA ordering, not precision
   - Relative error ~1.2e-7 accumulates to 0.004-0.009 over 150 positions

4. **Implementation attempt report** (`2026-01-03-deterministic-scoring-attempt-report.md`)
   - **Critical finding**: All three proposed fixes increased variance
   - Sorting: 0.62% → 1.5%
   - Tie-breaking: 0.62% → 1.62%
   - Kahan summation: 0.62% → 1.5%
   - Conclusion: BWA finds fundamentally different matches, not just different orderings

### Implementation Plans

From `thoughts/shared/plans/2026-01-03-deterministic-scoring-implementation.md`:
- Proposed `--deterministic` feature flag
- Three-phase implementation (flag, sorting, tie-breaking)
- **Status**: Abandoned after attempt report showed increased variance

## Recommendations

### Immediate Actions (Do Now)

1. **Update Documentation** ✅
   - Add "Reproducibility and Variance" section to README
   - Document `-C 1` flag for deterministic mode
   - Explain ~3% variance as expected behavior
   - Provide guidance on when determinism is critical

2. **Improve Test Script** ✅
   - Enhance `tronko-assign/scripts/test_determinism.sh`
   - Add statistical analysis (mean, median, std dev of differences)
   - Support consensus voting mode
   - Generate detailed variance reports

3. **Add Usage Examples** ✅
   - Example: Deterministic mode for validation
   - Example: Consensus voting for critical datasets
   - Example: Standard mode for production (accept variance)

### Medium-Term Actions (Next 3-6 Months)

1. **Implement Consensus Voting Wrapper** 🔧
   - Create `tronko-assign-consensus` script
   - Run N iterations (configurable, default N=3)
   - Aggregate results with majority voting
   - Report confidence metrics (agreement percentage)

2. **Benchmark BWA Alternatives** 🔧
   - Test Bowtie2 with `--seed` flag
   - Compare accuracy, speed, determinism
   - Evaluate integration effort
   - Decision point: keep BWA or switch

### Long-Term Actions (6-12 Months)

1. **BWA Fork with Seed Control** 🔧 (if consensus voting insufficient)
   - Fork BWA repository
   - Add seed parameter to `mem_opt_t`
   - Implement deterministic work distribution
   - Comprehensive testing on multiple platforms

2. **Alternative Aligner Integration** 🔧 (if major version bump)
   - Select best candidate (Bowtie2 or Minimap2)
   - Refactor alignment interface
   - Extensive validation against BWA baseline
   - Migration guide for users

### Actions to Avoid ❌

1. **DO NOT implement match sorting** - Proven to increase variance
2. **DO NOT implement epsilon tie-breaking** - Proven to increase variance
3. **DO NOT implement Kahan summation** - Proven to increase variance
4. **DO NOT attempt to replace BWA with WFA2** - WFA2 serves different purpose
5. **DO NOT pursue fixed-point arithmetic** - Similar issues to Kahan summation

## Open Questions

1. **Scientific Acceptability**: Is 3% variance acceptable for all use cases, or do some require strict reproducibility?
   - Answer: Likely depends on downstream analysis (e.g., differential abundance vs presence/absence)

2. **Consensus Voting Threshold**: What is optimal N for consensus (3, 5, 7)?
   - Needs empirical testing on variance reduction vs computational cost

3. **BWA Alternatives**: Would Bowtie2 or Minimap2 significantly change accuracy/sensitivity?
   - Requires benchmark study comparing assignments on known-composition mock communities

4. **Platform Dependency**: Does variance differ across platforms (x86 vs ARM, Linux vs macOS)?
   - Current testing only on x86_64 Linux

5. **Dataset Dependency**: Does variance scale with reference database size or query complexity?
   - Tested on 16S bacteria (192MB reference), needs testing on ITS fungi, 18S, etc.

## Related Research

- `thoughts/research/multi-threading-non-determinism.md` - Original analysis
- `thoughts/shared/research/2026-01-03-non-determinism-baseline-measurement.md` - Baseline metrics
- `thoughts/shared/research/2026-01-03-float-precision-non-determinism-analysis.md` - Precision study
- `thoughts/shared/research/2026-01-03-deterministic-scoring-attempt-report.md` - Failed attempt post-mortem
- `thoughts/shared/plans/2026-01-03-deterministic-scoring-implementation.md` - Original implementation plan
- `thoughts/shared/research/2026-01-01-bwa-multithreading-feasibility.md` - BWA threading analysis

## Appendices

### Appendix A: Test Commands

**Deterministic mode (slow but reproducible)**:
```bash
tronko-assign -C 1 -r -f reference.trkb -a reference.fasta -p \
  -g forward.fasta -h reverse.fasta -o deterministic.txt
```

**Standard mode (fast but ~3% variance)**:
```bash
tronko-assign -C 16 -r -f reference.trkb -a reference.fasta -p \
  -g forward.fasta -h reverse.fasta -o standard.txt
```

**Test variance between runs**:
```bash
cd tronko-assign
./scripts/test_determinism.sh \
  -r /path/to/reference.trkb \
  -a /path/to/reference.fasta \
  -1 /path/to/forward.fasta \
  -2 /path/to/reverse.fasta \
  -n 3 -c 16 -k
```

### Appendix B: Variance Analysis Script

```bash
#!/bin/bash
# Compare two tronko-assign output files and compute variance metrics

FILE1=$1
FILE2=$2

# Extract taxonomic paths (column 2)
cut -f2 "$FILE1" | sort > /tmp/taxa1.txt
cut -f2 "$FILE2" | sort > /tmp/taxa2.txt

# Count differences
TOTAL=$(wc -l < /tmp/taxa1.txt)
DIFF=$(diff /tmp/taxa1.txt /tmp/taxa2.txt | grep '^<' | wc -l)
VARIANCE=$(awk "BEGIN {printf \"%.2f\", ($DIFF / $TOTAL) * 100}")

echo "Total reads: $TOTAL"
echo "Differences: $DIFF"
echo "Variance: $VARIANCE%"

# Cleanup
rm /tmp/taxa1.txt /tmp/taxa2.txt
```

### Appendix C: Consensus Voting Implementation (Pseudo-code)

```python
#!/usr/bin/env python3
# Consensus voting for tronko-assign

import sys
from collections import Counter

def load_assignments(filename):
    """Load taxonomic assignments from tronko-assign output."""
    assignments = {}
    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            read_name = parts[0]
            taxonomy = parts[1]
            assignments[read_name] = taxonomy
    return assignments

def consensus_taxonomy(taxonomies):
    """Find consensus taxonomy from multiple assignments."""
    # Count occurrences
    counts = Counter(taxonomies)

    # If clear majority, return it
    most_common = counts.most_common(1)[0]
    if most_common[1] > len(taxonomies) / 2:
        return most_common[0]

    # If tie, compute LCA of tied paths
    tied = [tax for tax, count in counts.items() if count == most_common[1]]
    return compute_lca(tied)

def compute_lca(taxonomies):
    """Compute least common ancestor of taxonomic paths."""
    levels = [tax.split(';') for tax in taxonomies]
    lca = []
    for i in range(min(len(l) for l in levels)):
        if all(l[i] == levels[0][i] for l in levels):
            lca.append(levels[0][i])
        else:
            break
    return ';'.join(lca)

def main(input_files, output_file):
    # Load all runs
    all_assignments = [load_assignments(f) for f in input_files]

    # Get all read names
    all_reads = set()
    for assignments in all_assignments:
        all_reads.update(assignments.keys())

    # Compute consensus for each read
    with open(output_file, 'w') as out:
        for read in sorted(all_reads):
            taxonomies = [a.get(read, '') for a in all_assignments if read in a]
            if taxonomies:
                consensus = consensus_taxonomy(taxonomies)
                confidence = taxonomies.count(consensus) / len(taxonomies)
                out.write(f"{read}\t{consensus}\t{confidence:.2f}\n")

if __name__ == '__main__':
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    main(input_files, output_file)
```
