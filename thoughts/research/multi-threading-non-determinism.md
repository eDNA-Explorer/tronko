# Multi-Threading Non-Determinism in tronko-assign

## Summary

tronko-assign produces ~0.28% different taxonomic assignments between runs when using multiple threads, even with identical inputs. This is due to floating-point arithmetic order sensitivity combined with work-stealing parallelism.

## Evidence

Tested with 16S_Bacteria database (17,868 trees) and 100,000 paired-end reads:

| Comparison | Differences | Agreement |
|------------|-------------|-----------|
| Zstd Run 1 vs Uncompressed | 277 | 99.72% |
| Zstd Run 1 vs Zstd Run 2 | 278 | 99.72% |

The variability is identical whether comparing different database formats or the same database twice, confirming it's run-to-run threading variability, not a data issue.

### Direction of Differences

| Direction | Count | Percentage |
|-----------|-------|------------|
| Run A more specific | 127 | 45.8% |
| Run B more specific | 121 | 43.7% |
| Same depth, different taxon | 29 | 10.5% |

The roughly 50/50 split confirms random variability, not systematic bias.

## Root Causes

### 1. BWA Work Stealing (`bwa_source_files/kthread.c:25-32`)

```c
static inline long steal_work(kt_for_t *t) {
    int i, min_i = -1;
    long k, min = LONG_MAX;
    for (i = 0; i < t->n_threads; ++i)
        if (min > t->w[i].i) min = t->w[i].i, min_i = i;
    k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
    return k >= t->n? -1 : k;
}
```

BWA uses a work-stealing scheduler where threads can complete alignment tasks in different orders each run. This means the order of alignment results can vary.

### 2. Tie-Breaking with Cinterval (`placement.c:930`)

```c
if (nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) &&
    nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval)) {
    voteRoot[leaf_coordinates[i][0]][k]=1;
    index++;
}
```

When multiple nodes have scores within `Cinterval` of the maximum, they're all considered "tied" and an LCA is computed. Small floating-point rounding differences can cause nodes to fall in/out of this interval between runs.

### 3. Score Accumulation Order (`assignment.c:207-224`)

```c
score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
// ... repeated for each position
```

Log-probability additions are accumulated in a loop. Floating-point addition is not associative:
- `(a + b) + c ≠ a + (b + c)` due to rounding

If BWA returns alignments in different orders, the comparison between candidate nodes may involve slightly different accumulated scores.

## Potential Solutions

### Option 1: Sort BWA Results (Low Complexity)

Sort alignment results by (tree_id, position) before scoring to ensure consistent processing order regardless of thread scheduling.

**Pros:**
- Minimal code change
- Preserves all parallelism
- No accuracy impact

**Cons:**
- Small overhead for sorting
- Doesn't fix floating-point accumulation issue entirely

### Option 2: Stable Tie-Breaking (Low Complexity)

When scores are within floating-point epsilon (~1e-10), use deterministic tie-breaker:
1. Prefer lower tree_id
2. Prefer lower node_id

**Pros:**
- Very simple to implement
- No performance impact
- Handles the Cinterval edge cases

**Cons:**
- Doesn't address root cause (just masks it)

### Option 3: Kahan Summation (Medium Complexity)

Replace simple accumulation with compensated summation:

```c
// Current
score += value;

// Kahan summation
y = value - c;
t = score + y;
c = (t - score) - y;
score = t;
```

**Pros:**
- Reduces floating-point error by orders of magnitude
- Makes order matter much less

**Cons:**
- More complex code
- Small performance overhead
- Doesn't guarantee bit-exact results

### Option 4: Fixed-Point Scoring (High Complexity)

Convert log-probabilities to fixed-point integers for accumulation:

```c
// Convert to fixed-point (e.g., multiply by 2^20)
int64_t score_fixed = (int64_t)(log_prob * 1048576.0);
// Accumulate in integer
total_fixed += score_fixed;
// Convert back for final comparison
double final_score = total_fixed / 1048576.0;
```

**Pros:**
- Completely deterministic
- Integer addition is associative

**Cons:**
- Significant code changes
- Potential precision loss at extremes
- Need to handle overflow

### Option 5: Reproducible Math Libraries

Use libraries like Intel MKL with Conditional Numerical Reproducibility (CNR) mode.

**Pros:**
- Industry-standard solution
- Well-tested

**Cons:**
- External dependency
- May not cover all operations
- Platform-specific

## Recommendation

**Short-term (Quick Fix):**
1. Sort BWA results before processing
2. Add stable tie-breaking with tree_id/node_id

**Long-term (If Determinism Required):**
- Implement Kahan summation for score accumulation
- Consider fixed-point for critical path if needed

## Impact Assessment

The 0.28% variability affects:
- **CI/CD testing**: Cannot use exact output comparison
- **Reproducibility**: Same input may give slightly different results
- **Benchmarking**: Need statistical comparison, not exact match

However, the variability is:
- **Small**: 99.72% agreement
- **Unbiased**: 50/50 split between more/less specific
- **Scientifically acceptable**: Within normal biological uncertainty

## References

- `tronko-assign/bwa_source_files/kthread.c` - Work stealing implementation
- `tronko-assign/placement.c:930` - Cinterval tie-breaking
- `tronko-assign/assignment.c:172-238` - Score accumulation
- [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
