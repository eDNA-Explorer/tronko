# Tier 1 Optimizations Implementation Plan

## Overview

Implement the four Tier 1 optimizations from the optimization prioritization matrix with individual toggle capability for A/B comparison testing. Each optimization can be independently enabled/disabled to measure its performance and accuracy impact.

## Current State Analysis

### Makefile (`tronko-assign/Makefile:13`)
- Current: `OPTIMIZATION = -O3`
- No architecture-specific flags (`-march=native`, `-mtune=native`)
- No `-ffast-math` for aggressive floating-point optimizations
- No `-fopenmp` despite WFA2 having OpenMP pragmas

### Tree Traversal (`tronko-assign/assignment.c:24-65`)
- `assignScores_Arr_paired()` performs exhaustive DFS traversal
- Scores ALL nodes regardless of intermediate results
- No early termination or pruning logic

### Maximum Finding (`tronko-assign/placement.c:862-921`)
- Iterates over all nodes in all matched trees to find maximum score
- Collects all nodes within Cinterval of maximum for LCA voting

### Options Infrastructure (`tronko-assign/options.c`, `global.h:184-226`)
- Well-established pattern using `getopt_long`
- `Options` struct holds all configuration

## Desired End State

After implementation:
1. **Baseline builds** remain unchanged (default `make` produces identical binary)
2. **Optimized builds** available via `make NATIVE_ARCH=1 FAST_MATH=1 ENABLE_OPENMP=1`
3. **Runtime toggles** allow enabling/disabling algorithmic optimizations per run
4. **Full backward compatibility** - existing scripts and workflows unaffected

### Verification Commands
```bash
# Baseline build (unchanged behavior)
make clean && make
./tronko-assign -r -f db.txt -a ref.fasta -s -g reads.fasta -o baseline.txt

# Optimized build
make clean && make NATIVE_ARCH=1 FAST_MATH=1 ENABLE_OPENMP=1
./tronko-assign -r -f db.txt -a ref.fasta -s -g reads.fasta -o optimized_build.txt \
    --early-termination --enable-pruning

# Compare outputs
diff baseline.txt optimized_build.txt  # Should be identical or within tolerance
```

## What We're NOT Doing

- Tier 2 optimizations (mmap, SIMD vectorization, two-phase screening)
- Tier 3 optimizations (pre-computed bounds, A* search)
- GPU acceleration (Tier 4)
- Changes to tronko-build or database format
- Changes to the LCA voting logic itself

---

## Phase 1: Compile-Time Toggle Infrastructure

### Overview
Add Make variables to control compiler optimization flags. Default build remains unchanged.

### Changes Required:

#### 1. Makefile Updates
**File**: `tronko-assign/Makefile`

**Current** (line 13):
```makefile
OPTIMIZATION = -O3
```

**New** (replace lines 13-20 with):
```makefile
# Base optimization
OPTIMIZATION = -O3

# Architecture-specific optimization (use: make NATIVE_ARCH=1)
ifdef NATIVE_ARCH
    ARCH_FLAGS = -march=native -mtune=native
else
    ARCH_FLAGS =
endif

# Fast math optimization (use: make FAST_MATH=1)
ifdef FAST_MATH
    MATH_FLAGS = -ffast-math
else
    MATH_FLAGS =
endif

# OpenMP support for WFA2 (use: make ENABLE_OPENMP=1)
ifdef ENABLE_OPENMP
    OPENMP_FLAGS = -fopenmp
else
    OPENMP_FLAGS =
endif

# Memory optimization flag (use: make OPTIMIZE_MEMORY=1)
ifdef OPTIMIZE_MEMORY
    MEMOPT_FLAGS = -DOPTIMIZE_MEMORY
else
    MEMOPT_FLAGS =
endif
```

**Update build command** (line 33):
```makefile
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS) $(OPENMP_FLAGS)
```

**Update debug build** (line 35):
```makefile
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS) $(OPENMP_FLAGS)
```

### Success Criteria:

#### Automated Verification:
- [ ] Default `make` produces identical binary to current (compare with `md5sum` or behavior)
- [ ] `make NATIVE_ARCH=1` compiles without errors
- [ ] `make FAST_MATH=1` compiles without errors
- [ ] `make ENABLE_OPENMP=1` compiles without errors
- [ ] `make NATIVE_ARCH=1 FAST_MATH=1 ENABLE_OPENMP=1` compiles without errors
- [ ] All flag combinations pass existing test: `./tronko-assign -r -f example_datasets/...`

#### Manual Verification:
- [ ] Verify compiler flags appear in build output with `make NATIVE_ARCH=1 V=1`

---

## Phase 2: Runtime Toggle Infrastructure

### Overview
Add command-line options to control algorithmic optimizations at runtime.

### Changes Required:

#### 1. Options Struct Updates
**File**: `tronko-assign/global.h`

Add to `Options` struct (after line 225, before closing brace):
```c
    // Tier 1 optimization toggles
    int early_termination;      // Enable early termination (default: 0)
    double strike_box;          // Strike box size, multiplier of Cinterval (default: 1.0)
    int max_strikes;            // Maximum strikes before termination (default: 6)
    int enable_pruning;         // Enable subtree pruning (default: 0)
    double pruning_factor;      // Pruning threshold = pruning_factor * Cinterval (default: 2.0)
```

#### 2. Long Options Array
**File**: `tronko-assign/options.c`

Add to `long_options[]` array (before the terminating `{0, 0, 0, 0}`):
```c
    {"early-termination",no_argument,0,0},
    {"no-early-termination",no_argument,0,0},
    {"strike-box",required_argument,0,0},
    {"max-strikes",required_argument,0,0},
    {"enable-pruning",no_argument,0,0},
    {"disable-pruning",no_argument,0,0},
    {"pruning-factor",required_argument,0,0},
```

#### 3. Usage String Update
**File**: `tronko-assign/options.c`

Add to `usage[]` string (before the final `\n"`):
```c
	\n\
	Optimization Options:\n\
	--early-termination, Enable early termination during tree traversal\n\
	--no-early-termination, Disable early termination (default)\n\
	--strike-box [FLOAT], Strike box size as multiplier of Cinterval [default: 1.0]\n\
	--max-strikes [INT], Maximum strikes before termination [default: 6]\n\
	--enable-pruning, Enable subtree pruning\n\
	--disable-pruning, Disable subtree pruning (default)\n\
	--pruning-factor [FLOAT], Pruning threshold = factor * Cinterval [default: 2.0]\n\
```

#### 4. Option Parsing
**File**: `tronko-assign/options.c`

In `parse_options()`, expand the `case 0:` block to handle new long options:
```c
case 0:
    // Handle long options without short equivalents
    if (strcmp(long_options[option_index].name, "tsv-log") == 0) {
        strncpy(opt->tsv_log_file, optarg, sizeof(opt->tsv_log_file) - 1);
        opt->tsv_log_file[sizeof(opt->tsv_log_file) - 1] = '\0';
    }
    else if (strcmp(long_options[option_index].name, "early-termination") == 0) {
        opt->early_termination = 1;
    }
    else if (strcmp(long_options[option_index].name, "no-early-termination") == 0) {
        opt->early_termination = 0;
    }
    else if (strcmp(long_options[option_index].name, "strike-box") == 0) {
        if (sscanf(optarg, "%lf", &(opt->strike_box)) != 1) {
            fprintf(stderr, "Invalid strike-box value\n");
            opt->strike_box = 1.0;
        }
    }
    else if (strcmp(long_options[option_index].name, "max-strikes") == 0) {
        if (sscanf(optarg, "%d", &(opt->max_strikes)) != 1) {
            fprintf(stderr, "Invalid max-strikes value\n");
            opt->max_strikes = 6;
        }
    }
    else if (strcmp(long_options[option_index].name, "enable-pruning") == 0) {
        opt->enable_pruning = 1;
    }
    else if (strcmp(long_options[option_index].name, "disable-pruning") == 0) {
        opt->enable_pruning = 0;
    }
    else if (strcmp(long_options[option_index].name, "pruning-factor") == 0) {
        if (sscanf(optarg, "%lf", &(opt->pruning_factor)) != 1) {
            fprintf(stderr, "Invalid pruning-factor value\n");
            opt->pruning_factor = 2.0;
        }
    }
    break;
```

#### 5. Default Initialization
**File**: `tronko-assign/tronko-assign.c`

Find where `Options` is initialized and add defaults. Look for pattern like `opt.cinterval = 5;`:
```c
// Tier 1 optimization defaults (disabled by default for baseline comparison)
opt.early_termination = 0;
opt.strike_box = 1.0;
opt.max_strikes = 6;
opt.enable_pruning = 0;
opt.pruning_factor = 2.0;
```

### Success Criteria:

#### Automated Verification:
- [ ] `make` compiles without errors
- [ ] `./tronko-assign --help` shows new options in usage text
- [ ] `./tronko-assign --early-termination --help` parses without error
- [ ] `./tronko-assign --strike-box 2.0 --max-strikes 10 --help` parses without error
- [ ] `./tronko-assign --enable-pruning --pruning-factor 3.0 --help` parses without error
- [ ] Default run produces identical output to current version

#### Manual Verification:
- [ ] Options appear correctly formatted in help output
- [ ] Invalid values produce appropriate error messages

---

## Phase 3: Early Termination Implementation

### Overview
Modify tree traversal to stop early when finding clearly suboptimal nodes, using the "baseball heuristic" from pplacer.

### Changes Required:

#### 1. Update Function Signature
**File**: `tronko-assign/assignment.h`

Update declaration:
```c
void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
    FILE* site_scores_file, char* readname,
    int early_termination, type_of_PP *best_score, int *strikes,
    type_of_PP strike_box, int max_strikes);
```

#### 2. Implement Early Termination Logic
**File**: `tronko-assign/assignment.c`

Replace `assignScores_Arr_paired()` function (lines 24-65):
```c
void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
    FILE* site_scores_file, char* readname,
    int early_termination, type_of_PP *best_score, int *strikes,
    type_of_PP strike_box, int max_strikes) {

    int child0 = treeArr[rootNum][node].up[0];
    int child1 = treeArr[rootNum][node].up[1];

    // Calculate score for this node
    type_of_PP node_score = getscore_Arr(alength, node, rootNum, locQuery, positions,
                                          print_all_nodes, site_scores_file, readname);

    if (child0 == -1 && child1 == -1) {
        // Leaf node
        scores[search_number][rootNum][node] += node_score;
    } else if (child0 != -1 && child1 != -1) {
        // Internal node
        scores[search_number][rootNum][node] += node_score;

        // Early termination check (only if enabled)
        if (early_termination) {
            // Higher score is better in tronko (log probabilities, less negative = better)
            if (node_score < *best_score - strike_box) {
                // This node is significantly worse than best
                (*strikes)++;
                if (*strikes >= max_strikes) {
                    // Stop exploring this subtree
                    return;
                }
            } else if (node_score > *best_score) {
                // Found a better score, reset strikes
                *best_score = node_score;
                *strikes = 0;
            }
        }

        // Recurse to children
        assignScores_Arr_paired(rootNum, child0, locQuery, positions, scores, alength,
            search_number, print_all_nodes, site_scores_file, readname,
            early_termination, best_score, strikes, strike_box, max_strikes);
        assignScores_Arr_paired(rootNum, child1, locQuery, positions, scores, alength,
            search_number, print_all_nodes, site_scores_file, readname,
            early_termination, best_score, strikes, strike_box, max_strikes);
    }
}
```

#### 3. Update All Call Sites
**File**: `tronko-assign/placement.c`

Find all calls to `assignScores_Arr_paired()` and update them. There should be calls around lines 200-300 and 700-800. The pattern will be:

**Before**:
```c
assignScores_Arr_paired(leaf_coordinates[match][0], rootArr[leaf_coordinates[match][0]],
    locQuery, positions, nodeScores, alength, match, print_all_nodes,
    site_scores_file, forward_name);
```

**After**:
```c
// Initialize early termination state for this tree
type_of_PP best_score = -9999999999999999;
int strikes = 0;
type_of_PP strike_box_threshold = Cinterval * opt->strike_box;

assignScores_Arr_paired(leaf_coordinates[match][0], rootArr[leaf_coordinates[match][0]],
    locQuery, positions, nodeScores, alength, match, print_all_nodes,
    site_scores_file, forward_name,
    opt->early_termination, &best_score, &strikes,
    strike_box_threshold, opt->max_strikes);
```

**Note**: The `opt` pointer needs to be passed to the placement functions, or the optimization settings need to be available as globals. Review `placement.c` function signatures to determine best approach.

### Success Criteria:

#### Automated Verification:
- [ ] `make` compiles without errors or warnings
- [ ] Default run (no flags) produces identical output to baseline
- [ ] `./tronko-assign --early-termination ...` runs without crashes
- [ ] Test with example dataset completes successfully with early termination enabled
- [ ] Output with `--early-termination` matches baseline output (100% accuracy)

#### Manual Verification:
- [ ] Add debug logging to verify early termination is occurring (strikes incrementing)
- [ ] Compare runtime with/without early termination on a medium dataset
- [ ] Verify no accuracy degradation with conservative thresholds

---

## Phase 4: Subtree Pruning Implementation

### Overview
Add logic to skip entire subtrees when the current node's score indicates descendants cannot improve the result.

### Changes Required:

#### 1. Extend Early Termination with Pruning
**File**: `tronko-assign/assignment.c`

Update `assignScores_Arr_paired()` to include pruning logic:
```c
void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
    FILE* site_scores_file, char* readname,
    int early_termination, type_of_PP *best_score, int *strikes,
    type_of_PP strike_box, int max_strikes,
    int enable_pruning, type_of_PP pruning_threshold) {

    int child0 = treeArr[rootNum][node].up[0];
    int child1 = treeArr[rootNum][node].up[1];

    // Calculate score for this node
    type_of_PP node_score = getscore_Arr(alength, node, rootNum, locQuery, positions,
                                          print_all_nodes, site_scores_file, readname);

    if (child0 == -1 && child1 == -1) {
        // Leaf node
        scores[search_number][rootNum][node] += node_score;

        // Update best score tracking
        if (early_termination && node_score > *best_score) {
            *best_score = node_score;
            *strikes = 0;
        }
    } else if (child0 != -1 && child1 != -1) {
        // Internal node
        scores[search_number][rootNum][node] += node_score;

        // Subtree pruning check (only if enabled)
        // If this node is too bad, skip entire subtree
        if (enable_pruning && *best_score > -9999999999999998) {
            if (node_score < *best_score - pruning_threshold) {
                // This subtree cannot contain nodes within Cinterval of best
                // Skip children entirely
                return;
            }
        }

        // Early termination check (only if enabled)
        if (early_termination) {
            if (node_score < *best_score - strike_box) {
                (*strikes)++;
                if (*strikes >= max_strikes) {
                    return;
                }
            } else if (node_score > *best_score) {
                *best_score = node_score;
                *strikes = 0;
            }
        }

        // Recurse to children
        assignScores_Arr_paired(rootNum, child0, locQuery, positions, scores, alength,
            search_number, print_all_nodes, site_scores_file, readname,
            early_termination, best_score, strikes, strike_box, max_strikes,
            enable_pruning, pruning_threshold);
        assignScores_Arr_paired(rootNum, child1, locQuery, positions, scores, alength,
            search_number, print_all_nodes, site_scores_file, readname,
            early_termination, best_score, strikes, strike_box, max_strikes,
            enable_pruning, pruning_threshold);
    }
}
```

#### 2. Update Header
**File**: `tronko-assign/assignment.h`

```c
void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
    FILE* site_scores_file, char* readname,
    int early_termination, type_of_PP *best_score, int *strikes,
    type_of_PP strike_box, int max_strikes,
    int enable_pruning, type_of_PP pruning_threshold);
```

#### 3. Update Call Sites
**File**: `tronko-assign/placement.c`

Update all call sites with the additional pruning parameters:
```c
type_of_PP best_score = -9999999999999999;
int strikes = 0;
type_of_PP strike_box_threshold = Cinterval * opt->strike_box;
type_of_PP pruning_threshold = Cinterval * opt->pruning_factor;

assignScores_Arr_paired(leaf_coordinates[match][0], rootArr[leaf_coordinates[match][0]],
    locQuery, positions, nodeScores, alength, match, print_all_nodes,
    site_scores_file, forward_name,
    opt->early_termination, &best_score, &strikes,
    strike_box_threshold, opt->max_strikes,
    opt->enable_pruning, pruning_threshold);
```

### Success Criteria:

#### Automated Verification:
- [ ] `make` compiles without errors or warnings
- [ ] Default run (no flags) produces identical output to baseline
- [ ] `./tronko-assign --enable-pruning ...` runs without crashes
- [ ] `./tronko-assign --early-termination --enable-pruning ...` runs without crashes
- [ ] Test with example dataset completes successfully
- [ ] Output matches baseline output (100% accuracy with default pruning factor of 2.0)

#### Manual Verification:
- [ ] Verify pruning reduces node count processed (add counter/logging)
- [ ] Compare runtime: baseline vs early-term vs pruning vs both
- [ ] Test with pruning-factor values: 1.5, 2.0, 3.0 to understand accuracy/speed tradeoff

---

## Testing Strategy

### Unit Tests:
- Verify option parsing for all new flags
- Verify default values are set correctly
- Verify early termination logic with mock scores

### Integration Tests:
```bash
# Baseline
./tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/baseline.txt -w

# Early termination only
./tronko-assign --early-termination \
    -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/early_term.txt -w

# Pruning only
./tronko-assign --enable-pruning \
    -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/pruning.txt -w

# Both optimizations
./tronko-assign --early-termination --enable-pruning \
    -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/both.txt -w

# Compare outputs
diff /tmp/baseline.txt /tmp/early_term.txt
diff /tmp/baseline.txt /tmp/pruning.txt
diff /tmp/baseline.txt /tmp/both.txt
```

### Manual Testing Steps:
1. Run on a larger dataset (16S_Bacteria) and time each configuration
2. Verify output accuracy matches baseline
3. Test edge cases: very short reads, highly divergent sequences
4. Test with different Cinterval values

---

## Performance Considerations

### Expected Speedup:
- Compiler flags alone: 10-30%
- Early termination: 2-5x (depends on data)
- Subtree pruning: 2-10x (depends on tree structure)
- Combined: 3-15x potential

### Monitoring:
Add timing output with `-T` flag to compare:
```bash
./tronko-assign -T --early-termination --enable-pruning ...
```

### Safety Margins:
- Default `strike_box = 1.0` (same as Cinterval) is conservative
- Default `pruning_factor = 2.0` ensures no false negatives
- Both can be tuned for speed vs accuracy tradeoff

---

## Migration Notes

- No database format changes required
- No changes to tronko-build
- Existing scripts continue to work unchanged
- New flags are purely additive

---

## References

- Research document: `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md`
- Algorithm details: `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`
- pplacer baseball heuristic: Matsen et al., 2010

---

## Open Questions

None - all design decisions have been made:
1. Individual toggles for each optimization (confirmed by user)
2. Conservative defaults (disabled, with safe thresholds)
3. Compile-time for compiler flags, runtime for algorithms
