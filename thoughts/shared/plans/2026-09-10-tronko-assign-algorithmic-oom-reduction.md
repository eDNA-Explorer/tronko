# Tronko Assign Algorithmic OOM Reduction Implementation Plan

## Overview

Reduce `tronko-assign` peak memory on the `high-perf` branch without changing assignment parameters, numeric precision, BWA behavior, or infrastructure sizing. The first implementation series will make the BWA index and accession lookup process-owned, then replace reference-sized per-worker placement structures with bounded candidate-local workspaces.

The work is intentionally incremental. Each phase must preserve assignment output and produce a measurable memory reduction before the next phase begins. After the simple exact changes are measured, production suitability will be evaluated separately against the deployed Kubernetes resources. Larger reference-format and lazy-loading changes will be planned only if the simple changes do not resolve the operational OOM problem.

## Current State Analysis

`tronko-assign` already shares the loaded Tronko tree database across assignment workers, but it does not share all of the other reference-sized state:

- The process creates an outer pthread for each assignment chunk in `tronko-assign.c`.
- Every worker calls `run_bwa()`, which calls `main_mem()` independently ([`tronko-assign/tronko-assign.c:173`](../../../tronko-assign/tronko-assign.c#L173)).
- Every invocation of `main_mem()` loads and later destroys `BWA_IDX_ALL` ([`tronko-assign/bwa_source_files/fastmap.c:787`](../../../tronko-assign/bwa_source_files/fastmap.c#L787)). With 16 assignment workers, this can make 16 complete BWA index copies resident at once, and the work repeats for subsequent batches.
- Each `main_mem()` invocation also builds its own leaf-accession hashmap. Duplicate accessions can leak the rejected value because `hashmap_put()` failure is ignored.
- Each assignment worker allocates `nodeScores[MAX_NUM_BWA_MATCHES][numberOfTrees][nodesInTree]`, even though a read can use at most `MAX_NUM_BWA_MATCHES` candidate trees ([`tronko-assign/allocateMemoryForResults.c:12`](../../../tronko-assign/allocateMemoryForResults.c#L12)).
- Each worker additionally allocates full-reference `voteRoot`, `minNodes`, and `LCAnames` structures plus `leaf_coordinates[numberOfTrees]` ([`tronko-assign/allocateMemoryForResults.c:22`](../../../tronko-assign/allocateMemoryForResults.c#L22), [`tronko-assign/allocateMemoryForResults.c:75`](../../../tronko-assign/allocateMemoryForResults.c#L75)).
- Every assigned ASV clears or scans structures covering all reference nodes and trees even though only its retained BWA candidates can contribute ([`tronko-assign/tronko-assign.c:602`](../../../tronko-assign/tronko-assign.c#L602), [`tronko-assign/tronko-assign.c:626`](../../../tronko-assign/tronko-assign.c#L626)).
- The placement interface carries the dense three-dimensional score array and global vote array through both the WFA and Needleman-Wunsch paths ([`tronko-assign/placement.h:13`](../../../tronko-assign/placement.h#L13)).
- The current automated tests exercise basic workflows, but do not comprehensively lock down all seven assignment-output columns, BWA candidate identity, ties, or memory ownership. The Valgrind CI command also ignores failures.

For a historical reference containing 1,265,572 nodes, 16 workers, 25 retained BWA matches, and double scores, the dense `nodeScores` payload alone is approximately 3.77 GiB. That excludes its pointer/allocation topology, the other per-worker arrays, the shared posterior database, and the repeated BWA indexes. Repeated BWA index residency is therefore the first and highest-leverage target.

## Desired End State

After this plan:

1. One immutable BWA index is loaded for the process and shared by all assignment workers and batches.
2. One immutable accession-to-Tronko-coordinate map is built for the process and shared by all workers.
3. Every worker owns only mutable query/chunk state and a reusable placement workspace.
4. Placement score storage grows with the nodes in the ASV's retained candidate trees, not with all nodes in the reference multiplied by the match cap.
5. Winner counting and LCA calculation scan only active candidate trees and require no full-reference vote, name, or temporary-node arrays.
6. The same inputs and parameters produce byte-identical single-thread assignment rows across all seven columns.
7. The production-like FWH workload completes with no batch-to-batch RSS growth and a measured peak lower than the Phase 1 baseline.

Code-level acceptance does not set or assume a Kubernetes memory limit. After these changes pass their parity and memory-behavior gates, validate them under the deployed workload and resource configuration. If that operational validation still OOMs, collect a new memory breakdown and create a separate plan for reference packing and candidate-tree lazy loading.

### Key Discoveries

- BWA's `mem_process_seqs()` accepts the BWT, sequence annotations, and packed reference as read-only inputs; query-specific work structures are local to the call. Sharing a fully initialized index is compatible with the existing BWA execution model.
- The current `ignore_alt` path is disabled. Any one-time index mutation must be completed before the context is published to workers; workers must never mutate the shared index.
- `global_bns` in [`tronko-assign/bwa_source_files/bwamem.c:45`](../../../tronko-assign/bwa_source_files/bwamem.c#L45) is debugging state and becomes a formal data race when multiple workers share the index.
- `bwa_verbose` is assigned inside `main_mem()` ([`tronko-assign/bwa_source_files/fastmap.c:590`](../../../tronko-assign/bwa_source_files/fastmap.c#L590)); it must be configured once before workers start.
- `LCAnames` is populated during result reduction but is not consumed when making the assignment decision ([`tronko-assign/tronko-assign.c:670`](../../../tronko-assign/tronko-assign.c#L670)). It can be removed after parity coverage exists.
- Existing experimental `OPTIMIZE_MEMORY`/float builds are not an exact replacement: prior project measurements recorded assignment differences. Numeric representation must remain unchanged in this plan.
- The current `.trkb` representation is zstd-compressed and posterior values are transformed in place after loading. Direct mmap is therefore not a simple follow-on to this work.

## What We're NOT Doing

- Replacing BWA with minimap2 or adding any minimap2 execution path.
- Reducing thread/core count.
- Reducing the FASTA batch-size limit.
- Reducing `MAX_NUM_BWA_MATCHES` or changing the downstream production value.
- Enabling `OPTIMIZE_MEMORY`, converting scores from double to float, or otherwise changing numeric precision.
- Changing confidence intervals, alignment scores, pruning, early termination, or assignment thresholds.
- Increasing or reducing Kubernetes memory requests/limits as part of the implementation.
- Replacing the current outer pthread/chunk topology with BWA internal whole-batch threading. Chunk boundaries influence paired-end insert-size inference and can change tie identifiers.
- Bypassing SAM serialization and parsing with a new typed BWA-result API in this series.
- Flattening the entire Tronko reference, changing `.trkb`, memory-mapping it, or implementing lazy candidate-tree loading.
- Reordering outputs or changing diagnostics such as `site_scores.txt` and `scores_all_nodes.txt`.

## Implementation Approach

Use four independently reviewable implementation phases:

1. Establish parity, sanitizer, and memory baselines.
2. Make the BWA index and leaf map immutable process-owned state.
3. Replace reference-sized worker allocations with a bounded active-candidate score arena.
4. Reduce winners and LCAs directly from candidate score slices, remove dead storage, and run the production-like stop gate.

Preserve the current worker topology and observable behavior while changing ownership first, representation second, and reduction logic last. Do not keep permanent legacy/new runtime switches: use the baseline executable and golden fixtures for comparison, then delete superseded allocation paths after each phase passes.

## Phase 1: Lock Down Semantics and Memory Measurement

### Overview

Create a reliable safety net before changing BWA ownership or placement representation. This phase must be merged first because the existing tests do not sufficiently detect assignment drift.

### Changes Required

#### 1. Add deterministic assignment fixtures and golden outputs

**Files**:

- `tests/data/assignment/` (new focused fixtures)
- `tests/integration/test_assignment_parity.sh` (new)
- `tests/Makefile.simple`
- `run_tests.sh`

**Changes**:

- Add small checked-in references and queries covering:
  - Single-end and paired-end reads.
  - A single candidate tree and multiple candidate trees.
  - Forward and reverse-complement alignments.
  - Concordant and discordant paired hits.
  - No BWA hit.
  - More BWA hits than the retained cap.
  - Duplicate reference accessions.
  - A global score tie.
  - Scores exactly on both inclusive `Cinterval` boundaries.
  - Multiple winning nodes whose result is an internal-node LCA.
  - Equal winner counts in multiple trees, proving the lowest numeric tree ID remains the tie winner.
- Archive an executable from the unmodified `high-perf` commit, then add only the diagnostic instrumentation below. Verify that the instrumented executable is output-identical to the archived executable before adopting it as the Phase 1 baseline.
- Generate expected outputs from the instrumented Phase 1 baseline at `-C 1`.
- Compare every output row and all seven columns byte-for-byte; do not limit comparison to taxonomy columns.
- Run text, binary, and zstd reference variants through the same cases. Give each format its own baseline output; require old-versus-new byte parity within a format. Across formats, preserve the existing critical-column comparison and any explicitly documented numeric tolerance rather than assuming serialized floating-point scores are byte-identical.
- Keep the baseline executable outside version control and record its commit SHA in the test run output.

#### 2. Add candidate-level parity diagnostics

**Files**:

- `tronko-assign/bwa_source_files/fastmap.c`
- `tronko-assign/global.h`
- `tests/integration/test_assignment_parity.sh`

**Changes**:

- In the existing verbose/TSV diagnostic path, emit a stable candidate fingerprint containing:
  - Read identifier.
  - Concordant/discordant status.
  - Tree ID and node ID.
  - Forward/reverse start position.
  - Forward/reverse CIGAR.
  - Retained/dropped reason.
- Keep this diagnostic disabled during normal runs.
- Compare candidate fingerprints separately from final assignments. This distinguishes BWA/context drift from placement/reduction drift.

#### 3. Make memory-safety builds real and gating

**Files**:

- `tronko-assign/Makefile`
- `tests/Makefile.simple`
- `.github/workflows/test.yml`

**Changes**:

- Allow the build to accept additive `CPPFLAGS`, `CFLAGS`, and `LDFLAGS` instead of overwriting caller-provided sanitizer settings.
- Add AddressSanitizer/UndefinedBehaviorSanitizer and LeakSanitizer targets.
- Add a ThreadSanitizer target for the later shared-context phase.
- Remove `|| true` or equivalent status suppression from the Valgrind job.
- Keep sanitizer fixtures small enough for CI while preserving the larger local/profiling scenario.

#### 4. Record reproducible peak-memory baselines

**Files**:

- `tronko-assign/scripts/compare_memlogs.sh`
- `tests/integration/test_assignment_parity.sh`
- `docs/performance-logging.md`

**Changes**:

- Extend existing memory logging with phase labels for:
  - Reference loaded.
  - BWA context/index loaded.
  - Leaf map ready.
  - Worker workspaces allocated.
  - Each batch aligned.
  - Each batch placed and freed.
  - Final teardown.
- Log BWA index/map construction and destruction counts.
- Capture both process peak RSS and post-batch RSS.
- Record the exact executable SHA, reference/query hashes, cores, batch size, match cap, and scoring parameters with every comparison.
- Establish baselines for the small fixture and the production-like FWH workload before Phase 2.

### Success Criteria

#### Automated Verification

- [ ] `make -C tronko-assign clean && make -C tronko-assign` succeeds on the unmodified algorithm.
- [ ] `./run_tests.sh --unit-only` passes.
- [ ] `./run_tests.sh --integration-only` passes.
- [ ] The new parity test checks all output columns at `-C 1`.
- [ ] Text, binary, and zstd fixtures each match their own Phase 1 golden output, and the existing cross-format critical-column comparison passes.
- [ ] ASan/UBSan/LSan failures fail the test command.
- [ ] Valgrind findings fail CI rather than being ignored.
- [ ] The baseline candidate fingerprints and output goldens are recorded with the baseline commit SHA.

#### Manual Verification

- [ ] Baseline peak RSS and post-batch RSS are recorded for the production-like FWH dataset.
- [ ] The dataset manifest includes immutable hashes for the reference and ASV inputs.
- [ ] The baseline run uses the current production core count, batch size, match cap, and scoring parameters.

---

## Phase 2: Share the BWA Index and Leaf Map

### Overview

Remove the largest expected memory multiplier without changing BWA's alignment algorithm, chunking, or candidate ordering.

### Changes Required

#### 1. Introduce an opaque process-lifetime BWA context

**Files**:

- `tronko-assign/bwa_context.h` (new)
- `tronko-assign/bwa_context.c` (new)
- `tronko-assign/bwa_source_files_include.h`
- `tronko-assign/Makefile`

**Changes**:

- Define an opaque context API:

```c
typedef struct tronko_bwa_context tronko_bwa_context;

int tronko_bwa_context_create(
    const char *database_file,
    int number_of_trees,
    tronko_bwa_context **out_context);

int tronko_bwa_align_chunk(
    const tronko_bwa_context *context,
    int start,
    int end,
    bwaMatches *results,
    const tronko_bwa_chunk_options *options);

void tronko_bwa_context_destroy(tronko_bwa_context *context);
```

- Make the context own the `bwaidx_t` and the leaf-coordinate lookup data.
- Load `BWA_IDX_ALL` once in `tronko_bwa_context_create()`.
- Complete any initialization before returning the context, then expose only `const` access to workers.
- Keep BWA `mem_opt_t`, sequence buffers, insert-size inference, output buffers, and other query-specific state local to each alignment call.
- Use one cleanup path so partial initialization and worker failures cannot leak the index or map.

#### 2. Separate BWA alignment from index lifetime

**Files**:

- `tronko-assign/bwa_source_files/fastmap.c`
- `tronko-assign/bwa_source_files/bwa.c`
- `tronko-assign/bwa_source_files/bwa.h`
- `tronko-assign/bwa_source_files/bwamem.c`

**Changes**:

- Refactor `main_mem()` so alignment accepts a previously loaded immutable index instead of loading and destroying it per call.
- Preserve the existing number and boundaries of read chunks.
- Preserve current single-end and paired-end `mem_process_seqs()` calls and ordering.
- Remove the unused `global_bns` assignment or replace it with call-local debug data; do not introduce a lock solely for unused debug state.
- Set `bwa_verbose` during single-threaded process initialization and never write it from a worker.
- Keep the disabled alternate-contig behavior unchanged. If alternate-contig flags ever require mutation, apply it once before publishing the context.
- Return explicit status codes for initialization, alignment, parsing, and teardown errors.

#### 3. Build one immutable leaf-coordinate map

**Files**:

- `tronko-assign/bwa_context.c`
- `tronko-assign/hashmap.c`
- `tronko-assign/hashmap.h`

**Changes**:

- Build the accession map once after the Tronko reference metadata is available.
- Store `{tree_id, node_id}` values in one contiguous coordinate array with stable addresses.
- Borrow immutable leaf-name strings as keys rather than allocating a second copy.
- Reserve sufficient hashmap capacity before insertion.
- Preserve the existing first-seen-wins behavior for duplicate accessions.
- On a duplicate, record the collision and do not leak or replace the retained coordinate.
- Publish the completed map as read-only. Workers may perform lookups but never inserts, deletes, or resizes.

#### 4. Wire the context into the process and worker lifecycle

**Files**:

- `tronko-assign/global.h`
- `tronko-assign/tronko-assign.c`
- `tronko-assign/bwa_source_files_include.h`

**Changes**:

- Add `const tronko_bwa_context *bwa_context` to the worker argument structure.
- Create the context after the tree/reference metadata and BWA database path are ready, before the first assignment batch starts.
- Pass the same pointer to every worker in every batch.
- Replace `run_bwa()`'s database-file ownership with a call to `tronko_bwa_align_chunk()`.
- Join all workers before destroying the context.
- Destroy the context before the tree leaf-name storage borrowed by the leaf map is freed.
- Aggregate worker statuses and abort the batch cleanly on the first failed alignment without returning partial output as success.

### Correctness Invariants

- BWA options and scoring are unchanged.
- Existing read-chunk boundaries are unchanged.
- Candidate order and primary/secondary alignment behavior are unchanged.
- Paired insert-size inference sees exactly the same reads per BWA call as before.
- Shared index and map data are immutable after worker startup.
- A duplicate accession maps to the same first coordinate as the baseline.

### Success Criteria

#### Automated Verification

- [ ] All Phase 1 golden outputs remain byte-identical at `-C 1`.
- [ ] Candidate fingerprints remain identical to the Phase 1 baseline.
- [ ] Single-end, paired-end, duplicate-accession, and multi-batch fixtures pass.
- [ ] Instrumentation reports exactly one BWA index construction and one destruction per process.
- [ ] Instrumentation reports exactly one leaf-map construction and one destruction per process.
- [ ] ASan/UBSan/LSan report no index, duplicate-map, failure-path, or teardown leaks.
- [ ] TSAN reports no index, leaf-map, `global_bns`, or verbosity races.
- [ ] Existing unit and integration tests pass.

#### Manual Verification

- [ ] Peak and post-batch RSS are compared against the Phase 1 baseline using identical inputs and parameters.
- [ ] RSS no longer steps upward when a new batch starts because another set of BWA indexes/maps was created.
- [ ] The production-like FWH run completes and its outputs match the baseline.

---

## Phase 3: Allocate Placement Memory Only for Active Candidates

### Overview

Remove dense `MAX_NUM_BWA_MATCHES × all reference nodes` scoring storage and reference-sized per-worker candidate arrays while preserving placement arithmetic and scan order.

### Changes Required

#### 1. Introduce a reusable candidate workspace

**Files**:

- `tronko-assign/candidate_workspace.h` (new)
- `tronko-assign/candidate_workspace.c` (new)
- `tronko-assign/Makefile`
- `tronko-assign/global.h`

**Changes**:

- Define a fixed-capacity candidate array and a grow-only score arena:

```c
typedef struct {
    int tree_id;
    int leaf_node;
    size_t score_offset;
    size_t node_count;
    int forward_start;
    int reverse_start;
    char *forward_cigar;
    char *reverse_cigar;
} CandidatePlacement;

typedef struct {
    CandidatePlacement candidates[MAX_NUM_BWA_MATCHES];
    size_t candidate_count;
    type_of_PP *score_arena;
    size_t score_capacity;
    size_t score_length;
} CandidateWorkspace;
```

- Allocate one workspace per assignment worker.
- Retain no more than `MAX_NUM_BWA_MATCHES` candidates.
- Extract candidate selection/deduplication from the worker loop into a testable function.
- Preserve first-retained-hit-per-tree behavior, candidate cap behavior, and all overflow counters.
- Use a bounded linear lookup for deduplication. With at most 25 candidates, this avoids a new per-ASV hashmap and keeps behavior explicit.
- Grow `score_arena` only when `sum(nodes in retained candidate trees)` exceeds its capacity; reuse the high-water allocation for later ASVs.
- Assign each candidate an offset and node count into the arena.
- Clear only `score_length` elements, never the unused capacity or the full reference.
- Check all node-count sums and byte multiplications for overflow before allocation.

#### 2. Replace the dense result allocation

**Files**:

- `tronko-assign/allocateMemoryForResults.c`
- `tronko-assign/allocateMemoryForResults.h`
- `tronko-assign/global.h`

**Changes**:

- Remove allocation and teardown of:
  - `nodeScores[MAX_NUM_BWA_MATCHES][numberOfTrees][nodesInTree]`.
  - `leaf_coordinates[numberOfTrees]`.
- Store the candidate workspace in the per-worker result/scratch structure.
- Keep unrelated alignment scratch allocations unchanged in this phase.
- Ensure all partial-allocation failure paths free the workspace exactly once.

#### 3. Change placement and scoring to accept candidate slices

**Files**:

- `tronko-assign/assignment.c`
- `tronko-assign/assignment.h`
- `tronko-assign/placement.c`
- `tronko-assign/placement.h`

**Changes**:

- Change `assignScores_Arr_paired()` from three-dimensional indexing to a single candidate score slice:

```c
void assignScores_Arr_paired(
    int tree_id,
    int root_node,
    char *query,
    int *positions,
    type_of_PP *candidate_scores,
    int alignment_length,
    ...);
```

- Pass a `CandidatePlacement` and its `score_arena + score_offset` slice through both WFA and Needleman-Wunsch paths.
- Preserve the exact recursive child traversal and node write order.
- Preserve forward then reverse `+=` accumulation order for paired reads.
- Preserve all optional diagnostic output and its tree/node ordering.
- Determine the global maximum by scanning candidates in retained order and nodes in ascending numeric order, using the current strict `>` comparison.

#### 4. Remove full-tree candidate arrays and fix bounded loops

**Files**:

- `tronko-assign/tronko-assign.c`

**Changes**:

- Replace `trees_search[numberOfTrees]` with `workspace.candidate_count` and candidate-local lookups.
- Remove the two unused `max_query_length` worker stack arrays and other unused leaf-map/query locals.
- Audit every BWA-result and retained-candidate loop. Arrays sized to `MAX_NUM_BWA_MATCHES` must never be indexed by `numberOfTrees`.
- Keep result slots indexed by input position so worker completion order cannot reorder the output.

### Correctness Invariants

- Candidate retention and deduplication match the Phase 1 baseline.
- Scores use the same `type_of_PP` as the current build.
- Forward and reverse contributions are accumulated in the same order.
- Global maximum uses strict `>` and the same candidate/node traversal order.
- Diagnostic node/site score order remains unchanged.
- Mutable arenas are worker-private; no workspace pointer is shared between pthreads.

### Success Criteria

#### Automated Verification

- [ ] All Phase 1 golden outputs remain byte-identical at `-C 1`.
- [ ] Candidate fingerprints remain identical.
- [ ] Candidate-workspace unit tests cover zero, one, exactly `M`, and more than `M` input hits.
- [ ] Duplicate-tree tests prove the current first-retained-hit policy.
- [ ] Allocation-overflow and failed-growth paths return errors without corrupting the previous arena.
- [ ] A test-only allocation counter proves score payload equals the sum of active candidate-tree nodes, not `M × total_reference_nodes`.
- [ ] ASan/UBSan/LSan are clean at candidate counts `0`, `1`, `M - 1`, and `M`.
- [ ] Existing unit and integration tests pass.

#### Manual Verification

- [ ] Phase 3 peak RSS is compared with the Phase 1 and Phase 2 measurements using identical inputs and parameters.
- [ ] Repeated small-tree and single-large-tree candidate cases confirm the expected high-water behavior.
- [ ] No post-warm-up RSS growth occurs across repeated batches.

---

## Phase 4: Reduce Winners and LCAs Without Full-Reference Storage

### Overview

Remove the remaining full-reference worker arrays and `O(all reference nodes)` work performed for every ASV. Calculate candidate vote counts and LCAs directly while scanning active score slices.

### Changes Required

#### 1. Add candidate-local winner summaries

**Files**:

- `tronko-assign/candidate_workspace.h`
- `tronko-assign/candidate_workspace.c`
- `tronko-assign/placement.c`

**Changes**:

- Extend each candidate with:

```c
typedef struct {
    /* candidate identity and score slice */
    int winning_count;
    int lca_node;
    bool has_winner;
} CandidatePlacement;
```

- After the global maximum is known, scan only each candidate's score slice.
- A node wins only under the current inclusive predicate:

```c
score >= maximum - Cinterval &&
score <= maximum + Cinterval
```

- Initialize `lca_node` with the first winning node and fold subsequent winners using `getLCA_Arr()` in ascending node order.
- Count every unique winning node once.
- Build one winner summary for every tree with a positive winning-node count, enumerated in ascending numeric tree-ID order.
- Preserve the existing `maxRoot` rule used by the output fields: higher winner count wins; equal positive counts retain the lower numeric tree ID because the scan uses strict `>`.
- For a single positive tree, fold all of that tree's winning nodes into its LCA.
- For multiple positive trees, compute one LCA per tree, then preserve the current taxonomic reconciliation across every positive tree in ascending numeric tree-ID order. Do not replace this behavior with a winner-takes-all tree selection.
- Preserve the current output-field behavior for `maxRoot`, `LCA`, `taxRoot`, and `taxNode`, including cases where the multi-tree taxonomy is reconciled separately from the reported maximum-vote tree.

#### 2. Remove global voting and dead result storage

**Files**:

- `tronko-assign/allocateMemoryForResults.c`
- `tronko-assign/allocateMemoryForResults.h`
- `tronko-assign/global.h`
- `tronko-assign/tronko-assign.c`

**Changes**:

- Remove `voteRoot` allocation, full-reference clearing, and full-tree vote scans.
- Remove `minNodes` allocation and per-ASV clearing.
- Remove `LCAnames`; it is not part of the decision or final output.
- Replace `countVotes[numberOfTrees]`, `LCAs[count]`, and `maxRoots[count]` VLAs with bounded per-tree summaries derived from retained candidates.
- Delete `getLCAofArray_Arr_Multiple()` and `LCA_of_nodes()` only after equivalence tests pass and no diagnostic path calls them.
- Preserve existing no-hit and unassigned output behavior.

#### 3. Remove the multi-tree output-string leak

**Files**:

- `tronko-assign/tronko-assign.c`

**Changes**:

- Replace the repeated `asprintf()` score/root/LCA fragments with a single size-checked `snprintf()` into the existing result buffer or a growable helper with one owner.
- Check truncation and return an explicit worker error rather than emitting a partial row.
- Apply the same ownership rule to both single-tree and multi-tree branches.

#### 4. Add old/new LCA equivalence tests

**Files**:

- `tests/unit/test_assignment.c` (new)
- `tests/Makefile.simple`

**Changes**:

- Before deleting the old array-based LCA code, compare it with the incremental fold over:
  - One winning node.
  - Sibling leaves.
  - Ancestor/descendant winners.
  - Winners spanning the root.
  - Repeated/duplicate inputs.
  - Generated valid binary trees and randomized winner subsets with a fixed seed.
- Test the complete tree-selection reduction separately from LCA folding.
- Assert the current tree-ID tie behavior explicitly.

### Correctness Invariants

- Global maximum selection is unchanged.
- Confidence interval endpoints remain inclusive.
- Winning-node set semantics remain a boolean union.
- Incremental LCA matches the current array-based algorithm for every tested winner set.
- Equal tree vote counts retain the lowest numeric `maxRoot`, while multi-tree taxonomy reconciliation still considers every tree with positive votes.
- No-hit, unassigned, single-tree, and multi-tree rows are byte-identical.

### Success Criteria

#### Automated Verification

- [ ] All Phase 1 golden outputs remain byte-identical at `-C 1`.
- [ ] Generated-tree LCA equivalence tests pass for the fixed random seed.
- [ ] No references to `voteRoot`, `minNodes`, or `LCAnames` remain in production code.
- [ ] No per-ASV VLA or heap allocation is sized by total tree or node count.
- [ ] ASan/UBSan/LSan report no multi-tree output leak.
- [ ] TSAN remains clean with the shared BWA context.
- [ ] `./run_tests.sh --unit-only` passes.
- [ ] `./run_tests.sh --integration-only` passes.
- [ ] `./run_tests.sh --valgrind` passes and returns a nonzero status on injected leak failures.

#### Manual Verification

- [ ] Small-fixture memory logs show candidate-local score and reduction storage.
- [ ] The production-like FWH run completes with the same runtime parameters and byte-identical `-C 1` output.
- [ ] Multicore results are no more variable than the frozen Phase 1 baseline.
- [ ] Post-warm-up RSS remains stable across batches.
- [ ] Peak RSS is lower than the Phase 1 baseline on the production-like workload, with the absolute and relative reduction recorded.

---

## Testing Strategy

### Unit Tests

- BWA-context partial initialization and cleanup.
- Duplicate leaf accessions and first-seen coordinate retention.
- Candidate insertion, deduplication, cap handling, and overflow reporting.
- Arena sizing, growth, reuse, clearing, allocation failure, and integer-overflow rejection.
- Exact score accumulation for single-end and paired-end placement.
- Global maximum traversal and strict comparison behavior.
- Inclusive confidence-interval endpoints.
- Incremental versus legacy LCA equivalence.
- Tree vote-count and numeric-ID tie behavior.
- Safe result-row construction and truncation handling.

### Integration Tests

- Current executable versus phase executable at `-C 1`, comparing all output bytes and candidate fingerprints.
- Single-end and paired-end workflows.
- WFA and Needleman-Wunsch placement paths.
- Text, binary, and zstd Tronko references.
- Single and multiple batches, proving one BWA context per process.
- `-C 1`, `-C 2`, and the production core count.
- Duplicate accession, no-hit, cap-overflow, reverse-complement, and multi-tree fixtures.
- Existing diagnostic modes, including node/site score outputs.
- Sanitizer and Valgrind execution through real alignment and placement paths.

### Production-Like Verification Steps

1. Build the baseline from the recorded `high-perf` commit and the optimized executable with equivalent compiler options.
2. Verify SHA256 hashes of the FWH reference and ASV inputs against the Phase 1 manifest.
3. Run both executables with identical production parameters, environment, CPU allocation, and memory measurement method.
4. First run at `-C 1`; compare candidate fingerprints and all assignment output bytes.
5. Run the repeated-batch workload at the production core count.
6. Compare phase-labeled RSS, peak RSS, runtime, batch throughput, index/map creation counts, candidate counts, and output row counts.
7. Confirm post-warm-up RSS is flat rather than increasing with each batch.
8. Record the optimized absolute peak and percentage change from the Phase 1 baseline. Evaluate Kubernetes capacity separately during deployment validation.

## Performance and Memory Gates

The following gates determine whether the simple algorithmic series is complete:

| Gate | Required result |
|---|---|
| Assignment parity | Byte-identical `-C 1` output across all seven columns |
| BWA parity | Identical retained candidate fingerprints |
| BWA residency | Exactly one index and leaf map per process |
| Placement allocation | Score arena sized to active candidate nodes only |
| Batch stability | No continuing RSS growth after the first workspace high-water mark |
| Memory safety | Clean ASan, UBSan, LSan, Valgrind, and TSAN runs |
| Production-like peak | Lower than the Phase 1 baseline, with absolute and relative change recorded |
| Throughput | No more than 10% regression versus the Phase 1 baseline unless separately approved |

If all gates pass, stop. Do not implement reference packing, mmap, lazy loading, different precision, or parameter tuning.

If semantic gates pass but separate deployment validation still encounters OOM:

1. Capture a post-refactor memory breakdown separating shared Tronko posteriors, the single BWA index, query batches, worker arenas, and allocator overhead.
2. Confirm there is no remaining duplicate process- or worker-owned reference data.
3. Create a new plan choosing between contiguous reference slabs and a versioned per-tree block format with a bounded lazy-load cache.
4. Do not begin that redesign as an unreviewed extension of this plan.

## Rollout and Rollback

- Implement each phase as an independently reviewable commit or PR.
- Preserve the Phase 1 baseline executable and logs until the production-like gate passes.
- Do not deploy a phase that fails output or candidate parity, even if memory improves.
- If a phase fails parity, revert only that phase and use the separate candidate fingerprints to determine whether divergence begins in BWA, candidate retention, scoring, or winner reduction.
- If the shared BWA phase shows a race, restore process isolation temporarily and fix the shared mutable field; do not serialize all alignment behind a global lock as the final design.
- No production resource reductions are part of rollout. Any later memory-tier adjustment requires its own operational decision after sustained measurements.

## Migration Notes

- There is no user-facing data migration or CLI change.
- Existing `.trkb`, text, zstd, and BWA index files remain compatible.
- Existing runtime parameters retain their meanings and defaults.
- New context/workspace APIs are internal and can land without coordinating a database rebuild.
- Diagnostic log schemas should be versioned or extended compatibly so existing memory-log comparison tooling continues to work.

## References

- [`tronko-assign/tronko-assign.c`](../../../tronko-assign/tronko-assign.c)
- [`tronko-assign/allocateMemoryForResults.c`](../../../tronko-assign/allocateMemoryForResults.c)
- [`tronko-assign/global.h`](../../../tronko-assign/global.h)
- [`tronko-assign/placement.c`](../../../tronko-assign/placement.c)
- [`tronko-assign/assignment.c`](../../../tronko-assign/assignment.c)
- [`tronko-assign/bwa_source_files/fastmap.c`](../../../tronko-assign/bwa_source_files/fastmap.c)
- [`tronko-assign/bwa_source_files/bwamem.c`](../../../tronko-assign/bwa_source_files/bwamem.c)
- [`tests/README.md`](../../../tests/README.md)
- [`docs/performance-logging.md`](../../../docs/performance-logging.md)
- [`thoughts/shared/research/2026-01-01-bwa-multithreading-feasibility.md`](../research/2026-01-01-bwa-multithreading-feasibility.md)
- [`thoughts/shared/research/2026-01-01-memory-access-pattern-optimization.md`](../research/2026-01-01-memory-access-pattern-optimization.md)
- [`thoughts/shared/plans/2025-12-30-memory-optimizations.md`](2025-12-30-memory-optimizations.md)
