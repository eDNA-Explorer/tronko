---
date: 2026-03-11T12:00:00-07:00
researcher: Claude
git_commit: c0a65381d5a30c9806bd411b1c835edfeea0208e
branch: optimize-tronko-build
repository: tronko-fork
topic: "Tronko-assign amplicon assignment pipeline: current architecture and potential improvements"
tags: [research, codebase, tronko-assign, assignment, alignment, voting, LCA, accuracy]
status: complete
last_updated: 2026-03-11
last_updated_by: Claude
last_updated_note: "Added deep-dive on BWA alternatives for initial leaf matching (minimap2, VSEARCH, k-mer, MMseqs2, etc.)"
---

# Research: Tronko-assign Amplicon Assignment — Architecture & Potential Improvements

**Date**: 2026-03-11
**Researcher**: Claude
**Git Commit**: c0a6538
**Branch**: optimize-tronko-build
**Repository**: tronko-fork

## Research Question

How does tronko-assign's amplicon assignment pipeline work end-to-end, and what are the potential improvements to accuracy and efficiency?

## Summary

Tronko-assign processes amplicon reads through a 5-stage pipeline: (1) BWA alignment to leaf nodes, (2) sequence alignment to matched leaves via WFA or NW, (3) posterior probability scoring at every tree node, (4) confidence-interval voting, and (5) LCA consensus across trees. The system already has several accuracy-tuning knobs (soft voting, best-leaf override, consensus threshold) but several areas remain ripe for improvement: the multi-tree consensus algorithm, gap/missing-data scoring, early termination completeness, and alignment quality filtering.

---

## Detailed Findings

### 1. Pipeline Architecture (End-to-End)

```
Input: Query reads (FASTA/FASTQ, single or paired-end)
   ↓
Stage 1: BWA Alignment → leaf_coordinates[match][tree, node]
   ↓
Stage 2: WFA/NW Alignment → positions[], locQuery[]
   ↓
Stage 3: Posterior Scoring → nodeScores[match][tree][node]
   ↓
Stage 4: Confidence Voting → voteRoot[tree][node]
   ↓
Stage 5: LCA Consensus → taxonomic assignment
   ↓
Output: TSV with readname, taxonomy path, score, mismatch counts
```

### 2. Stage 1: BWA Leaf Matching

**Files**: `bwa_source_files/fastmap.c`, `tronko-assign.c:176-184, 384-574`

BWA MEM aligns query reads against a concatenated reference FASTA of all leaf node sequences. The SAM output is parsed to identify which reference leaf each read matches, then translated to (tree_id, node_id) coordinates via a hashmap.

**Current behavior**:
- Default `max_bwa_matches = 10` (configurable via `--max-bwa-matches`)
- Deduplication at leaf level in fastmap.c, then at **tree level** in tronko-assign.c (one match per tree)
- Concordant matches (both paired reads map to same leaf) prioritized over discordant
- Tree-level dedup uses O(n²) linear search of `trees_search[]` array

**Potential improvements**:

- **Alignment quality filtering**: BWA matches are accepted without any MAPQ or alignment score threshold. Adding a minimum MAPQ filter would eliminate ambiguous mappings that dilute the voting signal. Currently there is no filtering on alignment quality — every match that BWA returns is used.

- **Tree dedup with hash set**: The current O(n²) linear scan of `trees_search[]` (lines 426-444) becomes slow with many matches. A simple bitmap indexed by tree_id would be O(1).

- **Increase default max_bwa_matches**: The default of 10 can truncate matches for reads in conserved regions. The overflow logging (lines 577-598) shows this happens frequently. Raising to 20-50 with smarter filtering could improve accuracy.

### 3. Stage 2: Sequence Alignment (WFA / Needleman-Wunsch)

**Files**: `placement.c:1-860` (WFA path), `placement.c:1031+` (NW path), `needleman_wunsch.c`

After BWA identifies leaf matches, each query read is fully aligned against the matched leaf's reference sequence using either WFA (default) or NW.

**WFA configuration** (`placement.c:71-78`):
```c
attributes.affine_penalties.mismatch = 4;
attributes.affine_penalties.gap_opening = 6;
attributes.affine_penalties.gap_extension = 2;
attributes.alignment_form.span = alignment_endsfree;
```
- End-free alignment allows partial overlap (appropriate for amplicons)
- Single WFA aligner reused across all matches per thread (efficient)

**Position mapping** (`placement.c:190-217`):
- Alignment result is converted to `positions[]` (reference coordinate per aligned base) and `locQuery[]` (query nucleotide at each position)
- Gaps in reference → position = -1 (penalized as log(0.01) in scoring)
- Gaps in query → kept with reference position (penalized as log(0.25))

**Leaf portion mode** (`-e` flag):
- Instead of aligning to full leaf sequence, uses BWA CIGAR to extract a window around the alignment region
- Configurable padding (`-n` flag)
- Significantly faster for long reference sequences

**Potential improvements**:

- **Configurable alignment penalties**: WFA penalties (mismatch=4, gap_open=6, gap_ext=2) are hardcoded. Different amplicon markers (COI vs 16S vs ITS) have very different substitution profiles. Making these configurable or marker-aware could improve alignment quality.

- **Alignment score-based filtering**: After WFA alignment, there is no quality check on the alignment itself. A severely misaligned read gets scored against the tree just the same as a perfect match. Adding a minimum alignment identity or maximum mismatch fraction filter before scoring would prevent noise from propagating.

- **Mismatch counting only on first match**: Forward/reverse mismatch is only computed for `match==0` (lines 197-199, 596-598). This means the reported mismatch rate reflects only the first BWA match, not the best alignment. Computing it for the best-scoring match would be more informative.

### 4. Stage 3: Posterior Probability Scoring

**Files**: `assignment.c:168-236` (getscore_Arr), `tronko-assign.c:135-157` (store_PPs_Arr)

This is the core innovation of tronko — scoring query reads against **all nodes** (not just leaves) using pre-computed posterior probabilities.

**PP transformation** (`store_PPs_Arr`, `tronko-assign.c:135-157`):
```c
double d = 1 - c;          // c = score_constant (default 0.01)
double e = c / 3;          // uniform error prior
double f = d * PP[pos][nuc];     // correct nucleotide likelihood
double g = e * (1 - PP[pos][nuc]); // error nucleotide likelihood
PP[pos][nuc] = log(f + g);       // log-space for accumulation
```
- Blends phylogenetic posterior with sequencing error model
- Missing data (PP == -1) → set to 1.0 (treated specially)
- NaN/Inf → set to 1.0

**Scoring function** (`getscore_Arr`, `assignment.c:168-236`):
- For each aligned position: look up the log-transformed PP for the query nucleotide at that node
- Gap in query at aligned position: score += log(0.25)
- Unaligned position (pos == -1): score += log(0.01)
- Missing data at reference position: non-gap query nucleotide → log(0.01), gap → +0
- Uses LUT for nucleotide indexing (branchless, cache-friendly)

**Tree traversal** (`assignScores_Arr_paired`, `assignment.c:24-69`):
- Iterative DFS with explicit stack (replaces recursion)
- Scores all nodes: leaves AND internal nodes
- Early termination and pruning optionally skip low-scoring subtrees

**Potential improvements**:

- **Score normalization by alignment length**: Currently the score is a raw sum of log-likelihoods across all aligned positions. Longer alignments get lower (more negative) scores. Normalizing by alignment length would make scores comparable across reads of different lengths, improving voting fairness.

- **Position-specific confidence weighting**: All positions contribute equally to the score. Positions with high posterior probability certainty (PP near 1.0 or 0.0) are more informative than ambiguous positions. A confidence-weighted scoring could upweight discriminative sites.

- **Gap penalty model**: The current gap handling is simplistic — insertion in query gets log(0.01), gap in query at reference position gets log(0.25). A more nuanced model could use position-specific gap probabilities or reduce penalty for gaps in known variable regions.

- **Complete early termination**: The `strikes` tracking mechanism exists but is **not fully implemented** — strikes are reset on improvement but never checked to trigger termination. Completing this feature could give 2-5x speedup for clear-match reads. (This was independently confirmed in `2026-01-01-tier1-benchmark-results.md`.)

### 5. Stage 4: Confidence-Interval Voting

**File**: `placement.c:862-1026`

After scoring every node in every matched tree, the system determines which nodes "vote" for the final assignment.

**Algorithm** (`placement.c:911-925`):
1. Find maximum score across all matches and all nodes
2. All nodes within `±Cinterval` of maximum receive votes
3. Hard voting (default): weight = 1.0 for all voted nodes
4. Soft voting (optional): weight = `exp((score - max) / temperature)`

**Best-leaf override** (`placement.c:931-957`):
- If the best-scoring **leaf** exceeds a threshold AND total votes are below a limit, override all votes with just that leaf
- Disabled by default (`best_leaf_threshold=0, best_leaf_max_votes=0`)

**Potential improvements**:

- **Cinterval auto-tuning**: Cinterval is currently a single user-specified constant (default 5). A too-small Cinterval overfits to the single best node; too-large dilutes the signal. Auto-calibrating based on the score distribution (e.g., using the gap between best and second-best scores) could adapt per-read.

- **Separate leaf vs internal node voting**: Internal nodes always have higher scores than their descendants (more positions match the generic ancestral sequence). This means internal nodes near the root tend to accumulate votes even when a specific leaf is the correct assignment. Weighting leaf votes higher or scoring leaves and internals on separate scales could reduce over-conservative assignments.

- **Score-distribution-aware voting**: Instead of a fixed ±Cinterval window, use a statistical test (e.g., the score gap between ranked nodes) to determine which nodes are "significantly good." This would automatically handle varying degrees of certainty.

- **Enable soft voting by default**: The soft voting mechanism properly weights nodes by likelihood, which is more principled than hard binary voting. The current default (hard voting) can create ties between clearly-different nodes.

### 6. Stage 5: LCA Consensus Across Trees

**File**: `tronko-assign.c:660-787`

When reads match leaves in multiple trees, the system computes a per-tree LCA then finds consensus across trees.

**Single-tree path** (lines 705-721):
- Compute LCA of all voted nodes within one tree using `LCA_of_nodes()`
- Direct taxonomic assignment from the LCA node

**Multi-tree path** (lines 722-787):
1. Compute per-tree LCA for each tree with votes
2. Start at the most specific taxonomic level (`minLevel`)
3. Check if all tree pairs agree at that level
4. If `correctTax >= required_pairs`, assign at that level
5. Otherwise, move up to the next broader level (species → genus → family → ... → domain)

**Consensus threshold** (`--consensus-threshold`, default 1.0):
- 1.0 = unanimity required (all tree pairs must agree)
- < 1.0 = majority/plurality allowed
- `required_pairs = ceil(threshold * total_pairs)`

**Potential improvements**:

- **Weighted tree consensus**: Currently all trees contribute equally to consensus. Trees with higher-scoring matches should have more influence. A score-weighted consensus would prevent a single weak match in an unrelated tree from pushing the assignment to a higher taxonomic level.

- **Per-tree vote count weighting**: A tree with 50 voted nodes (strong signal) is treated the same as a tree with 1 voted node (weak signal). Incorporating the number/weight of votes per tree into the consensus would be more robust.

- **Progressive consensus**: The current algorithm tries each taxonomic level independently. A progressive approach that first identifies which trees are "reliable" (high scoring, many votes) and excludes outlier trees before computing consensus would prevent one bad tree from dragging the assignment up.

- **Break first-tree-wins tiebreaker**: When multiple trees have identical LCA taxonomy, `taxRoot=maxRoots[0]` is always used (line 772). This creates a systematic bias toward lower-numbered trees. Random or score-weighted tiebreaking would be fairer.

- **"NA" string comparison bug**: Line 764 compares taxonomy string to the string literal `"NA"` using `!=` (pointer comparison, not `strcmp`). This may not correctly filter out NA entries depending on string interning.

### 7. Cross-Cutting Improvement Opportunities

**Alignment quality as a confidence signal**:
The mismatch rate from WFA/NW alignment is recorded (`forward_mismatch`, `reverse_mismatch`) but only reported in output — it doesn't influence the assignment decision. Using alignment quality to:
- Weight votes (poorly-aligned reads get lower vote weights)
- Filter before scoring (skip reads with >X% mismatch)
- Adjust Cinterval (widen for poor alignments, tighten for good ones)

**Top-N results output**:
The infrastructure already supports this (nodeScores 3D array stores all scores). Outputting the top-N candidate assignments with scores would help downstream analysis. This was confirmed feasible in prior research (`2026-01-01-tie-handling-top-n-results.md`).

**Chimera detection**:
Chimeric amplicons (composed of sequences from two different organisms) will match leaves from divergent clades and push the LCA to a very high level. Detecting chimeras by checking whether the forward and reverse reads place in different clades, or whether the score profile has a bimodal distribution across the tree, could flag these reads.

**Read-length-dependent scoring**:
Short amplicons have fewer informative sites and tend to match many nodes with similar scores. The Cinterval should arguably be smaller for longer reads (more discriminative power) and larger for shorter reads. An adaptive Cinterval based on alignment length would be more principled.

---

## Code References

- `tronko-assign/tronko-assign.c:135-157` — PP transformation (store_PPs_Arr)
- `tronko-assign/tronko-assign.c:268-900` — Main per-read processing loop (runAssignmentOnChunk_WithBWA)
- `tronko-assign/tronko-assign.c:660-787` — Multi-tree LCA consensus
- `tronko-assign/placement.c:61-860` — WFA alignment and vote accumulation (place_paired)
- `tronko-assign/placement.c:862-1026` — Score aggregation, voting, best-leaf override
- `tronko-assign/placement.c:1031+` — NW alignment path (place_paired_with_nw)
- `tronko-assign/assignment.c:24-69` — Tree traversal with pruning (assignScores_Arr_paired)
- `tronko-assign/assignment.c:168-236` — Core scoring function (getscore_Arr)
- `tronko-assign/global.h:86-94` — Node struct with posterior probabilities
- `tronko-assign/global.h:29` — PP_IDX macro for 1D PP array access
- `tronko-assign/options.c` — All CLI options including accuracy tuning parameters
- `tronko-assign/bwa_source_files/fastmap.c` — BWA integration and match deduplication

## Architecture Documentation

### Key Data Structures

| Structure | Location | Purpose |
|-----------|----------|---------|
| `node` | global.h:86 | Tree node: children, parent, depth, posteriors, taxonomy |
| `bwaMatches` | global.h:197 | BWA match results per read: concordant/discordant hits |
| `resultsStruct` | global.h:116 | Per-thread working memory: scores, votes, alignment buffers |
| `mystruct` | global.h:152 | Thread arguments: all parameters for assignment chunk |
| `scoresStruct` | global.h:213 | Score pair with tree/node coordinates |

### Key Constants

| Constant | Default | Purpose |
|----------|---------|---------|
| `MAX_NUM_BWA_MATCHES` | 10 | Max leaf matches per read |
| `NUMCAT` | 1 | Gamma rate categories for substitution model |
| `SP_SCORE_MIN` | 0.8 | Minimum sum-of-pairs score |
| `MAXQUERYLENGTH` | 30000 | Maximum query read length |

### Tuning Parameters

| Parameter | Default | Effect on Accuracy |
|-----------|---------|-------------------|
| `-c` (Cinterval) | 5 | Wider → more conservative (higher taxa); narrower → more specific |
| `-u` (score_constant) | 0.01 | Controls PP transformation confidence blend |
| `--consensus-threshold` | 1.0 | Lower → allows more disagreement between trees |
| `--soft-voting` | off | When on, weights votes by score proximity |
| `--vote-temperature` | 1.0 | Lower → sharper voting weights (best scores dominate) |
| `--best-leaf-threshold` | 0 | Higher → only very good leaves can override LCA |
| `--best-leaf-max-votes` | 0 | When set, allows confident-leaf shortcutting |
| `--max-bwa-matches` | 10 | More matches → better coverage but slower |

## Historical Context (from prior research)

- `thoughts/shared/research/2026-01-01-tier1-benchmark-results.md` — Confirmed early termination is broken (0% accuracy); compiler flags NATIVE_ARCH and FAST_MATH are safe
- `thoughts/shared/research/2026-01-01-tie-handling-top-n-results.md` — Top-N results feasible with <1% overhead; data structures support it
- `thoughts/shared/research/2026-01-03-non-determinism-mitigation-strategies.md` — 3% variance from BWA threading is unavoidable without single-threaded mode; match sorting actually increases variance
- `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md` — Phased optimization plan: Tier 1 (3-15x), cumulative potential 100-500x
- `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md` — Algorithm details for branch-and-bound early termination
- Root-level `ABLATION.md` — Ablation system for database pruning studies
- Root-level `confidence_analysis.py` / `confidence_analysis2.py` — Tested 16+ candidate confidence signals for shallow assignments

## Prioritized Improvement Recommendations

### High Impact, Low Effort
1. **Complete early termination** — strikes tracking exists but isn't checked. Finish it for 2-5x speedup.
2. **Alignment quality filtering** — Add MAPQ or mismatch threshold before scoring. Prevents noisy matches from diluting votes.
3. **Fix "NA" string comparison** — Line 764 uses pointer comparison instead of strcmp.
4. **Enable soft voting by default** — More principled than hard binary voting.

### High Impact, Medium Effort
5. **Score normalization by alignment length** — Divide total score by number of aligned positions. Makes cross-read comparisons fair.
6. **Weighted tree consensus** — Use per-tree best scores to weight consensus, not just binary agreement.
7. **Cinterval auto-tuning** — Adapt per-read based on score gap distribution.
8. **Separate leaf vs internal node vote weighting** — Reduce over-conservative internal node dominance.

### Medium Impact, Higher Effort
9. **Configurable WFA alignment penalties** — Different markers have different substitution profiles.
10. **Top-N results output** — Infrastructure exists; surface multiple candidate assignments.
11. **Chimera detection** — Flag reads with bimodal clade placement.
12. **Progressive tree consensus** — Exclude outlier trees before computing consensus.

## Open Questions

1. What is the empirical accuracy impact of soft voting vs hard voting across different marker genes?
2. What is the optimal Cinterval for different amplicon lengths and marker genes?
3. How much accuracy is lost from the max_bwa_matches=10 cap on reads in conserved regions?
4. Would alignment-quality-weighted voting improve species-level accuracy without degrading higher-level assignments?
5. Is the score_constant (0.01) optimal, or should it vary by marker gene / tree size?

---

## Follow-up Research: Deep-Dive on 5 Key Improvement Areas

### A. Alignment Quality Filtering — Where MAPQ Gets Discarded

**BWA outputs MAPQ in SAM, but tronko-assign throws it away.**

In `bwa_source_files/fastmap.c:179`, the SAM parsing format string is:
```c
sscanf(token, "%s %d %s %d %*d %s %s %*d %*d %*s %*s %*s %*s %*s %*s %*[^.,;]",
       readname, &decimal, read1, &start_position, cigar, read2);
```

The 5th SAM field (MAPQ) is consumed by `%*d` — the `*` means "read and discard." BWA computes MAPQ in `bwamem.c:1036-1046` based on alignment score difference, seed coverage, and identity, but this information never reaches the bwaMatches struct.

**The bwaMatches struct (`global.h:197-211`) has no quality fields at all** — only tree/node IDs, CIGARs, and start positions. There is no MAPQ, no alignment score, no identity percentage.

**After WFA alignment, the alignment score is also discarded.** In `placement.c:109`:
```c
wavefront_align(wf_aligner, leaf_sequence, leaf_length, query_1, query_length);
```
The `wf_aligner` struct contains the alignment score, but the code only extracts the CIGAR via `wf_aligner->cigar`. The score is never read.

**Implementation path for MAPQ filtering:**

The cleanest approach is to filter in `fastmap.c` during SAM parsing:
1. Change `%*d` to `%d` and add an `int mapq` variable
2. Add a MAPQ threshold parameter (passed through `aux`)
3. `if (mapq < threshold) continue;` — skip low-quality matches before they enter bwaMatches

For WFA score filtering, capture `wf_aligner->score` after the `wavefront_align()` call in `placement.c:109` and skip matches where the alignment score indicates poor fit (e.g., more than X% gaps or mismatches).

### B. Internal Node vs Leaf Node Scoring Bias — The Real Mechanism

**The bias is NOT in getscore_Arr arithmetic. It's in tree geometry × Cinterval voting.**

**How posteriors differ between leaf and internal nodes** (`tronko-build/likelihood.c:1267-1327`):

For **leaf nodes** (line 1325-1327):
```c
if (i == b) posteriornc[s][i] = 1.0;  // Observed nucleotide
else        posteriornc[s][i] = 0.0;  // All others
```
Leaves are deterministic: PP=1.0 for the observed base, PP=0.0 for everything else. A mismatch at a leaf position contributes `log(c/3 × 1.0)` ≈ log(0.0033) ≈ **-5.7** per position.

For **internal nodes** (line 1290-1297):
```c
posteriornc[s][i] = likenc[s][i] * posteriornc[s][i] * pi[i];
posteriornc[s][i] = posteriornc[s][i] / sum;  // Normalize
```
Internal PPs are the product of conditional likelihood × parent posterior × nucleotide prior, normalized. For a clade where 80% of descendants have A and 20% have G, the internal PP for A might be ~0.75, G might be ~0.20, others ~0.025 each. A mismatch at an internal node only costs log(0.99×0.20 + 0.0033×0.80) ≈ log(0.20) ≈ **-1.6** per position.

**The scoring math (worked example, 100-bp read, score_constant=0.01):**

Scenario: Query matches species X perfectly. Species X is in a clade with species Y (95% similar).

| Node | Matching positions (95) | Mismatching positions (5) | Total Score |
|------|------------------------|--------------------------|-------------|
| Leaf X (correct) | 95 × log(0.99×1.0 + 0.0033×0.0) = 95 × (-0.010) = **-0.95** | 5 × log(0.0033) = 5 × (-5.71) = **-28.6** | **-29.5** |
| Leaf Y (wrong) | ~90 × (-0.010) = **-0.90** | ~10 × (-5.71) = **-57.1** | **-58.0** |
| Internal (parent of X,Y) | 95 × log(0.99×0.88 + 0.0033×0.12) ≈ 95 × (-0.13) = **-12.4** | 5 × log(0.99×0.12 + 0.0033×0.88) ≈ 5 × (-2.09) = **-10.5** | **-22.9** |

**The internal node scores -22.9 vs the correct leaf at -29.5.** The internal node wins by ~6.6 log-likelihood units because it doesn't pay the full -5.7 penalty per mismatch position.

With Cinterval=5, the voting window is [-22.9 - 5, -22.9 + 5] = [-27.9, -17.9]. The correct leaf (-29.5) **falls outside this window** and gets no vote! Only the internal node and possibly other internal ancestors vote.

**This is the fundamental problem:** When the best-scoring node is internal, the Cinterval window can exclude the correct leaf entirely, and the LCA of voting nodes climbs to genus/family level.

**Depth information is available for weighting** (`global.h:90`):
```c
typedef struct node {
    ...
    int depth;  // Computed during init via assignDepthArr()
    ...
} node;
```
Depth is set in `tronko-assign.c:158-163` and also loaded from binary reference files. It is **never used during voting** (placement.c:911-925).

**Concrete mitigation approaches:**

1. **Depth-weighted voting**: Multiply vote by `depth_weight = pow(depth / max_depth, alpha)` where alpha > 0 penalizes shallow nodes. Leaves (high depth) get full weight; root (depth 0) gets near-zero weight.

2. **Leaf-only maximum**: Find the best-scoring **leaf** instead of the best-scoring node overall, then set the Cinterval window around that leaf score. This ensures the voting window is centered on leaf-level scores, not inflated by internal node scores.

3. **Split scoring tracks**: Score leaves and internal nodes separately. Use leaf scores for voting/LCA, use internal scores only for confidence estimation. This is the most principled approach but requires restructuring placement.c.

4. **Internal node score penalty**: Apply a depth-based correction: `adjusted_score = raw_score - penalty × (max_depth - depth)`. This compensates for the inherent advantage of averaging.

### C. Score-Distribution-Aware Voting — What the Data Actually Shows

**Existing confidence analysis results dampen expectations for distribution-based approaches.**

The `confidence_analysis.py` and `confidence_analysis2.py` scripts tested 16+ candidate confidence signals on real data:

| Signal | Predictive Power | Notes |
|--------|-----------------|-------|
| best_leaf_score | **High** | Best single predictor of shallow correctness |
| total_votes | **High** | Low votes + high leaf score → reliable |
| gap_12 (1st-2nd leaf gap) | **Medium-High** | Measures leaf separation quality |
| voted_leaves | **High** | Usually 1-5 for shallow reads |
| voted_internal | **Low** | Mostly noise |
| max_score (overall) | **Low-Medium** | Only useful in extremes |

**Key finding from confidence_analysis2.py**: Cinterval-relative scaling (e.g., rules parameterized as `bl > -2×Cinterval`) did NOT improve over absolute thresholds. This suggests the score distribution shape doesn't vary enough across reads to benefit from adaptive scaling.

**What the trace output already reveals** (placement.c:961-1022):
- `voted_leaves` vs `voted_internal` — already tracked per read
- `best_leaf_score` vs `second_leaf_score` — already tracked with gap
- `maximum - best_leaf_score` — shows whether best node is a leaf

**What's NOT tracked but could be:**
- Score variance/stdev within the voting window (no existing code computes this)
- Score histogram shape (unimodal vs bimodal)
- Phylogenetic depth distribution of voted nodes

**Practical recommendation:**

Rather than building a full distribution-aware voting system, the most impactful changes are:

1. **Center the Cinterval window on best-leaf score, not best-node score** — this single change addresses both the internal-node bias and the score-distribution problem simultaneously, without needing to analyze the full distribution.

2. **Use the best-leaf-to-second-leaf gap as a confidence indicator** — reads where gap_12 is large (clear winner) should get species-level assignment; reads where gap_12 is small (ambiguous) should get genus-level. This is already computable from the trace output signals.

3. **Keep soft voting** — it already provides distribution-aware weighting via the temperature parameter. Tuning `vote_temperature` to be proportional to the read's score range would give adaptive behavior cheaply.

### D. Mismatch Counting — Only Match==0 Gets Tracked

**Current behavior** (`placement.c:197-199` for forward, `596-598` for reverse):
```c
if (pattern_alg[i] != text_alg[i] && match==0){
    forward_mismatch++;
}
```

Only the first BWA match (match==0) has its alignment quality tracked. For all subsequent matches (match=1,2,...), mismatch counting is skipped entirely. The reported mismatch values in the output (`minimum_score[1]` and `minimum_score[2]`) reflect only this first match.

**Problem**: Match==0 is the first concordant BWA hit, which may not be the best phylogenetic match. The best-scoring leaf (as determined by posterior probability scoring in stage 3) could be match==3 or match==7. The output mismatch rate is for a potentially irrelevant alignment.

**Impact**: Anyone using the mismatch columns for downstream quality filtering is filtering based on the wrong alignment's quality.

**Fix**: Compute mismatch for all matches and report the one corresponding to the best-scoring leaf (which is known after the scoring loop). Or at minimum, compute it for the match corresponding to `match_number` (the best match, found at placement.c:878-879).

### E. Combined Improvement: Leaf-Centered Voting

The most impactful single change would combine findings from sections B, C, and D into a "leaf-centered voting" approach:

**Current flow:**
```
1. Score ALL nodes (leaves + internal)
2. Find global maximum (often an internal node)
3. Vote within ±Cinterval of global max
4. Compute LCA of voted nodes → assignment
```

**Proposed leaf-centered flow:**
```
1. Score ALL nodes (leaves + internal) — unchanged
2. Find best-scoring LEAF (not global max)
3. Vote within ±Cinterval of best LEAF score
4. Optional: apply depth weight to internal node votes
5. Compute LCA of voted nodes → assignment
```

**Why this helps:**
- The Cinterval window is now centered where it matters — at leaf-level scores
- Internal nodes that scored higher than any leaf no longer set the window; they can still vote if they fall within the leaf-centered window, but they don't dominate it
- The correct leaf (and its close relatives) are much more likely to fall within the window
- Species-level assignment becomes possible even when internal nodes would have been the global maximum

**Code change magnitude**: ~15 lines in `placement.c:862-925` — change the max-finding loop to only consider leaf nodes (check `treeArr[t][k].up[0] == -1`), then keep the voting logic identical.

**Risk**: This could reduce assignment at higher taxonomic levels for genuinely ambiguous reads. The Cinterval parameter would need retuning (probably slightly wider to compensate).

**Testing approach**: Run the existing test datasets with `--trace-read '*'` to capture before/after vote distributions, then compare species/genus/family accuracy.

---

## Deep-Dive: BWA Alternatives for Initial Leaf Matching

### What BWA Is Actually Doing Here

BWA-MEM is **embedded as source code** (18 .c files in `bwa_source_files/`) and compiled directly into the tronko-assign binary. It is *not* called as an external process. The `main_mem()` function is called directly from `tronko-assign.c:180-182`.

The task BWA performs is narrow: given a query amplicon read (~150bp), find which **leaf reference sequences** (~250bp amplicon barcodes, typically 1,000–100,000+ sequences) it aligns to, and return the top matches with CIGAR strings. This is then used to identify candidate trees/nodes for full phylogenetic scoring.

This is fundamentally a **short-read-against-short-reference similarity search**, not genome mapping. BWA-MEM was designed for the latter — mapping reads against large genomes (Gb scale) using a BWT index optimized for long, contiguous references. Using it to search a collection of short amplicon sequences is valid but arguably a mismatch between tool design and use case.

### Why Consider Alternatives

1. **Sensitivity at high divergence**: BWA-MEM's seeding strategy (minimum seed length 17bp by default) can miss matches when the query diverges >10% from any reference. For novel species or fast-evolving markers, BWA may fail to find the correct candidate leaf, and no downstream scoring can fix a miss at this stage.

2. **Multi-mapping handling**: BWA-MEM reports secondary alignments, but its MAPQ calculation is tuned for genome mapping where one location is "correct." For amplicon matching where a read legitimately aligns to multiple similar leaves, BWA's duplicate suppression and MAPQ heuristics may discard informative secondary hits.

3. **Overhead for the task**: BWA-MEM's full Smith-Waterman extension, paired-end insert size modeling, and chimeric read detection are unnecessary here. The reads and references are similar length — we don't need splice-aware or large-gap-tolerant alignment. A simpler tool could be faster.

4. **Index size**: BWA's FM-index (BWT + suffix array) is space-efficient for genomes but overkill for ~100K short sequences totaling ~25MB of sequence. Simpler index structures could load faster and use less memory.

5. **CIGAR is partially redundant**: BWA produces a CIGAR string, but tronko-assign then re-aligns the read against the leaf sequence using WFA or Needleman-Wunsch (`placement.c`). The BWA CIGAR is only used for the `use_portion` optimization (extracting a window around the alignment). The actual alignment used for scoring is independent of BWA's alignment.

### Candidate Alternatives

#### 1. minimap2

**What it is**: General-purpose aligner by Heng Li (same author as BWA). Uses minimizer-based indexing instead of BWT.

**Relevant mode**: `-x sr` (short reads) or `-x map-ont` with tuned parameters. For amplicon-vs-amplicon, the `--eqx` flag gives detailed CIGAR with =/X operators.

**Pros**:
- 3-5x faster than BWA-MEM for short reads
- Better sensitivity for divergent matches — minimizer seeding tolerates higher error rates than BWA's exact-match seeding
- Reports all secondary alignments naturally (`-N` flag controls how many)
- Produces PAF format (simpler than SAM) or SAM — either works
- Lower memory footprint for small reference sets
- Can be embedded as a library (`libminimap2`) — has a clean C API (`mm_idx_build`, `mm_map`)

**Cons**:
- Designed for long reads primarily; short-read mode is functional but less battle-tested than BWA-MEM for ~150bp reads
- Minimizer seeding has a minimum window size that may cause missed seeds on very short overlaps
- Still designed for mapping reads to references, not sequence similarity search per se

**Integration effort**: Medium. minimap2 has a library API (`minimap2.h`) that could replace the embedded BWA source. The index build, query, and result extraction would need rewriting, but the interface is cleaner than BWA's. Index build: `mm_idx_build()`. Mapping: `mm_map()`. Result: array of `mm_reg1_t` structs with coordinates, CIGAR, MAPQ, alignment score.

**Best for**: Drop-in replacement that improves speed and sensitivity without major architectural changes.

#### 2. VSEARCH

**What it is**: Open-source reimplementation of USEARCH, specifically designed for amplicon/metabarcoding analysis.

**Relevant mode**: `--usearch_global` (global alignment search) or `--search_exact` for exact matching.

**Pros**:
- **Purpose-built for this exact use case** — amplicon reads against amplicon reference databases
- Global alignment (not local like BWA-MEM), which is more appropriate when query and reference are similar length
- Identity-based filtering (`--id 0.8`) directly answers "which references are >80% similar"
- Reports all hits above threshold (not just best), with `--maxaccepts` and `--maxrejects` controlling thoroughness
- Handles reverse-complement automatically
- Widely used in amplicon bioinformatics (QIIME2, mothur pipelines use it)
- Well-maintained, BSD-licensed, pure C

**Cons**:
- External process — no library API. Would need to be called via `system()` or `popen()`, adding I/O overhead (write queries to temp FASTA, read results back). Could pipe via stdin/stdout to reduce disk I/O.
- Heuristic search: uses k-mer pre-filtering then full alignment. For very large databases (>100K sequences), may be slower than indexed approaches like BWA
- No built-in paired-end handling — would need to process R1/R2 separately
- Output formats (BLAST6, SAM, UC) would require different parsing than current BWA SAM

**Integration effort**: High. Requires either embedding VSEARCH source (large codebase, ~50K lines) or calling as external process. The paired-end handling would need to be reimplemented. However, the output would be more directly useful — global identity percentage is exactly what the downstream scoring needs to know about.

**Best for**: If the goal is to align the initial matching step with what the amplicon community uses and expects. Natural fit for the data type.

#### 3. MMseqs2

**What it is**: Ultra-fast sequence search and clustering tool. Orders of magnitude faster than BLAST for large-scale searches.

**Relevant mode**: `mmseqs search` or `mmseqs map` for prefiltering + alignment.

**Pros**:
- Extremely fast for large reference databases (>100K sequences) — k-mer prefiltering + vectorized Smith-Waterman
- Returns alignment scores, identity, E-values — rich quality metrics
- Has a library API (`libmmseqs`) that can be embedded
- Handles both nucleotide and protein searches
- Can report all hits above a threshold, not just top N

**Cons**:
- Designed for protein-scale searches — nucleotide mode works but is less optimized
- Heavy dependency (C++ with SIMD intrinsics, OpenMP) — harder to embed than pure C tools
- Index format is different and more complex
- Overkill for databases of ~1,000-10,000 sequences — the speedup over BWA only materializes at much larger scales
- No paired-end support

**Integration effort**: High. C++ library with complex build system. Would need careful integration with the existing C codebase.

**Best for**: Very large reference databases (millions of sequences). Not ideal for typical tronko databases.

#### 4. Bowtie2

**What it is**: Short-read aligner using FM-index (like BWA) but with different seeding strategy.

**Relevant mode**: `--very-sensitive` mode with `--all` (report all alignments).

**Pros**:
- Excellent sensitivity for short reads with low divergence (<10%)
- Well-understood MAPQ semantics
- Can report all valid alignments (`-a` flag), unlike BWA which caps secondaries
- Mature, stable, widely benchmarked

**Cons**:
- External process only — no embeddable library API
- FM-index has same overhead issues as BWA for small reference sets
- Less sensitive than BWA-MEM at higher divergence (>10%)
- C++ codebase, harder to embed than C

**Integration effort**: High (external process). Not recommended as a BWA replacement here since it doesn't offer clear advantages for this use case.

#### 5. DIAMOND (nucleotide mode via blastx/blastn emulation)

Not appropriate — designed for protein search. Including for completeness.

#### 6. K-mer based approaches (Kraken2-style, sourmash)

**What it is**: Instead of alignment, use exact k-mer matching to identify candidate leaves.

**Concept**: Build a k-mer index of all leaf sequences. For each query read, extract k-mers and find which leaves share the most k-mers. Return the top N leaves as candidates.

**Pros**:
- Extremely fast — O(L) per query where L is read length
- No alignment needed at this stage (alignment happens downstream in WFA/NW anyway)
- Simple to implement from scratch (~200 lines of C for a hash-based k-mer index)
- Natural multi-match: returns all leaves above a k-mer overlap threshold
- No seeding sensitivity issues — if k-mers match, they're found

**Cons**:
- No positional information — can't produce CIGAR or start position. The `use_portion` optimization in placement.c would need an alternative approach (or just align to full reference, accepting slight overhead)
- K-mer size selection matters: k=15 works for ~97% identity, k=11 works for ~90% identity. May need adaptive k based on expected divergence
- Hash collisions at low k values could produce false positive candidates
- No MAPQ equivalent, though a "k-mer coverage fraction" metric serves a similar purpose

**Integration effort**: Low-medium. Could be implemented directly in C as a replacement module without external dependencies. The hash table infrastructure already exists in tronko-assign's codebase (the leaf name hashmap in fastmap.c).

**Best for**: Maximum speed with simplest possible implementation. Makes the most sense if BWA's alignment output is truly redundant (given that WFA re-aligns everything anyway).

#### 7. Parasail / edlib / KSW2 (library-level aligners)

**What it is**: Lightweight SIMD-accelerated alignment libraries. Not full mapping tools — they align a query against a single reference sequence.

**Concept**: Replace BWA's index-based search with a two-stage approach: (a) k-mer prefilter to find candidate leaves, (b) Parasail/edlib for fast pairwise alignment to compute identity and CIGAR.

**Pros**:
- Extremely fast pairwise alignment (edlib: O(n*m) with bit-parallel, ~10x faster than NW)
- edlib is MIT-licensed, single-header C library, trivial to embed
- Parasail has SIMD-accelerated NW/SW, reports score + CIGAR
- Full control over scoring parameters

**Cons**:
- Need a prefilter step (k-mer index) to avoid O(N) all-vs-all comparisons per read
- Two components to maintain (prefilter + aligner) instead of one (BWA)

**Integration effort**: Low for edlib (single .c/.h file). Medium for the prefilter.

### Comparison Matrix

| Tool | Speed vs BWA | Sensitivity | Paired-end | Embeddable | Integration effort | Best database size |
|------|-------------|-------------|------------|------------|-------------------|-------------------|
| BWA-MEM (current) | baseline | good <5% div | yes | yes (embedded) | — | any |
| minimap2 | 3-5x faster | better >5% div | yes | yes (library API) | medium | any |
| VSEARCH | ~1x | excellent (global) | no | no (external) | high | <500K seqs |
| MMseqs2 | 10-100x for large DB | excellent | no | difficult (C++) | high | >100K seqs |
| Bowtie2 | ~1x | good <10% div | yes | no (external) | high | any |
| k-mer index | 10-50x faster | good, k-dependent | trivial | yes (custom) | low-medium | any |
| edlib + k-mer | 5-20x faster | excellent | trivial | yes (single header) | low-medium | <100K seqs |

### Recommendation

**Short term (minimal disruption)**: **minimap2** is the strongest drop-in replacement. It has a C library API, handles paired-end reads, is faster than BWA-MEM, and has better sensitivity at higher divergence levels. The integration pattern would mirror the current BWA embedding: include minimap2 source, call `mm_idx_build()` to index leaf sequences, call `mm_map()` per read, extract hits from `mm_reg1_t` results. The SAM-parsing infrastructure in `fastmap.c` could be largely preserved since minimap2 can output SAM.

**Medium term (architectural improvement)**: **k-mer prefilter + WFA** eliminates the redundant alignment step entirely. Currently, BWA aligns the read to find candidates, then WFA re-aligns the read for scoring. With a k-mer prefilter, the first step becomes a simple hash lookup (which candidates?), and WFA does the only alignment. This removes ~40% of computation (BWA's Smith-Waterman extension) for the same result. The `use_portion` optimization would need to be replaced with "align to full reference" (slight cost increase per match, but far fewer total alignment operations if the k-mer filter is selective).

**Long term (if scaling to very large databases)**: **MMseqs2** or a custom minimizer index. At >1M reference sequences, BWA's index load time and memory become bottlenecks. But typical amplicon databases are 1K-100K sequences, so this is unlikely to matter soon.

### Integration Considerations

Regardless of which alternative is chosen, the integration needs to preserve:

1. **The `bwaMatches` output struct** — downstream code expects `(tree_id, node_id)` coordinate pairs. The leaf-name-to-tree-coordinate hashmap translation must be preserved or reimplemented.

2. **Paired-end concordance checking** — the current code separates concordant (both reads map to same tree) from discordant matches. Any alternative needs this concept, even if the tool itself doesn't natively support paired ends (can be done post-hoc by checking if R1 and R2 candidate sets overlap).

3. **CIGAR strings** (only if `use_portion` is active) — the start position and CIGAR from the initial alignment are used to extract a window from the reference sequence before WFA re-alignment. If the alternative doesn't produce CIGARs (e.g., k-mer approach), the `use_portion` path would need to be removed or replaced with a length-based heuristic.

4. **The max_bwa_matches cap** — the pipeline has a configurable limit on matches per read. The alternative should support returning the top N matches or all matches above a threshold.

5. **Thread safety** — BWA uses `kthread` for parallel processing. The alternative must be thread-safe or have its own parallelism strategy that integrates with tronko-assign's existing threading model.
