# Leaf-Centered Voting & Alignment Quality Filtering

## Overview

Two improvements to tronko-assign's amplicon assignment pipeline that address the same root problem: **over-conservative taxonomic assignments caused by internal tree nodes dominating the voting system.**

The first (leaf-centered voting) changes *which* nodes set the reference point for voting. The second (alignment quality filtering) prevents low-confidence matches from entering the pipeline at all. They are independent and can be implemented separately, but they reinforce each other.

---

## Current State: How Scoring and Voting Work Today

### Notation

Let's define the terms precisely:

- **Tree T** has N leaf nodes and N-1 internal nodes (2N-1 total).
- **Node n** has a posterior probability vector **PP(n, s, b)** at each alignment site s for each nucleotide b in {A, C, G, T}.
- **Query Q** is a read of length L, aligned to the tree's reference MSA, producing aligned positions **pos(i)** and query bases **q(i)** for i = 1..L.
- **Cinterval C** is a user-set constant (default 5) controlling the voting window.
- **score_constant c** (default 0.01) controls the error model in the PP transformation.

### Step 1: PP Transformation (store_PPs_Arr)

Before any reads are processed, all stored posteriors are converted from probability space to log-likelihood space with an error model:

```
PP'(n, s, b) = log( (1-c) * PP(n, s, b) + (c/3) * (1 - PP(n, s, b)) )
```

This blends the phylogenetic posterior with a uniform error rate. The transformed value is always negative (since the argument to log is < 1).

**Key insight**: This transformation is monotone. Higher PP yields a less-negative (closer to 0) transformed score. The ordering of nodes by PP at any site is preserved after transformation.

### Step 2: Per-Node Scoring (getscore_Arr)

For each node n in each matched tree, the score is:

```
S(n) = sum over i=1..L of:
  - PP'(n, pos(i), q(i))   if q(i) is A/C/G/T and pos(i) has data
  - log(0.25)               if q(i) is a gap character
  - log(0.01)               if pos(i) is missing data or unaligned
```

This is a sum of L log-likelihoods. **All values are negative**, so S(n) is a large negative number. A *better* match means a *less negative* (closer to zero) score.

### Step 3: Find Global Maximum

```
S_max = max over all nodes n in all matched trees of S(n)
```

Currently this considers **all** nodes — leaves and internal.

### Step 4: Voting

Every node within the Cinterval window of S_max receives a vote:

```
For each node n:
  if S_max - C <= S(n) <= S_max + C:
    vote(n) = 1    (hard voting, default)
```

### Step 5: LCA of Voted Nodes

The LCA (Lowest Common Ancestor) of all voted nodes determines the taxonomic assignment. More voted nodes spread across the tree push the LCA toward the root, yielding higher-level (less specific) taxonomy.

---

## The Problem: Why Internal Nodes Outscore Leaves

### Posterior Probability Differences

In tronko-build (`likelihood.c:1285-1298`), posteriors are computed differently for leaves vs internal nodes:

**Leaves** (observed data, line 1325-1327):
```
PP(leaf, s, b) = 1.0  if b == observed_base
PP(leaf, s, b) = 0.0  otherwise
```

**Internal nodes** (Bayesian reconstruction, line 1290-1291):
```
PP(internal, s, b) = L(n, s, b) * PP(parent, s, b) * pi(b) / Z
```
where L is the conditional likelihood (Felsenstein pruning), pi is the stationary frequency, and Z normalizes.

The critical consequence: **at any site where the query mismatches the leaf's observed base, the leaf gets PP = 0.0 for the correct query nucleotide**. After transformation:

```
PP'(leaf, s, wrong_base) = log(c/3) = log(0.0033) ≈ -5.7
```

But the internal node's PP for that same "wrong" base might be 0.15 (15% of descendants had it), giving:

```
PP'(internal, s, wrong_base) = log(0.99 * 0.15 + 0.0033 * 0.85) ≈ log(0.151) ≈ -1.89
```

**Per-mismatch cost: leaf pays -5.7, internal node pays only -1.9.** The internal node saves ~3.8 log-units per mismatch site.

### Worked Example

Consider a 150-bp read that matches species X with 5 mismatches (96.7% identity). The parent internal node of X represents a genus with 10 species.

| | Leaf X | Internal (genus) |
|---|---|---|
| 145 matching sites | 145 * (-0.010) = **-1.45** | 145 * (-0.14) = **-20.3** |
| 5 mismatching sites | 5 * (-5.71) = **-28.55** | 5 * (-1.89) = **-9.45** |
| **Total score** | **-30.0** | **-29.8** |

The internal node wins by 0.2 log-units. With Cinterval = 5:

- Window: [-29.8 - 5, -29.8 + 5] = **[-34.8, -24.8]**
- Leaf X at -30.0: **inside the window** (gets a vote)
- Internal node at -29.8: **inside the window** (gets a vote)
- All other internal nodes on the path to root: many also inside the window

The LCA of {leaf X, genus node, family node, ...} is pushed up to wherever the deepest voted ancestor sits. If the family-level node also falls within the window, the assignment becomes **family-level instead of species-level**, even though the read clearly belongs to species X.

### When it Gets Worse

The problem intensifies with:
- **More mismatches** — each one costs the leaf ~3.8 more than its parent
- **More species in the clade** — internal PPs become more "smeared," creating a wider range of internal nodes near the top
- **Larger Cinterval** — wider windows admit more internal ancestors
- **Short reads** — fewer sites mean each mismatch has proportionally larger impact

---

## Proposed Improvement 1: Leaf-Centered Voting

### The Idea

Instead of centering the voting window on the global maximum (which may be an internal node), **center it on the best-scoring leaf**.

### Algorithm

Define:

```
S_leaf_max = max over leaf nodes only of S(n)
```

Then the voting window becomes:

```
For each node n:
  if S_leaf_max - C <= S(n) <= S_leaf_max + C:
    vote(n) = ...
```

### Why This Works

Using our example: S_leaf_max = -30.0 (leaf X). New window: [-35.0, -25.0].

The internal genus node at -29.8 still falls inside (that's fine — it's the LCA of the correct species, so it should vote). But the key change is that **internal nodes that scored higher than all leaves no longer set the window**. The window is anchored to leaf-scale scores, so only nodes that are competitive with actual species-level matches get votes.

When the best-scoring node overall was internal and significantly above the best leaf, the old window was shifted upward, admitting many more shallow internal nodes. The new window stays grounded.

### The Three Scenarios That Matter

Leaf-centered voting needs to do the right thing in three distinct cases. Let's trace each one through the math.

#### Scenario A: Known species in DB (read IS species X)

The query matches leaf X at ~97% identity (5 mismatches per 150bp). As shown in the worked example above, leaf X scores -30.0 and the genus-level internal node scores -29.8.

**With leaf-centered voting**: Window centered on -30.0 (leaf X). The genus node at -29.8 is inside the window and votes. Sibling leaf Y (same genus, ~90% identity, scoring ~-58) is far outside and doesn't vote. The LCA of {leaf X, genus node} is the genus node, but since leaf X is the *only* leaf voting, the LCA lands at leaf X's species-level taxonomy.

Actually, the LCA function computes the LCA of all voted nodes. If only leaf X and its direct ancestors are voting, the LCA is leaf X itself (since it's the deepest voted node and its ancestors are all on the path to root). **Result: species-level assignment. Correct.**

#### Scenario B: Novel species, genus is in DB

This is the critical case you raised. The read comes from a species not in the database, but its genus contains 3 known species (leaves A, B, C) in the tree. The novel species is, say, 92% similar to A, 90% to B, 88% to C.

Let's trace the scores for a 150bp read:

**Leaf A** (closest known relative, 12 mismatches):
```
138 matches: 138 * (-0.010) = -1.38
 12 mismatches: 12 * (-5.71) = -68.5
Total: -69.9
```

**Leaf B** (next closest, 15 mismatches):
```
135 matches: 135 * (-0.010) = -1.35
 15 mismatches: 15 * (-5.71) = -85.7
Total: -87.0
```

**Genus-level internal node** (MRCA of A, B, C; PP smoothed across all three):
At matching positions the genus node has PP ~0.80 (not 1.0 because descendants vary).
At mismatching positions it has PP ~0.25 for each base (evenly distributed).
```
138 sites, query matches dominant base (PP~0.80):
  138 * log(0.99*0.80 + 0.0033*0.20) = 138 * log(0.793) = 138 * (-0.232) = -32.0

12 sites where query mismatches leaf A but might partially match genus:
  At these sites, some genus members have the query base. Say PP~0.20 for query base.
  12 * log(0.99*0.20 + 0.0033*0.80) = 12 * log(0.201) = 12 * (-1.60) = -19.2

Total genus node: ≈ -51.2
```

So:
- Leaf A: **-69.9** (best leaf)
- Leaf B: **-87.0** (second-best leaf)
- Genus node: **-51.2** (internal, scores much higher)

**delta_12 = -69.9 - (-87.0) = 17.1** — this is a large gap, but both leaves score *poorly* in absolute terms. The genus node beats both leaves by ~19 units.

**With leaf-centered voting**: Window centered on -69.9 (leaf A), width ±C (say ±5): [-74.9, -64.9].
- Leaf A at -69.9: **inside** (votes)
- Leaf B at -87.0: **outside** (doesn't vote)
- Genus node at -51.2: **outside** (above the window!)

**Problem**: Only leaf A votes, so the LCA is leaf A, giving a **species-level assignment to species A — which is wrong** (the read is from a novel species). The correct answer is genus-level.

This is the failure mode: leaf-centered voting *blindly* narrows to the nearest leaf even when that leaf isn't a good match. The genus node — which represents the right level of confidence — gets excluded because it scores above the window.

#### How to Fix Scenario B: Use delta_12 and Absolute Score Together

The key observation: **in scenario B, delta_12 is large (17.1), but the best leaf score is poor (-69.9).** In scenario A, delta_12 is also large, but the best leaf score is good (-30.0).

The distinction between "read IS species X" and "read is a novel species near X" shows up in the **absolute quality of the best leaf match**, not just the gap.

We need a rule like:

```
If best_leaf_score is good (close to 0, meaning few mismatches):
    → Species-level confidence. Use leaf-centered voting.
If best_leaf_score is poor AND delta_12 is large:
    → One leaf is closest but still not a great match.
      The read is probably a novel species in that leaf's genus.
      → Use standard voting (global max centered),
        which naturally lands on the genus-level internal node.
If best_leaf_score is poor AND delta_12 is small:
    → Multiple leaves match equally poorly.
      → Standard voting, lands on their shared higher-level ancestor.
```

More precisely, define a **leaf quality threshold** T_leaf. This can be derived from read length and expected divergence. For a 150bp read at 97% identity (typical intraspecific), the best leaf score is approximately:

```
T_leaf ≈ L * log(0.99*1.0 + 0.0033*0.0) + (L*0.03) * log(0.0033)
       ≈ 150 * (-0.01) + 4.5 * (-5.71)
       ≈ -1.5 - 25.7
       ≈ -27
```

So T_leaf ≈ -27 (or parametrized as a function of L and expected mismatch rate). Reads with best_leaf_score better than T_leaf are strong species-level matches; reads below T_leaf are likely novel or divergent.

The combined algorithm becomes:

```
1. Score all nodes
2. Find S_leaf_max and S_leaf_2nd among leaves only
3. delta_12 = S_leaf_max - S_leaf_2nd

4. Decide voting mode:
   if S_leaf_max > T_leaf:
       # Good leaf match → leaf-centered voting
       anchor = S_leaf_max
       effective_C = C * narrow_factor   if delta_12 > C
                   = C                   otherwise
   else:
       # Poor leaf match → standard voting (let internal nodes set level)
       anchor = max over ALL nodes of S(n)
       effective_C = C

5. Vote all nodes within [anchor - effective_C, anchor + effective_C]
6. LCA of voted nodes → assignment
```

#### Scenario B Revisited With This Rule

T_leaf ≈ -27. Best leaf score = -69.9 < -27 → poor match. Use standard voting.

anchor = -51.2 (genus node, global max). Window: [-56.2, -46.2].
- Genus node at -51.2: **inside** (votes)
- Leaf A at -69.9: **outside** (too far below)
- Leaf B at -87.0: **outside**

LCA of {genus node} = genus node. **Result: genus-level assignment. Correct.**

#### Scenario A Revisited With This Rule

Best leaf score = -30.0 > -27 → good match. Use leaf-centered voting.
delta_12 = -30.0 - (-58.0) = 28 > C=5 → use tight window.

anchor = -30.0. effective_C = 2.5 (0.5 * 5). Window: [-32.5, -27.5].
- Leaf X at -30.0: **inside**
- Internal genus node at -29.8: **inside**
- Higher internal nodes: likely outside the tight window

LCA of {leaf X, genus node} — since leaf X is a descendant of the genus node, the LCA is the genus node. But the genus node's taxonomy at the lowest level (species) is the same as leaf X's (because the genus node's taxonomic assignment is taken from the voted subtree). Actually, looking at the LCA code: `LCA_of_nodes` returns the deepest common ancestor. With one leaf and its parent, the LCA is the parent (genus). Hmm.

Actually, wait. Let me reconsider: in tronko's LCA, the voted nodes include **both** leaf X and the genus node. The LCA of those two is the genus node (since the leaf is a descendant of the genus node). The taxonomy is then taken from the genus node, which has `taxIndex[1]` = 1 (genus level).

This means even in the known-species case, if the genus-level internal node votes alongside the leaf, the assignment would be genus, not species. That's the existing behavior too — the *current* system has this same issue. The species-level assignment only happens when the only voted node within a tree is a single leaf (no internal ancestors within the window).

So the real question for species-level assignment is: **can the tight window exclude the genus-level internal node while keeping the leaf?**

In scenario A: leaf X = -30.0, genus node = -29.8. Difference = 0.2 units. Even with tight C = 2.5, the genus node is 0.2 above the leaf, so it's at -29.8 which is within [-32.5, -27.5]. **Both vote, LCA = genus.**

To get species-level assignment, the leaf must be the ONLY voted node. That requires the window to be tight enough to exclude even the nearest internal ancestor. With a 0.2-unit gap, that means C < 0.2 — impractically tight.

**This reveals that species-level assignment is structurally difficult under the current LCA approach** regardless of where the window is centered. As long as any internal ancestor falls in the window, the LCA climbs up.

However, there's a subtlety in the existing LCA code I should clarify: at `tronko-assign.c:707`, when `count == 1` (votes in only one tree), the LCA is computed from the voted nodes. If the voted set is {leaf X, genus node}, the LCA is the genus node, whose `taxIndex[1]` gives genus level (1). But if the voted set is {leaf X} alone, the LCA is leaf X, whose `taxIndex[1]` is 0 (species level).

So **species-level assignment requires that exactly one leaf and no internal nodes vote in that tree.**

This reframes the improvement: the goal of leaf-centered voting + tight window is to *shrink the window enough that internal nodes don't vote alongside the leaf*. That requires the gap between the best leaf and its nearest internal ancestor to exceed effective_C.

In practice, internal ancestors that are very close parents will almost always be within even a tight window (gap < 1 unit). The nodes we *can* exclude are the **distant** internal ancestors — family-level, order-level nodes. Excluding those prevents the LCA from climbing too high, even if it still lands at genus instead of species.

#### Revised Expectations

Given the structural constraint, leaf-centered + tight window achieves:

- **Prevents over-conservatism**: Stops LCA from climbing to family/order when genus is correct
- **Does NOT generally achieve species-level from genus-level**: The immediate parent internal node is almost always too close in score
- **For novel species**: Standard voting correctly lands at genus when best leaf score is poor

### Potential Concern: Truly Ambiguous Reads

If a read is equally similar to two species in different genera, the two best leaves will be from different genera, and the Cinterval window will include both. The LCA of those two leaves is still the correct higher-level taxon. Leaf-centering doesn't hurt this case — the window is the same whether centered on an internal node or the best leaf, because the best leaf score is close to the best internal score when multiple clades match equally well.

### Variant: Leaf-Only Voting

A more aggressive version: **only allow leaf nodes to vote, never internal nodes**. The LCA of voted leaves is then the assignment.

This would completely eliminate the internal-node bias. But it carries a risk: if the correct leaf doesn't fall within the window but its parent does, we'd lose that signal. I'd recommend starting with leaf-centered (shift the window) rather than leaf-only (exclude internal nodes), and measuring accuracy of both.

### Variant: Depth-Weighted Voting

A gentler version: keep the global maximum as the window center, but weight votes by depth:

```
vote(n) = (depth(n) / max_depth) ^ alpha
```

where alpha > 0 is a tuning parameter. Leaves (high depth) get full weight. Shallow internal nodes get reduced weight. This preserves the vote of every node but downranks the shallow ones that are causing over-conservatism.

**Pro**: No hard cutoff, graceful degradation.
**Con**: Introduces another tuning parameter (alpha). The right value of alpha depends on tree topology and divergence levels, so it's harder to set universally.

---

## Proposed Improvement 2: Score-Distribution-Aware Voting

### The Problem with Fixed Cinterval

Cinterval = 5 is a constant applied identically to every read, regardless of:
- How many nodes scored well (1 clear winner vs. 50 near-ties)
- How spread out the scores are
- How long the read is (150bp vs 300bp — longer reads have wider score ranges)

A read with one dominant species and a clear gap to the next competitor gets the same ±5 window as a read that matches 10 species nearly equally.

### What Existing Analysis Found

Previous confidence analysis (`confidence_analysis.py`, `confidence_analysis2.py`) tested Cinterval-relative rules (e.g., thresholds scaled by Cinterval) and found **no benefit over absolute thresholds**. This suggests the score *distribution shape* doesn't vary enough across reads for full distributional modeling to pay off.

However, one signal is clearly useful and currently ignored: **the gap between the best and second-best leaves** (delta_12).

### Revised Approach: Two-Signal Adaptive Voting

The novel-species analysis from Improvement 1 showed that delta_12 alone is insufficient — we also need the **absolute quality of the best leaf match** to distinguish "known species" from "novel species in a known genus."

The combined decision uses two signals:

```
S_leaf_max  = best leaf score (how well the read matches its best reference leaf)
delta_12    = S_leaf_max - S_leaf_2nd (how much better one leaf is than the next)
```

And one threshold:

```
T_leaf = quality threshold below which the best leaf is "not a good species match"
```

#### Setting T_leaf

T_leaf represents the expected score of a read that genuinely belongs to the best-matching species. For a read of length L with expected intraspecific mismatch rate m (typically 1-5% for amplicons):

```
T_leaf ≈ L * (1-m) * log(0.99) + L * m * log(0.0033)
       ≈ L * (1-m) * (-0.01) + L * m * (-5.71)
```

For L=150, m=0.05 (5% mismatch, conservative):
```
T_leaf ≈ 142.5 * (-0.01) + 7.5 * (-5.71)
       ≈ -1.4 - 42.8
       ≈ -44
```

For L=150, m=0.03 (3% mismatch, typical):
```
T_leaf ≈ 145.5 * (-0.01) + 4.5 * (-5.71)
       ≈ -1.5 - 25.7
       ≈ -27
```

In practice, T_leaf can be a CLI parameter (e.g., `--leaf-quality-threshold`) with a sensible default like -40 for typical amplicon lengths, or derived automatically as `T_leaf = L * threshold_per_site` where `threshold_per_site` is around -0.25.

#### The Decision Algorithm

```
1. Score all nodes
2. Find S_leaf_max and S_leaf_2nd among leaf nodes only
3. delta_12 = S_leaf_max - S_leaf_2nd
4. S_global_max = max over ALL nodes

5. Decide anchor and window width:

   Case 1: S_leaf_max > T_leaf AND delta_12 > C
     → Strong species-level match with clear winner
     → anchor = S_leaf_max, effective_C = C * narrow_factor

   Case 2: S_leaf_max > T_leaf AND delta_12 <= C
     → Good match but multiple species are close
     → anchor = S_leaf_max, effective_C = C
     (Leaf-centered, but standard width lets competing leaves vote.
      Their LCA gives genus-level, which is correct.)

   Case 3: S_leaf_max <= T_leaf
     → Poor leaf match — likely novel species or high divergence
     → anchor = S_global_max, effective_C = C
     (Standard global-max voting. Internal genus/family nodes
      that score well will set the level correctly.)

6. Vote all nodes n where: anchor - effective_C <= S(n) <= anchor + effective_C
7. LCA of voted nodes → assignment
```

#### Why This Handles the Novel-Species Case

In Case 3 (novel species), the best leaf score is poor (e.g., -69.9) — well below T_leaf. The system falls back to standard voting centered on the global max, which is typically the genus-level internal node. That node's posterior probabilities represent the ancestral consensus of the clade, which is the right level of generality for a novel species.

The key insight: **a novel species in a known genus will always have a poor best-leaf score** because no leaf in the DB shares enough sequence identity with it. The poor score is the signal that says "don't trust species-level, use the genus."

#### Why This Doesn't Over-Narrow Known Species

In Case 1 (known species, clear winner), the window is tight but centered on a leaf that genuinely scores well. The read actually *is* that species, so narrowing is correct. The immediate parent internal node (genus) typically scores within ~1-2 units of the leaf and still votes, so the LCA lands at genus or species depending on how tight the window is.

In Case 2 (known genus, multiple candidate species), delta_12 is small, so the standard-width window lets competing leaves from the same genus all vote. Their LCA is the genus node, yielding genus-level assignment — correct behavior when the read could belong to any of several congeners.

#### Remaining Assumptions

**Assumption 1**: T_leaf can be set to a value that separates "intraspecific matches" from "interspecific-but-same-genus matches" for a given marker gene and read length. This is plausible because the mismatch rates are qualitatively different: intraspecific divergence for amplicon markers is typically < 5%, while interspecific divergence within a genus is typically 5-15%.

**Assumption 2**: The genus-level internal node has higher score than any leaf when the read is a novel species. This follows from the math: internal nodes have smoothed PPs that don't pay the full -5.7 penalty per mismatch, while all leaves (being wrong species) accumulate many -5.7 penalties.

**Assumption 3**: For reads from truly novel genera (no genus-level representation in DB), all leaf scores will be very poor and the global-max node will be at family level or higher. This is correct behavior — the system should not claim genus-level confidence when the genus isn't in the DB.

#### What This Does NOT Do

- Does not attempt to distinguish novel species from sequencing error / chimeras (both produce poor leaf scores)
- Does not use the full score distribution — just two summary statistics (S_leaf_max and delta_12)
- Does not require re-tuning when switching marker genes, as long as T_leaf is set appropriately

---

## Proposed Improvement 3: Alignment Quality Filtering (Pre-Scoring)

### The Problem

Every BWA match that enters the pipeline gets a full alignment and full tree scoring, regardless of how good the initial BWA alignment was. BWA computes MAPQ (mapping quality) and the WFA/NW aligner computes an alignment score, but **both are currently discarded**.

A read that aligns ambiguously to 5 different clades (MAPQ = 0) is treated identically to one that aligns perfectly to one leaf (MAPQ = 60). The ambiguous matches dilute the voting signal.

### Where the Information Gets Lost

**BWA MAPQ** — `fastmap.c:179`:
```
sscanf(token, "... %*d ...")
                   ^^^^
                   MAPQ field, discarded by %*d
```

**WFA alignment score** — `placement.c:109`:
```
wavefront_align(wf_aligner, leaf_sequence, leaf_length, query, query_length);
// wf_aligner contains the alignment score, but only the CIGAR is extracted
```

**Mismatch counting** — `placement.c:197`:
```
if (pattern_alg[i] != text_alg[i] && match==0)  // Only first match!
    forward_mismatch++;
```

### Proposed Filtering Points

**Filter 1: MAPQ threshold (in fastmap.c, before match storage)**

```
Extract MAPQ from SAM field 5.
If MAPQ < threshold (e.g., 10), skip this match entirely.
```

This prevents truly ambiguous BWA alignments from entering the pipeline. Reads that map to many references equally (MAPQ = 0) are not informative and only add noise.

**Assumption**: Reads with low MAPQ are uninformative. In standard short-read bioinformatics, MAPQ < 10 is considered unreliable. Since tronko's references are leaf sequences from a phylogeny (often highly similar within genera), some MAPQ depression is expected. A conservative threshold (e.g., 5 or 10) would only remove the most ambiguous cases.

**Expected impact**: Modest. Most reads that match the database at all will have reasonable MAPQ to at least one leaf. But reads in ultra-conserved regions (e.g., rRNA stems) that map equally to dozens of taxa would be filtered, preventing them from inflating vote counts across many unrelated trees.

**Filter 2: Post-alignment mismatch threshold (in placement.c, after WFA)**

After WFA alignment, count mismatches (which is already partially done for match==0). If the mismatch rate exceeds a threshold:

```
mismatch_rate = mismatches / alignment_length
if mismatch_rate > max_mismatch (e.g., 0.15):
    skip this match, don't score the tree
```

**Assumption**: A 15%+ mismatch rate means the read is poorly placed on this particular leaf. The tree-scoring step will assign a very low score anyway, so skipping it saves computation and avoids contributing noise to the voting.

**Risk**: If ALL leaves have >15% mismatch (novel organism not in DB), all matches would be filtered and the read goes unassigned. This is arguably the correct behavior — a read that doesn't match any reference well should not be force-assigned. The `--print-unassigned` flag already handles outputting these.

### Relationship to the Other Improvements

Alignment quality filtering is complementary:
- **Without filtering**: Noisy matches enter voting, inflate internal node votes, push LCA up
- **With filtering**: Only high-quality matches enter, fewer noisy votes, leaf-centered window is even more effective
- **Combined effect**: The matches that survive filtering are confident; the leaf-centered window among confident matches yields specific assignments

---

## What We're NOT Doing

1. **Full Bayesian model** of assignment confidence — too complex, existing analysis shows simple thresholds work nearly as well
2. **SIMD/GPU acceleration** — separate concern, doesn't affect algorithmic accuracy
3. **Changing the PP computation in tronko-build** — the posteriors are mathematically correct; the problem is how voting *uses* them
4. **Replacing BWA** with another aligner — the BWA integration works; we're just extracting more information from it
5. **Adaptive per-position scoring** — weighting sites by informativeness is interesting but adds a lot of complexity for uncertain gain

---

## Implementation Phases

### Phase 1: Two-Signal Adaptive Voting

This is the core change. In `placement.c:862-925`, replace the flat max-finding and fixed-Cinterval voting with the three-case adaptive algorithm.

**What changes:**
1. The max-finding loop (lines 873-888) scans all nodes as before but **also** tracks S_leaf_max, S_leaf_2nd, and S_global_max separately
2. After the loop, apply the decision logic (Cases 1/2/3) to choose anchor and effective_C
3. The voting loop (lines 911-925) uses `anchor` and `effective_C` instead of `maximum` and `Cinterval`

**New CLI parameters:**
- `--leaf-quality-threshold` (default: 0 = disabled / auto). When 0, the system uses standard voting only (backward compatible). When set to e.g. -40, enables the three-case logic.
- `--narrow-factor` (default: 0.5). Multiplier for Cinterval in the "clear winner" case.

**What doesn't change:**
- Score computation (getscore_Arr) — untouched
- Tree traversal (assignScores_Arr_paired) — untouched
- LCA computation — untouched
- Multi-tree consensus — untouched
- All existing behavior preserved when leaf-quality-threshold = 0

**Verification:**

1. With `--leaf-quality-threshold 0` (default): output should be **identical** to current code.

2. With `--leaf-quality-threshold -40 --trace-read '*'`: inspect trace output.
   - Known-species reads (good leaf match): should see "Case 1" or "Case 2" in trace, tighter window, genus/species assignment.
   - Novel-species reads (poor leaf match): should see "Case 3" in trace, standard window, genus assignment.

3. Run on the example dataset simulated reads. Compare:
   - Species-level accuracy (should increase for known species)
   - Genus-level accuracy (should remain high or improve for novel species)
   - Over-specificity rate (reads wrongly assigned to species — should NOT increase)

### Phase 2: MAPQ Extraction and Filtering

Modify `fastmap.c` to extract MAPQ. Add `--min-mapq` CLI parameter (default 0, for backward compatibility). Add mismatch-rate tracking for all matches, not just match==0.

**Verification**: Run with `--min-mapq 10` and compare assignment accuracy. Unassigned count may increase slightly for ambiguous reads. Species-level accuracy should improve.

---

## Testing Strategy

### Accuracy measurement

Use the existing example datasets with known-species simulated reads (`missingreads_singleend_150bp_2error.fasta`). Compare the output taxonomy against the known ground truth at each taxonomic level.

### Key metrics

1. **Species-level accuracy**: % of reads assigned to correct species (should increase)
2. **Genus-level accuracy**: % of reads assigned to correct genus (should remain high)
3. **Over-specificity rate**: % of reads assigned to wrong species (should not increase significantly)
4. **Unassigned rate**: % of reads with no assignment (small increase is acceptable)

### Trace-based debugging

Use `--trace-read '*'` to inspect individual reads that change assignment. The existing trace output already reports all the signals we need:
- `SCORING: max_node is_leaf`
- `BEST_LEAF: score gap`
- `2ND_LEAF: score gap`
- `VOTES: leaves internal total`

---

## Performance Considerations

### Leaf-centered voting: Zero overhead

The max-finding loop already iterates all nodes. Adding an `if (up[0] == -1)` check is O(1) per node, negligible.

### Gap-adaptive Cinterval: Near-zero overhead

Finding second-best leaf requires tracking two scores instead of one. This is already done in the trace code, just needs to be moved outside the trace-only block.

### MAPQ filtering: Slight improvement

Filtering matches before WFA alignment and tree scoring saves computation for filtered matches. Net effect: equal or slightly faster.

## References

- Prior research: `thoughts/shared/research/2026-03-11-tronko-assign-amplicon-assignment-improvements.md`
- Confidence analysis results: `confidence_analysis.py`, `confidence_analysis2.py`
- PP computation: `tronko-build/likelihood.c:1267-1330`
- Scoring function: `tronko-assign/assignment.c:168-236`
- Voting: `tronko-assign/placement.c:862-925`
- LCA consensus: `tronko-assign/tronko-assign.c:660-787`
- BWA SAM parsing: `tronko-assign/bwa_source_files/fastmap.c:168-217`
