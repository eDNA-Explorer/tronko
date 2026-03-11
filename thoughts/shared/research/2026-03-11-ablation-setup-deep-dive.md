---
date: 2026-03-11T12:00:00-08:00
researcher: Claude
git_commit: 1ad46ba3d0c78145269c7ca5729826d406e905cd
branch: optimize-tronko-build
repository: tronko-fork
topic: "Deep dive into the tronko ablation system: logic, accuracy, and correctness analysis"
tags: [research, ablation, tronko-build, posteriors, phylogenetics, nw_prune]
status: complete
last_updated: 2026-03-11
last_updated_by: Claude
---

# Research: Tronko Ablation System — Logic, Accuracy, and Correctness Analysis

**Date**: 2026-03-11
**Researcher**: Claude
**Git Commit**: 1ad46ba
**Branch**: optimize-tronko-build
**Repository**: tronko-fork

## Research Question
Extensive deep dive into the ablation setup in tronko — assess its logic and accuracy for ablating sequences, with reference to the core tronko algorithmic concepts.

## Summary

The tronko ablation system (`ablate-tronko-db.sh`) removes sequences from a reference database without full rebuild. It prunes trees with `nw_prune`, filters MSAs/taxonomy with awk, and then invokes `tronko-build -y -v -f 999999 -a` to recompute posteriors on the pruned subtrees while skipping SP-score repartitioning. The approach is **mathematically sound** — all critical quantities (GTR model parameters, branch lengths, conditional likelihoods, posteriors) are fully re-estimated from scratch on the reduced data. The main trade-offs are (a) the MSA alignment is inherited rather than recomputed, and (b) partition boundaries are frozen from the original build.

---

## Core Tronko Concepts (Background)

### Posterior Computation Pipeline

Tronko-build's core computation for each subtree:

1. **Estimate GTR model parameters** (`estimatenucparameters_Arr` at `likelihood.c:1037`): Iterative MLE of 6 rate parameters + 4 base frequencies + branch lengths via multiple rounds of `estimatebranchlengths_Arr` and `maximizelikelihoodnc_globals_Arr`.

2. **Compute conditional likelihoods** (`makeconnc_Arr`): Bottom-up Felsenstein pruning. At each internal node, for each site and each possible ancestral nucleotide:
   ```
   like[s][i] = (Σ_j P(i→j|t_left) * left_child_like[s][j]) * (Σ_j P(i→j|t_right) * right_child_like[s][j])
   ```

3. **Compute posteriors** (`getposterior_nc_Arr` at `likelihood.c:1267`): Two phases:
   - **Top-down pass** (`makeposterior_nc_Arr` at `likelihood.c:1140`): Propagates information from parent+sibling down to each node
   - **Final combination**: For internal nodes: `posterior[s][i] = like[s][i] * posterior_from_above[s][i] * pi[i]`, normalized. For leaf nodes with known base: `posterior[s][base] = 1.0`, others 0.0. For missing data (base=4): posterior derived from parent's posterior through transition matrix.

4. **Missing data adjustment** (`changePP_Arr` + `changePP_parents_Arr` at `tronko-build.c:2125-2134`): Marks gaps with -1 posteriors and propagates upward.

5. **Log-posterior transform** (`set_posteriors` at `likelihood.c:1257`):
   ```
   stored[node][site][base] = log(0.01/4 + 0.99 * posterior[node][site][base])
   ```

### How tronko-assign Uses Posteriors

During assignment (`getscore_Arr` in `assignment.c`), for each query:
- Score at each node = Σ over aligned sites of `stored_posterior[node][site][query_base]`
- Find maximum score across all nodes/trees
- All nodes within `Cinterval` (default=5) of maximum "vote"
- LCA of all voting nodes determines taxonomic assignment

### SP-Score Partitioning

Sum-of-pairs score measures MSA column quality: +3 match, -2 mismatch, -1 gap. Normalized. If below threshold, the subtree is split, sequences re-aligned, new tree estimated. This ensures alignment quality backing each subtree.

---

## Ablation System Architecture

### Files

| File | Role |
|---|---|
| `ablate-tronko-db.sh` | Main ablation script |
| `ABLATION.md` | Documentation |
| `sweep_tronko_build.py` | Parameter sweep tool (not ablation-specific but shares infrastructure) |

### Workflow

```
Master Build (with -E flag)
    ↓
exported_subtrees/ directory
    ├── {N}_MSA.fasta (aligned)
    ├── RAxML_bestTree.{N}.reroot
    └── {N}_taxonomy.txt
    ↓
ablate-tronko-db.sh
    ├── Step 1: Prune each cluster
    │   ├── nw_prune (remove leaves from Newick tree)
    │   ├── awk (filter MSA)
    │   └── awk (filter taxonomy)
    ├── Step 2: tronko-build -y -v -f 999999 -a
    │   └── Full posterior recomputation (no repartitioning)
    └── Step 3: Finalize
        ├── Concatenate MSAs → marker.fasta (gaps stripped)
        ├── BWA index
        └── gzip reference_tree.txt
```

---

## Detailed Logic Analysis

### Step 1: Tree Pruning (`ablate-tronko-db.sh:82-140`)

For each cluster:
1. Extract accessions from MSA, sort them (`ablate-tronko-db.sh:94-96`)
2. Find intersection with remove list using `comm -12` (sorted set intersection)
3. Three outcomes:
   - **No removals**: Copy files as-is with renumbered index
   - **< 3 remaining**: Drop the entire cluster (too few sequences for meaningful phylogenetics)
   - **Otherwise**: Prune tree, filter MSA, filter taxonomy

**Tree pruning** (`ablate-tronko-db.sh:122-123`):
```bash
nw_prune -f "$tree_file" "$cluster_removes" > pruned_tree
```
`nw_prune -f` reads leaf names from a file and removes them from the Newick tree. When removing a leaf from a binary tree, the parent internal node becomes unnecessary — `nw_prune` collapses it by connecting the sibling directly to the grandparent, adding the collapsed node's branch length to the sibling's. The resulting tree is still binary.

**MSA filtering** (`ablate-tronko-db.sh:126-130`):
```awk
/^>/ { name=substr($0,2); split(name,a," "); skip=rm[a[1]] ? 1 : 0 }
!skip { print }
```
Removes entire FASTA entries (header + sequence) for matching accessions. Remaining sequences keep their original alignment (gap positions unchanged).

**Taxonomy filtering** (`ablate-tronko-db.sh:133-136`):
Removes matching rows from the tab-separated taxonomy file.

**Cluster renumbering**: Valid clusters are renumbered sequentially (0, 1, 2, ...) via `valid_cluster_count`, which is what tronko-build expects when reading from `-e` directory (`tronko-build.c:1654` calls `readFilesInDir` expecting `{i}_MSA.fasta`, `RAxML_bestTree.{i}.reroot`, `{i}_taxonomy.txt`).

### Step 2: Posterior Recomputation (`ablate-tronko-db.sh:167-173`)

```bash
tronko-build -y -e "$CACHE_DIR" -n "$valid_cluster_count" -d "$OUTPUT_DIR" -v -f 999999 -a -c "$THREADS"
```

**Flag analysis:**

| Flag | Value | Effect |
|---|---|---|
| `-y` | — | Partition directory mode (`use_partitions=1`) |
| `-e` | CACHE_DIR | Read pre-existing cluster files from this directory |
| `-n` | count | Number of clusters to read |
| `-v` | — | Enable min-leaves partitioning mode (`use_min_leaves=1`) |
| `-f` | 999999 | Min leaves threshold = 999999 |
| `-a` | — | Use VeryFastTree/polytomy resolution path (`fasttree=1`) |
| `-c` | THREADS | Thread count |

**The partition-skip trick**: At `tronko-build.c:1796`:
```c
if ( opt.use_partitions==1 && m->numspec > opt.min_leaves ){
    SPscoreArr[i]=0;
    createNewRoots(i,opt,max_nodename,max_lineTaxonomy,m);
}
```
Since no cluster will have >999,999 sequences, `createNewRoots` is never called. The clusters pass through unpartitioned.

**What then happens** (`tronko-build.c:2087-2105`):
```c
#pragma omp parallel for schedule(dynamic) private(i)
for(i=0; i<numberOfTrees; i++){
    double local_params[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    estimatenucparameters_Arr(local_params,i);
    getposterior_nc_Arr(local_params,i);
}
```

This is the **complete posterior pipeline**: GTR parameters are estimated from scratch (starting from uniform priors), branch lengths are re-optimized, conditional likelihoods are recomputed bottom-up, and posteriors are computed top-down. This runs on the pruned tree with the reduced MSA.

### Step 3: Finalization (`ablate-tronko-db.sh:184-192`)

- Concatenates all cluster MSAs, stripping alignment gaps (`gsub(/-/, "")`) to produce the unaligned FASTA for BWA indexing
- Builds BWA index for tronko-assign alignment
- Gzips `reference_tree.txt`

---

## Correctness Assessment

### What is Fully Recomputed (Correct)

1. **GTR model parameters**: `estimatenucparameters_Arr` re-estimates all substitution model parameters from the reduced data. The model is not inherited from the original build — it's estimated fresh with uniform starting values `{0.0, 1.0, 1.0, ..., 1.0}`.

2. **Branch lengths**: Multiple rounds of `estimatebranchlengths_Arr` re-optimize all branch lengths on the pruned tree topology. Branch lengths from `nw_prune` serve as starting points, but are re-optimized.

3. **Conditional likelihoods**: Fully recomputed bottom-up from the reduced MSA on the pruned tree (via `getlike_gamma_Arr` called inside `getposterior_nc_Arr`).

4. **Posteriors**: Fully recomputed via the complete top-down pass (`makeposterior_nc_Arr`) followed by normalization and the `set_posteriors` log-transform.

5. **Missing data handling**: `changePP_Arr` + `changePP_parents_Arr` run on the recomputed posteriors.

6. **BWA index**: Rebuilt from the remaining sequences.

### What is Inherited (Trade-offs)

1. **MSA alignment**: The alignment of the remaining sequences is the same alignment computed when all sequences were present. When you remove a sequence, the remaining sequences' alignment positions don't change, but:
   - The removed sequence may have influenced gap placement in the original MSA (e.g., FAMSA aligns all sequences simultaneously)
   - Some columns may become all-gaps or near-all-gaps
   - For **small removals** (leave-one-out, removing a few sequences), this effect is negligible
   - For **large removals** (removing an entire clade), the MSA quality impact could be more significant but alignment artifacts tend to be local

2. **Partition boundaries**: The SP-score partitioning that determined which sequences belong to which subtree is frozen. After removal:
   - Some subtrees might become small enough that they were unnecessary partitions
   - The SP-score of remaining sequences might differ from the original
   - The script mitigates this by dropping clusters with <3 sequences

3. **Tree topology**: The pruned tree inherits the topology from the original RAxML/VeryFastTree estimation. With sequences removed:
   - The tree may not be the ML tree for the remaining sequences
   - However, the pruned tree is still a valid, topologically correct subtree of the original
   - Branch length re-optimization partially compensates for topology imprecision

### Tree Pruning Correctness

`nw_prune` on a binary tree produces a binary tree:
- Removing a leaf → parent becomes unary → `nw_prune` collapses it: sibling connects to grandparent with branch length = sibling_bl + collapsed_parent_bl
- This preserves the total evolutionary distance between remaining taxa
- The `-a` flag ensures any accidental polytomies are resolved via `parseNewick` + `makeBinary`

### Edge Cases Handled

| Scenario | Handling |
|---|---|
| Sequence not in any cluster | Warning printed, count tracked (`ablate-tronko-db.sh:143-147`) |
| Cluster emptied completely | Dropped (< 3 remaining check at line 112) |
| All clusters emptied | Fatal error at line 157-159 |
| Missing tree/taxonomy for cluster | Skipped with warning (line 87-90) |
| Cluster with no removals | Copied as-is (line 102-109) |

### Potential Issues

1. **The `-a` flag (fasttree path) rewrites the tree file**: At `tronko-build.c:1761`, `exportTreeToNewick(m, buffer)` overwrites the tree file in the cache directory. This is fine for the ablation workflow since the cache is a copy, but means the original exported subtrees are not modified.

2. **Accession format sensitivity**: The ablation script uses exact string matching (sorted `comm -12`) between the remove list and FASTA headers. The remove list must use the post-sanitization format (`:` → `_`). The `ABLATION.md` documents this requirement.

3. **Alignment column degeneration**: If many sequences from one region of the tree are removed, some alignment columns may become uninformative (all same base or all gaps). This doesn't cause errors but may slightly dilute posterior signal. The posteriors at those columns will approach uniform (0.25 each), contributing minimal discriminative power during assignment.

---

## Accuracy Assessment

### When Ablation is Most Accurate

- **Leave-one-out studies**: Removing a single sequence has minimal impact on MSA quality and tree topology. The recomputed posteriors will be very close to what a full rebuild would produce.
- **Sparse removal across clusters**: Removing a few sequences distributed across many clusters — each cluster is minimally affected.
- **Removing redundant sequences**: If the removed sequence has close relatives remaining in the cluster, the tree topology and alignment are robust.

### When Ablation May Diverge from Full Rebuild

- **Removing an entire clade**: If all representatives of a genus/family are removed from a cluster, the tree topology in that region becomes meaningless (internal nodes that were ancestral to the clade are now collapsed). The posteriors at remaining nodes may differ from a fresh build.
- **Removing many sequences from one cluster**: Large removals may meaningfully change alignment quality and optimal tree topology.
- **Very small clusters post-removal**: Clusters with 3-5 remaining sequences have minimal phylogenetic signal; posteriors may be noisy.

### Theoretical Guarantee

The ablation produces a **mathematically valid** reference database: the posteriors are proper Bayesian posteriors conditioned on the reduced data and the (fixed) tree topology. The only question is whether the fixed topology is close to the optimal topology for the reduced data — and for small removals, it is.

---

## Code References

- `ablate-tronko-db.sh:82-140` — Step 1: tree pruning, MSA/taxonomy filtering
- `ablate-tronko-db.sh:167-173` — Step 2: tronko-build invocation with partition-skip trick
- `ablate-tronko-db.sh:184-192` — Step 3: FASTA concatenation + BWA indexing
- `tronko-build/tronko-build.c:1796` — Min-leaves check that enables partition-skip
- `tronko-build/tronko-build.c:2087-2105` — Posterior computation loop
- `tronko-build/likelihood.c:1037-1058` — `estimatenucparameters_Arr` (full model re-estimation)
- `tronko-build/likelihood.c:1140-1189` — `makeposterior_nc_Arr` (top-down posterior pass)
- `tronko-build/likelihood.c:1267-1338` — `getposterior_nc_Arr` (complete posterior computation)
- `tronko-build/likelihood.c:1257-1266` — `set_posteriors` (log-pseudocount transform)
- `tronko-build/tronko-build.c:2006-2071` — `-E` export subtrees implementation
- `tronko-build/tronko-build.c:1746-1775` — VeryFastTree polytomy resolution path

## Open Questions

1. **Empirical validation**: How do ablated database assignments compare to full-rebuild assignments on real datasets? Leave-one-out experiments would quantify the accuracy gap.
2. **MSA column trimming**: Would trimming all-gap or near-all-gap columns from the reduced MSA improve accuracy? (Currently not done.)
3. **Threshold for "too many removals"**: At what fraction of removed sequences does the ablation quality degrade significantly compared to a full rebuild?
