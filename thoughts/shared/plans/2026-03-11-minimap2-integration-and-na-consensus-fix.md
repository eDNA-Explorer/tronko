# minimap2 Integration & "NA" Consensus Fix

## Overview

Two independent improvements to tronko-assign:

1. **Replace embedded BWA-MEM with minimap2** for initial leaf matching — better sensitivity at high divergence, faster, cleaner library API
2. **Fix two bugs in the multi-tree consensus algorithm** — a pointer comparison against "NA" that never works, and a counter that accumulates across taxonomy levels instead of resetting

---

## Part A: The "NA" Pointer Comparison Bug (and Its Neighbor)

### What the consensus algorithm does

When a read matches trees in multiple partitions, tronko-assign needs a consensus: do the different trees agree on taxonomy?

The algorithm (`tronko-assign.c:753-780`) works like this:

```
1. Start at species level (minLevel = 0)
2. For every pair of trees (i, j) that produced an LCA:
     if tree_i and tree_j agree on the taxonomy name at minLevel:
       correctTax++
3. If correctTax >= required_pairs:
     → accept this level, stop
   else:
     → move up one level (genus, family, ...)
4. Repeat until we find consensus or exhaust all 7 levels
```

The idea: at species level, maybe only 2 of 5 trees agree. At genus level, maybe 4 of 5 agree. We want the most specific level where enough trees agree.

### Bug 1: The "NA" check never works

Line 764:
```c
taxonomyArr[...][minLevel] != "NA"
```

This compares a `char*` pointer against the string literal `"NA"` using `!=`. In C, this compares **memory addresses**, not string content. The taxonomy strings are dynamically allocated (via `calloc` + `strcpy` in `readreference.c:628-631`), so they live in heap memory. The string literal `"NA"` lives in the read-only data segment. These addresses will **never** be equal, so the `!= "NA"` check **always returns true**.

**Consequence**: When both tree_i and tree_j have `"NA"` at order level (common — many taxa have missing intermediate ranks), `strcmp(...)` returns 0 (they match), and the broken `!= "NA"` check passes, so `correctTax++` fires. Two trees with `"NA;NA;NA;Vahlkampfiidae"` at order level "agree" on order = "NA", which is treated as a genuine taxonomic consensus.

**Real-world impact**: For taxa with sparse taxonomy (e.g., protists, many marine invertebrates), "NA" values are common at class/order/family levels. The bug inflates consensus counts by counting "NA" = "NA" as agreement. This can cause the algorithm to stop at a **more specific level than warranted** — for example, accepting order-level consensus when the only "agreement" was that both trees had no order-level annotation.

**Fix**: Replace pointer comparison with `strcmp`:
```
strcmp(taxonomyArr[...][minLevel], "NA") != 0
```

### Bug 2: correctTax accumulates across levels

`correctTax` is initialized to 0 at line 753 (before the while loop) but **never reset** inside the loop body. Each iteration of `while(stop==0 && minLevel<=6)` adds to the running total.

**Worked example** with 3 trees (3 pairs, required_pairs = 2):

```
Species level (minLevel=0):
  pair(0,1): disagree → correctTax stays 0
  pair(0,2): disagree → correctTax stays 0
  pair(1,2): agree    → correctTax = 1
  Check: 1 >= 2? No → move to genus

Genus level (minLevel=1):
  pair(0,1): agree    → correctTax = 2   ← SHOULD be 1, but accumulated from species
  pair(0,2): disagree → correctTax stays 2
  pair(1,2): agree    → correctTax = 3
  Check: 3 >= 2? Yes → accept genus
```

But the correct count at genus level is 2 (only pairs 0,1 and 1,2 agree at genus). It happens to still pass, but with different numbers: `correctTax` reached 3 when it should have been 2. For edge cases where the accumulated species-level count pushes genus-level past the threshold early, the algorithm accepts genus **before finishing all pair comparisons**, which could give incorrect `taxRoot`/`taxNode` selection.

Wait — actually re-reading the code: the for-loops run to completion before the check at line 771. So correctTax accumulates ALL agreements across ALL levels tried so far. At genus level, correctTax = (species agreements) + (genus agreements). This means the threshold is easier to meet at higher levels not just because more pairs agree, but because previous levels' agreements carry over. This is incorrect — each level should be evaluated independently.

**Fix**: Reset `correctTax = 0` at the top of the while loop body, right after line 759.

### Interaction between the two bugs

They partially cancel each other out, which may explain why this code hasn't been flagged before:

- Bug 1 (NA always passes) **inflates** correctTax by counting NA=NA as agreement
- Bug 2 (no reset) **inflates** correctTax by accumulating across levels
- Combined: the algorithm over-counts agreements, making it stop at more specific levels than it should

After fixing both, the consensus algorithm becomes stricter. Some reads that previously got genus-level assignment might now get family-level or unassigned. This is **correct** behavior — the consensus was being artificially inflated.

### Assumptions

**Assumption 1**: "NA" in taxonomy files means "unknown at this level" and should NOT count as agreement between trees. Two trees both having "NA" at order level means neither knows the order — that's absence of information, not consensus. **This is the standard interpretation in biodiversity databases.**

**Assumption 2**: The consensus threshold is meant to evaluate each taxonomic level independently. Agreement at species level should not contribute to the genus-level count. Each level asks a fresh question: "do enough trees agree at THIS level?"

### Expected impact

- Taxa with complete taxonomy (all 7 levels filled): **no change** — no "NA" values, and if correctTax was passing at species it still will after reset
- Taxa with sparse taxonomy (NA at intermediate levels): **assignments may become less specific** — the correct behavior, since there genuinely isn't consensus at those levels
- Multi-tree databases with many partitions: **largest impact**, since more tree pairs means more opportunities for NA=NA inflation

---

## Part B: minimap2 as BWA Replacement

### What BWA does today (and what we need from a replacement)

BWA-MEM serves as a **rapid candidate finder**: given a query amplicon read (~150bp), find which leaf reference sequences (~250bp barcodes) it might match. The output is a list of `(tree_id, node_id)` pairs — candidate leaves for full phylogenetic scoring.

**The full pipeline**:
```
Query read (150bp)
    ↓
[BWA-MEM: seed → extend → align → SAM]     ← we're replacing this
    ↓
SAM parsing → leaf name → hashmap → (tree_id, node_id)
    ↓
[WFA/NW: re-align read to leaf reference]   ← NOT changing this
    ↓
Phylogenetic scoring → voting → LCA → taxonomy
```

**Critical observation**: BWA's alignment is partially redundant. WFA re-aligns the read from scratch against each candidate leaf. The only thing BWA's CIGAR is used for is the `use_portion` optimization (extracting a window around the hit in the reference sequence, via `getStartPosition`/`getEndPosition`). Otherwise, BWA is purely a candidate finder.

### Why BWA-MEM is suboptimal for this task

BWA-MEM was designed for mapping short reads to large genomes (human genome, 3Gb). Our use case is different in every dimension:

| Property | Genome mapping | Amplicon candidate finding |
|---|---|---|
| Reference size | 3 Gb | 250 KB – 25 MB |
| Reference sequences | 1-24 chromosomes | 1,000 – 100,000 short sequences |
| Expected matches per read | 1 (unique location) | 1-50 (similar species) |
| Query length | 150 bp | 150 bp |
| Reference sequence length | millions of bp | 226-315 bp |
| Goal | precise genomic coordinate | candidate list |

BWA's BWT index, Smith-Waterman extension, insert-size modeling, and chimeric read detection are overhead for our task. We need something optimized for many-short-sequences rather than few-long-sequences.

### Why minimap2

minimap2 is the strongest candidate because it:

1. **Uses minimizer indexing instead of BWT** — faster to build, faster to query, better sensitivity at >5% divergence (minimizers tolerate mutations in non-sampled positions, while BWA's exact-match seeds miss matches when mutations land in seed regions)

2. **Has a C library API** (`minimap2.h`) — can be embedded as source like BWA is today. No external process calls, no temp files, no SAM parsing overhead.

3. **Handles paired-end natively** — `mm_map()` with `mm_mapopt_t` supports paired mapping

4. **Returns structured results** — `mm_reg1_t` structs with CIGAR, MAPQ, alignment score, identity. No need to parse SAM text.

5. **Is 3-5x faster** for short reads in practice — less overhead per query, vectorized alignment extension.

### How minimizer indexing works (vs BWA's BWT)

**BWA's BWT approach**:
```
1. Build suffix array of entire reference concatenation
2. Construct Burrows-Wheeler Transform (BWT) from suffix array
3. For each query: find Maximal Exact Matches (MEMs) via backward search in BWT
4. Chain MEMs into candidate alignments
5. Smith-Waterman extension for each candidate
```

The key limitation: MEMs require **exact matches** of length ≥ seed_length (default 17bp). If a query has mutations every 15bp (which happens at ~7% divergence), no 17bp exact match exists, and BWA finds nothing.

**minimap2's minimizer approach**:
```
1. For each reference sequence: extract (w, k)-minimizers
   - Slide a window of w positions
   - In each window, hash every k-mer and keep the minimum hash
   - Minimizers are a representative subset of k-mers (~2/(w+1) of all k-mers)
2. Build hash table: minimizer → list of (sequence_id, position)
3. For each query: extract the same minimizers, look up in hash table
4. Chain matching minimizers into candidate alignments
5. Base-level alignment (KSW2) for each candidate
```

The advantage: minimizers are **sampled** k-mers, not **maximal exact matches**. Even if some k-mers are mutated, the query and reference share enough minimizers to find the match. With k=15 and w=10, minimap2 finds matches at up to ~15% divergence where BWA would fail.

**For our amplicon data** (references ~250bp, queries ~150bp, expected divergence 0-15%):
- BWA seeds need ≥17bp exact matches — fails at >~7% divergence in worst case
- minimap2 minimizers with k=15, w=5 — works up to ~15% divergence
- This matters most for **novel species** where the nearest leaf may be 8-12% divergent

### The integration design

#### What changes

```
BEFORE:
tronko-assign.c → run_bwa() → main_mem() [bwa_source_files/fastmap.c]
                                  ↓
                              SAM text output
                                  ↓
                              SAM parsing → leaf hashmap → bwaMatches struct

AFTER:
tronko-assign.c → run_minimap2() → mm_idx_build() + mm_map() [minimap2 lib]
                                        ↓
                                   mm_reg1_t array (structured)
                                        ↓
                                   leaf coordinate extraction → bwaMatches struct
```

The `bwaMatches` struct stays the same. Everything downstream of it is untouched.

#### Index building

**Current (BWA)**: Pre-built index loaded from disk (`.bwt`, `.sa`, `.pac`, `.ann`, `.amb` files). Built externally or via tronko-assign's `-6` flag.

**With minimap2**: Build index at runtime from the reference FASTA.

```
mm_idx_t *mi = mm_idx_build(reference_fasta_path, w, k, 0);
```

Where:
- `w` = minimizer window size (recommend 5 for short amplicons)
- `k` = k-mer size (recommend 15 for amplicon-level divergence)
- The `0` flag = build from file

This takes <1 second for typical amplicon databases (1K-100K sequences, ~25MB total). Since tronko-assign already loads the full database into memory at startup, adding ~0.5s for index building is negligible.

Alternatively, minimap2 supports pre-built `.mmi` index files for instant loading. The index could be built once during tronko-build and stored alongside `reference_tree.txt`.

#### Per-read mapping

```
For each batch of reads:
  mm_reg1_t *regs = mm_map(mi, query_length, query_sequence, &n_regs, tbuf, &mopt, NULL);

  for each reg in regs[0..n_regs-1]:
    leaf_name = mi->seq[reg->rid].name    // reference sequence name
    leaf_map  = hashmap_get(&map, leaf_name)  // same hashmap as today
    tree_id   = leaf_map->root
    node_id   = leaf_map->node
    mapq      = reg->mapq                 // NEW: MAPQ available for free
    score     = reg->score                // NEW: alignment score available
    cigar     = reg->p->cigar             // CIGAR for use_portion
    start_pos = reg->rs                   // start position on reference
```

**Key difference from BWA**: no SAM text intermediary. The results are structured C arrays. No `sscanf`, no `strtok_r`, no string parsing. This is both faster and less error-prone.

#### Paired-end handling

minimap2 has `mm_map_frag()` for paired/fragment mapping that handles:
- Insert size estimation
- Concordant/discordant classification
- Proper pair orientation

Alternatively, we can map R1 and R2 independently (both are amplicons of the same marker) and determine concordance ourselves: if both map to leaves in the same tree, it's concordant. This matches the current logic in `fastmap.c:196-271`.

#### What the leaf hashmap sees

No change. The hashmap maps leaf names → (tree_id, node_id). Whether those names come from parsing SAM text (BWA) or from `mi->seq[reg->rid].name` (minimap2), the hashmap lookup is identical.

#### CIGAR and start position (use_portion mode)

minimap2's `mm_reg1_t` provides:
- `reg->rs`: start position on reference (0-based, unlike BWA's 1-based SAM POS)
- `reg->p->cigar`: CIGAR as integer array (not string), each element encoding op + length

For the `use_portion` path, we'd either:
1. Convert minimap2's integer CIGAR to a string and feed it to the existing `getEndPosition()` parser — minimal change
2. Rewrite `getStartPosition()`/`getEndPosition()` to consume integer CIGARs directly — cleaner but more code to change

Option 1 is recommended for initial integration. The CIGAR string conversion is ~5 lines.

### Sensitivity improvement: why this matters for assignment

The biggest practical impact isn't speed — it's **finding candidates that BWA misses**.

Consider a novel species with 10% divergence from its nearest leaf. In a 150bp read:
- ~15 mismatched positions, randomly distributed
- Average spacing between mismatches: ~10bp
- BWA needs a 17bp exact match: probability of finding one = (1 - 0.10)^17 ≈ 0.167 per position. Over ~134 possible positions, expected MEMs ≈ 22. But this is the *average* — for unlucky mutation distributions (clustered mismatches), BWA may find 0 seeds and report no match.
- minimap2 with k=15 needs a 15bp exact match, but only samples one per w=5 window: more tolerant of clustered mutations because it looks at more positions.

When BWA misses a candidate leaf entirely, the read either goes unassigned or gets assigned based on a less-related leaf. minimap2's better sensitivity means more correct candidates enter the scoring pipeline, which improves downstream assignment accuracy — especially for novel species in known genera (the scenario from the leaf-centered voting plan).

### What we're NOT changing

- **WFA/NW alignment** — untouched, still re-aligns from scratch
- **Posterior probability scoring** — untouched
- **Voting/LCA** — untouched
- **The `bwaMatches` struct** — stays as the interface; we just populate it differently
- **Thread model** — minimap2 wrapper runs single-threaded (matching current BWA behavior)
- **All downstream output** — identical format

### Assumptions

**Assumption 1**: minimap2's minimizer indexing is well-suited for short amplicon references. This is well-established — minimap2 is routinely used for amplicon-length sequences in tools like medaka and longshot. The `-x sr` preset or custom `w=5, k=15` parameters are appropriate.

**Assumption 2**: The sensitivity improvement at >5% divergence will translate to better assignment accuracy for divergent reads. This follows from the pipeline structure: more correct candidates → better scoring → better assignment. A candidate that BWA misses can never be rescued downstream.

**Assumption 3**: Building the minimap2 index at runtime (~0.5s) is acceptable. If not, pre-built `.mmi` files can be generated by tronko-build. This adds a build step but eliminates runtime cost.

**Assumption 4**: minimap2 source can be embedded in the tronko-assign build tree the same way BWA is today. minimap2 is pure C (with some SSE/NEON intrinsics for SIMD), MIT-licensed, and has ~30 source files. The Makefile change mirrors the current BWA file list.

### Risks

**Risk 1**: Different default behavior. minimap2 may return more or fewer matches than BWA for the same query. The `max_bwa_matches` cap still applies, but the *set* of candidates may differ. This is intentional (better candidates), but means results won't be bit-identical to BWA-based runs. **Mitigation**: Keep BWA as a compile-time or runtime option for the transition period.

**Risk 2**: CIGAR format differences. minimap2 uses `=/X` (sequence match/mismatch) operators by default with `--eqx`, while BWA uses `M` (match-or-mismatch). The `getEndPosition()` parser handles `M`, `I`, `D`, `H`, `S` — would need `=` and `X` added. **Mitigation**: Either add `=/X` to the parser or configure minimap2 to output `M`-style CIGARs.

**Risk 3**: Memory usage. minimap2's minimizer index for ~25MB of reference sequence is ~50-100MB. BWA's BWT index for the same data is ~30-60MB. Slightly higher, but both are small relative to the posterior probability arrays that dominate tronko-assign's memory. **Mitigation**: Monitor with the existing resource logging.

### Testing strategy

1. **Candidate comparison**: Run both BWA and minimap2 on the same reads, compare which leaf candidates each returns. Expect minimap2 to find a superset (same + additional divergent matches).

2. **Assignment accuracy**: Compare final taxonomy assignment accuracy on the example dataset (`missingreads_singleend_150bp_2error.fasta`). Species/genus accuracy should be equal or better.

3. **Novel species simulation**: Take reads from species X, remove species X from the database, rebuild, and assign. BWA may fail to find candidate leaves at >8% divergence; minimap2 should still find them.

4. **Performance**: Time both approaches on the same dataset. Expect 2-5x speedup from minimap2.

5. **Regression**: Verify that for reads with clear, low-divergence matches (the common case), BWA and minimap2 produce identical assignments. The downstream scoring should be deterministic given the same candidates.

---

## Implementation Phases

### Phase 1: Fix consensus bugs (small, safe, immediate)

Two changes in `tronko-assign.c:753-780`:

1. **Reset correctTax** at the start of each while-loop iteration (after line 759)
2. **Fix "NA" comparison** — replace `!= "NA"` with `strcmp(..., "NA") != 0`

**Verification**:
- With all-complete taxonomy (no NA values): output should be identical
- With sparse taxonomy (NA at intermediate levels): some assignments may become less specific
- The `--trace-read` output will show different `correctTax` and `minLevel` values

### Phase 2: Add minimap2 source and dual-aligner infrastructure

Both BWA and minimap2 must remain fully functional. A runtime `--aligner` flag selects which one runs.

#### 2a. Add minimap2 source

1. Add minimap2 source files to `tronko-assign/minimap2_source_files/`
2. Update Makefile to compile **both** BWA and minimap2 source files into the binary
   - BWA files stay in `bwa_source_files/` — no changes
   - minimap2 files added as additional object files
   - Both are always compiled and linked

#### 2b. CLI flag: `--aligner bwa|minimap2`

Add to `options.c`:
- `--aligner` with string argument, values `"bwa"` (default) or `"minimap2"`
- Store in `Options.aligner` as an enum or string
- Default to `"bwa"` so existing behavior is unchanged — users opt into minimap2 explicitly
- Print selected aligner in startup log

#### 2c. Write `run_minimap2()` wrapper

New function in a new file (e.g., `minimap2_wrapper.c`) with the same contract as `run_bwa()`:

```
int run_minimap2(int start, int end, bwaMatches *results, int concordant,
                 int ntree, char *databaseFile, int paired,
                 int max_query_length, int max_readname_length,
                 int max_acc_name, int max_bwa_matches,
                 int mm2_kmer, int mm2_window);
```

This function:
1. Builds minimizer index from reference FASTA via `mm_idx_build()`
   - Index is built once and cached for the lifetime of the process (stored in a static or global pointer)
   - Subsequent calls reuse the cached index
2. Maps each read via `mm_map()`, extracting:
   - `mi->seq[reg->rid].name` → leaf name → hashmap → `(tree_id, node_id)`
   - `reg->rs` → start position (convert from 0-based to 1-based to match BWA convention)
   - `reg->p->cigar` → convert integer CIGAR to string for `use_portion` compatibility
   - `reg->mapq` → MAPQ (stored in bwaMatches if/when MAPQ field is added)
3. Populates `bwaMatches` struct identically to how `main_mem()` does
4. Handles concordant/discordant paired-end logic:
   - Map R1 and R2 independently
   - If both map to leaves in the same tree → concordant
   - Otherwise → discordant
   - Same logic as `fastmap.c:196-271`, just applied to minimap2 results

The leaf name hashmap (`leafMap`) is currently built inside `main_mem()` (`fastmap.c:112-123`). For the dual-aligner design, this hashmap should be **extracted** to a shared location:
- Build the hashmap once in `tronko-assign.c` before the aligner call
- Pass it as a parameter to both `run_bwa()` and `run_minimap2()`
- This avoids duplicating the hashmap construction code

#### 2d. Dispatch in tronko-assign.c

At the current BWA call site (`tronko-assign.c:379`), add a branch:

```
if aligner == "minimap2":
    run_minimap2(start, end, bwa_results, concordant, ...)
else:
    run_bwa(start, end, bwa_results, concordant, ...)
```

Everything after this point — the `bwa_results` → `leaf_coordinates` extraction, WFA alignment, scoring, voting, LCA — is identical regardless of aligner. No changes needed downstream.

#### 2e. Index handling

**BWA path** (unchanged): Expects pre-built `.bwt`/`.sa`/`.pac`/`.ann`/`.amb` index files on disk. Loaded via `bwa_idx_load()`.

**minimap2 path**: Two options, both supported:
- If a `.mmi` file exists alongside the reference FASTA → load it via `mm_idx_reader_open()` (instant)
- If no `.mmi` file → build index at runtime from reference FASTA via `mm_idx_build()` (~0.5s for typical databases)
- Optionally add `--save-minimap2-index` flag to write the `.mmi` file after building, for faster subsequent runs

**Verification**:
- With `--aligner bwa`: output **bit-identical** to current code (no BWA code paths changed)
- With `--aligner minimap2`: compare candidates and final assignments against BWA baseline
- Both paths tested on single-tree and multi-tree example datasets

### Phase 3: Tune minimap2 parameters for amplicon data

1. Test k-mer sizes (k=11, 13, 15) and window sizes (w=3, 5, 10)
2. Evaluate sensitivity vs specificity tradeoff on example datasets
3. Set optimal defaults for typical amplicon databases
4. Add `--minimap2-kmer` and `--minimap2-window` CLI options for power users
5. Both BWA and minimap2 remain available — tuning only affects minimap2 path

---

## References

- Prior research: `thoughts/shared/research/2026-03-11-tronko-assign-amplicon-assignment-improvements.md` (BWA alternatives section)
- Related plan: `thoughts/shared/plans/2026-03-11-leaf-centered-voting-and-alignment-filtering.md`
- NA bug location: `tronko-assign/tronko-assign.c:764`
- correctTax bug location: `tronko-assign/tronko-assign.c:753-780`
- BWA integration: `tronko-assign/bwa_source_files/fastmap.c:577` (`main_mem`)
- BWA call site: `tronko-assign/tronko-assign.c:176-183`
- Leaf hashmap: `tronko-assign/bwa_source_files/fastmap.c:112-123`
- bwaMatches struct: `tronko-assign/global.h:197-211`
- use_portion CIGAR parsing: `tronko-assign/getSequenceinRoot.c:36-99`
- minimap2 library API: `minimap2.h` (`mm_idx_build`, `mm_map`, `mm_reg1_t`)
