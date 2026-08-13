# Best-Leaf Override + Cross-Tree Score Normalization Implementation Plan

## Overview

Add two **opt-in, flag-gated, default-OFF** accuracy features to `tronko-assign` on the `high-perf-aligner` branch so they can be grid-searched for optimal settings on AncestralClust (AC) databases:

1. **Best-leaf override** — ported from `optimize-tronko-build`. After normal Cinterval voting, if a single leaf clearly dominates and the vote spread is small, collapse the LCA to that one leaf (forces a confident species call). Flags: `--best-leaf-threshold [FLOAT]`, `--best-leaf-max-votes [INT]`.
2. **Cross-tree score normalization** — new (exists on no branch). Normalize each candidate's accumulated log-posterior by its **informative-position count** before the cross-tree winner/vote selection, removing the MSA-length bias that makes short-MSA trees win. Flag: `--normalize-scores`.

Both default OFF, so the binary's default behavior stays byte-identical to today. The point is to enable a grid search over `--best-leaf-threshold × --best-leaf-max-votes × --normalize-scores × --Cinterval × -u(score_constant)`.

## Current State Analysis

- `high-perf-aligner` is a deliberately stripped-down reimplementation: it has the minimap2 aligner and the runtime `--max-leaf-matches` cap, but **no** best-leaf override, adaptive-cinterval, soft-voting, or score normalization (`grep -rn "best_leaf\|adaptive\|normaliz\|override" tronko-assign/*.c *.h` → no matches).
- Benchmark evidence (`~/assignment-tool-benchmarking`): on the vert12S grid, minimap2's accuracy edge over BWA (+0.042 F1) is **unlocked by the best-leaf override** (`bl-2_mv5`: 0.6139 vs no-override 0.5726 ≈ BWA). Parameter importance: **LCA Cutoff/Cinterval 0.237 (#1), score_constant 0.205 (#2), override 0.159 (#3)**, … `max_leaf_matches` 0.0035 (last). Score normalization is the recommended root-cause fix for AC over-classification (MSA-length bias `score ≈ −2.84·positions − 52.3`; displaces ~9.2% of winners overall, 27% on traced ASVs) but was never implemented.
- The two placement entry points exist on both branches with compatible structure: `place_paired()` (WFA/minimap2 path) and `place_paired_with_nw()` (NW path) in `tronko-assign/placement.c`.

### Key Discoveries

- **Voting loop (where override hooks)** — `placement.c:927-937` (WFA `place_paired`) and `placement.c:1783-1793` (NW `place_paired_with_nw`). Votes nodes within `maximum ± Cinterval`.
- **Cross-tree argmax (where normalization acts)** — `placement.c:889-904` (WFA) and `:1746-1760` (NW); `maximum` initialized at `:878`/`:1735`. The variable to normalize is `nodeScores[i][leaf_coordinates[i][0]][k]` (= `results->nodeScores`, `type_of_PP ***`).
- **`voteRoot` is `int **` on this branch** (it was `double **` on `optimize-tronko-build`) — the ported override block must use integer literals `0`/`1`/`>0`, not `0.0`/`1.0`/`>0.0`. `place_paired`/`place_paired_with_nw` signatures at `placement.c:60`/`:947` take `int **voteRoot`.
- **Per-node score accumulation + the maskable positions** — `getscore_Arr()` at `assignment.c:172-239`. The three "uninformative" cases the normalization denominator must exclude are already branched there: `positions[i]==-1` (missing, `:188-192`), reference gap-insertion column `posteriornc[PP_IDX(positions[i],0)]==1` (`:200-204`), and query gap `locQuery[i]=='-'` (`:226-230`); the informative match cases are `:206-225`. Scores accumulate (forward `+=` reverse) in `assignScores_Arr_paired()` at `assignment.c:24-94` (`scores[...][node] += node_score` at `:40`/`:49`).
- **Flag-wiring template** — the existing `max_leaf_matches`/`aligner`/`minimap2_kmer` options were added end-to-end the same way; mirror them. `Cinterval` (`--Cinterval`/`-c` → `opt.cinterval` → global `Cinterval` at `tronko-assign.c:1172`) and `score_constant` (`-u`) are already CLI-tunable.
- **No existing informative-position counter** exists; only raw `alength` (includes gaps) and `forward_mismatch`/`reverse_mismatch`. One must be added.

## Desired End State

`tronko-assign` accepts three new flags, all default-OFF:
- `--best-leaf-threshold [FLOAT]` (default `0.0` = disabled), `--best-leaf-max-votes [INT]` (default `0` = disabled).
- `--normalize-scores` (default disabled).

With none of them set, output is **byte-identical** to the current binary (verified via the single_tree golden regressions). With them set, they alter placement as designed and compose cleanly for a grid search. Verification: golden diffs at default; functional runs of each flag on the 12S MiFish AC db without crash; a demonstrated winner-change under `--normalize-scores` on a known length-biased case.

## What We're NOT Doing

- **Adaptive-cinterval** (`--adaptive-cinterval`/`--adaptive-gap-scale`) — out of scope (importance ~0.006; not requested). Noted only because it lives adjacent in the same functions.
- **Soft-voting** — never implemented anywhere; no evidence basis.
- **Changing any default** — `--max-leaf-matches` stays 100, aligner stays `bwa`, minimap2 k/w stay as set; the three new flags default to off.
- **Adding/altering `--Cinterval` or `-u` code** — already CLI flags; we only document co-sweeping them.
- **Auto-coupling Cinterval to normalization in code** — kept orthogonal so the grid search can explore the interaction.

## Implementation Approach

Two independent, default-off features, each behind flags, each verified to leave default behavior byte-identical. Phase 1 (override) is a faithful port with one `int`-vs-`double` adaptation. Phase 2 (normalization) is approach **A (faithful, per-node)**: `getscore_Arr()` also returns the per-call informative-position count; it is accumulated into a parallel `informativeCounts[match][root][node]` array (mirroring `nodeScores`); the cross-tree argmax and vote loops divide by it when `--normalize-scores` is on. Raw scores are left untouched (so early-termination/pruning, which assume raw sums, are unaffected). Phase 3 verifies and documents grid-search readiness.

---

## Phase 1: Best-Leaf Override

### Overview
Port `--best-leaf-threshold` / `--best-leaf-max-votes` and the override block into both placement functions. Default-off (`max_votes=0` gate) ⇒ no behavior change.

### Changes Required

#### 1. `tronko-assign/global.h`
**Changes**: Add fields to both structs (mirror `max_leaf_matches`).
- `Options` struct (after `int max_leaf_matches;`, ~`:263`):
```c
double best_leaf_threshold; // Best-leaf override score threshold (default: 0 = disabled)
int    best_leaf_max_votes; // Max total votes for best-leaf override (default: 0 = disabled)
```
- `mystruct` (after `int max_leaf_matches;`, ~`:182`):
```c
type_of_PP best_leaf_threshold;
int        best_leaf_max_votes;
```

#### 2. `tronko-assign/options.c`
- `long_options[]` (after `{"max-bwa-matches",...}`, ~`:52`):
```c
{"best-leaf-threshold",required_argument,0,0},
{"best-leaf-max-votes",required_argument,0,0},
```
- `case 0:` parse chain (after the `max-leaf-matches` branch, ~`:182`):
```c
else if (strcmp(long_options[option_index].name, "best-leaf-threshold") == 0) {
    if (sscanf(optarg, "%lf", &(opt->best_leaf_threshold)) != 1) {
        fprintf(stderr, "Invalid best-leaf-threshold value; disabling\n");
        opt->best_leaf_threshold = 0.0;
    }
}
else if (strcmp(long_options[option_index].name, "best-leaf-max-votes") == 0) {
    if (sscanf(optarg, "%d", &(opt->best_leaf_max_votes)) != 1 || opt->best_leaf_max_votes < 0) {
        fprintf(stderr, "Invalid best-leaf-max-votes value (must be >= 0); disabling\n");
        opt->best_leaf_max_votes = 0;
    }
}
```
- `usage[]` help (after the `--max-leaf-matches` line):
```
--best-leaf-threshold [FLOAT], Best-leaf override score threshold [default: 0 = disabled]\n\
--best-leaf-max-votes [INT], Max total votes for best-leaf override [default: 0 = disabled]\n\
```

#### 3. `tronko-assign/tronko-assign.c`
- Defaults block (after `opt.max_leaf_matches = MAX_NUM_LEAF_MATCHES;`, ~`:946`):
```c
opt.best_leaf_threshold = 0.0;
opt.best_leaf_max_votes = 0;
```
- Worker locals (after `int max_leaf_matches = mstr->max_leaf_matches;`, ~`:289`):
```c
type_of_PP best_leaf_threshold = mstr->best_leaf_threshold;
int        best_leaf_max_votes = mstr->best_leaf_max_votes;
```
- `mstr[i]` plumbing in BOTH the paired (~`:1310`) and single-end (~`:1638`) setup loops (after `mstr[i].max_leaf_matches = opt.max_leaf_matches;`):
```c
mstr[i].best_leaf_threshold = opt.best_leaf_threshold;
mstr[i].best_leaf_max_votes = opt.best_leaf_max_votes;
```
- The 4 `place_paired(...)` / `place_paired_with_nw(...)` call sites (~`:602, 604, 608, 610`): append `, best_leaf_threshold, best_leaf_max_votes` after the trailing `pruning_factor` argument.

#### 4. `tronko-assign/placement.c` (+ `placement.h` prototypes)
- Add `, type_of_PP best_leaf_threshold, int best_leaf_max_votes` to the signatures of `place_paired` (`:60`) and `place_paired_with_nw` (`:947`), and to their prototypes in `placement.h`.
- Insert the override block **after the WFA voting loop closes at `:937`** (before `minimum_score[0] = maximum;`) and **after the NW voting loop at `:1793`** (before `:1796`). Block (note **integer** literals for this branch's `int **voteRoot`):
```c
/* Best-leaf override: if the best-scoring leaf exceeds threshold and total
   votes are below max, collapse LCA to that single leaf. Default off (max_votes=0). */
if (best_leaf_max_votes > 0 && number_of_matches > 0) {
    type_of_PP bl_best = -9999999999999999;
    int bl_node = -1, bl_tree = -1, bl_total = 0;
    for (i = 0; i < number_of_matches; i++) {
        int t = leaf_coordinates[i][0];
        if (t == -1) continue;
        int nn = 2*numspecArr[t]-1;
        for (k = 0; k < nn; k++) {
            if (treeArr[t][k].up[0] == -1 && treeArr[t][k].up[1] == -1) {
                if (nodeScores[i][t][k] > bl_best) { bl_best = nodeScores[i][t][k]; bl_node = k; bl_tree = t; }
            }
            if (voteRoot[t][k] > 0) bl_total++;
        }
    }
    if (bl_node >= 0 && bl_best > best_leaf_threshold && bl_total < best_leaf_max_votes) {
        for (i = 0; i < number_of_matches; i++) {
            int t = leaf_coordinates[i][0];
            if (t == -1) continue;
            int nn = 2*numspecArr[t]-1;
            for (k = 0; k < nn; k++) voteRoot[t][k] = 0;
        }
        voteRoot[bl_tree][bl_node] = 1;
    }
}
```

### Success Criteria

#### Automated Verification
- [x] Clean build: `cd tronko-assign && make clean && make` (no new errors/warnings).
- [x] Debug build: `make debug`.
- [x] Help shows flags: `./tronko-assign -h 2>&1 | grep -q -- '--best-leaf-threshold'` and `--best-leaf-max-votes`.
- [ ] **Default-off regression (byte-identical):** single_tree single-end and paired golden runs with no new flags match the committed `example_datasets/single_tree/missingreads_*_results.txt` exactly (`diff` empty).

#### Manual Verification
- [ ] On the 12S MiFish AC db, `--best-leaf-threshold -3 --best-leaf-max-votes 10` runs to completion (no crash) and produces more species-level calls than default (override forces commitment).
- [x] Invalid input (`--best-leaf-max-votes -1`) prints the validation warning and disables the feature.

---

## Phase 2: Cross-Tree Score Normalization (approach A — faithful, per-node)

### Overview
Add `--normalize-scores`. Track an informative-position count per accumulated node score; when enabled, compare `nodeScores / informativeCounts` at the cross-tree argmax and vote loops. Raw scores unchanged ⇒ early-termination/pruning unaffected; default-off ⇒ no behavior change.

### Changes Required

#### 1. `tronko-assign/assignment.c` — emit the informative count
- `getscore_Arr()` (`:172-239`): add a trailing out-param `int *informative_count` (callers that don't need it pass `NULL`). Initialize `*informative_count = 0` if non-NULL; increment it **only in the informative match branches (`:206-225`)** — i.e. not for `positions[i]==-1`, not for the `posteriornc[...0]==1` gap-insertion column, not for query gaps (`locQuery[i]=='-'`), not for both-gap. Leave the returned `score` unchanged.
- `assignScores_Arr_paired()` (`:24-94`) and any single-end counterpart: pass an `int n_inf_node` to each `getscore_Arr(...)` call and accumulate it into a parallel array exactly where `scores[...][node] += node_score` happens (`:40` leaf, `:49` internal): `informativeCounts[search_number][rootNum][node] += n_inf_node;` (forward and reverse both `+=`, matching the score accumulation).

#### 2. `tronko-assign/global.h` — results field + option/struct flags
- `resultsStruct` (alongside `type_of_PP ***nodeScores;`): add `int ***informativeCounts;`.
- `Options` (tier-1 flag block): `int normalize_scores; // per-informative-position normalization (default: 0)`.
- `mystruct` (tier-1 block): `int normalize_scores;`.

#### 3. `tronko-assign/allocateMemoryForResults.{c,h}` — allocate/free the parallel array
- In `allocateMemForResults()`, mirror the `nodeScores` allocation exactly (it is sized `[max_leaf_matches][numberOfTrees][2*numspecArr[j]-1]`): allocate `results->informativeCounts` with the identical triple-loop, as `int`.
- In `freeMemForResults()`, free it with the matching triple-loop.
- (Signatures already carry `max_leaf_matches`, `numberOfTrees`; no new params needed.)

#### 4. `tronko-assign/tronko-assign.c` — reset + plumb the flag
- Wherever `nodeScores` is zeroed per read (~`:619`), zero `informativeCounts` the same way (so `+=` accumulation starts clean each read).
- Defaults (~`:946`): `opt.normalize_scores = 0;`.
- Worker local (~`:289`): `int normalize_scores = mstr->normalize_scores;`.
- `mstr[i]` plumbing in both setup loops: `mstr[i].normalize_scores = opt.normalize_scores;`.
- Append `, normalize_scores` to the 4 `place_paired*` call sites and to `assignScores_Arr_paired`'s call path so `informativeCounts` is available (it lives in `results`, already passed).

#### 5. `tronko-assign/options.c`
- `long_options[]`: `{"normalize-scores",no_argument,0,0},`.
- `case 0:` chain: `else if (strcmp(...,"normalize-scores")==0) { opt->normalize_scores = 1; }`.
- `usage[]`: `--normalize-scores, Normalize scores per informative position before LCA [default: off]\n\`.

#### 6. `tronko-assign/placement.c` (+ prototypes) — normalized comparison
- Add `int normalize_scores` to `place_paired`/`place_paired_with_nw` signatures and prototypes.
- In the WFA argmax (`:889-904`) and vote (`:927-937`) loops, and the NW argmax (`:1746-1760`) and vote (`:1783-1793`) loops, replace each read of `nodeScores[i][j][k]` used for *comparison* with a helper value:
```c
type_of_PP s = nodeScores[i][j][k];
if (normalize_scores) {
    int nc = informativeCounts[i][j][k];
    if (nc > 0) s = s / (type_of_PP)nc;   /* else: leave raw (guard div-by-zero) */
}
```
and compare/vote on `s` (and set `maximum` from `s`). `informativeCounts` reaches `place_paired*` via `results` (already a parameter) — confirm the results struct is in scope there; if not, pass `results->informativeCounts` as an argument alongside `nodeScores`.

### Success Criteria

#### Automated Verification
- [x] Clean build + debug build pass.
- [x] Help shows `--normalize-scores`.
- [ ] **Default-off regression (byte-identical):** single_tree single-end and paired golden runs (no `--normalize-scores`) match committed `*_results.txt` exactly.
- [x] No memory errors: `leaks --atExit -- ./tronko-assign ... --normalize-scores` reports 0 leaks on the single_tree example (the parallel array is freed).

#### Manual Verification
- [ ] On the 12S MiFish AC db, `--normalize-scores` runs without crash and **changes the assigned tree/taxonomy for at least some reads** vs default (the winner-change the report predicts).
- [ ] On a hand-checked length-biased read (short-MSA tree winning on raw score), `--normalize-scores` flips the winner toward the longer-MSA / correct tree, or collapses the LCA to a higher rank.
- [ ] `--normalize-scores --Cinterval <small>` (per grid-search guidance) produces sane assignments; a too-large Cinterval under normalization over-collapses to high ranks (expected, demonstrates the coupling).

---

## Phase 3: Verification & Grid-Search Readiness

### Overview
Confirm the features build on both arches, are default-byte-identical, compose for a sweep, and document the Cinterval coupling.

### Changes Required
- No code; verification + a short README/usage note (optional) documenting the new flags and the Cinterval coupling.

### Success Criteria

#### Automated Verification
- [ ] arm64 build (`make`) and Linux build via CI (push triggers `build`/`docker-tests`) both pass.
- [x] All three flags appear in `-h`.
- [ ] Full default-off golden regression: single_tree single-end + paired byte-identical to committed results.
- [x] All flags can be combined in one invocation without error (e.g. `--normalize-scores --best-leaf-threshold -3 --best-leaf-max-votes 10 --Cinterval 0.01 -u 0.0005`).

Implementation note (2026-05-26): a clean `HEAD` snapshot and this implementation produce byte-identical default outputs for the single_tree single-end and paired runs. Both differ from the committed `example_datasets/single_tree/missingreads_*_results.txt` files in the same way, so the committed-golden criteria remain unchecked pending golden refresh or separate investigation.

#### Manual Verification
- [ ] A small grid (e.g. `normalize∈{off,on} × Cinterval∈{0.5,0.05,0.01} × best_leaf_max_votes∈{0,5,10}`) runs end-to-end on the 12S MiFish AC db with a mock community, and scored results vary across the grid (confirming the knobs are live and independent).

---

## Grid-Search Guidance (the Cinterval coupling)

- **Scores change scale under `--normalize-scores`.** Raw = summed log-posteriors (~hundreds; cross-tree gaps ~100+). Normalized = per-informative-position averages (~single digits; cross-tree gaps ~0.01–0.1). The vote window `maximum − Cinterval` therefore needs **~÷(avg informative positions)** smaller Cinterval when normalization is on.
- **Suggested sweep ranges:**
  - `--normalize-scores`: {off, on}
  - `--Cinterval`: raw mode {0.5, 1, 3, 5}; normalized mode {0.005, 0.01, 0.05, 0.1} (the most important parameter — importance 0.237).
  - `-u` (score_constant): {0.0001, 0.0005, 0.001, 0.01} (importance 0.205).
  - `--best-leaf-threshold`: {0(off), −1, −2, −3}; `--best-leaf-max-votes`: {0(off), 5, 10}.
  - Aligner: {bwa, minimap2} (minimap2 k=11/w=3 already default).
- Keep `--max-leaf-matches` fixed (e.g. 100; it's the least-important knob, importance 0.0035).
- For AC over-classification, score against a **holdout DB** (species removed) with a **truth-depth metric** that penalizes over-classification — not in-DB per-read species match — or the cap/normalization benefits won't be visible (see `~/assignment-tool-benchmarking` harness).

## Testing Strategy

### Regression (most important)
Default-off must be byte-identical. Before changes, capture single_tree single-end + paired baselines (or use the committed golden `*_results.txt`); after each phase, re-run with no new flags and `diff` (empty).

### Functional (AC db)
On `12S_MiFish_U/lca/ac/default/sp0.05` (pulled locally) with the trimmed MiFish mock community: run each flag alone and combined; confirm no crash, sane assignments, and that `--normalize-scores` and the override visibly change calls.

### Edge cases
- `--best-leaf-max-votes 0` / no flags ⇒ override never triggers (gate).
- `informativeCounts == 0` for a match ⇒ normalization falls back to raw (no div-by-zero).
- Invalid flag values ⇒ validation warning + disable.

## Performance Considerations

- `informativeCounts` is an `int` array parallel to `nodeScores` (`type_of_PP`): ~+50% of the `nodeScores` footprint. On the 11k-tree AC db at `--max-leaf-matches 100` this is on the order of +100 MB (acceptable; the run already uses ~4.4 GB). Default-off does not avoid the allocation — consider only allocating/zeroing it when `opt.normalize_scores` is set, to avoid the memory cost when unused.
- Per-node normalization is one extra int load + a divide at the comparison loops only when the flag is on; negligible.
- Override is an O(matches × nodes) scan per read, only when `best_leaf_max_votes > 0`.

## Migration Notes

No database/format change; existing DBs and the default invocation are unaffected (all three flags default off). Callers opt in explicitly.

## References

- Best-leaf override reference implementation: `optimize-tronko-build:tronko-assign/placement.c:945-1004` (WFA) and `:1905-1922` (NW); flags/struct/plumbing in `options.c`, `global.h`, `tronko-assign.c` on that branch.
- Score-normalization design: `~/assignment-tool-benchmarking/projects/assignment_benchmarks/findings/tronko-db-investigation-report.html` (`:697-763` formula/masks/regression, `:894-898` 9.2% displacement, `:1533-1536` "~20 lines in assignment.c"); `tronko-real-data-investigation.md`; `tronko-db-investigation-analysis-notes.md`.
- Target hook points (high-perf-aligner): `placement.c:889-904`/`:1746-1760` (argmax), `:927-937`/`:1783-1793` (vote), `assignment.c:172-239` (`getscore_Arr`), `:24-94` (`assignScores_Arr_paired`).
- Parameter importance / benchmark context: this session's reconciliation of the vert12S grid (`notebooks/archive/vert12s_v1/tronko_grid_analysis.ipynb`) and the AC investigation.
- Flag-wiring template: the existing `max_leaf_matches`/`aligner`/`minimap2_kmer` options.
