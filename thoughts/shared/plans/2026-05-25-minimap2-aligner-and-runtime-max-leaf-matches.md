# minimap2 Aligner + Runtime max-leaf-matches Implementation Plan

## Overview

Add a minimap2 aligner option and a runtime `--max-leaf-matches` cap to `tronko-assign` on the `high-perf-aligner` branch, and **raise the default cap from 10 to 100**. Both are classify-time changes only — no database/tree rebuild required.

**Why a higher default:** we are locked to the AncestralClust (AC) database (~15k tiny trees) and cannot repartition it, so the aligner cap is the only lever against its over-classification failure mode. On the AC db, a cap of 10 maps reads to only 3–5 unique trees that often lack the correct species, producing confidently-wrong species calls; 10→100 cut wrong species calls 93% (analysis April 2026). The synthetic benchmark's preference for a cap of 10 is the known artifact that investigation exposed — its holdouts contain the species in-db, so the correct genus-level abstention that a higher cap produces *looks* like lower F1 but is the honest call on real data. 100 is the validated sweet spot; 500+ only adds diminishing returns and a pathological call. The work is **reimplemented** against this branch's current code, using the existing implementation on `optimize-tronko-build` (commits `c0a6538`, `4aa29fc`) only as a reference. We do **not** pull files from that branch, because it carries many unrelated changes (best-leaf override, soft-voting, adaptive-cinterval, the speedup refactor, ablation system, trace-read, etc.).

Score normalization is explicitly **deferred** to a later effort (see `thoughts/shared/research/2026-04-10-cross-tree-score-normalization.md`).

## Current State Analysis

- `tronko-assign` on `high-perf-aligner` uses **BWA only**. The candidate-leaf cap is a compile-time constant `MAX_NUM_BWA_MATCHES = 10` (`tronko-assign/global.h:43`), used in ~30 places: per-read buffer allocations (`tronko-assign.c:319-329`), init/free loops, ~12 `leaf_iter < MAX_NUM_BWA_MATCHES` voting guards, and the BWA seeding loops in `bwa_source_files/fastmap.c`.
- Prior analysis (`assignment-tool-benchmarking` findings, April 2026) showed the cap of 10 is the dominant driver of over-confident misassignment: BWA/minimap2 seeds map to only 3–5 unique trees, often not containing the correct species, so a wrong leaf wins its tree by default. Raising the cap (to ~50–100) lets LCA voting collapse ambiguous calls correctly.
- minimap2 integration on `optimize-tronko-build` is a drop-in alternative to BWA: `run_minimap2()` has the **same contract** as `run_bwa()` — it fills the same `bwaMatches` struct (`global.h:188-202`), so `placement.c`/`assignment.c` are aligner-agnostic and need no changes.
- `high-perf-aligner` is an ancestor of `optimize-tronko-build`; the integration boundary (`bwaMatches`) and the `leafMap` hashmap pattern (`fastmap.c:111-123`) are byte-identical across both branches, so the reference maps cleanly.

## Desired End State

`tronko-assign` accepts:
- `--max-leaf-matches [INT]` (**default 100**, alias `--max-bwa-matches`) — runtime cap on candidate leaf matches per read, honored by both aligners.
- `--aligner [bwa|minimap2]` (default `bwa`), `--minimap2-kmer [INT]` (default 15), `--minimap2-window [INT]` (default 5).

With `--aligner bwa` (default) and the cap pinned to the old value (`--max-leaf-matches 10`), output is **byte-identical** to the current binary (this is how we verify the refactor; see Success Criteria). At the new default of 100, the BWA path intentionally changes — more candidate trees are scored, which is the over-classification fix. With `--aligner minimap2`, reads are seeded by minimap2 (building its own in-memory index from the `-a` reference FASTA; BWA index files are not created), and produce sane assignments on the example datasets.

### Key Discoveries

- **Integration boundary** is `bwaMatches` (`global.h:188-202`): `concordant_matches_{roots,nodes}` and `discordant_matches_{roots,nodes}` are parallel, `-1`-terminated int arrays of (tree_index, leaf_node_index), deduped, capped; `cigars_forward/starts_forward` (+ reverse for the mate) filled only when `use_portion==1` with BWA-style CIGAR (`=`/`X`→`M`) and 1-based starts.
- **minimap2 index target** is the `-a` reference FASTA (`mstr->databasefile` = `opt.fasta_file`), the same file BWA indexes — *not* `reference_tree.txt`. `readreference.c` differs by only 2 lines between branches, so the DB format is unchanged and existing databases work as-is.
- **Dispatch is a single site** (`tronko-assign.c:368`, inside the one worker used for both paired and single-end), plus two BWA-index-skip guards (`:1193`, `:1514`) and `mstr` plumbing (after `:1300`, `:1623`).
- **Allocation/cap consistency**: the same `opt.max_leaf_matches` value sizes the mallocs *and* bounds every loop/guard, so there is no "request more than allocated" gap. `MAX_NUM_LEAF_MATCHES` survives only as the default value and invalid-input fallback.
- **minimap2 build is arch-sensitive**: arm64 (this Mac) needs the `sse2neon` path (`-DKSW_SSE2_ONLY -D__SSE2__`); x86_64 needs `-msse2` + `ksw2_dispatch.c`. minimap2 is built as a static `libminimap2.a` (excluding `main.c`/`kthread.c`/`misc.c`, whose symbols collide with BWA); `minimap2_shim.c` supplies the few symbols `misc.c` would have provided.
- **Index cache thread-safety**: the reference wrapper caches the index in file-scope statics built lazily on first `run_minimap2` call — a data race across worker threads. We fix this with a serialized lazy build (mutex/`pthread_once`-style).

## What We're NOT Doing

- **Score normalization** — deferred to a separate plan.
- **Best-leaf override** (`--best-leaf-threshold`, `--best-leaf-max-votes`).
- **Adaptive cinterval**, **trace-read** (`--trace-read`, `g_trace_read`).
- **Soft-voting** (`voteRoot` `int**`→`double**`, `soft_voting`/`vote_temperature`) and the **dead voting-param removal** cleanup.
- The **`getSequenceinRoot.c` crash guard** and other defensive null-checks bundled into `4aa29fc`/`ababc0e` (port separately if crash-hardening is wanted; not required here).
- Any **`placement.c`/`assignment.c` signature or scoring change** — the `bwaMatches` contract keeps placement aligner-agnostic.
- Changing the **default aligner**: stays `aligner=bwa`. (Only the cap default changes, from 10 to 100; minimap2 remains opt-in.)

## Implementation Approach

Two phases, each independently buildable and testable. Phase 1 (runtime cap) stands alone, benefits the existing BWA path, and is a prerequisite for minimap2's `mopt.best_n` cap. Phase 2 (minimap2) layers the new aligner on top. Naming follows the `optimize-tronko-build` tip: constant `MAX_NUM_LEAF_MATCHES` (now **100**), field/flag `max_leaf_matches` / `--max-leaf-matches`, with `--max-bwa-matches` as a back-compat alias. Refactor correctness is verified separately from the default change by pinning `--max-leaf-matches 10` and diffing against the pre-change baseline.

---

## Phase 1: Runtime `--max-leaf-matches`

### Overview

Convert the compile-time `MAX_NUM_BWA_MATCHES` cap into a runtime `opt.max_leaf_matches` (default **100**), threaded through the worker struct, `run_bwa`→`main_mem`, the allocation/loop sites, and every voting guard. ~9 files. The refactor itself is behavior-preserving at equal cap; the new default of 100 intentionally raises candidate coverage on the BWA path.

### Changes Required

#### 1. `global.h`
- `:43` — rename and raise: `#define MAX_NUM_BWA_MATCHES 10` → `#define MAX_NUM_LEAF_MATCHES 100` (this constant is now the default value + invalid-input fallback only; allocations/loops use the runtime value).
- `Options` struct — add `int max_leaf_matches;` (after `pruning_factor`).
- `mystruct` (per-thread worker, lines ~152-186) — add `int max_leaf_matches;` so the worker can read it.

#### 2. `bwa_source_files_include.h`
- Add `, int max_leaf_matches` as the last parameter of the `main_mem(...)` prototype.

#### 3. `bwa_source_files/fastmap.c`
- `ktp_aux_t` struct (~`:45`) — add `int max_leaf_matches;`.
- `:196` and `:269` — `for(k=0; k<MAX_NUM_BWA_MATCHES; k++)` → `... k<aux->max_leaf_matches ...` (concordant + discordant seeding loops).
- `main_mem(...)` signature (~`:569`) — add `, int max_leaf_matches`; body (~`:827`, by `aux.max_acc_name = ...`) — add `aux.max_leaf_matches = max_leaf_matches;`.
- (Recommended defensive cap: `&& k < aux->max_leaf_matches` on any post-loop cigar write.) **Exclude** the unrelated `if (!leaf_map) continue;` null-checks and `MEM_F_ALL` change.

#### 4. `allocateMemoryForResults.h`
- Add `, int max_leaf_matches` as last param of `allocateMemForResults(...)` and `freeMemForResults(...)`.

#### 5. `allocateMemoryForResults.c`
- Add the `int max_leaf_matches` param to both definitions; replace `MAX_NUM_BWA_MATCHES` with it at the sizing/loop sites: lines `12, 13, 30, 32, 34, 36, 38, 94, 110`. (Do **not** touch any `voteRoot` lines.)

#### 6. `tronko-assign.c`
- Defaults block (~`:934`, after `opt.pruning_factor = 2.0;`): `opt.max_leaf_matches = MAX_NUM_LEAF_MATCHES;`.
- `run_bwa` signature (`:173`) — add `, int max_leaf_matches`; its two `main_mem(...)` calls (`:177`, `:179`) — append `, max_leaf_matches`.
- Worker fn (owns `mstr`, ~`:267`) — add local `int max_leaf_matches = mstr->max_leaf_matches;`.
- `bwa_results[i]` alloc block (`:319-329`) — replace `MAX_NUM_BWA_MATCHES` with `max_leaf_matches` in all 8 mallocs; init loops (`:332`, `:344`) likewise.
- `run_bwa(...)` dispatch call (`:368`) — append `, max_leaf_matches`.
- 12 voting guards `leaf_iter < MAX_NUM_BWA_MATCHES` → `... < max_leaf_matches`: lines `397, 418, 423, 437, 457, 462, 479, 499, 504, 522, 542, 547`.
- Crash/warn sites (`:576`, `:581`) and free loop (`:850`) → `max_leaf_matches`.
- `allocateMemForResults(...)` calls (4 sites) and `freeMemForResults(...)` calls (4 sites) — append `, opt.max_leaf_matches`.
- `mstr[i].max_leaf_matches = opt.max_leaf_matches;` in both setup loops (single-end + paired).
- End-of-run overflow summary warnings (`:1445-1446`, `:1713-1714`) — point at `opt.max_leaf_matches`.

#### 7. `options.c`
- `long_options[]` (after `pruning-factor`): `{"max-leaf-matches",required_argument,0,0}` and `{"max-bwa-matches",required_argument,0,0}` (alias).
- `usage[]`: `--max-leaf-matches [INT], Maximum leaf matches per read [default: 10] (alias: --max-bwa-matches)`.
- `parse_options` `case 0:` chain — add a branch matching either name: `sscanf("%d", &opt->max_leaf_matches)`, require `>= 1`, else warn and reset to `MAX_NUM_LEAF_MATCHES`.

### Success Criteria

#### Automated Verification
- [x] Builds clean: `cd tronko-assign && make clean && make` (no warnings introduced beyond baseline).
- [x] Debug build works: `make debug`.
- [x] Help shows the flag: `./tronko-assign -h 2>&1 | grep -q -- '--max-leaf-matches'`.
- [x] BWA refactor regression — output identical at **equal cap**:
  - Save baseline from the current binary on the CLAUDE.md single-tree example **before** Phase 1; after Phase 1 run the same example with `--max-leaf-matches 10`; `diff` must be empty. (This isolates refactor correctness from the default change.)

#### Manual Verification
- [x] Default run (no cap flag) now uses 100, runs to completion without crash/leak on the single-tree example, and evaluates more candidate trees than the old binary (observed via the existing `g_max_potential_matches` warning / candidate count).
  - Verification note: local single-tree and bundled small multi-tree fixtures completed at cap 100; they did not emit dropped-match/candidate-count warnings or produce cap 10 vs 100 final-assignment differences, so the larger-candidate-tree effect was not externally observable on these fixtures.
- [x] Alias `--max-bwa-matches 50` behaves identically to `--max-leaf-matches 50`.
- [x] Invalid input (`--max-leaf-matches 0`) prints the validation warning and falls back to the default (100).
- [x] Memory at the higher default is acceptable: check peak RSS on a realistically-sized batch (the per-read `cigars_*` buffers scale as cap × `MAX_CIGAR`; confirm no blowup at 100, especially with `use_leaf_portion`).

---

## Phase 2: minimap2 Aligner

### Overview

Add `--aligner minimap2` as a drop-in seeding alternative. New vendored sources + a reimplemented wrapper that fills `bwaMatches` exactly like BWA, plus the Makefile, option fields/flags, dispatch, index-skip guards, and `mstr` plumbing.

### Changes Required

#### 1. Vendored sources (new files)
**Files**: `tronko-assign/minimap2_src/` (93 files), `tronko-assign/minimap2_shim.c`
**Changes**: These are upstream third-party / boilerplate, not "their changes" — copy verbatim from `optimize-tronko-build`:
```
git checkout optimize-tronko-build -- tronko-assign/minimap2_src tronko-assign/minimap2_shim.c
```
`minimap2_shim.c` provides `mm_verbose`, `mm_dbg_flag`, `mm_realtime0`, `peakrss`, `mm_err_*`, and radix-sort instantiations (would otherwise come from minimap2's `misc.c`, which collides with BWA's `utils.c`).

#### 2. `minimap2_wrapper.h` (new — reimplement)
**Changes**: Declare the BWA-compatible entry point:
```c
void run_minimap2(int start, int end, bwaMatches *bwa_results, int concordant,
                  int numberOfTrees, char *databaseFile, int paired,
                  int max_query_length, int max_readname_length,
                  int max_acc_name, int max_leaf_matches,
                  int mm2_kmer, int mm2_window);
```

#### 3. `minimap2_wrapper.c` (new — reimplement against high-perf)
**Changes**: Same contract as `run_bwa`; fills `bwaMatches`. Structure:
- **Serialized lazy index cache** (decision: mutex/`pthread_once`-style). File-scope `static mm_idx_t *g_mm2_idx`, `g_mm2_db_path`, `g_mm2_kmer`, `g_mm2_window`, guarded by a `static pthread_mutex_t` (double-checked): the first worker to enter builds the index from `databaseFile` (the `-a` FASTA) with `iopt.k=mm2_kmer; iopt.w=mm2_window;`; others reuse it. Rebuild only if path/k/w change.
- **Mapping opts**: `mm_set_opt(0,&iopt,&mopt); mm_set_opt("sr",&iopt,&mopt); mopt.flag |= MM_F_CIGAR; mopt.best_n = max_leaf_matches; mm_mapopt_update(&mopt, mi);`.
- **leafName→(root,node) hashmap**: build like `fastmap.c:111-123` over `treeArr`/`numspecArr` (leaf nodes `[numspecArr[t]-1, 2*numspecArr[t]-1)`), value = `leafMap{root,node}` (`global.h:96-100`). Requires the same `extern` decls (`treeArr`, `numspecArr`, query matrices).
- **Per hit**: `mi->seq[r->rid].name` → `hashmap_get` → `(root,node)`; `add_match()` finds first `-1` slot, dedupes on `(root,node)`, respects `max_leaf_matches`; when `use_portion==1`, decode minimap2's packed CIGAR (`(len<<4)|op`, mapping `=`/`X`→`M`) into `cigars_forward[k]`/`starts_forward[k]` (1-based) — reverse slot for the mate.
- **Concordant/discordant** (mirror `fastmap.c:196-339`): paired — R1 hit is concordant if any R2 hit lands in the same tree (`root` match), storing the mate CIGAR in the reverse slot; else discordant; unmatched R2 hits → discordant; honor `concordant` flag. Single-end — all hits → discordant.
- **Cleanup**: free per-read `mm_reg1_t` arrays + `.p`, destroy `mm_tbuf_t` per call; free hashmap (`leafMap` values then `hashmap_cleanup`) at end.

#### 4. `tronko-assign/Makefile`
**Changes**: Port the minimap2 build block from `optimize-tronko-build` (`git show 4aa29fc -- tronko-assign/Makefile`):
- `MINIMAP2_DIR = minimap2_src`; `MINIMAP2_CORE_SRCS` (kalloc, bseq, sketch, sdust, options, index, lchain, align, hit, seed, jump, map, format, pe, esterr, splitidx) and `MINIMAP2_KSW2_SRCS` (ksw2_ll_sse + extz2/extd2/exts2; +dispatch on x86).
- Arch dispatch via `UNAME_M`: `arm64` → `MINIMAP2_KSW2_FLAGS = ... -DKSW_SSE2_ONLY -D__SSE2__ -I$(MINIMAP2_DIR)/sse2neon`; else `-msse2` + `ksw2_dispatch.c`.
- `libminimap2.a: $(MINIMAP2_ALL_OBJS)` via `$(AR) -csru`; per-object rules using `MINIMAP2_COMMON_FLAGS` (`-DHAVE_KALLOC -I$(MINIMAP2_DIR) -w`) and the kseq symbol renames to avoid BWA collision.
- Add `minimap2_shim.c minimap2_wrapper.c` to `SOURCES`; add `libminimap2.a` + `-ldl` and `-I$(MINIMAP2_DIR) -DHAVE_KALLOC` to both the `all` and `debug` link lines; `all`/`debug` depend on `libminimap2.a`.
- Include the macOS `_XOPEN_SOURCE`/zstd-include bits as needed for Darwin builds.

#### 5. `global.h`
**Changes**: Add to **both** `Options` and `mystruct`:
```c
char aligner[16];     // "bwa" (default) or "minimap2"
int  minimap2_kmer;
int  minimap2_window;
```
(Do **not** add best-leaf/trace/adaptive fields.)

#### 6. `options.c`
**Changes**: Add `long_options[]` entries `{"aligner",...}`, `{"minimap2-kmer",...}`, `{"minimap2-window",...}`; parsing branches in `case 0:` — `aligner` must be `"bwa"`/`"minimap2"` (else error), kmer/window `>= 1`; help text.

#### 7. `tronko-assign.c`
**Changes**:
- `#include "minimap2_wrapper.h"` with the other includes.
- Defaults block (~`:935`): `strcpy(opt.aligner, "bwa"); opt.minimap2_kmer = 15; opt.minimap2_window = 5;`.
- Dispatch (`:368`): wrap the `run_bwa(...)` call:
  ```c
  if (strcmp(mstr->aligner, "minimap2") == 0) {
      run_minimap2(mstr->start, end, bwa_results, mstr->concordant, mstr->ntree,
                   mstr->databasefile, paired, max_query_length, max_readname_length,
                   max_acc_name, max_leaf_matches, mstr->minimap2_kmer, mstr->minimap2_window);
  } else {
      run_bwa(mstr->start, end, bwa_results, mstr->concordant, mstr->ntree,
              mstr->databasefile, paired, max_query_length, max_readname_length,
              max_acc_name, max_leaf_matches);
  }
  ```
- BWA-index-skip guards (`:1193`, `:1514`): add `&& strcmp(opt.aligner, "minimap2") != 0` to the `if (opt.skip_build==0 ...)` conditions wrapping `bwa_index(2, opt.fasta_file)`.
- `mstr` plumbing (after `:1300` single, `:1623` paired): copy `opt.aligner` (via `strncpy`), `opt.minimap2_kmer`, `opt.minimap2_window` into `mstr[i]`.

### Success Criteria

#### Automated Verification
- [x] minimap2 static lib builds: `cd tronko-assign && make clean && make` produces `libminimap2.a` and links `tronko-assign` (arm64 path on this Mac).
- [x] `make debug` builds.
- [x] Help shows aligner flags: `./tronko-assign -h 2>&1 | grep -q -- '--aligner'`.
- [x] Default-aligner regression: with `--aligner bwa --max-leaf-matches 10`, output still byte-identical to the original pre-Phase-1 baseline (adding minimap2 code leaves the BWA path inert unless selected).
- [x] minimap2 mode does not create BWA index files (`.bwt`/`.sa`/`.pac`/`.ann`/`.amb`) for the reference when run with `--aligner minimap2`.

#### Manual Verification
- [x] `--aligner minimap2` on the single-tree example completes and produces non-empty, sane assignments (spot-check against the BWA run — broadly concordant at higher ranks).
- [x] `--aligner minimap2 --max-leaf-matches 100` runs without crash; evaluates more trees than at cap 10.
  - Verification note: cap 10 and 100 both completed on local single-tree and small multi-tree examples; these fixtures produced identical final assignment files, so they did not expose a visible cap-dependent assignment change.
- [x] `--minimap2-kmer 11 --minimap2-window 3` changes seeding behavior without error.
  - Verification note: the option combination completed without error; final assignments were unchanged on the single-tree fixture.
- [x] Paired-end (`-p`) and single-end both work under minimap2.
- [x] No data race / crash under multithreading (`-T <ncores>`): index built once, shared safely.
- [x] `valgrind`/leak spot-check on a small run (per-read minimap2 buffers and the hashmap are freed).
  - Verification note: `valgrind` is not installed on this Mac; `/usr/bin/leaks --atExit` reported 0 leaks for the minimap2 single-tree run.

---

## Testing Strategy

### Golden BWA Regression (most important)
Capture a baseline from the current binary **before** any change, then re-run after each phase **at equal cap** (`--max-leaf-matches 10`, since the old binary hardcoded 10) and `diff`:
```bash
# baseline (current HEAD, before edits) — old binary, hardcoded cap 10
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/baseline.txt -w
# after Phase 1 / Phase 2 — pin cap to 10 to isolate refactor correctness
tronko-assign ... --max-leaf-matches 10 ... -o /tmp/after.txt -w
diff /tmp/baseline.txt /tmp/after.txt   # must be empty
```
Separately, confirm the **default** (no cap flag → 100) runs cleanly and produces *different* (more conservative) output — that change is the intended improvement, not a regression.

### minimap2 Functional
- Single-tree example with `--aligner minimap2` (merged/single and paired).
- `--max-leaf-matches 100 --aligner minimap2` on the multi-cluster example DB if available locally.

### Edge Cases
- Invalid `--aligner foo` → clean error.
- `--max-leaf-matches 0` → warning + fallback.
- Multithreaded run (`-T 4`) → no race on first minimap2 index build.
- Reads with zero minimap2 hits → handled like BWA (no candidates, no crash).

## Performance Considerations

- Raising `max_leaf_matches` increases per-read buffer allocation and the number of trees scored (~4× at 100 vs 10 per prior analysis) — acceptable, and the prior sweep showed diminishing returns past ~100.
- minimap2 index is built once and cached; build cost is amortized across all reads. Serializing the first build adds negligible overhead (one mutex acquisition per worker on first call).
- minimap2 is generally faster at seeding than BWA for short reads; net runtime should be comparable or better at equal cap.

## Migration Notes

- No database/format migration: existing tronko databases work unchanged (only the `-a` reference FASTA is re-used for the minimap2 index, built in-memory at runtime).
- Default aligner stays `bwa` (minimap2 is opt-in), but the **default cap changes 10→100**. Callers that relied on the implicit cap of 10 will see more conservative (more-correct on the AC db) assignments; pin `--max-leaf-matches 10` to restore the old behavior. Production configs that already set the cap explicitly are unaffected.

## References

- Reference implementation: `optimize-tronko-build` commits `c0a6538` (runtime accuracy params incl. max-leaf-matches) and `4aa29fc` (minimap2 aligner). Reimplemented here, not pulled.
- Integration boundary: `tronko-assign/global.h:188-202` (`bwaMatches`), `tronko-assign/bwa_source_files/fastmap.c:111-123` (leafMap hashmap) and `:196-339` (concordant/discordant fill).
- Deferred work: `thoughts/shared/research/2026-04-10-cross-tree-score-normalization.md` (score normalization), and the Cinterval-rescale analysis (normalized Cinterval ≈ raw_optimum / informative_positions, marker-dependent — to be swept, not assumed).
- Motivation: `assignment-tool-benchmarking/projects/assignment_benchmarks/findings/tronko-db-investigation-report.html` (the `max_leaf_matches` fix) and `tronko-real-data-investigation.md` (aligner/partition routing).
