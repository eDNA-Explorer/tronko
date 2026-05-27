# Add minimap2 aligner (new default), score normalization & best-leaf override, with benchmark-tuned assign defaults [edn-976]

## Summary

This PR reworks the seeding and LCA-voting path of `tronko-assign` and ships a new set of **benchmark-validated default settings**. The motivating problem is specific to AncestralClust (AC) reference databases — the ~11k–15k tiny per-cluster trees we run in production and cannot repartition. On those DBs the old BWA-only path, with a hard-coded candidate cap and *raw* (length-unnormalized) cross-tree scoring, systematically **over-classifies**: it emits confident species calls that are wrong, because short-MSA trees win the cross-tree argmax for the wrong reasons and a single mis-seeded leaf can carry its tree.

The work was driven end-to-end by a grid benchmark (the "AC-grid") run against held-out mock communities, and every default in this PR is justified by that benchmark. The headline changes:

1. **minimap2 aligner — now the default** (`--aligner minimap2`, BWA still available via `--aligner bwa`). Drop-in seeder that fills the same `bwaMatches` contract, so placement stays aligner-agnostic.
2. **Cross-tree score normalization** (`--normalize-scores`, **default on**) — divides each candidate's accumulated log-posterior by its informative-position count before the winner/vote selection, removing the MSA-length bias that was the root cause of AC over-classification.
3. **Best-leaf override** (`--best-leaf-threshold` / `--best-leaf-max-votes`, **default on** at `-0.1` / `10`) — collapses the LCA to a single dominant leaf when the vote spread is small.
4. **Runtime `--max-leaf-matches` cap** (alias `--max-bwa-matches`) — replaces the compile-time `MAX_NUM_BWA_MATCHES` constant; honored by both aligners.
5. **New default LCA knobs**: `--Cinterval 0.02` (was an integer `5`) and `--score-constant 0.0001` (was `0.01`), tuned for normalized scoring.
6. **BWA segfault fix** on real multi-tree DBs (mate-suffix read-name mismatch + unbounded sync loop).
7. **Clean, pinned minimap2 v2.30** vendored source + CodeQL/CI hardening.

> **Note on the title / scope:** the original title ("default 10→100") described an early hypothesis. The full benchmark reversed that decision — `--max-leaf-matches` is the *least* important knob and the default landed back at **10**. The real wins came from minimap2, score normalization, and the override. The PR has been retitled accordingly.

---

## Why — the experiments that justify these changes

### AC-grid benchmark (the basis for every default)

Source: `assignment-tool-benchmarking` @ `55e55a39e`, report `projects/assignment_benchmarks/findings/tronko_ac_grid_settings_report.md`.
Scope: **15,120 assignments = 1,008 `tronko-assign` configurations × 3 mock communities × 5 LCA markers.** Configs were scored against **holdout databases** (target species removed) with a truth-depth objective that *penalizes over-classification* — the only way the normalization/cap benefits become visible (in-DB per-read species match hides them).

Per-marker best objective and winning config:

| Marker | Best objective | Best config |
|---|---:|---|
| `12s_mifish_lca` | 0.5491 | `minimap2 --normalize-scores --Cinterval 0.05 --score-constant 0.00005 --max-leaf-matches 10 --best-leaf-threshold -0.1 --best-leaf-max-votes 10` |
| `16smamm_lca` | 0.5327 | `minimap2 --normalize-scores --Cinterval 0.05 --score-constant 0.00005 --max-leaf-matches 10 --best-leaf-threshold -0.1 --best-leaf-max-votes 5` |
| `its2_plants_lca` | 0.4913 | `minimap2 --normalize-scores --Cinterval 0.02 --score-constant 0.00005 --max-leaf-matches 10 --best-leaf-threshold -0.02 --best-leaf-max-votes 10` |
| `vert12s_lca` | 0.4512 | `minimap2 --normalize-scores --Cinterval 0.005 --score-constant 0.0001 --max-leaf-matches 10` |
| `18s_euk_lca` | 0.3282 | `bwa --normalize-scores --Cinterval 0.01 --score-constant 0.001 --max-leaf-matches 10 --best-leaf-threshold -0.1 --best-leaf-max-votes 10` |

Key takeaways that set the shipped defaults:

- **The single common default reaches 95.9% of each marker-specific best on average** — i.e. one config is nearly as good as per-marker tuning, so we ship it as the default.
- **minimap2 wins decisively on seeding.** Controlled paired comparisons gave minimap2 **504/504 paired wins on both 16Smamm and 12S_MiFish**, and it was the best aligner on 4 of 5 markers (18S is the only one where BWA wins). → minimap2 is the global default.
- **Score normalization helped every single marker.** → on by default.
- **`--max-leaf-matches 10` was faster *and* marginally more accurate than 100** across the grid. → default stays at 10 (see below).
- **Best-leaf override at `-0.1` / `10` votes was net-positive or neutral** for the common config. → on by default.

### Parameter-importance ranking (vert12S grid)

A feature-importance analysis over the vert12S grid ranked the knobs:

| Rank | Parameter | Importance |
|---:|---|---:|
| 1 | LCA cutoff / `--Cinterval` | **0.237** |
| 2 | `--score-constant` (`-u`) | **0.205** |
| 3 | best-leaf override | **0.159** |
| … | … | … |
| last | `--max-leaf-matches` | **0.0035** |

So tuning effort (and the new defaults) concentrate on Cinterval and score-constant; the candidate cap is effectively irrelevant — which is why the early "raise it to 100" plan was abandoned.

On the same grid, **minimap2's accuracy edge over BWA (+0.042 F1) is unlocked specifically by the best-leaf override**: the `bl-2_mv5` override config scored **0.6139 F1 vs 0.5726** for the no-override config (≈ BWA-level). Without the override, minimap2 ≈ BWA.

### Root cause that score normalization fixes

The AC over-classification was traced to an **MSA-length bias** in raw cross-tree scoring: a candidate's summed log-posterior scales with alignment length, fit by `score ≈ −2.84·positions − 52.3`. This lets short-MSA trees win the cross-tree argmax regardless of fit. Normalizing per informative position **changes the cross-tree winner for ~9.2% of reads overall, and ~27% of traced ASVs** — pushing those calls toward the correct (longer-MSA) tree or honestly collapsing them to a higher rank.

### Why `--max-leaf-matches` stayed at 10

An earlier isolated analysis (April 2026) found that at a cap of 10 reads map to only **3–5 unique trees**, and raising 10→100 cut wrong species calls by **93%** *in isolation*. That motivated the original PR title. But once score normalization and the best-leaf override were added, they addressed over-classification far more directly; the full AC-grid then showed the cap is the least-important parameter (0.0035) and that **10 is both faster and marginally more accurate than 100**. Net: the cap became a runtime knob (good for ad-hoc experiments) but the default reverted to 10.

---

## New default configuration

The shipped defaults are equivalent to invoking:

```bash
tronko-assign --aligner minimap2 --normalize-scores --Cinterval 0.02 \
  --score-constant 0.0001 --max-leaf-matches 10 \
  --best-leaf-threshold -0.1 --best-leaf-max-votes 10
```

(Documented in `README.md` under "Recommended defaults & evidence".)

## New / changed CLI flags

| Flag (alias) | Default | Notes |
|---|---|---|
| `--aligner [bwa\|minimap2]` | `minimap2` | seeding backend; BWA still fully supported |
| `--minimap2-kmer [INT]` | `11` | minimap2 k-mer size |
| `--minimap2-window [INT]` | `3` | minimap2 minimizer window |
| `--max-leaf-matches [INT]` (`--max-bwa-matches`) | `10` | runtime candidate-leaf cap (was compile-time) |
| `--normalize-scores` / `--no-normalize-scores` | **on** | per-informative-position normalization before LCA; pass `--no-normalize-scores` to use raw summed scores |
| `--best-leaf-threshold [FLOAT]` | `-0.1` | best-leaf override score threshold |
| `--best-leaf-max-votes [INT]` | `10` | best-leaf override vote ceiling |
| `-c` / `--Cinterval [FLOAT]` | `0.02` | **was integer `5`** — now a float voting-window width |
| `-u` / `--score-constant [FLOAT]` | `0.0001` | **was `0.01`** |

All new options validate input and fall back to their default (or `exit(1)` for an unknown `--aligner`).

---

## ⚠️ Behavior changes & migration

- **The default invocation now produces different (more conservative, more correct) output than before.** This is intentional and is the point of the PR. With minimap2 + normalization + override on and the new Cinterval/score-constant, default assignments change versus the old binary, especially on AC databases where over-classification is curbed.
- **To reproduce legacy behavior**, pin the old aligner and knobs: `--aligner bwa` (and set `--Cinterval`/`-u` to your prior values; disable the new features with `--best-leaf-max-votes 0` and `--no-normalize-scores`). Production configs that already set these explicitly are only affected by the aligner default.
- **`-c` is now a float.** Callers that passed `-c 5` (the old integer LCA cut-off) should move to the float Cinterval scale (default `0.02`); a literal `5` is now an enormous window and will over-collapse to high ranks.
- **No database/format migration.** Existing tronko-build databases work unchanged. minimap2 builds its index in-memory from the `-a` reference FASTA at runtime — **no on-disk index files are created** in minimap2 mode (BWA `.bwt/.sa/.pac/.ann/.amb` are only written under `--aligner bwa`).
- **Paired + minimap2** auto-enables read-2 reverse-complement (equivalent to `-z`), since the minimap2 wrapper does not self-handle strand the way BWA does. Logged at startup.

---

## Component-level changes

**Diffstat:** 63 files, +13,754 / −161. The vast majority of the additions (~13k lines) are the **vendored, pinned minimap2 v2.30 upstream source** (`tronko-assign/minimap2_src/`, incl. `sse2neon` for arm64). The hand-written changes are concentrated in ~10 files:

- `minimap2_wrapper.{c,h}`, `minimap2_shim.c` — new wrapper that mirrors `run_bwa`'s contract (fills `bwaMatches`), with a serialized lazy index cache (thread-safe across workers) and BWA-compatible CIGAR/concordance handling.
- `tronko-assign/Makefile` — builds minimap2 as a static `libminimap2.a` with arch dispatch (sse2neon on arm64, `-msse2`/`ksw2_dispatch` on x86_64); symbol renames to avoid colliding with BWA.
- `placement.c` / `placement.h` — best-leaf override block in both `place_paired` (WFA) and `place_paired_with_nw` (NW); normalized comparison in the cross-tree argmax and vote loops, gated on `--normalize-scores`.
- `assignment.c`, `allocateMemoryForResults.{c,h}` — emit/accumulate a parallel `informativeCounts` array alongside `nodeScores`; allocate/free it.
- `bwa_source_files/fastmap.c`, `bwa_source_files_include.h` — runtime cap threaded through `main_mem`; **segfault fix** (bound the read-sync loops; `if (!leaf_map) continue;` null-guards).
- `readreference.c` — strip trailing `/<digit>` mate suffixes from query names in all four readers so they match BWA's `trim_readno` (the other half of the segfault fix); no-op for the legacy `_1`/`_2` convention.
- `options.c`, `global.h`, `tronko-assign.c` — new option parsing, `Options`/`mystruct` fields, defaults, and worker plumbing.
- `README.md` — documents the new flags and the "Recommended defaults & evidence" section.
- CI: CodeQL advanced workflow scanning C/C++ + actions, vendored trees excluded from analysis; fix for pre-existing CI failures (Unity submodule + Docker zstd).

---

## Verification

Build and functional checks run locally on this branch (Ubuntu 24.04, x86_64):

- [x] **Clean build** — `cd tronko-assign && make clean && make` succeeds (exit 0). Only pre-existing `crash_debug.c` format-truncation warnings; no new warnings from this work.
- [x] **Help advertises all new flags** — `--aligner`, `--minimap2-kmer`, `--minimap2-window`, `--max-leaf-matches`, `--best-leaf-threshold`, `--best-leaf-max-votes`, `--normalize-scores` all present in `-h`.
- [x] **Default (minimap2 + normalize + override) single-end run** on the `single_tree` example completes (exit 0); 164/164 reads assigned; scores are on the normalized per-position scale (e.g. `-0.34`, `-0.14`).
- [x] **`--aligner bwa` run** on the same example completes (exit 0); taxonomic paths are identical to the minimap2 run on this single-tree fixture (the two seeders agree where the DB is unambiguous).
- [ ] **AC-database functional check** (12S MiFish AC db with a mock community) — exercised during the benchmark, not re-run in this PR's CI; recommend a spot-check before merge. The single-tree example fixtures are too small to externally expose the cap/normalization/override effects (they only manifest on the many-tree AC DBs).
- [ ] **CI** (`build` + `docker-tests` + CodeQL) — will run on push; confirm green before merge.

> **Example outputs:** the committed `example_datasets/single_tree/missingreads_*_results.txt` goldens have been **refreshed** in this PR (commit `5dd1912`) to match the new defaults — taxonomic assignments are unchanged; only the `Score` (now normalized) and `Node_Number` columns differ. The `example_datasets/multiple_trees/*` goldens are left as-is: regenerating them requires rebuilding the partition / multi-MSA databases, which re-estimates trees non-deterministically (and one mode currently errors on a missing `list.txt`). They depend on `tronko-build` (untouched here) and no test/CI consumes these committed result files, so this is a documentation artifact, not a regression.

---

## Updates since first review

- **Added `--no-normalize-scores`** (commit `4ac2a1c`): normalization defaults on, so an explicit off-switch is required for the documented legacy-reproduction path and for raw-vs-normalized comparison. Verified: default output is byte-identical to the goldens; `--no-normalize-scores` reverts to raw summed scores.
- **CodeQL** (now passing, 0 alerts). Two root causes were addressed:
  - *Vendored code:* `paths-ignore` does not filter alerts for compiled C/C++, so the vendored minimap2 v2.30 sources are instead kept out of the CodeQL database entirely — `libminimap2.a` is pre-built *before* CodeQL tracing starts, so the traced `make` reuses it and never recompiles `minimap2_src/*.c`.
  - *Query suite:* the workflow had opted into `security-and-quality` (the noisiest suite), whose medium-precision `cpp/missing-check-scanf` flagged ~70 pre-existing unchecked `scanf`/`fscanf` sites across the whole module — including files this branch never touches — all counted as "new" because `main` has no baseline with that suite. Reverted to the default high-precision `code-scanning` suite, which still covers real high-confidence vulnerabilities.
  - As genuine hardening (kept regardless), the two reference-DB parses in `readreference.c` now check their `sscanf` return inline and bound parsed tree dimensions (`MAX_TREE_LEAVES`); no behavior change on valid databases (goldens byte-identical).
- Refreshed stale default-value comments in `global.h`. (Greptile's "freeMemForResults param order" note does not reproduce — header and implementation already agree.)

## Changelog

`tronko-assign`: add a minimap2 aligner (new default) alongside BWA; add cross-tree score normalization and a best-leaf override (both on by default); make the candidate-leaf cap a runtime flag (`--max-leaf-matches`); retune default `--Cinterval` (0.02) and `--score-constant` (0.0001); fix a BWA segfault on multi-tree databases. Defaults are tuned from a 15,120-assignment AC-grid benchmark and curb over-classification on AncestralClust databases. Existing databases work unchanged; default output changes intentionally (pin `--aligner bwa` + legacy knobs to restore old behavior).
