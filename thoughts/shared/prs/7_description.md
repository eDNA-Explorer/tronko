# Add minimap2 aligner and runtime `--max-leaf-matches` cap (default 10→100)

Linear: **edn-976**

## Summary

Two classify-time changes to `tronko-assign`, both **opt-in or behavior-preserving at equal settings**, and **neither requires a database/tree rebuild**:

1. **Runtime `--max-leaf-matches` cap** (default raised **10 → 100**). The candidate-leaf cap was a compile-time constant `MAX_NUM_BWA_MATCHES = 10` used in ~30 sites; it is now a runtime option threaded through the worker, both aligners, and every allocation/voting site. Alias `--max-bwa-matches` is accepted for back-compat.
2. **`--aligner minimap2`** — a drop-in seeding alternative to BWA (default stays `bwa`). minimap2 fills the exact same `bwaMatches` contract, so placement/assignment are unchanged. Tunable via `--minimap2-kmer` (default 15) and `--minimap2-window` (default 5).

## Why (motivation)

We are locked to the AncestralClust (AC) reference database (~15k tiny trees) and cannot repartition it, so the aligner cap is the only lever against its over-classification failure mode. With a cap of 10, BWA/minimap2 seeds map each read to only **3–5 unique trees** that often do not contain the correct species, so a wrong leaf wins its tree by default and produces a confidently-wrong species call. Raising the cap lets LCA voting collapse ambiguous calls correctly: **10→100 cut wrong species calls 93%** (investigation, April 2026). 100 is the validated sweet spot; 500+ only adds diminishing returns and a pathological call.

minimap2 is added as a faster short-read seeding alternative that uses the same integration boundary, enabling head-to-head aligner comparison without touching the scoring/placement code.

Score normalization is explicitly **deferred** to a later effort.

## What changed

### Phase 1 — runtime `--max-leaf-matches`

- **`global.h`**: `#define MAX_NUM_BWA_MATCHES 10` → `#define MAX_NUM_LEAF_MATCHES 100` (now only the default/invalid-input fallback). Added `max_leaf_matches` (+ `aligner`, `minimap2_kmer`, `minimap2_window`) to both the `Options` and per-thread `mystruct` structs.
- **`options.c`**: new `--max-leaf-matches` (+ `--max-bwa-matches` alias) parsing with `>= 1` validation and fallback-to-default-with-warning; help text.
- **`bwa_source_files_include.h` / `bwa_source_files/fastmap.c`**: `main_mem(...)` and the `ktp_aux_t` struct gain `max_leaf_matches`; the concordant/discordant seeding loops and the second-slot writes are bounded by it instead of the constant.
- **`allocateMemoryForResults.{c,h}`**: `allocateMemForResults`/`freeMemForResults` take `max_leaf_matches`; all per-read buffer sizing and free loops use it.
- **`tronko-assign.c`**: defaults set; `run_bwa` plumbs the cap into `main_mem`; per-read `bwaMatches` allocations/init/free use it; ~12 `leaf_iter < MAX_NUM_BWA_MATCHES` voting guards switch to `< max_leaf_matches`; the overflow-summary warnings now report `opt.max_leaf_matches` and point at `--max-leaf-matches`.
- **Latent sizing mismatch fixed**: `results->leaf_coordinates` (in `allocateMemoryForResults.c`) and the stack array `trees_search` (in `tronko-assign.c`) were previously sized to `numberOfTrees` / `mstr->ntree` but indexed by `leaf_iter`/loop counters bounded by the match cap. Both are now consistently sized and iterated by `max_leaf_matches`, removing a size/index mismatch.

### Phase 2 — minimap2 aligner

- **Vendored third-party sources** (`tronko-assign/minimap2_src/`, ~93 files): upstream minimap2, built as a static `libminimap2.a` (excludes `main.c`/`kthread.c`/`misc.c`, whose symbols collide with BWA's `utils.c`).
- **`minimap2_shim.c`** (new): supplies the few symbols `misc.c` would have provided (`mm_verbose`, `mm_dbg_flag`, `mm_realtime0`, `peakrss`, `mm_err_*`, radix-sort instantiations) without the BWA collisions.
- **`minimap2_wrapper.{c,h}`** (new): `run_minimap2()` mirrors `run_bwa()`'s contract and fills `bwaMatches` (roots/nodes arrays, plus CIGAR/start slots when `use_portion==1`, mapping minimap2's `=`/`X` ops → `M`). Builds an in-memory minimap2 index from the `-a` reference FASTA. The index is cached in file-scope statics behind a **`pthread_mutex_t` (serialized lazy build)** — every worker takes the lock, the first builds the index and the rest reuse it; the cache-hit hold time is a strcmp + a few int compares, so contention is negligible. It is rebuilt only if the path/k/w change, and released via an `atexit` finalizer.
- **`tronko-assign.c`**: `#include "minimap2_wrapper.h"`; the single dispatch site selects `run_minimap2` vs `run_bwa` on `mstr->aligner`; both BWA-index-build sites are skipped under `--aligner minimap2`; `aligner`/k/w copied into each worker's `mstr`.
- **`tronko-assign/Makefile`**: builds `libminimap2.a` with arch-aware flags (`arm64` → `sse2neon` via `-DKSW_SSE2_ONLY -D__SSE2__`; `x86_64` → `-msse2` + `ksw2_dispatch.c`), renames `kseq_*` symbols to avoid BWA collision, adds the shim/wrapper to `SOURCES`, links `libminimap2.a -ldl`, and adds macOS/Linux platform-lib and Homebrew-zstd include/lib handling (`_XOPEN_SOURCE`/`_DARWIN_C_SOURCE` for Darwin).

### Docs
- `thoughts/shared/plans/2026-05-25-minimap2-aligner-and-runtime-max-leaf-matches.md` — full implementation plan, rationale, and per-phase success criteria.

## Diff at a glance

106 files / +27,621 −103, but **~93 of those files are the vendored `minimap2_src/` third-party tree**. The hand-written changes are ~13 files: `Makefile`, `global.h`, `options.c`, `tronko-assign.c`, `allocateMemoryForResults.{c,h}`, `bwa_source_files/fastmap.c`, `bwa_source_files_include.h`, `crash_debug.h`, plus the new `minimap2_shim.c` and `minimap2_wrapper.{c,h}`.

## Behavior & compatibility

- **Default aligner is unchanged** (`bwa`); minimap2 is strictly opt-in.
- **Default cap changes 10 → 100.** Callers relying on the implicit cap of 10 will see **more conservative (more-correct on the AC db)** assignments. Pin `--max-leaf-matches 10` to restore the exact old behavior. Configs that already set the cap explicitly are unaffected.
- At equal settings (`--aligner bwa --max-leaf-matches 10`) the BWA path is behavior-preserving relative to the previous binary.
- **No database/format migration.** Existing tronko databases work unchanged. minimap2 builds its index **in-memory** from the `-a` reference FASTA; it does **not** create BWA index files (`.bwt`/`.sa`/`.pac`/`.ann`/`.amb`).

## Testing & verification

### Automated
- [x] **Clean build**: `cd tronko-assign && make clean && make` succeeds and produces `libminimap2.a` + `tronko-assign` (arm64 / Apple Silicon path). Only pre-existing warnings (WFA2 typedef-redefinition, `logger.c` format) — none introduced by this change.
- [x] **Help shows the flags**: `--max-leaf-matches` (alias `--max-bwa-matches`), `--aligner`, `--minimap2-kmer`, `--minimap2-window` all appear in `-h`.

### Functional regression vs golden results (single-tree example)
Ran the documented CLAUDE.md example — DB from `tronko-build/example_datasets/single_tree`, query reads from the repo-root `example_datasets/single_tree` (`missingreads_*`), diffed against the committed golden `*_results.txt`:

- [x] **BWA `--max-leaf-matches 10`, single-end (164 reads)** — exit 0; output **byte-identical** to `missingreads_singleend_150bp_2error_results.txt`. Confirms the runtime-cap refactor is behavior-preserving.
- [x] **BWA `--max-leaf-matches 10`, paired-end (9 reads)** — exit 0; output **byte-identical** to `missingreads_pairedend_150bp_2error_results.txt`.
- [x] **BWA default cap (100), single-end** — exit 0; identical to the golden/cap-10 output on this single-tree DB (expected: the cap effect only manifests on the multi-tree AC database).
- [x] **`--aligner minimap2`, single-end & paired-end** — exit 0; **identical taxonomic path on every read** vs the BWA run (0/164 single-end, 0/9 paired-end differ in assignment). Differences are confined to the likelihood score and a ±1 reverse-mate mismatch count — i.e. fully concordant assignments, different seeder arithmetic.

### Not verified here (needs reviewer / CI attention)
- [ ] **The 10→100 accuracy improvement itself** — only observable on the full AncestralClust database, which is not in this repo. The bundled single-tree fixture is in-db, so cap-10 and cap-100 produce identical output here. The 93%-reduction figure comes from the April 2026 investigation referenced in the plan, not from this PR's local run.
- [ ] **Linux build path** (`-msse2` + `ksw2_dispatch.c`, `-lrt -ldl`) — built and tested on macOS arm64 only; recommend confirming the x86_64/Linux Makefile path in CI.
- [ ] **Multithreaded minimap2 index race** (`-T <n>`) and a leak/`valgrind` pass on a larger batch — the plan reports these were checked previously (`leaks --atExit` clean); not re-run here.
- [ ] **Multi-tree example** (`example_datasets/multiple_trees`, partition + multiple_MSA golden results exist) — not run in this pass; would exercise the cap/aligner across >1 tree.

## Migration notes

No action required for existing databases or for callers that already pass an explicit `--max-leaf-matches`. Callers that depended on the implicit cap of 10 and want byte-identical legacy output should add `--max-leaf-matches 10`.

## References

- Plan: `thoughts/shared/plans/2026-05-25-minimap2-aligner-and-runtime-max-leaf-matches.md`
- Reference implementation (reimplemented here, not pulled): `optimize-tronko-build` commits `c0a6538`, `4aa29fc`
- Integration boundary: `tronko-assign/global.h` (`bwaMatches`), `tronko-assign/bwa_source_files/fastmap.c` (leafMap hashmap + concordant/discordant fill)
- Deferred: cross-tree score normalization

🤖 Generated with [Claude Code](https://claude.com/claude-code)
