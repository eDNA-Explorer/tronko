# Tronko minimap2 release analysis

Release-readiness analysis for PR #8 (minimap2 strand reorientation +
sub-optimal hit filtering), run against **real production project data**.

## Context

Synthetic ground-truth benchmarking of PR #8 across four markers found that of the
assignments it moved, **71% moved toward the correct answer**, with truth-depth F1
improving on every marker. Review feedback was that this is not sufficient on its
own to call the update release-ready, and asked for:

1. Real project outputs reviewed in depth — old vs new assignments, assignments
   versus unassigned reads, and how placements move in the reference trees.
2. Bias assessed across markers and datasets, not only the initially affected
   marker or freshwater invertebrates.
3. BLAST-based checks on the reads whose assignments changed, automated.
4. Structured, repeatable release gates rather than discovering issues ad hoc
   through user feedback.

Alongside that, a colleague reported that invertebrate recovery differs between the
November 2025 (BWA) and August 2026 (fixed minimap2) production runs of the BioSCape
CO1 data. That question is folded in here as Strand E.

**The central question, answerable entirely on real data with no synthetic input:**
when an ASV's classification changed, does the new call agree better with its top
BLAST hit than the old one did? And specifically — when the call goes *deeper on the
same lineage*, is that extra depth on the lineage the top hit actually supports, or
is it unearned?

**Constraint: this analysis does not modify tronko source.** It varies only the
binary under test and the inputs fed to it.

## Machine requirements

**tronko-assign + minimap2 against the CO1_Metazoa reference peaks at 24.6 GB RSS**
(measured with `/usr/bin/time -l`: 24,641,882,520 bytes). This is the single most
important operational fact in this document — running two such jobs concurrently
OOM'd a 48 GB machine.

- Minimum **64 GB RAM** to run CO1 cells with headroom.
- ~60 GB disk for databases, reads, and outputs.
- **One CO1-scale job at a time. Never two.**
- Gate the next job on memory actually being released, *not* on an output file
  appearing — a finished-looking file does not mean the process freed its pages.
  That mistake is what caused the OOM.
- vert12S is ~2 GB and may overlap with itself, but never with a CO1 job.

## Bootstrap on a fresh instance

### 1. Build the three binaries

`OPTIMIZE_MEMORY=1` is required, not optional: it makes `type_of_PP` float rather
than double (`global.h:14-19`), which changes which nodes fall inside the
`Cinterval` voting window. Production builds with it. `make clean` is mandatory —
the Makefile lists only `$(TARGET).c libminimap2.a` as prerequisites, so changes to
`minimap2_wrapper.c` do not otherwise trigger a relink.

| binary | commit | role |
|---|---|---|
| `tronko-assign.pre` | `f3bdfac` | pre-fix, the A/B baseline |
| `tronko-assign.post` | `530caf9` | post-fix, what production runs |
| `tronko-assign.bwanov` | `71f6ec3` | BWA-era build for the November arm |

```bash
for spec in pre:f3bdfac post:530caf9 bwanov:71f6ec3; do
  name=${spec%%:*}; sha=${spec##*:}
  git worktree add "wt-$name" "$sha"
  make -C "wt-$name/tronko-assign" clean
  make -C "wt-$name/tronko-assign" OPTIMIZE_MEMORY=1 \
    CC="gcc -O3 -fcommon -Wno-error -Wno-implicit-function-declaration \
        -Wno-incompatible-pointer-types -Wno-int-conversion"
  cp "wt-$name/tronko-assign/tronko-assign" "bin/tronko-assign.$name"
done
```

Provenance assertion — `post` carries a guard string that `pre` does not
(`tronko-assign.c:975`):

```bash
strings bin/tronko-assign.post | grep -c "leaf-portion mode) is not supported"  # 1
strings bin/tronko-assign.pre  | grep -c "leaf-portion mode) is not supported"  # 0
```

### 2. Data manifest — everything to pull from GCS

Everything this analysis uses lives in GCS. Nothing is local-only except the built
binaries (§1) and the scripts in this branch, so a fresh instance can reproduce the
whole substrate from the URIs below. Use `gcloud storage`, **not** `gsutil` —
`gsutil -m` stalled with 51 workers at 0% CPU and zero bytes transferred.

Two buckets are involved, and the distinction matters:

- `gs://edna-reference-databases` — the canonical reference databases production
  loads (`REFERENCE_BUCKET`). Note `gs://edna-project-files-ca/CruxV2/` **also**
  holds a `CruxV2` tree with the same version labels but *different content* and no
  `.trkb`; it is not what production uses.
- `gs://edna-project-files-ca` — project reads and production assignment outputs.

#### 2a. Reference databases

Production loads `reference_tree.trkb` directly and raises `FileNotFoundError` if it
is missing (`assign_core.py:377-415`) — there is no conversion at assign time. **Only
the builds below carry a `.trkb`**, which is how you can tell which build production
actually uses; `20250421` and `20250709` have none and would fail. Confirmed
empirically: a post-fix run on Jalama vert12S reproduced production at **100.00%
exact path match on 5,000 reads**.

Path: `gs://edna-reference-databases/CruxV2/<marker>/<build>/tronko/`

| marker | build | `.trkb` | full dir |
|---|---|---:|---:|
| vert12S | `2023-04-07` | 0.23 GiB | 1.29 GiB |
| CO1_Metazoa | `2023-04-07` | 1.11 GiB | 6.48 GiB |
| 18S_Euk | `2023-04-07` | 0.41 GiB | — |
| ITS2_Plants | `2023-04-07` | 0.54 GiB | — |
| ITS1_Fungi | `2023-04-07` | 0.94 GiB | — |
| 16S_Bacteria | `2023-04-07` | 0.19 GiB | — |
| BF3_BR2 | `20260216` | 1.16 GiB | — |

Minimum per marker for the minimap2 arms — `reference_tree.trkb`, `<marker>.fasta`
(also the BLAST subject), `<marker>_taxonomy.txt` (the BLAST lineage join):

```bash
REF=gs://edna-reference-databases/CruxV2
for M in vert12S CO1_Metazoa 18S_Euk ITS2_Plants ITS1_Fungi; do
  mkdir -p "db/$M"
  gcloud storage cp \
    "$REF/$M/2023-04-07/tronko/reference_tree.trkb" \
    "$REF/$M/2023-04-07/tronko/$M.fasta" \
    "$REF/$M/2023-04-07/tronko/${M}_taxonomy.txt" "db/$M/"
done
# Strand E arm A runs BWA, so CO1 additionally needs the BWA index (~752 MB):
gcloud storage cp "$REF/CO1_Metazoa/2023-04-07/tronko/CO1_Metazoa.fasta.*" db/CO1_Metazoa/
```

#### 2b. Project reads and production assignments

Each production assign run stores the **exact FASTAs it was fed** next to its output,
so the A/B needs no read regeneration:

```
gs://edna-project-files-ca/projects/<projectId>/assign/<qcRunId>/<markerId>/<arm>/
    assignments.parquet                       # what production produced (G3 target)
    <project>-<primer>-<arm>.fasta.zst        # the exact input reads
```
where `<arm>` is `paired`, `unpaired_F`, or `unpaired_R`.

The `<markerId>` values are opaque and cannot be guessed — this table is the mapping,
resolved from the object listings:

| project | projectId | qcRunId | marker | markerId |
|---|---|---|---|---|
| CAT | `cmrphey44000004l4q7tbq4q3` | `cmrwlyw7d000k04jwjhtysee3` | BF3_BR2 | `cwnhxk377w3a0ipcdk1ukbbo` |
| | | | CO1_Metazoa | `jdj9w90a01cg264h82tcb9r4` |
| | | | CO1_fwhF2_EPTDr2n | `v8m0qyrbqi4npawnfy3s97ou` |
| Jalama | `cmrv09djy000004l2dnoqgsih` | `cmrwlzarp00000bjfpqktu3k4` | vert12S | `iw4ajr4phwz0z1gu5mdmj7v3` |
| | | | BF3_BR2 | `quw7dv2bnnray4geulac2et6` |
| | | | CO1_Metazoa | `ub1yracacte1k1j8zcz9k1kh` |
| | | | fwhF2_fwhR2n | `xtg8offb0reuny1efk7oy07v` |
| BioSCape / Plate TI | `cm4n7f2bb0001ssovx4mh0dt4` | `cmsqafuvp000004ld45fuy60z` | vert12S | `83be2cf9-38b4-42ed-8cfb-a47d166a5fc8` |
| | | | 18S_Euk | `9eb33690-3c96-4daa-a2c8-89416fdcb4b3` |
| | | | ITS1_Fungi | `a41ba16c-9d78-41d7-aee9-38bd1d63c570` |
| | | | ITS2_Plants | `d51abee8-6e69-486d-a726-356bd60f952e` |
| | | | CO1_Metazoa | `bfe214ab-fd47-4e47-92be-3a73ad1d98b1` |
| | | | 16S_Bacteria | `6b2b89a5-eaf6-4721-91ed-e9b8e071f7c5` (unpaired arms only) |

Read volumes already measured, so cells can be budgeted before fetching:

| cell | paired pairs | compressed FASTA |
|---|---:|---:|
| Jalama vert12S | 210,316 | 2.4 MB |
| BioSCape vert12S | 1,199,348 | 16.3 MB |
| BioSCape 18S_Euk | — | 149.2 MB |
| BioSCape ITS1_Fungi | — | 150.9 MB |
| BioSCape CO1_Metazoa | — | 166.8 MB |
| BioSCape ITS2_Plants | — | 73.8 MB |
| CAT CO1_Metazoa | — | 233.6 MB |
| CAT BF3_BR2 | 4,296,021 | 206.3 MB |
| Jalama CO1_Metazoa | — | 48.4 MB |

```bash
P=gs://edna-project-files-ca/projects
CELL=$P/cmrv09djy000004l2dnoqgsih/assign/cmrwlzarp00000bjfpqktu3k4/iw4ajr4phwz0z1gu5mdmj7v3/paired
mkdir -p reads/jalama_vert12S
gcloud storage cp "$CELL/*.fasta.zst" "$CELL/assignments.parquet" reads/jalama_vert12S/
```

Skip `paired.manifest.json.zst` — it is tens of MB of per-sequence hashes and is not
used here.

#### 2c. BioSCape CO1 run history (Strand E)

All under `gs://edna-project-files-ca/projects/cm4n7f2bb0001ssovx4mh0dt4/assign/`:

| run | qcRunId | note |
|---|---|---|
| Nov 2025 (BWA) | `cmhwjvd4z000ejs042nriini2` | **v2 FASTA is corrupted — do not use** (see Strand E) |
| Nov 2025 legacy | `CO1_Metazoa/paired/` | **the intact copy — use this** |
| Mar 2026 (BWA) | `cmmjfc1ph0003js040lgonrcr` | internally consistent |
| Jul 2026 (minimap2) | `n1b2grbmyiri37t0mdxrouud` | carries the reference-indexing bug |
| Aug 2026 (minimap2) | `cmsqafuvp000004ld45fuy60z` | current |

The legacy November prefix holds `*-CO1_Metazoa-paired.txt.zst` (the TSV assignments)
alongside `*-paired_F.fasta.zst` / `*-paired_R.fasta.zst`.

## Strand A/B/C — controlled pre/post A/B

Same database, same reads, same parameters; only the binary differs. Production
parameters are also the compiled-in defaults at both commits
(`tronko-assign.c:939-956`):

```
-r -f <trkb> -a <marker.fasta> -p -z -w -1 F.fasta -2 R.fasta
-6 --number-of-cores 1 --Cinterval 0.02 --aligner minimap2
-u 0.0001 --max-leaf-matches 10 --best-leaf-threshold -0.1 --best-leaf-max-votes 10
```

`--number-of-cores 1` is required: tronko is ~3.2% non-deterministic above `-C 1`,
which would sit under the signal being measured.

Note tronko's option struct uses fixed 200-byte buffers for `-f`/`-a`/`-o`
(`global.h:227`). Longer paths silently overrun the NUL terminator into the adjacent
field and fail with exit 255 and an empty error. Keep working paths short.

Cells, cheapest first so the pipeline is proven before the expensive runs:

| dataset | marker | pairs |
|---|---|---|
| Jalama | vert12S | 210,316 (full) |
| BioSCape | vert12S | 300,000 (seeded subsample of 1,199,348) |
| BioSCape | 18S_Euk, ITS2_Plants, ITS1_Fungi, CO1_Metazoa | 300,000 each |
| CAT, Jalama | CO1_Metazoa | 300,000 each |

Subsamples are drawn once with a fixed seed and reused across both arms so the diff
stays paired.

### Strand A — change anatomy

Per cell: old vs new path with the rank of first divergence; the full
assigned↔unassigned transition matrix; **placement movement** from `Tree_Number` /
`Node_Number` (same tree same node / same tree different node / different tree, with
the mismatch delta per bucket); and taxon winners and losers by read count.

Placement movement needs no new instrumentation — tronko already emits both columns
(`parquet_writer.c:26-62`). That same signal is what identified the July 2026
reference-indexing regression, where not one read landed in the first 3,339 trees.

### Strand B — bias audit

Two directions, so marker effects and dataset effects separate: multi-marker within
BioSCape (same samples, same prep) and multi-dataset within CO1 (CAT / Jalama /
BioSCape). Per clade, and against covariates that would indicate systematic rather
than uniform improvement: reference representation for the clade, ASV abundance, and
top-hit identity band.

### Strand C — BLAST arbiter (the core)

For each changed ASV, blastn both mates against the tronko database's **own**
`marker.fasta`, sum bitscores per subject across mates, and map the winning
accession to a lineage via `<marker>_taxonomy.txt`. Because both the tronko paths
and the BLAST lineages come from the same reference, no taxonomy remapping is
needed. That also scopes the claim: this adjudicates **placement**, not reference
coverage — a taxon absent from the reference is invisible to both tools alike.

| change type | condition | adjudication |
|---|---|---|
| **extension** | pre is a prefix of post | do the *added* ranks match the BLAST lineage? → confirmed / contradicted / unverifiable |
| **retraction** | post is a prefix of pre | did the dropped ranks disagree with BLAST? then the retraction was right |
| **divergence** | paths split | which side agrees with BLAST further down |
| **assigned → unassigned** | post gave nothing | was the lost call BLAST-supported, or BLAST-contradicted? |
| **unassigned → assigned** | pre gave nothing | is the gained call BLAST-supported? |

Deepening is credited **only** when the added ranks sit on the lineage the top hit
supports. Deeper-but-off-lineage is `extension_contradicted` and counts against the
change, not for it.

Everything is stratified by top-hit identity band (≥99, 97–99, 90–97, 80–90, 50–80,
<50) because an 85%-identity hit cannot license a species-level call; species-rank
agreement is only treated as meaningful at ≥97%. Results are reported both ASV-wise
and read-weighted, since a few abundant ASVs can carry the read-level story in the
opposite direction from the ASV-level one.

BLAST top hit is an independent second opinion, not ground truth. Exact bitscore
ties are reported as their own bucket rather than broken arbitrarily.

## Strand E — CO1 invertebrate case study

Both the November 2025 BWA run and the August 2026 minimap2 run indexed the **same
2023 CO1 reference**, so the reference is already held constant in production and no
2×2 reconstruction is needed.

- **Input**: November 2025 CO1 sequences, 2,557 dereplicated paired ASVs
- **Reference**: CO1_Metazoa 2023-04-07
- **Arm A** — BWA, `71f6ec3`
- **Arm B** — fixed minimap2, `530caf9`
- **Arm C** — minimap2 with `normalize_scores` / `best_leaf_threshold` disabled

Arm C is what makes this diagnostic rather than merely descriptive. The
BWA→minimap2 swap also carried placement-scoring changes, which are flag-gated and
default-off since `cd534e5`. If the invertebrates return with those flags off, the
aligner is exonerated and the scoring change is the cause — a distinction worth
having before concluding the aligner is at fault.

Changed ASVs are then adjudicated with the Strand C arbiter, and filtered to the
invertebrate clades (Insecta, Diptera, Chironomidae, Coleoptera, Arthropoda) for the
case study view.

### Blocking prerequisite — use the legacy November files

The November `assignments.parquet` and the November FASTA under the same **v2** GCS
prefix **do not describe each other**:

| pairing | `forward_length` agreement |
|---|---:|
| March parquet ↔ March FASTA | 100.0% (2557/2557) |
| **November parquet ↔ November FASTA** | **22.2%** (567/2557) |
| November parquet ↔ March FASTA | 20.9% |

22.2% is below the 31.7% chance baseline, so the two are uncorrelated rather than
merely reordered — unrecoverable in that layout. A retroactive March 2026 backfill
overwrote the November FASTA. The **legacy pre-backfill copy is intact** and must be
used instead:

```
gs://edna-project-files-ca/projects/cm4n7f2bb0001ssovx4mh0dt4/assign/CO1_Metazoa/paired/
```

It is 100% self-consistent, and its TSV is confirmed as November by 100%
forward-mismatch and 98.1% score/tree agreement with the November parquet, against
~1% for March and July. (98.1% rather than 100% is consistent with November having
run at `-c 20`, where tronko is nondeterministic.)

**This is a production data-integrity bug independent of PR #8** and is reported
separately rather than folded into the release verdict.

## Release gates

| gate | check |
|---|---|
| G1 determinism | post-vs-post null run byte-identical at `cores=1` |
| G2 provenance | recorded binary SHA256s + the post-only guard string |
| G3 reproduces prod | local post-fix vs production `assignments.parquet` ≥99% path agreement |
| G4 BLAST concordance | post-agrees-better share exceeds pre-agrees-better, per marker |
| G5 earned depth | `extension_contradicted` share below threshold in the ≥97% band |
| G6 no clade bias | no well-powered clade inverting the pooled direction |
| G7 unassigned | no marker's unassigned rate rises by more than 2 pp |

G3 is the gate that would have caught the July reference-indexing regression before
users did. Thresholds for G4–G7 are calibrated from the first full run and then
ratified; until then they are not meaningful, and reports say so.

## Verification

1. **Memory** — every cell records peak RSS. A CO1 cell starting while another
   tronko process still holds >8 GB is a queue bug, not something to retry.
2. **Null run** — post vs post byte-identical on one real cell; otherwise every
   delta needs a noise band.
3. **Reproduce production** — G3 per cell. Low agreement means the database build or
   parameters are mismatched and that cell's numbers are void.
4. **Read identity** — assert the pre and post arms cover the same readname set
   before diffing; a mangled header set reads as churn rather than as an error.
5. **Input pairing** — for Strand E, re-assert the legacy November FASTA↔TSV
   agreement is 100% before running.
6. **Dereplication round-trip** — ASV read weights must sum back to the input count.
7. **BLAST join** — unresolved top-hit accessions are counted and reported, never
   silently dropped.
8. **Hand spot-check** — ~20 changed ASVs spanning extension-confirmed,
   extension-contradicted, retraction and divergence, verified by eye.

## Scripts

Under `scripts/release_analysis/`:

| script | role |
|---|---|
| `tronko_blast_arbiter.py` | Strand C — dereplicate, BLAST, classify, adjudicate, band-stratify |
| `tronko_change_anatomy.py` | Strand A — path diff, transition matrix, placement movement, taxon deltas |
| `tronko_bias_audit.py` | Strand B — per-clade and covariate breakdowns |
| `run_cell.sh` | one `{pre,post}` × cell run at production parameters |
| `subsample.py` | seeded, pair-preserving subsample |

`run_cell.sh` validates that the output row count matches the query count before
treating a cell as complete; a killed run leaves a partial TSV that would otherwise
be silently accepted as finished.

## Out of scope

- Modifying tronko source.
- BLAST against `core_nt` (which would separate "misplaced" from "reference lacks
  it") — a follow-up if the marker-reference arbiter proves ambiguous.
- Wiring the gates into CI; the heavy arms cannot run there on time or storage.
- Pooling Strand E's November-vs-August numbers with the controlled pre/post
  numbers. That contrast spans two pipeline versions, so the two remain separate
  evidence lines.
