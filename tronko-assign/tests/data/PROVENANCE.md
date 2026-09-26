# Test data: every committed file and its source

Every data file the tests commit is either derived from the repository's own example dataset,
`tronko-build/example_datasets/single_tree/` (the Charadriiformes reference: one tree, 1,466
leaves, `Charadriiformes.fasta` and `reference_tree.txt`), or, under `real/`, public reads from the
European Nucleotide Archive with production's goldens for them
([Real-reference read sets](#tronko-assigntestsdatareal-real-reference-read-sets)). The
references of the real-reference tests are not committed; [`../manifests/`](../manifests/) lists
them by URL with their SHA-256.

Fixture goldens (`expected_*.tsv`) are the TSV output of production `tronko-assign` (`high-perf`, commit
`71f6ec3`) at `--number-of-cores 1`, built as committed (`MAX_NUM_BWA_MATCHES` 10: `goldens/`,
`multibatch/`) and with production's cap patch (25: `goldens-cap25/`, `multibatch-cap25/`), with
the argument lists of the script that reads them. The fixture goldens under `tronko-assign/tests/data/`
(all but `real/`) were recorded by the scripts' record mode (`TRONKO_TESTS_RECORD=1`) from a plain
`make` build on aarch64 (Ubuntu 24.04, gcc 13.3); the goldens of `tests/data/assignment/` are
byte-identical on aarch64 and x86-64. Until the SAM SEQ fix (commit "bwa: convert read bases before
indexing the SAM SEQ literals"), `mt_single_F` and `mt_unpaired_R` depend on the build: aarch64
builds and production's Debian build reproduce them, Ubuntu gcc 13.3 x86-64 builds do not
([README](../README.md)). They hold read names, taxonomy paths and scores,
no sequence. The exception is `gaps/mate-rescue/`: its goldens were regenerated from the fixed
code by commit "tests: regenerate the mate-rescue goldens from the fixed code", because production's
output on that fixture depends on the `-a` path string ([README](../README.md)).

The fixture generators are in `tronko-assign/tests/tools/` ([README](../tools/README.md)):
`sim_reads.py`, `multitree_fixture.py` and `gap_fixtures.py`, and `make_fixture_inputs.sh`, which
regenerates every fixture input below from the example dataset and compares each file with the
committed one. Each fixture below names the tool that wrote it.

## `tests/data/assignment/` (repository root)

| File | Source |
|---|---|
| `reference_tree.trkb` | `tronko-convert` of the example `reference_tree.txt` (one tree, 1,466 leaves) |
| `single_4000.fasta` | 4,000 single-end reads simulated from example leaves of at least 250 bases (a random 150-base window, 2 substitutions each; names `r<i>_<accession>`), `tools/sim_reads.py`, seed 7 |
| `paired_2000_1.fasta`, `paired_2000_2.fasta` | 2,000 read pairs from the same simulation: the first 150 bases of a leaf and the reverse complement of its last 150 bases, 2 substitutions each (names `p<i>_<accession>`), `tools/sim_reads.py`, seed 7 |
| `single.fasta`, `paired_1.fasta`, `paired_2.fasta` | four reads and four pairs from example leaf `GU572157.1` (the original parity test's inputs; their generator is not in the repository) |
| `expected_single_nw_trkb.tsv`, `expected_unpaired_r_nw_trkb.tsv`, `expected_paired_nw_trkb.tsv` | goldens of `tests/integration/test_assignment_production_parity.sh` (unpaired forward, unpaired reverse, paired) |
| `expected_single_nw.tsv`, `expected_single_wfa.tsv`, `expected_paired_wfa.tsv` | goldens of `tests/integration/test_assignment_parity.sh` on the four-read inputs |

Added by the parity test of pull request #10 (commit "Add assignment parity test with fixtures and
goldens") and by the commit "tests: assignment parity with the production invocation".

## `tronko-assign/tests/data/multibatch/`, `multibatch-cap25/`

| File | Source |
|---|---|
| `expected_paired_nw_trkb_L400.tsv`, `expected_unpaired_r_nw_trkb_L400.tsv` | goldens of `integration/test_multibatch_parity.sh`: the `tests/data/assignment/` pairs at `-L 400` (ten batches) |

## `tronko-assign/tests/data/multitree/` (`tools/multitree_fixture.py`, seed 7)

| File | Source |
|---|---|
| `multitree.fasta` | every leaf of example tree copy 0, every second leaf of copy 1, every leaf of copy 2 (leaf names suffixed `_1`, `_2`) |
| `multitree.fasta.amb`, `.ann`, `.bwt`, `.pac`, `.sa` | the BWA index of `multitree.fasta`, built by `tronko-assign` in a run without `-6` (`tools/make_fixture_inputs.sh`) |
| `reference_tree.trkb` | `tronko-convert` of three copies of the example tree (copy 2 with its taxonomy lines rotated by one) |
| `single_F.fasta`, `single_R.fasta` | 2,000 single-end reads each from example leaves: plain (150 bases), and chimeras of two or three example segments |
| `goldens/`, `goldens-cap25/`: `expected_mt_paired.tsv`, `expected_mt_single_F.tsv`, `expected_mt_unpaired_R.tsv` | goldens of `integration/test_multitree_parity.sh` (the paired case reads `tests/data/assignment/paired_2000_{1,2}.fasta`) |

## `tronko-assign/tests/data/gaps/` (`tools/gap_fixtures.py`, `tools/make_fixture_inputs.sh`)

| Fixture | Files | Source |
|---|---|---|
| `mate-rescue/` | `mr_1.fasta`, `mr_2.fasta` | 620 pairs on example leaves in which one mate cannot be seeded: a substitution every 12th base, or random sequence (uniform A, C, G, T); 40 with a 260-base mate; `gap_fixtures.py mate-rescue`, seed 11 |
| `prod-cmdline/` | `paired_2000_1.fasta.zst`, `paired_2000_2.fasta.zst`, `single_4000.fasta.zst` | the `tests/data/assignment/` reads, compressed as the pipeline compresses them (`zstd -3 --no-check`) |
| `read-content/` | `rc_1.fasta`, `rc_2.fasta`, `rc_single.fasta` | example-leaf reads with one N, random reads that align nowhere, and reads of 30 to 60 bases; `gap_fixtures.py read-content`, seed 13 |
| `prod-names/` | `pn_paired_F.fasta`, `pn_paired_R.fasta`, `pn_unpaired_F.fasta`, `pn_unpaired_R.fasta` | the `tests/data/assignment/` reads renamed as the pipeline names them (`12S_MiFish_U_paired_F_<idx>`); the sequences are the example-derived reads |
| each | `goldens/`, `goldens-cap25/`: `expected_<case>.tsv` | goldens of `integration/test_gap_fixtures.sh` |

## `tronko-assign/tests/data/ksw_fixture_calls.kswd.gz` (`tools/record_ksw_calls.sh`)

2,115 records, the input of `unit/test_ksw_extend.c`: every hundredth `ksw_extend2()` call of each
of the three cases of `tests/integration/test_assignment_production_parity.sh` at one thread
(paired, single, unpaired reverse, in that order), with its inputs and the outputs of BWA's original
scalar code. The record format is described in `unit/test_ksw_extend.c`.
`tools/record_ksw_calls.sh <file>` records it again from the commit before the vectorised
`ksw_extend2`; the decompressed bytes equal the committed file's
(`cmp <(gzip -dc a) <(gzip -dc b)`).

## `tronko-assign/tests/data/real/`: real-reference read sets

Public Illumina MiSeq amplicon reads of the two production markers, from the European Nucleotide
Archive (INSDC data, available without restriction; cite the study accession). Each set's
`reads.json` names every source run with its sample, the ENA FASTQ URLs and their MD5, the primers
and the trimming rules, and the SHA-256 of the committed files. The reads were primer-trimmed (the
last 14 bases of each primer, IUPAC-aware, everything after the match kept), pairs with an N or
shorter than 100 bases dropped; `dev` and `heldout` also drop adapter dimers, remove duplicate pairs
and interleave ten or twelve runs; `dev5k` is a random sample of 5,000 pairs of one run. Read names
are `<run>.<spot>` (`dev5k`) or `<run>_<spot>`.

| Set | Pairs | Study | Runs |
|---|---|---|---|
| `mifish/dev5k` | 5,000 | [PRJNA985779](https://www.ebi.ac.uk/ena/browser/view/PRJNA985779) (SRP444959), NORCE: 18S and MiFish Hywind Scotland water samples metbarcoding data | SRR25108797 |
| `mifish/dev` | 30,000 | [PRJNA985779](https://www.ebi.ac.uk/ena/browser/view/PRJNA985779) (SRP444959), NORCE: 18S and MiFish Hywind Scotland water samples metbarcoding data | 10 runs, SRR25108806 … (`reads.json`) |
| `mifish/heldout` | 30,000 | [PRJNA894161](https://www.ebi.ac.uk/ena/browser/view/PRJNA894161) (SRP405273), AZTI: Comparison of teleo and miFish barcodes | 10 runs, SRR22098887 … (`reads.json`) |
| `fwh/dev5k` | 5,000 | [PRJNA664693](https://www.ebi.ac.uk/ena/browser/view/PRJNA664693) (SRP284209), University of Duisburg-Essen: Improved freshwater macroinvertebrate detection from environmental DNA through minimized non-target amplification | SRR12681741 |
| `fwh/dev` | 30,000 | [PRJNA664693](https://www.ebi.ac.uk/ena/browser/view/PRJNA664693) (SRP284209), University of Duisburg-Essen: Improved freshwater macroinvertebrate detection from environmental DNA through minimized non-target amplification | 10 runs, SRR12681741 … (`reads.json`) |
| `fwh/heldout` | 30,000 | [PRJDB15936](https://www.ebi.ac.uk/ena/browser/view/PRJDB15936) (DRP010099), Kanagawa Environmental Research Center: Environmental DNA samples in Kanagawa rivers | 12 runs, DRR495858 … (`reads.json`) |

MiFish is `12S_MiFish_U` (primers GTCGGTAAAACTCGTGCCAGC / CATAGTGGGGTATCTAATCCCAGTTTG), FWH is
`CO1_fwhF2_EPTDr2n` (GGDACWGGWTGAACWGTWTAYCCHCC / CAAACAAATARDGGTATTCGDTY). Beside the reads:
the goldens, production `tronko-assign` (`71f6ec3`) at one thread built as the pipeline builds it,
packed with each run's command line and `time -v` record in the set's `goldens.tar.zst` (every file listed with its SHA-256 in `goldens.sha256`), and the set's
`provenance.txt` ([README](../README.md#real-reference-goldens)).
