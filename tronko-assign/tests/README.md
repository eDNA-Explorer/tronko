# tronko-assign tests

These tests pin what production `tronko-assign` writes. Every integration test runs the binary
with the flags the eDNA Explorer pipeline uses and compares its output, byte for byte, with a
golden: the output of production `tronko-assign` (`high-perf`, commit `71f6ec3`) at
`--number-of-cores 1`. A change that alters any assignment, score digit or row order fails.

## Running them

```
make -C tronko-assign clean && make -C tronko-assign      # the binary under test
make -C tronko-assign/tests test                          # every test; same as: make -C tronko-assign test
```

`tests/run_tests.sh` prints one `==` line per test: PASS, FAIL or SKIP (a missing fixture is a
skip, exit 77). It exits 1 if any test failed. The whole run takes a few minutes at one thread,
most of it the three-tree and gap fixtures.

Environment, read by every integration script (`integration/golden_lib.sh`):

| Variable | Default | Meaning |
|---|---|---|
| `TRONKO_ASSIGN_BIN` | `tronko-assign/tronko-assign` | binary under test |
| `TRONKO_ASSIGN_CORES` | `1` | thread counts; each must reproduce every golden |
| `TRONKO_MATCH_CAP` | read from `global.h` | `MAX_NUM_BWA_MATCHES` of the binary: 10 as committed, 25 as production builds it |
| `TRONKO_TESTS_XFAIL` | empty | case names whose mismatch is reported as XFAIL, not a failure |
| `TRONKO_TESTS_RECORD` | `0` | `1` writes each output as the golden instead of comparing (one thread count only) |

**Production's match cap.** Production patches `global.h` before it builds
(`sed -i 's/#define MAX_NUM_BWA_MATCHES 10/#define MAX_NUM_BWA_MATCHES 25/' global.h`). Build a
tree patched that way and run the tests as above: the scripts read the cap from `global.h` and
compare with the cap-25 goldens (`<dir>-cap25/`). On these fixtures the cap-25 goldens are
byte-identical to the cap-10 ones, because no read here has more than three candidate trees.

**Thread counts.** Production's output depends on the thread count (each thread runs BWA on its
own slice of the batch), so the default is one thread. `TRONKO_ASSIGN_CORES="1 4"` shows the
dependence: at four threads most cases differ from their goldens.

## What each test pins

| Test | Fixture | What it reaches |
|---|---|---|
| `../../tests/integration/test_assignment_production_parity.sh` | `tests/data/assignment/` at the repository root: the example Charadriiformes reference (one tree, 1,466 leaves) as `.trkb`, 4,000 single-end reads, 2,000 pairs | the three read modes production runs: paired (`-p -z`), unpaired forward (`-s`), unpaired reverse (`-s -v`) |
| `integration/test_multibatch_parity.sh` | the same pairs at `-L 400` | ten batches per run: what is loaded, estimated and reset per batch |
| `integration/test_multitree_parity.sh` | `data/multitree/`: three copies of the example tree in one reference | reads with two or three candidate trees: candidate gathering, votes in several trees, the multi-tree LCA |
| `integration/test_gap_fixtures.sh` | `data/gaps/` | code production runs that the other fixtures never reach (below) |

The four gap fixtures, with the example reference:

| Fixture, cases | Reads | What it reaches |
|---|---|---|
| `mate-rescue`: `mr_paired` | 620 pairs in which one mate cannot be seeded (a substitution every 12th base, or random sequence), 40 of them with a 260-base mate | BWA's mate rescue: `ksw_align2` about 4,000 times per run, including the 16-bit kernel |
| `prod-cmdline`: `pc_single`, `pc_unpaired_r`, `pc_paired` | the repository fixture's reads as `.fasta.zst`, written as the pipeline writes them (one frame, no checksum) | production's full command line, `-R -T --tsv-log <file> -V2`: the zstd reader, logging, the resource monitor, the crash-handler set-up; the goldens equal the repository fixture's |
| `read-content`: `rc_single`, `rc_unpaired_r`, `rc_paired` | reads with one N, random reads that align nowhere, reads of 30 to 60 bases | reads production's QC keeps (`bbduk maxns=1`): an N in scoring, single-end reads with no hit, unaligned mates |
| `prod-names`: `pn_single`, `pn_unpaired_r`, `pn_paired` | the repository fixture's reads named as the pipeline names them, `12S_MiFish_U_paired_F_<idx>` | the read-name handling of the paired readers for names that start with a digit |

## Cases whose golden depends on the build

Two defects in the BWA code that `tronko-assign` carries make some outputs depend on bytes read
out of bounds. The goldens record what production's build writes; a build that lays out memory
differently can differ on these cases only:

- **`mt_single_F`, `mt_unpaired_R`.** BWA's SAM writer formats the SEQ column from bases that were
  never converted to 2-bit codes, reading past two string literals; where the byte read is 0 the
  SAM text ends there and the read's later candidate trees are lost. Production's Debian build and
  aarch64 gcc builds read no 0 for these reads, which align to the reverse strand. Ubuntu gcc 13.3
  x86-64 builds read a 0 for a reverse-strand G and differ in about 620 rows of each case.
- **`mr_paired`.** Mate rescue passes the same unconverted bases to `ksw_align2`, whose query
  profile reads past the 5 × 5 score matrix into a freed heap copy of the `-a` path. The rescue
  scores of the 260-base mates therefore depend on that path string, so the script reads every
  input of this case through one fixed directory, `/tmp/tronko-assign-tests/mate-rescue`
  (it copies the files there). With the same strings the output is the same on every run.

CI on Ubuntu x86-64 lists these three cases in `TRONKO_TESTS_XFAIL`. Both defects are fixed, and
these goldens regenerated, by a later pull request of this series (branch `pr/03-ub-fixes`).

## Data layout

```
tronko-assign/tests/
  run_tests.sh, Makefile
  integration/golden_lib.sh                  shared by the scripts: cap, threads, compare, record
  integration/test_*.sh
  data/multibatch[-cap25]/expected_<case>_L400.tsv
  data/multitree/                            reference_tree.trkb, multitree.fasta + BWA index, single_F.fasta, single_R.fasta
  data/multitree/goldens[-cap25]/expected_<case>.tsv
  data/gaps/<fixture>/                       reads
  data/gaps/<fixture>/goldens[-cap25]/expected_<case>.tsv
```

The pairs of the three-tree fixture are `tests/data/assignment/paired_2000_{1,2}.fasta`; the
gap fixtures use `tronko-build/example_datasets/single_tree/Charadriiformes.fasta` (with its BWA
index) and `tests/data/assignment/reference_tree.trkb`.

## Regenerating the goldens

Only from the commit whose behaviour the goldens pin, at one thread, once per cap:

```
TRONKO_ASSIGN_BIN=<binary built from that commit> TRONKO_MATCH_CAP=10 TRONKO_TESTS_RECORD=1 \
    bash tronko-assign/tests/integration/test_gap_fixtures.sh      # and the other scripts
```
