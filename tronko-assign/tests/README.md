# tronko-assign tests

These tests pin what `tronko-assign` writes, as the eDNA Explorer pipeline runs it. Every
integration test runs the binary with the flags the pipeline uses and compares its output, byte
for byte, with a golden: the output at `--number-of-cores 1` of the commit whose behaviour the
golden pins: production `tronko-assign` (`high-perf`, commit `71f6ec3`) with the determinism
change, the undefined-behaviour fixes, F1 removed and the zero insert-size spread guarded
(commit `ef420cf`), made on x86-64 ([below](#goldens-regenerated-by-removing-f1)). A change that
alters any assignment, score digit or row order fails.

## Running them

```
make -C tronko-assign clean && make -C tronko-assign      # the binary under test
make -C tronko-assign/tests test                          # every test; same as: make -C tronko-assign test
```

`make test` also builds the unit tests (`tests/unit/`) with the binary's compiler flags.
`tests/run_tests.sh` prints one `==` line per test: PASS, FAIL or SKIP (a missing fixture is a
skip, exit 77). It exits 1 if any test failed. The whole run takes a few minutes, most of it the
three-tree and gap fixtures at one thread.

Environment, read by every integration script (`integration/golden_lib.sh`):

| Variable | Default | Meaning |
|---|---|---|
| `TRONKO_ASSIGN_BIN` | `tronko-assign/tronko-assign` | binary under test |
| `TRONKO_ASSIGN_CORES` | `1 4 16` | thread counts; each must reproduce every golden |
| `TRONKO_MATCH_CAP` | read from `global.h` | `MAX_NUM_BWA_MATCHES` of the binary: 10 as committed, 25 as production builds it |
| `TRONKO_TESTS_XFAIL` | empty | case names whose mismatch is reported as XFAIL, not a failure |
| `TRONKO_TESTS_RECORD` | `0` | `1` writes each output as the golden instead of comparing (one thread count only) |

**Production's match cap.** Production patches `global.h` before it builds
(`sed -i 's/#define MAX_NUM_BWA_MATCHES 10/#define MAX_NUM_BWA_MATCHES 25/' global.h`). Build a
tree patched that way and run the tests as above: the scripts read the cap from `global.h` and
compare with the cap-25 goldens (`<dir>-cap25/`). On these fixtures the cap-25 goldens are
byte-identical to the cap-10 ones, because no read here has more than three candidate trees.

**Thread counts.** Every thread count must reproduce the one-thread goldens: BWA runs once per
batch over the whole batch, and tie-breaks are keyed to a read's position in the batch, not to
the thread that places it. The default checks 1, 4 and 16 threads.

## What each test pins

| Test | Fixture | What it reaches |
|---|---|---|
| `../../tests/integration/test_assignment_production_parity.sh` | `tests/data/assignment/` at the repository root: the example Charadriiformes reference (one tree, 1,466 leaves) as `.trkb`, 4,000 single-end reads, 2,000 pairs | the three read modes production runs: paired (`-p -z`), unpaired forward (`-s`), unpaired reverse (`-s -v`) |
| `integration/test_multibatch_parity.sh` | the same pairs at `-L 400` | ten batches per run: what is loaded, estimated and reset per batch |
| `integration/test_multitree_parity.sh` | `data/multitree/`: three copies of the example tree in one reference | reads with two or three candidate trees: candidate gathering, votes in several trees, the multi-tree LCA |
| `integration/test_gap_fixtures.sh` | `data/gaps/` | code production runs that the other fixtures never reach (below) |
| `unit/test_sam_seq` | none | BWA's SAM writer: for every byte a read can hold, on both strands, the SEQ column is the read as `A C G T N` and the SAM text has no NUL before its end, so the parse finds every record (7,710 checks) |
| `unit/test_mate_rescue` | none | mate rescue in all four orientations, 8-bit and 16-bit kernels, with `mem_opt_t` placed against a page with no access rights: every index into the score matrix is in bounds and the mate is found where it was taken from |
| `integration/test_path_options.sh` | none | each of the thirteen path options at the longest length its buffer holds (accepted) and one longer (refused, exit 1, message naming the option) |
| `integration/test_slot_cap.sh` | chimeric reads made from the example reference | a copy built with `MAX_NUM_BWA_MATCHES 2` and AddressSanitizer on reads with three candidate leaves: the SAM parse stops at the last slot and placement reads no further |

The four gap fixtures, with the example reference:

| Fixture, cases | Reads | What it reaches |
|---|---|---|
| `mate-rescue`: `mr_paired` | 620 pairs in which one mate cannot be seeded (a substitution every 12th base, or random sequence), 40 of them with a 260-base mate | BWA's mate rescue: `ksw_align2` about 4,000 times per run, including the 16-bit kernel |
| `prod-cmdline`: `pc_single`, `pc_unpaired_r`, `pc_paired` | the repository fixture's reads as `.fasta.zst`, written as the pipeline writes them (one frame, no checksum) | production's full command line, `-R -T --tsv-log <file> -V2`: the zstd reader, logging, the resource monitor, the crash-handler set-up; the goldens equal the repository fixture's |
| `read-content`: `rc_single`, `rc_unpaired_r`, `rc_paired` | reads with one N, random reads that align nowhere, reads of 30 to 60 bases | reads production's QC keeps (`bbduk maxns=1`): an N in scoring, single-end reads with no hit, unaligned mates |
| `prod-names`: `pn_single`, `pn_unpaired_r`, `pn_paired` | the repository fixture's reads named as the pipeline names them, `12S_MiFish_U_paired_F_<idx>` | the read-name handling of the paired readers for names that start with a digit |

## Goldens regenerated by the undefined-behaviour fixes

Before the undefined-behaviour fixes (tested by `unit/test_sam_seq` and `unit/test_mate_rescue`),
two out-of-bounds reads in the BWA code that `tronko-assign` carries made some outputs depend on
memory layout. Both are fixed, and the goldens they affected were regenerated from the fixed code
at one thread:

- **`mr_paired` (regenerated).** Mate rescue passed unconverted ASCII bases to `ksw_align2`,
  whose query profile read past the 5 × 5 score matrix into a freed heap copy of the `-a` path,
  so rescue scores depended on that string. With the mate converted to 2-bit codes, 37 of the
  620 rows changed: 22 pairs whose unseedable mate has 150 bases are now rescued by the 8-bit
  kernel (before, the garbage scores saturated it and every such rescue was rejected), and 15
  pairs with a 260-base mate get the rescue score upstream BWA computes. One of the 37 is
  assigned two ranks deeper (family to species); the others change only in score. The old goldens held only for one
  `-a` string; the new ones hold for any. The script still reads the inputs of this case through
  `/tmp/tronko-assign-tests/mate-rescue`, so that the goldens of earlier commits can be checked
  with the same arguments.
- **`mt_single_F`, `mt_unpaired_R` (unchanged).** BWA's SAM writer formatted the SEQ column from
  unconverted bases, reading past two string literals; where the byte read was 0 the SAM text
  ended there and a read's later candidate trees were lost. Which letters read a 0 depended on
  the build: none in production's Debian build, a reverse-strand G in Ubuntu gcc 13.3 x86-64
  builds (about 620 rows of each case differed there), a forward-strand T in aarch64 builds.
  The goldens were production's output all along; with the bases converted, every build
  reproduces them. Production's own build writes the same output before and after the fix.

## Goldens regenerated by removing F1

The BWA code `tronko-assign` carries is a fork of BWA 0.7.17 (r1194). Its `worker1` called
`mem_align1()` where upstream `bwa mem` calls `mem_align1_core()`, and the fork's `mem_align1()`
marked a read's primary hit an extra time (`mem_mark_primary_se()` keyed by an `lrand48()` draw,
which `pr/02-determinism` replayed so that every thread count drew the same value). That extra
marking, F1, reorders a read's hits before `worker2` marks them as upstream does, so among hits with
equal alignment score the fork took a different primary hit, and with it a different candidate
tree, from upstream BWA on the same reads. With F1 removed, the four SAM fields `tronko-assign`
reads (QNAME, FLAG, RNAME, RNEXT) equal those of upstream BWA 0.7.17-r1194 run once per batch.

Every golden whose reads have equally scored hits on different trees changed: 18 of the 21 fixture
goldens (all but the three four-read goldens of `test_assignment_parity.sh`) and every
real-reference golden of MiFish and FWH (`dev5k`, `dev`, `heldout`; TSV and parquet). The commit that regenerates them lists the rows per
golden; every changed row belongs to a read whose SAM fields differ between the code before and
after the change, and the primary hits that differ have equal alignment scores.

The commit "bwa: take a zero insert-size spread as no deviation, as BWA-MEM3 does", after it, removes
an undefined behaviour BWA-MEM3 has already fixed: when every pair of a batch has the same spacing,
the insert-size estimate has a standard deviation of 0 and the pairing score was computed from 0/0.
Among the goldens only the three paired fixtures that read the 2,000 simulated pairs as one batch
reach it (`tests/data/assignment` paired, `pc_paired`, `pn_paired`): 309 of their 2,000 rows change
again, and they were recorded once more with the guarded build, in the same way. With it the SAM
fields `tronko-assign` reads equal BWA-MEM3 v0.14.0 `--compat=bwa-mem`'s run once per batch on
every fixture.

How they were made, at one thread, on x86-64 (Intel Cascade Lake, Ubuntu 24.04, gcc 13.3):

- fixture goldens: the record mode of each script (`TRONKO_TESTS_RECORD=1`) from a plain `make`
  build at match cap 10 and at match cap 25 (the two give the same goldens), and the three cases of
  `../../tests/integration/test_assignment_production_parity.sh` with that script's arguments;
  recorded twice, byte-identical; an aarch64 build (Ubuntu 24.04, gcc 13.3) records the same
  bytes for every fixture golden;
- real-reference goldens (all three sets): the pipeline's recipe (below), production's full
  command line, as `make update-goldens` runs them, from the F1 commit; `provenance.txt` in each set
  gives commit, machine, dates and wall times. No batch of these read sets makes an insert-size
  estimate, so the zero-spread guard cannot act on them: the guarded build at 64 threads writes
  every one of them byte for byte.

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
  data/real/<marker>/<set>/                  real-reference read pairs and goldens (below)
  manifests/<marker>.sha256                  the real references, by URL with SHA-256
  real/                                      Makefile, lib.sh, test-real.sh, verify-branches.sh, make-goldens.sh
  tools/                                     the fixture generators: make_fixture_inputs.sh and its Python tools
```

The pairs of the three-tree fixture are `tests/data/assignment/paired_2000_{1,2}.fasta`; the
gap fixtures use `tronko-build/example_datasets/single_tree/Charadriiformes.fasta` (with its BWA
index) and `tests/data/assignment/reference_tree.trkb`.

## Real-reference goldens

The fixtures above are small. The real-reference tests run production's own references,
`12S_MiFish_U` (MiFish, 342 MiB) and `CO1_fwhF2_EPTDr2n` (FWH, 6.8 GiB), on public reads: `dev5k`
(5,000 pairs), `dev` and `heldout` (30,000 pairs each, from two other studies). The read pairs and
the goldens are in the repository; the references are fetched by URL and checked by checksum.

```
tronko-assign/tests/data/real/<marker>/<set>/
    <set>_F.fasta.zst, <set>_R.fasta.zst   the read pairs, compressed as the pipeline writes them
    reads.json                             their source: study, runs, ENA files with MD5, primers, trimming
    goldens.tar.zst                        the goldens and the record of each golden run, one archive
    goldens.sha256                         the SHA-256 of every file in goldens.tar.zst
    provenance.txt                         commit, build recipe, threads, machine, date
tronko-assign/tests/manifests/<marker>.sha256   the reference files, by URL with SHA-256
```

`goldens.tar.zst` holds `<marker>_<case>.tsv` (the TSV goldens, `-o`),
`parquet/<marker>_<case>.parquet` (`--parquet`) and
`provenance/<marker>_<case>.<format>[.cap10].cmd.txt` and `.time.txt` (the command line and `time -v`
record of each golden run). The scripts unpack it into the cache and check every file against
`goldens.sha256` before using it; to look at the goldens, `zstd -dc goldens.tar.zst | tar -x -C
<dir>`. One archive per set keeps a pull request's diff small: `goldens.sha256` shows which goldens
a change touches, and `update-goldens` prints the rows that moved.

The goldens are the output at one thread of the code with F1 removed
([above](#goldens-regenerated-by-removing-f1)), built as the pipeline builds (below, cap 25). The references are production's objects of 2023-04-07,
`gs://edna-reference-databases/CruxV2/<reference>/2023-04-07/tronko/`: `manifests/<marker>.sha256`
lists the seven files `tronko-assign` opens (`reference_tree.trkb`, the FASTA and its five BWA
index files) with their SHA-256. The checksums, not the URLs, identify the files: a copy from
anywhere passes if its bytes match.

Four targets in `real/Makefile` use them. Each fetches the reference files by URL into a cache
(`~/.cache/tronko-test-real`, or `TRONKO_REAL_CACHE`), checks each against the manifest when it
arrives and again before anything runs, and stops naming the file on any mismatch. Fetching uses
`gcloud storage cp`, else `gsutil cp`, else `curl` from `storage.googleapis.com`; it needs read
access to the bucket. `TRONKO_REAL_MIRROR=<dir>` reads a local copy laid out like the bucket
instead ([Running without eDNA Explorer's Google Cloud](#running-without-edna-explorers-google-cloud)).
Every target builds with the eDNA Explorer pipeline's recipe (`build_pipeline` in `real/lib.sh`):
the cap patched to 25, then
`make CC="gcc -O3 -fcommon -Wno-error -Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-int-conversion"`,
with `ENABLE_PARQUET=1` for the parquet build. Every run uses production's full command line:
`.fasta.zst` reads and `-w -6 --Cinterval 10 -R -T --tsv-log <file> -V2`, with `-o` (TSV) or
`--parquet` (parquet).

| Target | What it does |
|---|---|
| `test-real` | builds this checkout, runs every case at each of `THREADS`, compares every output with its golden with `cmp`; writes `results.tsv` in the work directory |
| `verify-branches` | `test-real` on each revision in `BRANCHES` (exported from this repository, built separately), against this checkout's goldens; one summary table, `summary.tsv` |
| `goldens` | builds `GOLDEN_COMMIT` (default `ef420cf`), makes every golden again at one thread, compares each with the committed one: `IDENTICAL`, `DIFFER` or `NEW`; fails unless all are identical |
| `update-goldens` | the same, and writes every golden that is new or differs into its set's `goldens.tar.zst`, with its run record, rewrites `goldens.sha256`, and prints the rows of each golden that differs; a set whose goldens are all identical is left untouched, so `git status` names exactly the sets that changed |

| Variable | Default | Meaning |
|---|---|---|
| `MARKERS` | `mifish` | `mifish`, `fwh` |
| `SETS` | `dev5k` | `dev5k`, `dev`, `heldout` |
| `CASES` | `paired unpF unpR` | read modes: paired (`-p -z`), unpaired forward (`-s`), unpaired reverse (`-s -v`) |
| `FORMATS` | `tsv parquet` | parquet needs `tronko-assign/carquet` (`git submodule update --init`) |
| `THREADS` | `16` | `test-real`, `verify-branches`: thread counts, each must reproduce every golden |
| `BRANCHES` | none | `verify-branches`: branches or commits |
| `GOLDEN_COMMIT` | `ef420cf` | `goldens`, `update-goldens`: the commit the goldens are made from |
| `TRONKO_REAL_WORK` | a new `mktemp -d` | builds, runs and outputs |

**Thread counts.** Every thread count must reproduce the one-thread goldens; the default is 16,
the thread count production runs.

Memory: MiFish needs about 14 GiB at one thread and 16 GiB at 16 threads; FWH about 124 GiB at
one thread and 146 GiB at 16 threads (one copy of the BWA index), so FWH needs a machine with at
least 192 GiB. One case at one thread takes, by the `time -v` records of the golden runs
(`provenance/` in each set's `goldens.tar.zst`, 32-vCPU Intel Cascade Lake): MiFish `dev5k` 14 to 29 minutes, `dev`
and `heldout` 1.3 to 3.8 hours; FWH `dev5k` 8 to 12 minutes, `dev` and `heldout` 20 to 62 minutes.
Making every golden again, (a) below, is about 38 hours of one-thread runs; `MARKERS`, `SETS` and
`CASES` select a part.

### (a) Reproduce every golden at one thread

```
git fetch origin pr2/04-f1-removed                 # the goldens' commit ef420cf must be present
git submodule update --init tronko-assign/carquet
make -C tronko-assign/tests/real goldens MARKERS="mifish fwh" SETS="dev5k dev heldout"
```

This exports `ef420cf` (with its carquet commit) from the repository, builds it plain and parquet
with the pipeline's recipe, runs every marker, set, case and format at `--number-of-cores 1`, and
compares each output with the committed golden. It prints one line per golden and fails if any
differs.

### (b) Verify a branch, or every branch of a series, against the goldens

```
git checkout <the merged branch>
git submodule update --init tronko-assign/carquet
make -C tronko-assign/tests/real test-real MARKERS="mifish fwh" SETS="dev5k dev heldout" THREADS="1 16"
make -C tronko-assign/tests/real verify-branches BRANCHES="<branch> <branch> ..." THREADS=16
```

`make -C tronko-assign/tests/real test-real` alone runs MiFish `dev5k` at 16 threads, the smallest
check (six runs, a few minutes in all).

### (c) Add goldens, or replace them after an intended change

A new read set needs its pairs, `data/real/<marker>/<set>/<set>_{F,R}.fasta.zst` (`zstd -3
--no-check`), and a `reads.json` that names their source; a new marker also needs
`manifests/<marker>.sha256` with its seven reference files. Then

```
make -C tronko-assign/tests/real update-goldens MARKERS=<marker> SETS=<set>
git status tronko-assign/tests/data/real        # the sets with new goldens, and nothing else
```

When a change is meant to alter the output, make the goldens from the commit that defines the new
behaviour (`GOLDEN_COMMIT=<commit>`), review the rows that moved (`update-goldens` prints them; `git diff` of `goldens.sha256` names the goldens), and commit them
with the change.

### Running without eDNA Explorer's Google Cloud

Nothing in these tests needs eDNA Explorer's infrastructure except the location of the reference
files. For anyone outside it:

- **References.** `gs://edna-reference-databases/` is not publicly readable. With a copy of the
  seven files of a marker, place them under a directory laid out like the bucket,
  `<dir>/CruxV2/<reference>/2023-04-07/tronko/<file>`, and set `TRONKO_REAL_MIRROR=<dir>`; each
  file is checked against the SHA-256 in `manifests/<marker>.sha256`, so any copy with the same
  bytes works. Files that differ (another build of the reference) fail that check by design.
- **The public reference libraries.** eDNA Explorer publishes its reference libraries on Zenodo,
  [doi:10.5281/zenodo.15353120](https://doi.org/10.5281/zenodo.15353120) (CC BY 4.0), as
  `<marker>.fasta` and `<marker>.tax.tsv`, for example `12S_MiFish_U.fasta`. They are not the files
  these goldens were made with: the record holds no `reference_tree.trkb` or BWA index, and its
  sequences are another build. Its `12S_MiFish_U.fasta` has 153,684 sequences against the 153,675
  of production's 2023-04-07 FASTA (56 accessions only on Zenodo, 47 only in production, 225 with a
  different sequence); its `CO1_fwhF2_EPTDr2n.fasta` is 436,339,047 bytes against production's
  1,885,411,333 (alignment rows). To test against a public reference instead: build a Tronko
  reference from the library with `tronko-build`, convert its `reference_tree.txt` to
  `reference_tree.trkb` with `tronko-convert`, and index its FASTA (a name ending in `.fasta`) with
  `bwa index`. List the seven files in a new `manifests/<marker>.sha256`, one
  `<sha256>  gs://<bucket>/<path>` line each; any bucket name will do, because with
  `TRONKO_REAL_MIRROR=<dir>` each file is read from `<dir>/<path>`. Copy a set's read pairs and
  `reads.json` to `data/real/<marker>/<set>/`, make goldens for it from `ef420cf` with
  `make -C tronko-assign/tests/real update-goldens MARKERS=<marker> SETS=<set>`, and verify changes
  against those.
- **Build recipe.** The pipeline's recipe is written out in full above and in `real/lib.sh`
  (`build_pipeline`); nothing outside this repository is needed to build.
- **Machine.** Any Linux machine with gcc, make, zstd, zlib and enough memory (above). The goldens
  were made on x86-64 (Intel Cascade Lake, Ubuntu 24.04, gcc 13.3). Floating-point results can
  depend on the compiler and CPU architecture; if `make goldens` reports `DIFFER` on another
  machine, make goldens from `ef420cf` there (`update-goldens` on a scratch branch) and verify
  changes against those.
- **Reads.** The read pairs are committed; `reads.json` names the public ENA runs they come from.

## Regenerating the goldens

Only from the commit whose behaviour the goldens pin, at one thread, once per cap:

```
TRONKO_ASSIGN_BIN=<binary built from that commit> TRONKO_MATCH_CAP=10 TRONKO_TESTS_RECORD=1 \
    bash tronko-assign/tests/integration/test_gap_fixtures.sh      # and the other scripts
```

The real-reference goldens are made with `make update-goldens` ([above](#c-add-goldens-or-replace-them-after-an-intended-change)).

The fixture inputs themselves (reads, references, BWA index) are regenerated from the example
dataset by `tools/make_fixture_inputs.sh <out-dir>`, which also compares every file with the
committed one ([tools/README.md](tools/README.md)); `data/PROVENANCE.md` names the tool that wrote
each fixture.
