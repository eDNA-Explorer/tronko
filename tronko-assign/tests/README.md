# tronko-assign tests

These tests pin what `tronko-assign` writes in production, as the eDNA Explorer pipeline runs it.
Every integration test runs the binary with the flags the pipeline uses and compares its output,
byte for byte, with a golden: the output of production `tronko-assign` (`high-perf`, commit
`71f6ec3`) at `--number-of-cores 1`, except the mate-rescue golden, which is the output of the fixed
code (below). A change that alters any assignment, score digit or row order fails.

## Running them

```
make -C tronko-assign clean && make -C tronko-assign      # the binary under test
make -C tronko-assign/tests test                          # every test
make -C tronko-assign/tests test-san                      # the unit tests of tests/Makefile under ASan and UBSan
make -C tronko-assign/tests test-shadow                   # tronko-assign built with -DNW_SHADOW on the repository fixtures
make -C tronko-assign/tests test-tsan                     # the threaded unit tests under ThreadSanitizer (x86-64)
make -C tronko-assign/tests clean
```

`make -C tronko-assign/tests test` builds and runs the unit tests of `tests/Makefile` with
production's flags (gcc `-O3 -fno-stack-protector`, no `-march`), then `tronko-assign`'s own
`make test`. Requirements: gcc and zlib (`build-essential zlib1g-dev libzstd-dev`), Linux on
x86-64 or aarch64. If ThreadSanitizer stops with "unexpected memory mapping" (kernels with high
ASLR entropy), run it as `setarch x86_64 -R make -C tronko-assign/tests test-tsan`.

`tronko-assign`'s `make test` builds the unit tests of the undefined-behaviour fixes with the
binary's compiler flags and runs `tests/run_tests.sh`.
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
| `unit/test_vote_tally` | random forests and reads | the vote tally over a read's candidate trees equals the pass over every tree (below) |
| `unit/test_ksw_extend` | `data/ksw_fixture_calls.kswd.gz` and random calls | the vectorised seed-extension kernels return what BWA's scalar `ksw_extend2` returns (below) |
| `unit/test_nodestore` | random trees and reads | the blocked node-store scoring equals production's per-node scoring bit for bit (below) |
| `integration/test_nodestore_paths.sh` | `tests/data/assignment` | every node-store variant (block sizes, budget, kernels, the legacy loop, reads without candidates at block edges) reproduces the goldens (below) |
| `unit/test_nw_fill` | `data/nw_fixture_alignments.txt.gz`, small and random alignments | the two-pass Needleman-Wunsch fill equals the original fill cell by cell (below) |
| `integration/test_nw_shadow.sh` (`make test-shadow`) | `tests/data/assignment` | a build that runs both fills on every alignment and aborts on a difference reproduces the goldens |
| `unit/test_dense_sa` | the example BWA index and synthetic indexes | the dense suffix-array sample returns what BWA's sampled lookup returns for every entry (below) |
| `unit/test_chain_flt` | `data/chain_flt_fixture.chf.gz`, random and constructed calls | the grouped chain filter equals BWA's `mem_chain_flt` (below) |
| `unit/test_pos2rid` | the example BWA index and random intervals | the bucketed accession lookup equals BWA's `bns_intv2rid`/`bns_pos2rid` (below) |
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
  data/ksw_fixture_calls.kswd.gz             ksw_extend2 calls recorded on the fixtures (tools/record_ksw_calls.sh)
  data/nw_fixture_alignments.txt.gz          needleman_wunsch_align() inputs recorded on the fixtures (tools/record_nw_alignments.sh)
  data/chain_flt_fixture.chf.gz              mem_chain_flt() calls recorded on the fixtures (tools/record_chain_flt.sh)
  data/real/<marker>/<set>/                  real-reference read pairs and goldens (below)
  manifests/<marker>.sha256                  the real references, by URL with SHA-256
  real/                                      Makefile, lib.sh, test-real.sh, verify-branches.sh, make-goldens.sh
  tools/                                     the fixture generators: make_fixture_inputs.sh and its Python tools
  tools/record_ksw_calls.sh, tools/ksw_record.patch   record data/ksw_fixture_calls.kswd.gz again
  tools/record_nw_alignments.sh              record data/nw_fixture_alignments.txt.gz again
  tools/record_chain_flt.sh                  record data/chain_flt_fixture.chf.gz again
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

The goldens are production's (`71f6ec3`) output at one thread, built as the pipeline builds it
(below, cap 25). The four FWH single-end parquet goldens of `dev` and `heldout` come from the cap-10
build, which writes the same output on these reads (their `provenance.txt` says how this was
checked). The references are production's objects of 2023-04-07,
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
| `goldens` | builds `GOLDEN_COMMIT` (default `71f6ec3`), makes every golden again at one thread, compares each with the committed one: `IDENTICAL`, `DIFFER` or `NEW`; fails unless all are identical |
| `update-goldens` | the same, and writes every golden that is new or differs into its set's `goldens.tar.zst`, with its run record, rewrites `goldens.sha256`, and prints the rows of each golden that differs; a set whose goldens are all identical is left untouched, so `git status` names exactly the sets that changed |

| Variable | Default | Meaning |
|---|---|---|
| `MARKERS` | `mifish` | `mifish`, `fwh` |
| `SETS` | `dev5k` | `dev5k`, `dev`, `heldout` |
| `CASES` | `paired unpF unpR` | read modes: paired (`-p -z`), unpaired forward (`-s`), unpaired reverse (`-s -v`) |
| `FORMATS` | `tsv parquet` | parquet needs `tronko-assign/carquet` (`git submodule update --init`) |
| `THREADS` | `16` | `test-real`, `verify-branches`: thread counts, each must reproduce every golden |
| `BRANCHES` | none | `verify-branches`: branches or commits |
| `GOLDEN_COMMIT` | `71f6ec3` | `goldens`, `update-goldens`: the commit the goldens are made from |
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

### (a) Reproduce every golden from production's single-threaded release

```
git fetch origin high-perf                         # production's commit 71f6ec3 must be present
git submodule update --init tronko-assign/carquet
make -C tronko-assign/tests/real goldens MARKERS="mifish fwh" SETS="dev5k dev heldout"
```

This exports `71f6ec3` (with its carquet commit) from the repository, builds it plain and parquet
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
  `reads.json` to `data/real/<marker>/<set>/`, make goldens for it from `71f6ec3` with
  `make -C tronko-assign/tests/real update-goldens MARKERS=<marker> SETS=<set>`, and verify changes
  against those.
- **Build recipe.** The pipeline's recipe is written out in full above and in `real/lib.sh`
  (`build_pipeline`); nothing outside this repository is needed to build.
- **Machine.** Any Linux machine with gcc, make, zstd, zlib and enough memory (above). The goldens
  were made on x86-64 (Intel Cascade Lake, Ubuntu 24.04, gcc 13.3). Floating-point results can
  depend on the compiler and CPU architecture; if `make goldens` reports `DIFFER` on another
  machine, make goldens from `71f6ec3` there (`update-goldens` on a scratch branch) and verify
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

## `unit/test_vote_tally.c`

After placement, `runAssignmentOnChunk_WithBWA` tallies a read's votes: how many nodes of each
tree lie inside the score window, the list of voted nodes (`minNodes`), the tree with the most
votes (`maxRoot`, lowest index on a tie) and the ascending list of trees with votes
(`maxRoots`). It used to pass over every node of every tree in the reference, 7,086,092 nodes on
CO1_fwhF2_EPTDr2n, for every read; `vote_tally.h` now passes only over the read's candidate trees
(`vote_hit_trees()`, `vote_tally_hit_trees()`, `vote_reset_hit_trees()`). The old pass is kept
there as `vote_tally_all_trees()`.

The test runs both on random forests (1 to 4,000 trees, one-leaf trees, one large tree) and
random reads, with the buffers reused from read to read as in `tronko-assign`, and requires the
same `numMinNodes`, `minNodes`, per-tree counts, count of voting trees, `maxRoot`, `max` and
`maxRoots` for every read, and `voteRoot` all zero after each reset. The candidate lists come out
of tree order, with duplicate trees, with -1 and out-of-range ids, longer than the tree count and
empty; the votes include none, one node, sparse, dense, every node and equal counts in several
trees. It prints one line, for example:

```
test_vote_tally: 140000 reads on 6 forests; 55613 with more than one hit tree, ... ; without the sort 12713 results would change, without the duplicate removal 31529; 0 differ
```

The last counts show that the inputs exercise the rules: that many reads would get a different
result if the candidate trees were taken in hit order, or with duplicates.

## `integration/test_multitree_parity.sh` with `-DVOTE_SHADOW_CHECK`

Built with `make ARCH_FLAGS=-DVOTE_SHADOW_CHECK`, `tronko-assign` also recomputes every read's
vote tally with the old full pass and aborts on any difference, and prints one line of counters at
exit; the three-tree script passes that line on:

```
PASS: mt_paired, --number-of-cores 4
  VOTE_SHADOW reads=2000 multi_hit=1247 multi_vote=1184 unsorted=630 neg_in_prefix=0 dup_in_prefix=0 cap_reached=0 written_past_prefix=0
```

## `unit/test_ksw_extend.c`

`bwa_source_files/ksw.c` keeps BWA's seed-extension kernel as `ksw_extend2_scalar()` and makes
`ksw_extend2()` dispatch at run time to a vectorised kernel: AVX2 when the CPU has it, SSE2
otherwise (on aarch64 through `sse2neon.h`). The vectorised kernels must return the same six values
(`score`, `qle`, `tle`, `gtle`, `gscore`, `max_off`) as the scalar code for every input, or the
alignments BWA reports, and with them Tronko's output, change.

The test compares `ksw_extend2()`, `ksw_extend2_sse2()` and, on a CPU with AVX2,
`ksw_extend2_avx2()` with the expected values on:

- the 2,115 calls in `data/ksw_fixture_calls.kswd.gz`, inputs and outputs recorded from the
  unmodified kernel while `tronko-assign` ran the repository's fixture cases (every hundredth call
  of the paired, single and unpaired-reverse runs on the Charadriiformes example reference), so
  the file holds only repository data;
- 67,000 random calls from a fixed seed, compared with `ksw_extend2_scalar()`: general inputs, and
  inputs aimed at the rules a vectorised kernel breaks most easily (cells beyond the band that are
  read again when it grows, the z-drop's gap-extension factors, tiny bands, empty targets, and
  scores at and beyond the 16-bit bound where the scalar code takes over).

It prints one line, for example on an x86-64 CPU with AVX2:

```
test_ksw_extend: 2115 recorded calls, 67000 random calls, 4 kernels (with avx2), 276460 comparisons, 0 differ
```

On a CPU without AVX2, and on aarch64, the AVX2 kernel is not run and the line says
`3 kernels (without avx2)`. Run the test on an x86-64 machine with AVX2 after any change to the
AVX2 code.

`tools/record_ksw_calls.sh <output.kswd.gz>` records the calls again: it applies
`tools/ksw_record.patch` to a copy of the commit before the vectorised kernel, builds
`tronko-assign`, runs the three fixture cases at one thread and writes every hundredth
`ksw_extend2` call of each, with its inputs and outputs. The decompressed result equals the
committed file byte for byte:

```
bash tronko-assign/tests/tools/record_ksw_calls.sh /tmp/calls.kswd.gz
cmp <(gzip -dc /tmp/calls.kswd.gz) <(gzip -dc tronko-assign/tests/data/ksw_fixture_calls.kswd.gz)
```

## `unit/test_nodestore.c` and `integration/test_nodestore_paths.sh`

`nodestore.c` lays each tree's posteriors out node-innermost at load and scores reads in blocks
against node tiles. `unit/test_nodestore.c` (built by `tronko-assign`'s `make test` with the
binary's flags, about 2 s) checks, on random trees whose node counts straddle the tile widths and
on reads with sentinel columns, NaN and tied posteriors, `-1` positions, gaps, lowercase and other
characters, one- and two-mate reads and empty leaf strings, that every (read, node) score of the
blocked path equals production's `assignScores_Arr_paired`/`getscore_Arr` (linked unchanged) bit
for bit, for blocks of 1, 7, 64 and 1,000 reads and every kernel the CPU supports; that every leaf
string from the store equals `getSequenceInNodeWithoutNs`; and that `ns_reduce` equals the
reduction of `place_paired_with_nw`, including scores on the vote window's edge. Nine deliberately
broken variants of `nodestore.c` are all caught.

`integration/test_nodestore_paths.sh` runs the repository fixtures against their goldens: the
three read modes, the store self-check (`TRONKO_NODESTORE_CHECK=1`), blocks of 1 and 7 reads, a
1 MiB block budget, each supported kernel (`TRONKO_NS_KERNEL`), the legacy loop
(`TRONKO_NODESTORE=0`), and reads without candidates at block starts and ends against the legacy
loop of the same binary.

## `unit/test_nw_fill.c`

`needleman_wunsch_align()` aligns every candidate leaf to its read (production runs `-w`). The
three integer score matrices it fills decide the traceback, and with it the alignment strings, the
positions and the mismatch counts that placement scores. `alignment.c` keeps seq-align's original
fill, `alignment_fill_matrices()`, and adds `alignment_fill_matrices_fast()`, which `aligner_align()`
runs for every scoring without `no_mismatches`, `no_gaps_in_a` or `no_gaps_in_b` (tronko's
included): substitution scores from a per-alignment profile built with `scoring_lookup()` itself,
and each row in two vectorisable passes (match and gap-in-A, which read only the previous row) and
one serial pass (gap-in-B, which reads the cell to its left). It must produce the same integer in
every cell of every matrix, or tronko's output can change.

The test includes `alignment.c` (both fills are `static`) and, for every alignment, runs
`aligner_align()`, copies the three matrices, runs the original fill on the same aligner, and
compares all `(len_a + 1) × (len_b + 1)` cells of each matrix with `memcmp`. Inputs:

- `data/nw_fixture_alignments.txt.gz` (71 KB, one `leaf TAB read` per line): every fifth
  `needleman_wunsch_align()` call of the three production-parity fixture cases (single, unpaired
  reverse, paired) on the Charadriiformes example reference, 2,000 alignments of leaves of about
  311 bases and reads of 150, so repository data only; `tools/record_nw_alignments.sh 71f6ec3
  <file>` records the same file again, byte for byte (`data/PROVENANCE.md`);
- every pair of lengths 0 to 8, with every scoring and also as Smith-Waterman (which must take the
  original fill);
- 200,000 random alignments from a fixed seed (`make test`; 20,000 in `make test-san`, set
  `NW_RANDOM_SAN` for more), lengths 0 to 300 and leaves up to 1,174 bases (the longest CO1 leaf):
  ACGT; both cases with `N` and `-`; reads derived from the leaf with substitutions, indels, lower
  case and `N`; all ASCII bytes 1 to 127; and reads with 129 to 255 distinct byte values drawn
  from 1 to 255 (so bytes 0x80 to 0xFC too), against DNA leaves or leaves of any byte;
- scorings: tronko's (match 2, mismatch −1, gap open −3, gap extend −1, no start or end gap
  penalty, case-insensitive), each end-gap and start-gap setting, case-sensitive, a wildcard `N`
  and a swap table, another penalty scale, and `no_mismatches` and `no_gaps_in_a`, which must take
  the original fill.

It also checks that `aligner_align()` ran the fast fill where it should (the profile was built)
and the original fill where it should not, and it hashes the score and both alignment strings of
`needleman_wunsch_align2()` for the fixture alignments and the ASCII random alignments. Those
hashes are pinned to the values the original seq-align code produces: `alignment.c`,
`alignment_scoring.c` and `needleman_wunsch.c` as at `71f6ec3`, unchanged up to the parent of the
commit that adds the two-pass fill. To recompute them, compile the test against that code, from
`tronko-assign/tests`:

```
mkdir -p /tmp/nw-reference
git archive 71f6ec3 tronko-assign | tar -x -C /tmp/nw-reference
gcc -O3 -DNW_TEST_REFERENCE -I/tmp/nw-reference/tronko-assign -o /tmp/nw-reference/ref unit/test_nw_fill.c \
    /tmp/nw-reference/tronko-assign/alignment_scoring.c /tmp/nw-reference/tronko-assign/needleman_wunsch.c -lz
/tmp/nw-reference/ref data/nw_fixture_alignments.txt.gz
```

It prints the values that `EXPECTED_FIXTURE_HASH` and `checkpoint_expected` in the test hold:

```
test_nw_fill reference: fixture 0xcb40bf124bff058dULL, random after 20000 0x7a4f15eddcf82f66ULL, after 200000 0x2517715214af8866ULL (2000 fixture, 200000 random alignments)
```

The test itself prints one line, for example:

```
test_nw_fill: 2000 fixture, 2916 small, 200000 random alignments (19820 with 129+ distinct bytes); 207832 fills compared, 0 differ; fast fill not used 0, fallback misused 0; end-to-end hash vs the original code: fixture equal, random equal at 20000 and 200000
```

Bytes above 0x7F: where `char` is signed (x86-64), the unchanged `scoring_lookup()` indexes its
bitsets and tables with negative values for such bytes. That is undefined behaviour of the
original code, reached only by non-DNA input, and identical for both fills. The test places the
scoring struct after 8 KiB of zeroed memory so those reads stay inside its own allocation, and
`make test-san` compiles `alignment_scoring.c` with AddressSanitizer only; everything else,
`alignment.c` included, runs under both sanitizers.

## `integration/test_nw_shadow.sh`

Built with `-DNW_SHADOW`, `tronko-assign` recomputes every alignment the two-pass fill computes
with the original fill, compares the three matrices and aborts on any difference; at exit it
prints `NW_SHADOW alignments_compared=<n> mismatches=0`. The script runs the three
production-parity fixture cases at 1 and 4 threads (`TRONKO_ASSIGN_CORES`) with such a binary
(`TRONKO_ASSIGN_BIN`) and requires exit 0, the line with `n > 0`, and output identical to the
production goldens. `make test-shadow` builds `tests/tronko-assign-nw-shadow` with the main build's
compiler line (a minute or two) and runs it; `tronko-assign/tronko-assign` is left alone. The same
build can be pointed at any other read set.

## `unit/test_dense_sa.c`

`bwa_source_files/bwt.c` can replace the BWA index's suffix-array sample (every 32nd index, from
the `.sa` file) with a denser one built in memory when the index is loaded
(`bwt_densify_sa()`, called once per process from `main_mem` in `fastmap.c`): a `uint32` value
every 4th index, so that `bwt_sa()` walks about 3 steps of the BWT instead of 31 to turn a
suffix-array index into a text position. BWA's seeds, chains and hits, and so Tronko's output,
stay the same only if `bwt_sa()` returns exactly the value the original walk over the file sample
returns, for every index. The test checks that:

- on the repository's Charadriiformes index (`tronko-build/example_datasets/single_tree/`), whose
  `primary % 4` is 2, so walks through `primary` step into index 0 and read the `$` sentinel: every
  dense entry against the original lookup of its index, the dense lookup against the original for
  every index (910,485) at 4 threads and for 100,000 random indexes at the other thread counts and
  intervals, and the checking lookup (`TRONKO_CHECK_SA`) on a sixteenth of them;
- on 289 synthetic indexes whose suffix array the test computes itself by sorting the suffixes,
  independently of BWA (random, two-letter, periodic and single-letter texts of 1 to 270,000
  bases, file intervals 32, 16 and 8, dense intervals 2 to 16, 1 to 9 threads): every entry and
  every index against the true suffix array and the original lookup. All four residues of
  `primary % 4` occur, and one index has more file samples than the build has segments;
- the build's self-check: a corrupted file sample, at a segment start and elsewhere, makes it give
  up and leave the index on the file sample; it skips a `seq_len` at or above `UINT32_MAX`, a
  missing file sample, an interval not below the file's, and a second build;
- the checking mode computes both lookups and aborts on a corrupted dense entry;
- without the checking mode the file sample is freed and its pointer cleared, the lookups still
  equal those of an untouched copy of the index, and `bwt_destroy()` releases everything (the
  sanitizer build runs with leak detection).

It prints, on success:

```
fixture: seq_len 910484, primary 349034 (% 4 = 2), 18 dense walks through the $ sentinel
fixture: file sample 1 corrupted: status 4 (mismatches 1, unwritten 1), file sample kept
fixture: file sample 28452 corrupted: status 4 (mismatches 1, unwritten 0), file sample kept
fixture: TRONKO_CHECK_SA lookup aborts on a corrupted dense entry
synthetic: 289 indexes, primary % 4 residues seen 0xf, 3739 dense walks through the $ sentinel
test_dense_sa: 6339561 entries and lookups compared, 0 failures
```

## `unit/test_chain_flt.c`

BWA keeps or drops each chain of seeds in `mem_chain_flt()` (`bwa_source_files/bwamem.c`) before it
extends the kept ones. `tronko-assign` marks the chains with `mem_chain_flt_grouped()`, which compares
a chain with one representative per group of kept chains that share their query interval and ALT
flag, instead of with every kept chain. The two must set the same `.kept` and `.first` on every
chain and return the same chains in the same order, or BWA extends different chains and the
assignments change.

The test includes `bwamem.c` and compares the current code with a verbatim copy of BWA's original
function on:

- `data/chain_flt_fixture.chf.gz`: 288 recorded `mem_chain_flt()` calls (26,769 chains) from the
  repository's three fixture cases (every 40th call, and every call with 250 chains or more), with
  the output the original binary produced; both functions must reproduce it;
  `tools/record_chain_flt.sh 71f6ec3 <file>` records the same calls again (`data/PROVENANCE.md`);
- 20,000 random chain sets through the whole function (weights from the seeds, the sort, the
  marking, the cap and the compaction), with options away from BWA's defaults, among them a
  negative and a NaN `drop_ratio`, which take the original loop;
- 50,000 random chain lists with weights set directly, through the marking step alone: ALT chains,
  weights of 2^24 and above with dense ties, non-default `mask_level`, `drop_ratio`, `min_seed_len`
  and `max_chain_gap`, more than 64 groups and more than 4,096 chains (the heap paths);
- six constructed lists (a break at the very first kept chain, a group opened after others were
  covered, a group refilled after all its members were covered, ALT chains on either side).

`./test_chain_flt data/chain_flt_fixture.chf.gz 10` runs ten times as many random calls. It prints
one line, for example:

```
test_chain_flt: 288 recorded calls (26769 chains), 20000 random whole-function calls, 50000 random marking calls, 6 constructed; 120588 comparisons, 0 differ
```

## `unit/test_pos2rid.c`

`bns_pos2rid()` (`bwa_source_files/bntseq.c`) maps a reference position to the accession that holds
it. `tronko-assign` builds a table of the answer at the start of every 256 bases when the index is
loaded and steps forward from it, instead of binary-searching every accession. The test compares
`bns_pos2rid()` and `bns_intv2rid()` with verbatim copies of BWA's originals on the repository's
fixture index (every position, every bucket, intervals at every accession boundary on both strands),
on 2,000 synthetic references with zero-length accessions and lengths around and at multiples of a
bucket (every position), on references where no table may be built, and on the fixture index
converted to BWA's shared-memory layout, which must keep BWA's own layout and use the binary search.
It takes the fixture index prefix as its argument (the Makefile passes it) and prints one line:

```
test_pos2rid: fixture index 455248 positions and its accession boundaries, 262186522 synthetic positions; 264669519 comparisons, 0 differ
```
