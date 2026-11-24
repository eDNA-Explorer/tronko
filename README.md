# tronko

A phylogeny-based method for assigning species to metabarcoding reads. It has two programs:

- `tronko-build` builds a reference database from trees, alignments and taxonomy.
- `tronko-assign` assigns reads to species using that database.

`tronko-build -h` and `tronko-assign -h` list every option.

## Compile

```
sudo apt-get install -y zlib1g-dev libzstd-dev    # macOS: brew install zlib zstd gcc
make -C tronko-build
make -C tronko-assign && make -C tronko-assign bwa-mem3
```

## Build a database from an rCRUX release

One command builds a database on a Jetstream2 VM and opens its release PR:

```
ops/release/launch.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0
```

It creates the VM, runs the release image there, prints the PR and deletes the VM. Merging the PR tags the release `<marker>_tronko_<version>`. You can close your laptop while it runs: `launch.sh list` shows the VMs, `launch.sh logs NAME` the build output, and `launch.sh wait NAME` picks one up again. A VM whose build failed is kept, powered off, until `launch.sh delete NAME`.

It takes every build setting below as an option, plus `--ref` (the tronko commit to build with, default `main`), `--flavor`, `--volume` and `--dry-run`. `launch.sh --help` lists what it needs on your machine.

The image is `Dockerfile.release`: one tronko commit and every tool the build needs, at fixed versions. CI pushes one for each commit on `main`; for another commit, run `gh workflow run build-push-tronko-release.yml -f ref=<commit>`.

To build without a VM (it lands in `./work/out`, nothing is published):

```
docker build -f Dockerfile.release --build-arg TRONKO_REF=<commit or tag> -t tronko-release .
docker run --rm -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -v "$PWD/work:/work" \
    tronko-release --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0 \
    --work /work --steps fetch,build --threads 32
```

What each step does:

| Step | What happens |
|---|---|
| fetch | Reads the curated LCA FASTA and taxonomy from rcrux-py's DVC remote at `--rcrux-tag`. |
| build | Runs `build-tronko-db.sh`: clustering, alignment, trees, `tronko-build`, `tronko-convert`, bwa and BWA-MEM3 indexes. |
| publish | Pushes the database to this repository's `tronko` DVC remote and commits `databases/<marker>/tronko.dvc` on a `publish/...` branch. |
| tag | On merge, the "Tag tronko releases" workflow tags the commit. The version is the next after the marker's existing tags, or `--tronko-version`. |

The image is `linux/amd64` only, so every release comes from one platform. An arm64 host runs it under emulation, which is only practical for a small test.

### Build settings

The defaults are our standard recipe, in `scripts/tronko_build.conf`:

| Setting | Default |
|---|---|
| `THREADS` | all CPUs |
| `CLUSTER_CUTOFF` | 25000 (cluster with AncestralClust above this many sequences) |
| `ANCESTRALCLUST_BIN_SIZE`, `ANCESTRALCLUST_MIN_DESCENDANTS` | 10000, 75 |
| `TREE_TOOL` | `veryfasttree` (or `fasttree`, `raxml`) |
| `PARTITION`, `SP_THRESHOLD`, `MIN_LEAF` | `sum_of_pairs`, 0.1, 500 |
| `SP_FORMULA` | `legacy` (or `corrected`) |
| `CLUSTER_JOBS`, `TRONKO_BUILD_JOBS` | 1, 1 |
| `ANCESTRALCLUST_ARGS`, `FAMSA_ARGS`, `TREE_ARGS`, `TRONKO_BUILD_ARGS`, `BWA_INDEX_ARGS` | empty (extra arguments for that tool) |

Change one for a run with `--set KEY=VALUE`, or several with `--config FILE`. `--print-config` shows what will be used. Every value used is recorded in the database's `metadata.json`.

## Build a database from your own trees

Input files:

- **Taxonomy:** `name<TAB>domain;phylum;class;order;family;genus;species`, one line per sequence, all 7 ranks filled (`NA` for an unknown one).
- **Alignment:** FASTA with each sequence on one line.
- **Names:** identical in the tree, the alignment and the taxonomy.

One tree:

```
tronko-build -l -t tree.nwk -m alignment.fasta -x taxonomy.txt -d outdir
```

Several clusters, numbered in one directory (`1_MSA.fasta`, `1_taxonomy.txt`, `RAxML_bestTree.1.reroot`, ...), split further by sum-of-pairs score or by leaf count:

```
tronko-build -y -e clusters_dir -n 2 -d outdir -s -u 0.1
tronko-build -y -e clusters_dir -n 2 -d outdir -v -f 500
```

Splitting needs `famsa` on PATH, plus `raxmlHPC-PTHREADS` and `nw_reroot`, or `FastTree` with `-a`. `-J N` splits N clusters at once.

The database is `outdir/reference_tree.txt`. `tronko-convert -i reference_tree.txt -o reference_tree.trkb` makes the binary form, which loads fastest.

## Assign reads

```
tronko-assign -r -6 -f reference_tree.trkb -a marker.fasta -p -1 forward.fasta -2 reverse.fasta -o results.txt
tronko-assign -r -6 -f reference_tree.trkb -a marker.fasta -s -g reads.fasta -o results.txt
```

- Reads can be FASTA or FASTQ (`-q`), plain, gzip or zstd.
- Orientation is not detected: `-z` reverse-complements the second read of a pair, `-v` a single read.
- `-6` uses the BWA-MEM3 index beside the reference FASTA, and builds it there if it is missing.

The output is tab-separated, one row per read:

```
Readname	Taxonomic_Path	Score	Forward_Mismatch	Reverse_Mismatch	Tree_Number	Node_Number
```

## Try it on the example data

```
tronko-build -l -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta \
    -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt \
    -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot -d outdir

tronko-assign -r -f outdir/reference_tree.txt \
    -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta -o results.txt -w
```

## Test

```
tests/test_tronko_build_golden.sh
./run_tests.sh
```

## More

- [docs/](docs/README.md): how the build and the assignment work, with step-by-step workflows in [docs/examples/](docs/examples/).
- [docs/performance-logging.md](docs/performance-logging.md) and [docs/crash-debugging.md](docs/crash-debugging.md): the `-V`, `-R`, `-T` and `-l` flags.
- [tests/README.md](tests/README.md): the test suites.
- Pipes L, Nielsen R (2024). A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. eLife. https://elifesciences.org/articles/85794
