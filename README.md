# tronko

A phylogeny-based method for assigning species to metabarcoding reads. It has two programs:

- `tronko-build` builds a reference database from trees, alignments and taxonomy.
- `tronko-assign` assigns reads to species using that database.

## Compile

```
sudo apt-get install -y zlib1g-dev libzstd-dev    # macOS: brew install zlib zstd gcc
make -C tronko-build
make -C tronko-assign && make -C tronko-assign bwa-mem3
```

## Build a database from an rCRUX release

Build the image once per tronko commit:

```
docker build -f Dockerfile.release --build-arg TRONKO_REF=<commit or tag> -t tronko-release:<ref> .
```

Build a database (it lands in `./work/out`):

```
docker run --rm -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -v "$PWD/work:/work" \
    tronko-release:<ref> --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0 \
    --work /work --steps fetch,build --threads 32
```

To publish it too, drop `--steps` and add a clean tronko checkout and a git identity:

```
    -e GIT_AUTHOR_NAME -e GIT_AUTHOR_EMAIL -e GIT_COMMITTER_NAME -e GIT_COMMITTER_EMAIL \
    -v "$HOME/tronko-release:/release"   ...   --tronko-checkout /release
```

Then push the `publish/...` branch from that checkout and open a PR. Merging it tags the release `<marker>_tronko_<version>`.

Settings and their defaults are in `scripts/tronko_build.conf`; change one with `--set KEY=VALUE`. Everything else is in `scripts/release_tronko.sh --help`.

## Build a database from your own tree

```
tronko-build -l -t tree.nwk -m alignment.fasta -x taxonomy.txt -d outdir
```

The taxonomy file is `name<TAB>domain;phylum;class;order;family;genus;species`, one line per sequence, with the same names as the tree and the alignment. `tronko-build -h` lists every option.

## Assign reads

```
tronko-assign -r -6 -f reference_tree.trkb -a marker.fasta -p -1 forward.fasta -2 reverse.fasta -o results.txt
tronko-assign -r -6 -f reference_tree.trkb -a marker.fasta -s -g reads.fasta -o results.txt
```

Add `-q` for FASTQ input and `-z` to reverse-complement the second read. `tronko-assign -h` lists every option.

## Test

```
tests/test_tronko_build_golden.sh
./run_tests.sh
```

## More

- [docs/usage.md](docs/usage.md): input file formats, building from several clusters, output columns and the example data.
- [docs/](docs/README.md): how the build and the assignment work.
- Cite: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg)](https://doi.org/10.5281/zenodo.7407318)
