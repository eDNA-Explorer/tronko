# tronko usage

`tronko-build -h` and `tronko-assign -h` list every option. This page covers what they don't.

## Input files for tronko-build

- **Taxonomy:** `name<TAB>domain;phylum;class;order;family;genus;species`, one line per sequence, all 7 ranks filled (`NA` for an unknown one).
- **Alignment:** FASTA with each sequence on one line, no line breaks.
- **Names:** identical in the tree, the alignment and the taxonomy.

For several clusters, put them in one directory, numbered:

```
1_MSA.fasta   1_taxonomy.txt   RAxML_bestTree.1.reroot
2_MSA.fasta   2_taxonomy.txt   RAxML_bestTree.2.reroot
```

## Building

One tree:

```
tronko-build -l -t tree.nwk -m alignment.fasta -x taxonomy.txt -d outdir
```

Several clusters, split further by sum-of-pairs score or by leaf count:

```
tronko-build -y -e clusters_dir -n 2 -d outdir -s -u 0.1
tronko-build -y -e clusters_dir -n 2 -d outdir -v -f 500
```

Splitting needs `famsa` on PATH, plus `raxmlHPC-PTHREADS` and `nw_reroot`, or `FastTree` with `-a`. `-J N` splits N clusters at once.

The database is `outdir/reference_tree.txt`. `tronko-convert -i reference_tree.txt -o reference_tree.trkb` makes the binary form tronko-assign loads fastest.

## Assigning

- Reads can be FASTA or FASTQ (`-q`), plain, gzip or zstd.
- tronko-assign does not detect read orientation: `-z` reverse-complements the second read of a pair, `-v` a single read.
- `-6` uses the BWA-MEM3 index beside the reference FASTA, and builds it there if it is missing.

The output is tab-separated, one row per read (the forward read's name for a pair):

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

## Further reading

- [examples/](examples/): the single-tree and multiple-tree workflows step by step.
- [performance-logging.md](performance-logging.md) and [crash-debugging.md](crash-debugging.md): the `-V`, `-R`, `-T` and `-l` flags.
- [../tests/README.md](../tests/README.md): the test suites.
- Pipes L, Nielsen R (2024). A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. eLife. https://elifesciences.org/articles/85794
