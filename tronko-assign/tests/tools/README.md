# Fixture generators

These tools write the input files of the `tronko-assign` test fixtures from the repository's
example dataset, `tronko-build/example_datasets/single_tree/` (the Charadriiformes reference:
`Charadriiformes.fasta` and `reference_tree.txt`). Each input file they write is byte for byte
the committed one.

They do not make goldens. The goldens are made by the test scripts' record mode
(`TRONKO_TESTS_RECORD=1`, see [Regenerating the goldens](../README.md#regenerating-the-goldens)),
with `tronko-assign` built from production commit `71f6ec3`, at `--number-of-cores 1`.

## Regenerating every input

```
tronko-assign/tests/tools/make_fixture_inputs.sh <out-dir>
```

The script builds `tronko-convert` and `tronko-assign` from the checkout in a temporary directory,
writes every input under `<out-dir>` with the repository's layout, and then compares each file
with the committed one, printing `IDENTICAL` or `DIFFERENT`. It exits 1 if any file differs. It
needs bash, gcc, make, zlib and libzstd headers, the `zstd` command, python3 and awk, and takes
about a minute.

## What each tool makes

| Tool | Makes | Seed |
|---|---|---|
| `sim_reads.py` | `tests/data/assignment/single_4000.fasta` (150 bases, 2 substitutions each), `paired_2000_1.fasta`, `paired_2000_2.fasta` (first and reverse-complemented last 150 bases of a leaf, 2 substitutions each) | 7 |
| `multitree_fixture.py` | `tronko-assign/tests/data/multitree/`: the three-tree `reference_tree.txt` (converted to `reference_tree.trkb`), `multitree.fasta`, `single_F.fasta`, `single_R.fasta` | 7 |
| `gap_fixtures.py mate-rescue` | `tronko-assign/tests/data/gaps/mate-rescue/mr_1.fasta`, `mr_2.fasta` | 11 |
| `gap_fixtures.py read-content` | `tronko-assign/tests/data/gaps/read-content/rc_single.fasta`, `rc_1.fasta`, `rc_2.fasta` | 13 |
| `make_fixture_inputs.sh` | all of the above, and: both `reference_tree.trkb` files (`tronko-convert`), `gaps/prod-names/` (the `tests/data/assignment/` reads renamed, `12S_MiFish_U_paired_F_<idx>` and so on) and `gaps/prod-cmdline/` (the same reads compressed with `zstd -3 --no-check` from standard input) | |

The Python tools use the standard library only. Each prints its usage when run without
arguments. `sim_reads.py` draws 20,000 single-end reads before the pairs, as the original run
did, and keeps the first 4,000 reads and the first 2,000 pairs.

## Files reproduced byte for byte

`make_fixture_inputs.sh` reproduces these 20 committed files byte for byte (checked with Python
3.12, gcc 13.3 and zstd 1.5.5 on Linux):

- `tests/data/assignment/`: `reference_tree.trkb`, `single_4000.fasta`, `paired_2000_1.fasta`,
  `paired_2000_2.fasta`
- `tronko-assign/tests/data/multitree/`: `reference_tree.trkb`, `multitree.fasta`,
  `single_F.fasta`, `single_R.fasta`
- `tronko-assign/tests/data/gaps/mate-rescue/`: `mr_1.fasta`, `mr_2.fasta`
- `tronko-assign/tests/data/gaps/read-content/`: `rc_single.fasta`, `rc_1.fasta`, `rc_2.fasta`
- `tronko-assign/tests/data/gaps/prod-names/`: `pn_paired_F.fasta`, `pn_paired_R.fasta`,
  `pn_unpaired_F.fasta`, `pn_unpaired_R.fasta`
- `tronko-assign/tests/data/gaps/prod-cmdline/`: `single_4000.fasta.zst`,
  `paired_2000_1.fasta.zst`, `paired_2000_2.fasta.zst`

The `.zst` files depend on the `zstd` version; another version can compress differently while
decompressing to the same reads. The Python tools use `random.Random(seed)` with `randrange`,
`choice`, `sample` and `shuffle`; Python guarantees only `random()` to stay the same across
versions, so a future Python could draw different reads from the same seed.

Not made by these tools: `tests/data/assignment/single.fasta`, `paired_1.fasta` and
`paired_2.fasta`, the four reads and four pairs of the original parity test (pull request #10).
They are reads of example leaf `GU572157.1`, but not exact slices of it, and their generator is
not in the repository.
