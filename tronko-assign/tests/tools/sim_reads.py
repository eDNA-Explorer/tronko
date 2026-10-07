#!/usr/bin/env python3
"""Simulated reads of tests/data/assignment/ (repository root), from the example Charadriiformes
reference.

Usage: sim_reads.py <Charadriiformes.fasta> <out_dir>

<Charadriiformes.fasta> is tronko-build/example_datasets/single_tree/Charadriiformes.fasta. Writes to
<out_dir>:

  single_4000.fasta        4,000 single-end reads: 150 bases from a random position of a random
                           leaf, 2 substitutions; names r<i>_<accession>
  paired_2000_1.fasta      2,000 pairs: mate 1 is the first 150 bases of a random leaf, mate 2 the
  paired_2000_2.fasta      reverse complement of its last 150 bases, 2 substitutions each; names
                           p<i>_<accession>

Leaves are used without alignment gaps ("-" removed), upper case, and only if at least 250 bases
long. One random stream, seed 7. The original run drew 20,000 single-end reads and then 10,000
pairs from that stream; the fixture keeps the first 4,000 reads and the first 2,000 pairs. This
script draws the same 20,000 reads, so that the pairs start at the same point of the stream, and
writes only the ones the fixture keeps.

Standard library only. Goldens are not made here: see tools/README.md.
"""
import os
import random
import sys

SEED = 7
N_SINGLE_DRAWN, N_SINGLE_KEPT = 20000, 4000
N_PAIRS_KEPT = 2000

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def read_leaves(path):
    seqs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs.append((name, "".join(buf).replace("-", "").upper()))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs.append((name, "".join(buf).replace("-", "").upper()))
    return [s for s in seqs if len(s[1]) >= 250]


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    fasta, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(SEED)
    seqs = read_leaves(fasta)

    def mutate(s, k=2):
        s = list(s)
        for _ in range(k):
            i = rng.randrange(len(s))
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
        return "".join(s)

    with open(os.path.join(out_dir, "single_4000.fasta"), "w") as f:
        for i in range(N_SINGLE_DRAWN):
            n, s = rng.choice(seqs)
            st = rng.randrange(0, len(s) - 150)
            read = mutate(s[st:st + 150])
            if i < N_SINGLE_KEPT:
                f.write(">r%d_%s\n%s\n" % (i, n, read))
    with open(os.path.join(out_dir, "paired_2000_1.fasta"), "w") as f1, \
            open(os.path.join(out_dir, "paired_2000_2.fasta"), "w") as f2:
        for i in range(N_PAIRS_KEPT):
            n, s = rng.choice(seqs)
            f1.write(">p%d_%s\n%s\n" % (i, n, mutate(s[:150])))
            f2.write(">p%d_%s\n%s\n" % (i, n, mutate(rc(s[-150:]))))


if __name__ == "__main__":
    main()
