#!/usr/bin/env python3
"""Inputs of the three-tree fixture, tronko-assign/tests/data/multitree/, from the repository's
example Charadriiformes reference.

Usage: multitree_fixture.py <example_dir> <out_dir>

<example_dir> is tronko-build/example_datasets/single_tree. Writes to <out_dir>:

  reference_tree.txt  three copies of the one example tree (convert with tronko-convert):
                      tree 0 is the original; tree 1 is the same tree whose leaves are only half
                      present in the BWA FASTA; tree 2 has its taxonomy lines rotated by one, so
                      the same leaf carries a different name in tree 2. Leaf names in trees 1 and
                      2 get the suffixes _1 and _2 (the header's name width, 13, still covers the
                      10-character accessions plus two).
  multitree.fasta     every tree-0 leaf, every second tree-1 leaf, every tree-2 leaf.
  single_F.fasta      2,000 single-end reads, forward orientation (run with -s -g)
  single_R.fasta      2,000 single-end reads, reverse-complemented (run with -s -v -g)

Read design. On this reference the same sequence sits in two or three trees, but BWA reports only
a primary and supplementary alignments (tronko-assign never sets -a and ignores XA), so a plain
single-end read yields one candidate tree. A read made of segments from different leaves, joined
out of reference order so that no single alignment spans them, yields one alignment per segment,
and each lands in whichever copy BWA's hash picks. The sets mix three kinds, drawn from one random
stream with seed 7:
  plain  (25 %) 150 bases of one leaf, 2 substitutions (as sim_reads.py): one candidate tree;
  chim2  (50 %) 75 + 75 bases of two leaves, the second segment upstream of the first in
                alignment coordinates, 1 substitution per segment: up to two candidate trees;
  chim3  (25 %) 50 + 50 + 50 bases of three leaves, each upstream of the one before: up to three,
                the most the three-tree reference allows (tronko-assign reads at most ntree hits).
Orientation. The BWA code in tronko-assign formats SAM SEQ from unconverted ASCII bases (it indexes
"ACGTN" with the character codes 65 to 84, reading past the literal), and the SAM text is then
parsed as a C string. In the builds this fixture was designed with (gcc on aarch64), the byte read
for a forward-strand T is NUL, so the text of a read aligned to the forward strand ends inside the
primary record's SEQ and its supplementary records are never seen; for the reverse strand only N
maps to NUL. Every chimeric read is therefore written so that it aligns to the reverse strand in
the mode it is run in (reverse-complemented in single_F.fasta, run with -s; forward in
single_R.fasta, which -v reverse-complements at load). Plain reads take either orientation at
random, so both SAM paths stay covered. Which strand is affected depends on the bytes the linker
places after the two literals, so the multi-candidate share is a property of the build, not only
of the reads (tronko-assign/tests/README.md, "Cases whose golden depends on the build").
Paired reads are tests/data/assignment/paired_2000_{1,2}.fasta, used as they are.

Standard library only. Goldens are not made here: see tools/README.md.
"""
import os
import random
import sys

SUFFIX = ["", "_1", "_2"]
K = 3


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_reference(example_dir, out_dir):
    lines = open(os.path.join(example_dir, "reference_tree.txt")).read().split("\n")
    if lines[-1] == "":
        lines = lines[:-1]
    assert int(lines[0]) == 1, "expected a one-tree reference"
    max_nodename, max_tax, max_line = int(lines[1]), lines[2], lines[3]
    numbase, root, numspec = map(int, lines[4].split("\t"))
    tax = lines[5:5 + numspec]
    body = lines[5 + numspec:]
    leaf_names = []
    with open(os.path.join(out_dir, "reference_tree.txt"), "w") as f:
        # widest leaf name after suffixing must fit the header's name width
        width = max_nodename
        f.write("%d\n%d\n%s\n%s\n" % (K, width, max_tax, max_line))
        for _ in range(K):
            f.write("%d\t%d\t%d\n" % (numbase, root, numspec))
        for t in range(K):
            for x in (tax if t != 2 else tax[1:] + tax[:1]):
                f.write(x + "\n")
        for t in range(K):
            i = 0
            while i < len(body):
                hdr = body[i].split("\t")
                assert hdr[0] == "0", hdr
                hdr[0] = str(t)
                if int(hdr[2]) == -1 and int(hdr[3]) == -1:
                    hdr[8] = hdr[8] + SUFFIX[t]
                    assert len(hdr[8]) <= width, (hdr[8], width)
                    if t == 0:
                        leaf_names.append(hdr[8])
                f.write("\t".join(hdr) + "\n")
                for r in body[i + 1:i + 1 + numbase]:
                    f.write(r + "\n")
                i += 1 + numbase
    return numbase, numspec, len(leaf_names)


def write_fasta(example_dir, out_dir):
    recs = read_fasta(os.path.join(example_dir, "Charadriiformes.fasta"))
    n = 0
    with open(os.path.join(out_dir, "multitree.fasta"), "w") as f:
        for t in range(K):
            for k, (name, seq) in enumerate(recs):
                if t == 1 and k % 2 == 1:
                    continue
                f.write(">%s%s\n%s\n" % (name, SUFFIX[t], seq))
                n += 1
    return recs, n


COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def mutate(rng, s, k):
    s = list(s)
    for _ in range(k):
        i = rng.randrange(len(s))
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def make_reads(rng, recs, n, prefix):
    """n reads as (name, forward sequence, kind, strand to align to); ungapped leaf coordinates."""
    seqs = [(name, seq.replace("-", "").upper()) for name, seq in recs]
    seqs = [(nm, s) for nm, s in seqs if len(s) >= 250 and set(s) <= set("ACGT")]
    out = []
    for i in range(n):
        u = rng.random()
        if u < 0.25:
            nm, s = rng.choice(seqs)
            st = rng.randrange(0, len(s) - 150)
            orient = rng.choice(("fwd", "rev"))
            out.append(("%s%d_plain_%s_%s" % (prefix, i, orient, nm), mutate(rng, s[st:st + 150], 2), "plain", orient))
        elif u < 0.75:
            (x, sx), (y, sy) = rng.sample(seqs, 2)
            # segment A from x at [p, p+75), segment B from y at [q, q+75) with q+75 <= p - 25
            p = rng.randrange(100, len(sx) - 75 + 1)
            q = rng.randrange(0, p - 100 + 1)
            a, b = sx[p:p + 75], sy[q:q + 75]
            read = mutate(rng, a, 1) + mutate(rng, b, 1)
            out.append(("%s%d_chim2_%s_%s" % (prefix, i, x, y), read, "chim2", "rev"))
        else:
            (x, sx), (y, sy), (z, sz) = rng.sample(seqs, 3)
            # three 50-base segments, each ending at least 20 bases before the previous one starts
            p = rng.randrange(140, len(sx) - 50 + 1)
            q = rng.randrange(70, p - 70 + 1)
            r = rng.randrange(0, q - 70 + 1)
            read = mutate(rng, sx[p:p + 50], 1) + mutate(rng, sy[q:q + 50], 1) + mutate(rng, sz[r:r + 50], 1)
            out.append(("%s%d_chim3_%s_%s_%s" % (prefix, i, x, y, z), read, "chim3", "rev"))
    return out


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    example_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    numbase, numspec, nleaves = write_reference(example_dir, out_dir)
    recs, nfasta = write_fasta(example_dir, out_dir)
    assert nleaves == len(recs) == numspec
    rng = random.Random(7)
    kinds = {}
    # single_F is run with -s (reads used as written); single_R with -s -v (reverse-complemented at
    # load). A read meant to align to the reverse strand is written reverse-complemented in
    # single_F and as is in single_R.
    for label, prefix, flipped_at_load in (("single_F", "sf", False), ("single_R", "sr", True)):
        reads = make_reads(rng, recs, 2000, prefix)
        with open(os.path.join(out_dir, label + ".fasta"), "w") as f:
            for name, seq, kind, strand in reads:
                write_rc = (strand == "rev") != flipped_at_load
                f.write(">%s\n%s\n" % (name, rc(seq) if write_rc else seq))
                key = (label, kind + ("" if kind != "plain" else "_" + strand))
                kinds[key] = kinds.get(key, 0) + 1
    print("reference: 3 trees x %d leaves, %d alignment columns; FASTA %d records" % (numspec, numbase, nfasta))
    for (label, kind), c in sorted(kinds.items()):
        print("%s %s %d" % (label, kind, c))


if __name__ == "__main__":
    main()
