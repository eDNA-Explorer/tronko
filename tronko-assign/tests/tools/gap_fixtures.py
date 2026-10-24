#!/usr/bin/env python3
"""Read sets of two coverage-gap fixtures, tronko-assign/tests/data/gaps/mate-rescue/ and
tronko-assign/tests/data/gaps/read-content/, from the repository's example Charadriiformes
reference.

  gap_fixtures.py mate-rescue <Charadriiformes.fasta> <out-dir>     writes mr_1.fasta, mr_2.fasta
  gap_fixtures.py read-content <Charadriiformes.fasta> <out-dir>    writes rc_single.fasta,
                                                                     rc_1.fasta, rc_2.fasta

<Charadriiformes.fasta> is tronko-build/example_datasets/single_tree/Charadriiformes.fasta.

mate-rescue: pairs on the example Charadriiformes reference, written for production's paired
mode (-p -z: read 2 is given as the reverse complement and tronko-assign reverse-complements it
back), so that BWA estimates an insert-size distribution and then attempts mate rescue
(mem_matesw -> ksw_align2) for the pairs in which only one mate aligns. Classes, shuffled
together with seed 11:
  norm      340 pairs, both mates 150 bases, 1 substitution each, ends jittered by 0-15 bases;
            240 of them from the 33 leaves whose first and last 150 bases differ from every
            other leaf in at least 8 positions, because mem_pestat uses only pairs whose mates
            have no second hit within 80 % of the best score, and most Charadriiformes leaves
            have a near-identical sibling (the insert-size estimate; its spread keeps std > 0,
            unlike tests/data/assignment/, whose pairs all have insert 163)
  r2unseed   80 pairs, mate 2 has a substitution every 12th base: no 19-base exact seed, so BWA
            cannot align it, but a correct Smith-Waterman in the rescue window would (score ~85)
  r1unseed   60 pairs, the same with mate 1 unseedable
  r2random   60 pairs, mate 2 is random sequence
  r2long     40 pairs, mate 2 is 260 bases and unseedable: l_ms * a >= 250, so ksw_align2 takes
            the 16-bit kernel (ksw_i16)
  r2fwd      40 pairs, mate 2 written in forward orientation (not reverse-complemented), giving
            BWA a second orientation to estimate and rescue in (the is_rev branch)

read-content: see read_content() below; seed 13.

Standard library only. Goldens are not made here: see tools/README.md.
"""
import random
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs.append((name, "".join(buf).upper()))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs.append((name, "".join(buf).upper()))
    return seqs


def mutate(rng, s, k):
    s = list(s)
    for i in rng.sample(range(len(s)), k):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def unseedable(rng, s, step=12):
    s = list(s)
    for i in range(rng.randrange(step), len(s), step):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def unique_leaves(leaves, window=150, min_mm=8):
    """Leaves whose first and last `window` bases differ from every other leaf's, compared at
    the same offset, in at least `min_mm` positions."""
    # Each window is packed into one integer with one bit per (position, base); two windows then
    # agree at popcount(a & b) positions.
    def pack(w):
        v = 0
        for c in w:
            v = (v << 4) | (1 << "ACGT".index(c))
        return v

    ok = [True] * len(leaves)
    for part in (lambda s: s[:window], lambda s: s[-window:]):
        m = [pack(part(s)) for n, s in leaves]
        for i in range(len(m)):
            if not ok[i]:
                continue
            for j in range(len(m)):
                if j != i and window - bin(m[i] & m[j]).count("1") < min_mm:
                    ok[i] = False
                    break
    return [leaves[i] for i in range(len(leaves)) if ok[i]]


def mate_rescue(fasta, out):
    rng = random.Random(11)
    leaves = [(n, s) for n, s in read_fasta(fasta) if len(s) >= 300 and set(s) <= set("ACGT")]
    pairs = []

    def pick():
        n, s = rng.choice(leaves)
        s1 = rng.randrange(0, 16)
        e = len(s) - rng.randrange(0, 16)
        return n, s, s1, e

    unique = unique_leaves(leaves)
    for i in range(340):
        if i < 240:
            n, s = rng.choice(unique)
            s1, e = rng.randrange(0, 16), len(s) - rng.randrange(0, 16)
        else:
            n, s, s1, e = pick()
        pairs.append(("norm", n, mutate(rng, s[s1:s1 + 150], 1), rc(mutate(rng, s[e - 150:e], 1))))
    for i in range(80):
        n, s, s1, e = pick()
        pairs.append(("r2unseed", n, mutate(rng, s[s1:s1 + 150], 1), rc(unseedable(rng, s[e - 150:e]))))
    for i in range(60):
        n, s, s1, e = pick()
        pairs.append(("r1unseed", n, unseedable(rng, s[s1:s1 + 150]), rc(mutate(rng, s[e - 150:e], 1))))
    for i in range(60):
        n, s, s1, e = pick()
        pairs.append(("r2random", n, mutate(rng, s[s1:s1 + 150], 1), "".join(rng.choice("ACGT") for _ in range(150))))
    for i in range(40):
        n, s, s1, e = pick()
        pairs.append(("r2long", n, mutate(rng, s[s1:s1 + 150], 1), rc(unseedable(rng, s[e - 260:e]))))
    for i in range(40):
        n, s, s1, e = pick()
        pairs.append(("r2fwd", n, mutate(rng, s[s1:s1 + 150], 1), mutate(rng, s[e - 150:e], 1)))
    rng.shuffle(pairs)
    kinds = {}
    with open(out + "/mr_1.fasta", "w") as f1, open(out + "/mr_2.fasta", "w") as f2:
        for i, (kind, n, a, b) in enumerate(pairs):
            name = "mr%d_%s_%s" % (i, kind, n)
            f1.write(">%s\n%s\n" % (name, a))
            f2.write(">%s\n%s\n" % (name, b))
            kinds[kind] = kinds.get(kind, 0) + 1
    for k in ["norm", "r2unseed", "r1unseed", "r2random", "r2long", "r2fwd"]:
        print("%-9s %d pairs" % (k, kinds[k]))


def read_content(fasta, out):
    """Reads whose content production's QC lets through but the other fixtures lack: one N
    (bbduk maxns=1 keeps them), no alignment at all (non-target sequence), and short reads
    (bbduk minlen is the k-mer length). Seed 13.
      rc_single.fasta  single-end, for -s -g and -s -v -g: 150 plain 150-base reads, 150 with one
                       N, 100 random, 60 short (30-60 bases), shuffled
      rc_1/rc_2.fasta  pairs for -p -z: 150 plain, 100 with one N in one mate, 60 with a random
                       mate 2, 40 with both mates random, 50 with a short mate 2"""
    rng = random.Random(13)
    leaves = [(n, s) for n, s in read_fasta(fasta) if len(s) >= 300 and set(s) <= set("ACGT")]

    def seg(length=150):
        n, s = rng.choice(leaves)
        st = rng.randrange(0, len(s) - length)
        return n, s, st, mutate(rng, s[st:st + length], 1)

    def with_n(r):
        i = rng.randrange(len(r))
        return r[:i] + "N" + r[i + 1:]

    def rand(k):
        return "".join(rng.choice("ACGT") for _ in range(k))

    single = []
    for i in range(150):
        n, s, st, r = seg(); single.append(("plain", n, r))
    for i in range(150):
        n, s, st, r = seg(); single.append(("oneN", n, with_n(r)))
    for i in range(100):
        single.append(("random", "none", rand(150)))
    for i in range(60):
        n, s, st, r = seg(rng.randrange(30, 61)); single.append(("short", n, r))
    rng.shuffle(single)
    kinds = {}
    with open(out + "/rc_single.fasta", "w") as f:
        for i, (k, n, r) in enumerate(single):
            f.write(">rc%d_%s_%s\n%s\n" % (i, k, n, r))
            kinds["single " + k] = kinds.get("single " + k, 0) + 1

    def pair():
        n, s = rng.choice(leaves)
        s1 = rng.randrange(0, 16)
        e = len(s) - rng.randrange(0, 16)
        return n, mutate(rng, s[s1:s1 + 150], 1), rc(mutate(rng, s[e - 150:e], 1)), s, e

    pairs = []
    for i in range(150):
        n, a, b, s, e = pair(); pairs.append(("plain", n, a, b))
    for i in range(100):
        n, a, b, s, e = pair()
        if rng.random() < 0.5:
            a = with_n(a)
        else:
            b = with_n(b)
        pairs.append(("oneN", n, a, b))
    for i in range(60):
        n, a, b, s, e = pair(); pairs.append(("r2random", n, a, rand(150)))
    for i in range(40):
        pairs.append(("random", "none", rand(150), rand(150)))
    for i in range(50):
        n, a, b, s, e = pair()
        k = rng.randrange(30, 61)
        pairs.append(("r2short", n, a, rc(mutate(rng, s[e - k:e], 1))))
    rng.shuffle(pairs)
    with open(out + "/rc_1.fasta", "w") as f1, open(out + "/rc_2.fasta", "w") as f2:
        for i, (k, n, a, b) in enumerate(pairs):
            name = "rcp%d_%s_%s" % (i, k, n)
            f1.write(">%s\n%s\n" % (name, a))
            f2.write(">%s\n%s\n" % (name, b))
            kinds["pairs " + k] = kinds.get("pairs " + k, 0) + 1
    for k in sorted(kinds):
        print("%-16s %d" % (k, kinds[k]))


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("mate-rescue", "read-content"):
        print(__doc__)
        sys.exit(1)
    {"mate-rescue": mate_rescue, "read-content": read_content}[sys.argv[1]](sys.argv[2], sys.argv[3])
