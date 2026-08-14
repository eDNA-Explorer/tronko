#!/usr/bin/env python3
"""Draw a seeded, pair-preserving subsample of a paired FASTA cell.

The same index set is applied to both mates so pairing is preserved, and the
draw is seeded so the pre and post arms see byte-identical input.  Writes
pF.fasta / pR.fasta in place of the full files (originals kept as *.full).
"""

import random
import sys
from pathlib import Path


def read_records(path):
    recs, header, chunks = [], None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(chunks)))
                header, chunks = line.rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if header is not None:
        recs.append((header, "".join(chunks)))
    return recs


def main():
    cell = Path(sys.argv[1])
    cap = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1234

    f = read_records(cell / "pF.fasta")
    r = read_records(cell / "pR.fasta")
    if len(f) != len(r):
        raise SystemExit(f"mate count mismatch: {len(f)} vs {len(r)}")
    if len(f) <= cap:
        print(f"{cell.name}: {len(f):,} pairs <= cap {cap:,}; leaving full")
        return 0

    idx = sorted(random.Random(seed).sample(range(len(f)), cap))
    for tag, recs in (("F", f), ("R", r)):
        src = cell / f"p{tag}.fasta"
        if not (cell / f"p{tag}.fasta.full").exists():
            src.rename(cell / f"p{tag}.fasta.full")
        with (cell / f"p{tag}.fasta").open("w") as fh:
            for i in idx:
                h, s = recs[i]
                fh.write(f"{h}\n{s}\n")
    print(f"{cell.name}: {len(f):,} -> {cap:,} pairs (seed {seed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
