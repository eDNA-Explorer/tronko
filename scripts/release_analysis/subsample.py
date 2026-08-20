#!/usr/bin/env python3
"""Draw a seeded, pair-preserving subsample of a paired FASTA cell.

The same index set is applied to both mates so pairing is preserved, and the
draw is seeded so the pre and post arms see byte-identical input.  Writes
pF.fasta / pR.fasta in place of the full files (originals kept as *.full).
"""

import random
import sys
from pathlib import Path


def same_mate_pair(hf: str, hr: str) -> bool:
    """Do a forward and reverse header name the same mate pair?

    Two conventions are seen in real data: byte-identical headers (e.g.
    bioscape_co1_fetch_nov25.py's output), and production's actual
    `<marker>_paired_F_<i>` / `<marker>_paired_R_<i>` -- identical except
    for a single F/R token. Accept either: headers must be equal length and
    differ at exactly one character, which must be 'F' paired with 'R'.
    """
    if hf == hr:
        return True
    if len(hf) != len(hr):
        return False
    diffs = [(a, b) for a, b in zip(hf, hr) if a != b]
    return len(diffs) == 1 and set(diffs[0]) == {"F", "R"}


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

    def source(tag):
        # Re-running with a different cap must resample from the original set,
        # not from a previous subsample -- read the .full backup once it exists.
        full = cell / f"p{tag}.fasta.full"
        return full if full.exists() else cell / f"p{tag}.fasta"

    f = read_records(source("F"))
    r = read_records(source("R"))
    if len(f) != len(r):
        raise SystemExit(f"mate count mismatch: {len(f)} vs {len(r)}")
    mismatches = [i for i, ((hf, _), (hr, _)) in enumerate(zip(f, r))
                  if not same_mate_pair(hf, hr)]
    if mismatches:
        i = mismatches[0]
        raise SystemExit(
            f"mate identity mismatch at {len(mismatches):,} position(s) -- e.g. "
            f"record {i} is {f[i][0]!r} in pF.fasta but {r[i][0]!r} in pR.fasta; "
            "the two files are not aligned by mate"
        )
    if len(f) <= cap:
        print(f"{cell.name}: {len(f):,} pairs <= cap {cap:,}; leaving full")
        return 0

    idx = sorted(random.Random(seed).sample(range(len(f)), cap))
    for tag, recs in (("F", f), ("R", r)):
        full = cell / f"p{tag}.fasta.full"
        if not full.exists():
            (cell / f"p{tag}.fasta").rename(full)
        with (cell / f"p{tag}.fasta").open("w") as fh:
            for i in idx:
                h, s = recs[i]
                fh.write(f"{h}\n{s}\n")
    print(f"{cell.name}: {len(f):,} -> {cap:,} pairs (seed {seed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
