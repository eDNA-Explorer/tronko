#!/usr/bin/env python3
"""The NCBI taxid of every lineage path tronko-assign can print: <marker>_lineage_taxids.txt.

Written by build-tronko-db.sh beside <marker>_taxonomy.txt and <marker>_taxids.txt,
from those two and the NCBI taxdump the lineages were named from. Standard
library only, so the build runs it with any python3.

Usage:
    python3 lineage_taxids.py --marker MARKER --build-dir DIR --taxdump DIR
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path

TAXONOMY_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
TAXDUMP_FILES = ("nodes.dmp", "names.dmp", "merged.dmp")
# The NCBI ranks each of the seven lineage slots is named from (rCRUX's main ranks);
# a sub-rank of the same name (subgenus Drosophila) loses to them.
SLOT_RANKS = (
    ("superkingdom", "domain"), ("phylum",), ("class",), ("order",), ("family",), ("genus",), ("species",)
)


def lineage_taxids_name(marker: str) -> str:
    return f"{marker}_lineage_taxids.txt"


def read_taxdump(taxdump_dir: Path) -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    """Each taxid's parent and rank, and its scientific name. A merged taxid's parent is its new taxid."""
    parents, ranks, names = {}, {}, {}
    with open(taxdump_dir / "nodes.dmp") as handle:
        for line in handle:
            taxid, parent, rank = line.split("\t|\t", 3)[:3]
            parents[int(taxid)], ranks[int(taxid)] = int(parent), rank
    with open(taxdump_dir / "names.dmp") as handle:
        for line in handle:
            fields = line.split("\t|\t")
            if fields[3].startswith("scientific name"):
                names[int(fields[0])] = fields[1]
    with open(taxdump_dir / "merged.dmp") as handle:
        for line in handle:
            old, new = line.split("\t|\t")[:2]
            parents[int(old)] = int(new.rstrip("\t|\n"))
            ranks[int(old)] = "merged"
    return parents, ranks, names


def lineage_taxids(taxonomy: Path, taxids: Path, taxdump_dir: Path) -> tuple[list[tuple[str, str, str]], int]:
    """Every lineage path tronko-assign can print, with the NCBI taxid of its deepest named rank.

    tronko-assign prints a read's lineage as the first ranks of one reference
    lineage, so every prefix of every lineage is a key. Each name is found on
    its record's own NCBI path (from its taxid, or every taxid of a record
    rCRUX merged), so no name is looked up on its own and a homonym cannot be
    picked. A path ending in NA takes the taxid of its deepest named rank, and
    the rank says which. A path no record resolves has an empty taxid. Returns
    the rows (path, rank, taxid) and how many paths had records disagreeing.
    """
    parents, ranks, names = read_taxdump(taxdump_dir)

    def ancestry(taxid: int) -> list[int]:
        path = []
        while taxid in parents and taxid != 1:
            if ranks[taxid] != "merged":
                path.append(taxid)
            taxid = parents[taxid]
        return path

    # rCRUX keeps a merged record's rank when every taxid naming it agrees and
    # skips those naming none, so a taxid whose path lacks the name is skipped too.
    def find(paths: list[list[int]], name: str, slot: int) -> int | None:
        found = set()
        for path in paths:
            if hits := [node for node in path if names.get(node) == name]:
                found.add(next((n for n in hits if ranks[n] in SLOT_RANKS[slot]), hits[0]))
        return found.pop() if len(found) == 1 else None

    votes: dict[str, Counter] = defaultdict(Counter)
    with open(taxonomy) as lineages, open(taxids) as ids:
        for row, id_row in zip_longest(lineages, ids, fillvalue=""):
            accession, lineage = (row.rstrip("\n").split("\t") + [""])[:2]
            id_accession, id_list = (id_row.rstrip("\n").split("\t") + [""])[:2]
            if accession != id_accession:
                raise ValueError(f"{taxonomy.name} and {taxids.name} disagree at {accession or id_accession}")
            paths = [ancestry(int(taxid)) for taxid in id_list.split(",") if taxid.isdigit()]
            names_in = lineage.split(";")
            deepest = None
            for slot, name in enumerate(names_in):
                if name not in ("", "NA"):
                    taxid = find(paths, name, slot)
                    deepest = (TAXONOMY_RANKS[slot], taxid) if taxid else None
                votes[";".join(names_in[: slot + 1])][deepest] += 1
    rows, conflicts = [], 0
    for path, counts in sorted(votes.items()):
        resolved = [(choice, n) for choice, n in counts.most_common() if choice]
        conflicts += len(resolved) > 1
        rank, taxid = resolved[0][0] if resolved else ("", "")
        rows.append((path, rank, str(taxid)))
    return rows, conflicts


def write_lineage_taxids(build_dir: Path, marker: str, taxdump_dir: Path) -> list[str]:
    """Write <marker>_lineage_taxids.txt beside the taxonomy; the lines to log about it."""
    rows, conflicts = lineage_taxids(
        build_dir / f"{marker}_taxonomy.txt", build_dir / f"{marker}_taxids.txt", taxdump_dir
    )
    with open(build_dir / lineage_taxids_name(marker), "w") as handle:
        handle.writelines("\t".join(row) + "\n" for row in rows)
    # A path of NA only names nothing, so it has no taxid to find.
    missing = [path for path, _, taxid in rows if not taxid and set(path.split(";")) != {"NA"}]
    lines = [f"lineage taxids: {len(rows):,} paths, {len(missing):,} without a taxid, {conflicts:,} with records disagreeing"]
    lines += [f"WARNING: no NCBI taxid for lineage {path}" for path in missing[:20]]
    if conflicts:
        lines.append(f"WARNING: {conflicts:,} lineage paths have records naming different taxids; the commonest is kept")
    return lines



def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--marker", required=True)
    parser.add_argument("--build-dir", required=True, type=Path, help="Holds <marker>_taxonomy.txt and <marker>_taxids.txt")
    parser.add_argument("--taxdump", required=True, type=Path, help="Directory with nodes.dmp, names.dmp, merged.dmp")
    args = parser.parse_args(argv)
    try:
        print("\n".join(write_lineage_taxids(args.build_dir, args.marker, args.taxdump)))
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
