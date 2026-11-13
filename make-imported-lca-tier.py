#!/usr/bin/env python3
"""Write an LCA tier for an imported reference, the input build-rcrux-marker-db.sh
--imported builds a tronko master from.

An imported reference (PLANiTS, UNITE, PR2, Greengenes2, SILVA, as the
eDNA Explorer assignment benchmarks import them) ships a species tier only:
<ref>/dedup/species.fasta plus species_taxonomy.txt, with one row per label, so
a sequence filed under several taxa has several rows. A tronko leaf holds
exactly one lineage, so this collapses each sequence's labels to their lowest
common ancestor -- the ranks every label agrees on, NA from the first rank they
differ -- as rCRUX's own lca tier does.

Writes <ref>/lca/lca.fasta (the species tier's records, unchanged),
<ref>/lca/lca_taxonomy.txt (accession, taxid 0, 7-rank lineage) and
<ref>/lca/lca_tier.json (input hashes, counts, the rank each collapse kept).

Usage: make-imported-lca-tier.py <imported reference dir>
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lca(lineages: list[list[str]]) -> list[str]:
    """Ranks every lineage agrees on, NA from the first rank they differ."""
    out: list[str] = []
    for ranks in zip(*lineages, strict=True):
        if any(rank != ranks[0] for rank in ranks):
            break
        out.append(ranks[0])
    return out + ["NA"] * (len(RANKS) - len(out))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    ref = Path(sys.argv[1]).resolve()
    fasta = ref / "dedup" / "species.fasta"
    taxonomy = ref / "dedup" / "species_taxonomy.txt"
    manifest = ref / "build_manifest.json"
    for path in (fasta, taxonomy, manifest):
        if not path.is_file():
            raise SystemExit(f"missing {path}")

    labels: dict[str, list[list[str]]] = {}
    with taxonomy.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise SystemExit(f"{taxonomy}:{number}: expected 2 columns")
            accession, lineage = fields
            ranks = lineage.split(";")
            if len(ranks) != len(RANKS) or any(not rank for rank in ranks):
                raise SystemExit(f"{taxonomy}:{number}: expected 7 non-empty ranks")
            labels.setdefault(accession, []).append(ranks)

    out_dir = ref / "lca"
    out_dir.mkdir(exist_ok=True)
    seen: set[str] = set()
    kept_rank: Counter[str] = Counter()
    n_multi = 0
    with (
        fasta.open() as source,
        (out_dir / "lca.fasta").open("w") as fasta_out,
        (out_dir / "lca_taxonomy.txt").open("w") as tax_out,
    ):
        for line in source:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                if accession in seen:
                    raise SystemExit(f"{fasta}: duplicate record {accession}")
                if accession not in labels:
                    raise SystemExit(f"{fasta}: {accession} has no taxonomy row")
                seen.add(accession)
                rows = labels[accession]
                lineage = lca(rows)
                if len(rows) > 1:
                    n_multi += 1
                    named = [rank for rank in lineage if rank != "NA"]
                    kept_rank[RANKS[len(named) - 1] if named else "none"] += 1
                tax_out.write(f"{accession}\t0\t{';'.join(lineage)}\n")
            fasta_out.write(line)
    missing = set(labels) - seen
    if missing:
        raise SystemExit(f"{len(missing)} taxonomy accessions have no sequence")

    record = {
        "schema": "imported-lca-tier-1",
        "reference": ref.name,
        "rule": "lca of every label of a sequence; NA from the first differing rank",
        "inputs": {
            "species.fasta": sha256(fasta),
            "species_taxonomy.txt": sha256(taxonomy),
            "build_manifest.json": sha256(manifest),
        },
        "sequences": len(seen),
        "multi_label_sequences": n_multi,
        "multi_label_kept_rank": dict(sorted(kept_rank.items())),
        "outputs": {
            name: sha256(out_dir / name) for name in ("lca.fasta", "lca_taxonomy.txt")
        },
    }
    (out_dir / "lca_tier.json").write_text(json.dumps(record, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: record[k]
                for k in (
                    "reference",
                    "sequences",
                    "multi_label_sequences",
                    "multi_label_kept_rank",
                )
            }
        )
    )


if __name__ == "__main__":
    main()
