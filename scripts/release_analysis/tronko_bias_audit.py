#!/usr/bin/env python3
"""Is the change an even improvement, or does it favour particular clades?

Rachel's concern is that a fix validated on one marker and one taxonomic group
may not behave the same elsewhere.  This reads the per-cell ASV verdict tables
produced by tronko_blast_arbiter.py and breaks the result down two ways:

  * multi-marker within one dataset  (e.g. BioSCape's markers, same samples)
  * multi-dataset within one marker  (e.g. CO1 across CAT / Jalama / BioSCape)

and then tests whether the direction of change tracks covariates that would
indicate systematic bias rather than a uniform improvement: how well the clade
is represented in the reference, how abundant the ASV is, and how good its top
BLAST hit was.

Usage:
    tronko_bias_audit.py --verdicts LABEL=path.csv [LABEL=path.csv ...] \
        [--taxonomy ref_taxonomy.txt] --out-dir OUT
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tronko_blast_arbiter import BANDS  # noqa: E402

MIN_N = 30  # below this a per-clade share is noise, and is reported as such


def share(rows) -> dict:
    pb = sum(1 for r in rows if r["verdict"] == "post_better")
    prb = sum(1 for r in rows if r["verdict"] == "pre_better")
    tie = sum(1 for r in rows if r["verdict"] == "tie")
    unv = sum(1 for r in rows if r["verdict"] == "unverifiable")
    ext = [r for r in rows if r["change"] == "extension"]
    ext_contra = sum(1 for r in ext if r["detail"] == "extension_contradicted")
    return {
        "n": len(rows),
        "post_better": pb,
        "pre_better": prb,
        "tie": tie,
        "unverifiable": unv,
        "post_better_share": round(pb / (pb + prb), 4) if (pb + prb) else None,
        "decided": pb + prb,
        "extensions": len(ext),
        "extension_contradicted": ext_contra,
        "extension_contradicted_share": (
            round(ext_contra / len(ext), 4) if ext else None
        ),
        "underpowered": (pb + prb) < MIN_N,
    }


def load(path: Path) -> list[dict]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["abundance"] = int(r["abundance"])
        r["is_control"] = r["is_control"] == "True"
        r["blast_pident"] = float(r["blast_pident"]) if r["blast_pident"] else None
    return [r for r in rows if not r["is_control"]]


def clade_of(row: dict, rank_index: int) -> str | None:
    """Clade name at a rank, preferring the BLAST lineage as the neutral label."""
    for key in ("blast_lineage", "post_path", "pre_path"):
        v = row.get(key)
        if not v:
            continue
        parts = v.split(";")
        if rank_index < len(parts) and parts[rank_index] not in ("NA", ""):
            return parts[rank_index]
    return None


def reference_representation(taxonomy: Path) -> dict[str, int]:
    """How many reference sequences each clade has, by name at any rank."""
    counts: Counter[str] = Counter()
    with taxonomy.open() as fh:
        for line in fh:
            _, _, lineage = line.partition("\t")
            for name in lineage.strip().split(";"):
                if name and name != "NA":
                    counts[name] += 1
    return dict(counts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verdicts", nargs="+", required=True,
                    help="LABEL=path.csv, LABEL formatted dataset__marker")
    ap.add_argument("--taxonomy", type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--clade-rank", type=int, default=2,
                    help="0=domain 1=phylum 2=class (default) 3=order")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    cells: dict[str, list[dict]] = {}
    for spec in args.verdicts:
        label, _, path = spec.partition("=")
        cells[label] = load(Path(path))

    report: dict = {"cells": {}, "by_marker": {}, "by_dataset": {}}

    for label, rows in cells.items():
        report["cells"][label] = share(rows)

    # multi-dataset within a marker, and multi-marker within a dataset
    by_marker: dict[str, list] = defaultdict(list)
    by_dataset: dict[str, list] = defaultdict(list)
    for label, rows in cells.items():
        dataset, _, marker = label.partition("__")
        by_marker[marker or label].extend(rows)
        by_dataset[dataset].extend(rows)
    report["by_marker"] = {k: share(v) for k, v in by_marker.items()}
    report["by_dataset"] = {k: share(v) for k, v in by_dataset.items()}

    # per-clade, per cell
    clade_rows: list[dict] = []
    for label, rows in cells.items():
        groups: dict[str, list] = defaultdict(list)
        for r in rows:
            c = clade_of(r, args.clade_rank)
            if c:
                groups[c].append(r)
        for clade, sub in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            s = share(sub)
            clade_rows.append({"cell": label, "clade": clade, **s})

    if clade_rows:
        with (out / "clade_breakdown.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(clade_rows[0].keys()))
            w.writeheader()
            w.writerows(clade_rows)

    # covariates
    pooled = [r for rows in cells.values() for r in rows]

    by_band = {}
    for name, _, _ in BANDS:
        sub = [r for r in pooled if r["band"] == name]
        if sub:
            by_band[name] = share(sub)

    def abundance_bin(a: int) -> str:
        if a == 1:
            return "1"
        if a <= 5:
            return "2-5"
        if a <= 50:
            return "6-50"
        if a <= 500:
            return "51-500"
        return ">500"

    by_abund = {}
    for b in ("1", "2-5", "6-50", "51-500", ">500"):
        sub = [r for r in pooled if abundance_bin(r["abundance"]) == b]
        if sub:
            by_abund[b] = share(sub)

    by_repr = {}
    if args.taxonomy and args.taxonomy.exists():
        rep = reference_representation(args.taxonomy)

        def rep_bin(r: dict) -> str | None:
            c = clade_of(r, args.clade_rank)
            if c is None:
                return None
            n = rep.get(c, 0)
            if n < 10:
                return "<10"
            if n < 100:
                return "10-99"
            if n < 1000:
                return "100-999"
            return ">=1000"

        for b in ("<10", "10-99", "100-999", ">=1000"):
            sub = [r for r in pooled if rep_bin(r) == b]
            if sub:
                by_repr[b] = share(sub)

    report["covariates"] = {
        "by_identity_band": by_band,
        "by_asv_abundance": by_abund,
        "by_reference_representation": by_repr,
    }

    # the bias flag: any well-powered clade whose direction inverts the pooled one
    pooled_share = share(pooled)["post_better_share"]
    inversions = [
        r for r in clade_rows
        if not r["underpowered"]
        and r["post_better_share"] is not None
        and pooled_share is not None
        and r["post_better_share"] < 0.5 <= pooled_share
    ]
    report["pooled"] = share(pooled)
    report["direction_inversions"] = sorted(
        inversions, key=lambda r: r["post_better_share"]
    )

    (out / "bias_audit.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(
        {
            "pooled": report["pooled"],
            "by_marker": report["by_marker"],
            "by_dataset": report["by_dataset"],
            "direction_inversions": [
                {k: r[k] for k in ("cell", "clade", "n", "post_better_share")}
                for r in report["direction_inversions"]
            ],
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
