#!/usr/bin/env python3
"""Seeding, filtering, or placement? Name the mechanism behind a lost call.

When an ASV that BWA assigned comes back worse (or unassigned) under minimap2,
there are three distinct causes with three different fixes:

  seeding    minimap2 never produced the leaf as a candidate at all
  filtering  the leaf was a candidate but was dropped by the score == best
             filter or the best_n cap
  placement  the leaf survived as a candidate and was simply not chosen

This repo already carries evidence pointing at seeding for CO1 -- f0a608a
measured minimap2 434/500 against BWA 475/500, and 66a545c's k-mer sweep put
the shipped k11/w3 default at 461 against BWA's 475 -- but nobody has checked
whether the deficit is concentrated in arthropods. This script answers that.

Inputs, and how to produce them:

  --leaf-map      tronko-assign -5 FILE ...    (accession, tree, node per leaf;
                  note -5 prints the map and then exits, so it is its own
                  one-off invocation and cannot share a pass with -P)
  --bwa-align     the -P alignment output of the BWA arm
  --minimap2-align  the -P alignment output of the minimap2 arm
  --joined        bioscape_co1_aligner_ab.py's *_joined.csv

Usage:
    bioscape_co1_hit_forensics.py --leaf-map leafmap.txt \
        --bwa-align bwa.aln --minimap2-align mm2.aln \
        --joined bioscape_co1_joined.csv --out-dir OUT
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

# tronko caps candidates at --max-leaf-matches; more than this from a parse
# means the parse is wrong, not that tronko misbehaved.
MAX_LEAF_MATCHES = 10


def load_leaf_map(path: Path) -> dict[str, tuple[int, int]]:
    """`-5` output: accession \t tree \t node, one line per leaf."""
    out: dict[str, tuple[int, int]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, tree, node = parts[0], parts[1], parts[2]
            try:
                out[acc] = (int(tree), int(node))
            except ValueError:
                continue
    return out


def load_candidates(path: Path) -> dict[str, set[str]]:
    """`-P` alignment output: header lines are `>readname\treference\tlength`."""
    cands: dict[str, set[str]] = defaultdict(set)
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\n").split("\t")
            if len(parts) >= 2:
                cands[parts[0]].add(parts[1])
    return dict(cands)


def classify(
    chosen_bwa: tuple[int, int] | None,
    mm2_leaves: set[tuple[int, int]],
    mm2_accs: set[str],
) -> tuple[str, str]:
    if not mm2_accs:
        return "seeding", "minimap2 produced no candidate at all"
    if chosen_bwa is None:
        return "unknown", "BWA's chosen leaf is not resolvable from the leaf map"
    if chosen_bwa in mm2_leaves:
        return "placement", "BWA's leaf was a minimap2 candidate but was not chosen"
    if any(t == chosen_bwa[0] for t, _ in mm2_leaves):
        return "filtering", (
            "BWA's leaf absent but its tree is represented -- dropped by the "
            "score == best filter or the best_n cap"
        )
    return "seeding", "BWA's leaf and its whole tree are absent from minimap2's candidates"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--leaf-map", required=True, type=Path)
    ap.add_argument("--minimap2-align", required=True, type=Path)
    ap.add_argument("--bwa-align", type=Path)
    ap.add_argument("--joined", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="bioscape_co1")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    leaf_map = load_leaf_map(args.leaf_map)
    by_coord: dict[tuple[int, int], str] = {v: k for k, v in leaf_map.items()}
    mm2 = load_candidates(args.minimap2_align)
    bwa = load_candidates(args.bwa_align) if args.bwa_align else {}
    print(f"leaf map: {len(leaf_map):,} leaves")
    print(f"minimap2 candidates: {len(mm2):,} reads")
    if bwa:
        print(f"BWA candidates:      {len(bwa):,} reads")

    # candidate-set sanity: a parse bug shows up here, not as a tronko problem
    over = [r for r, c in mm2.items() if len(c) > MAX_LEAF_MATCHES]
    if over:
        print(f"WARNING: {len(over):,} reads show more than {MAX_LEAF_MATCHES} "
              f"candidates; the -P parse is probably wrong", file=sys.stderr)

    with args.joined.open() as fh:
        rows = list(csv.DictReader(fh))

    def as_int(v):
        try:
            return int(float(v))
        except (TypeError, ValueError):
            return None

    verdicts = Counter()
    detail_rows = []
    for r in rows:
        # the ASVs of interest: November/BWA had a call, minimap2 lost or changed it
        lost = r.get("lost_by_minimap2") in ("True", "true", "1")
        if not lost:
            continue
        name = r["readname"]
        c_tree, c_node = as_int(r.get("C_tree")), as_int(r.get("C_node"))
        chosen = (c_tree, c_node) if c_tree is not None and c_node is not None and c_tree >= 0 else None
        # fall back to the November record when the local BWA arm is absent
        if chosen is None:
            a_tree, a_node = as_int(r.get("A_tree")), as_int(r.get("A_node"))
            chosen = (a_tree, a_node) if a_tree is not None and a_node is not None and a_tree >= 0 else None

        mm2_accs = mm2.get(name, set())
        mm2_leaves = {leaf_map[a] for a in mm2_accs if a in leaf_map}
        verdict, why = classify(chosen, mm2_leaves, mm2_accs)
        verdicts[verdict] += 1
        detail_rows.append({
            "readname": name,
            "verdict": verdict,
            "why": why,
            "is_arthropod": r.get("is_arthropod"),
            "bwa_chosen_leaf": by_coord.get(chosen) if chosen else None,
            "bwa_chosen_tree_node": f"{chosen[0]}:{chosen[1]}" if chosen else None,
            "minimap2_n_candidates": len(mm2_accs),
            "bwa_n_candidates": len(bwa.get(name, set())) if bwa else None,
            "A_path": r.get("A_path"),
            "C_path": r.get("C_path"),
            "B_path": r.get("B_path"),
        })

    arth = [d for d in detail_rows if d["is_arthropod"] in ("True", "true", "1")]
    non = [d for d in detail_rows if d not in arth]

    def tally(rs):
        c = Counter(d["verdict"] for d in rs)
        n = len(rs)
        return {k: {"n": v, "pct": round(100.0 * v / n, 1) if n else None}
                for k, v in c.most_common()}

    report = {
        "label": args.label,
        "lost_calls_examined": len(detail_rows),
        "verdicts": tally(detail_rows),
        "arthropoda": {"n": len(arth), **{"verdicts": tally(arth)}},
        "non_arthropoda": {"n": len(non), **{"verdicts": tally(non)}},
        "candidate_parse_warning": len(over),
        "interpretation": {
            "seeding": "fix lives in --minimap2-kmer / --minimap2-window or the "
                       "mm_set_opt(\"sr\") preset; sweep k and w before concluding",
            "filtering": "fix lives in the score == best filter or the best_n cap",
            "placement": "the aligner is exonerated; this is scoring",
        },
    }

    (out / f"{args.label}_forensics.json").write_text(json.dumps(report, indent=2))
    if detail_rows:
        with (out / f"{args.label}_forensics.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(detail_rows[0].keys()))
            w.writeheader()
            w.writerows(detail_rows)

    print(json.dumps(report, indent=2))
    if not detail_rows:
        print("\nNo lost calls found -- nothing to explain, which is itself a result.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
