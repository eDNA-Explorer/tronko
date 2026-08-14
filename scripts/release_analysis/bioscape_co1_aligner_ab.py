#!/usr/bin/env python3
"""BioSCape CO1: November BWA vs fixed minimap2, on identical November ASVs.

Three answers per ASV:
  A  November as recorded  (legacy pre-backfill TSV / BigQuery paths)
  C  local BWA             (high-perf @ 71f6ec3, the November-era build)
  B  local minimap2        (530caf9, post PR #8)

**A-vs-C is the fidelity gate and is not optional.** It establishes that the
harness reproduces November at all. Only once it passes does C-vs-B mean
"the aligner and placement delta on identical input" rather than "some
difference between two things we did not control".

The non-arthropod ASVs are a free internal control: if they move as much as the
arthropods, this is a global effect and not an invertebrate story at all.

Usage:
    bioscape_co1_aligner_ab.py --november A.tsv --bwa C.tsv --minimap2 B.tsv \
        --out-dir OUT
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tronko_blast_arbiter import depth_of, read_assignments, split_path  # noqa: E402

# The clades the colleague's report is about, plus the control.
CLADES = ("Arthropoda", "Insecta", "Diptera", "Chironomidae", "Coleoptera")


def path_of(rec: dict | None) -> str:
    return (rec or {}).get("path") or ""


def ranks(rec: dict | None) -> list[str]:
    return [r for r in path_of(rec).split(";") if r and r != "NA"]


def has_clade(rec: dict | None, clade: str) -> bool:
    return clade in ranks(rec)


def pct(n: int, d: int) -> float | None:
    return round(100.0 * n / d, 2) if d else None


def agreement(rows: list[dict], x: str, y: str) -> dict:
    n = len(rows)
    same_path = sum(1 for r in rows if path_of(r[x]) == path_of(r[y]))
    same_tree = sum(1 for r in rows if (r[x] or {}).get("tree") == (r[y] or {}).get("tree"))
    same_node = sum(
        1 for r in rows
        if (r[x] or {}).get("tree") == (r[y] or {}).get("tree")
        and (r[x] or {}).get("node") == (r[y] or {}).get("node")
    )
    return {
        "n": n,
        "path_pct": pct(same_path, n),
        "tree_pct": pct(same_tree, n),
        "tree_node_pct": pct(same_node, n),
    }


def arm_summary(rows: list[dict], key: str) -> dict:
    n = len(rows)
    d = [depth_of(split_path(path_of(r[key]))) for r in rows]
    return {
        "assigned_pct": pct(sum(1 for x in d if x > 0), n),
        "mean_depth": round(sum(d) / n, 3) if n else None,
        "species_level_pct": pct(sum(1 for x in d if x >= 7), n),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--november", required=True, type=Path, help="arm A, recorded")
    ap.add_argument("--bwa", required=True, type=Path, help="arm C, local BWA")
    ap.add_argument("--minimap2", required=True, type=Path, help="arm B, local minimap2")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="bioscape_co1")
    ap.add_argument(
        "--fidelity-floor",
        type=float,
        default=95.0,
        help="A-vs-C path agreement below this is reported as a FAILED gate",
    )
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    A = read_assignments(args.november)
    C = read_assignments(args.bwa)
    B = read_assignments(args.minimap2)

    common = sorted(set(A) & set(C) & set(B))
    rows = [{"readname": k, "A": A[k], "C": C[k], "B": B[k]} for k in common]
    if not rows:
        raise SystemExit("no overlapping readnames across the three arms")

    coverage = {
        "A": len(A), "C": len(C), "B": len(B), "joined": len(rows),
        "note": (
            "readname sets differ; investigate before trusting anything below"
            if not (len(A) == len(C) == len(B) == len(rows)) else "all arms cover the same set"
        ),
    }

    fidelity = agreement(rows, "A", "C")
    gate_pass = (fidelity["path_pct"] or 0) >= args.fidelity_floor

    report: dict = {
        "label": args.label,
        "coverage": coverage,
        "fidelity_gate_A_vs_C": {
            **fidelity,
            "floor_pct": args.fidelity_floor,
            "status": "PASS" if gate_pass else "FAIL",
            "meaning": (
                "harness reproduces November; C-vs-B is a controlled aligner comparison"
                if gate_pass else
                "harness does NOT reproduce November -- C-vs-B is indicative only. "
                "Escalate to 77ade9e + the text reference tree, or re-derive the "
                "November-era parameters, before quoting any delta below."
            ),
        },
        "comparisons": {
            "C_vs_B_aligner_delta": agreement(rows, "C", "B"),
            "A_vs_B_bundled": agreement(rows, "A", "B"),
        },
        "arms": {
            "A_november_recorded": arm_summary(rows, "A"),
            "C_local_bwa": arm_summary(rows, "C"),
            "B_local_minimap2": arm_summary(rows, "B"),
        },
    }

    # Retention of November's calls, per clade
    clade_rows = []
    for clade in CLADES:
        sub = [r for r in rows if has_clade(r["A"], clade)]
        if not sub:
            continue
        ck = sum(1 for r in sub if has_clade(r["C"], clade))
        bk = sum(1 for r in sub if has_clade(r["B"], clade))
        clade_rows.append({
            "clade": clade,
            "november_n": len(sub),
            "bwa_keeps": ck, "bwa_keeps_pct": pct(ck, len(sub)),
            "minimap2_keeps": bk, "minimap2_keeps_pct": pct(bk, len(sub)),
            "bwa_lost": len(sub) - ck,
            "minimap2_lost": len(sub) - bk,
        })

    control = [r for r in rows
               if not has_clade(r["A"], "Arthropoda") and depth_of(split_path(path_of(r["A"]))) > 0]
    if control:
        cs = sum(1 for r in control if path_of(r["C"]) == path_of(r["A"]))
        bs = sum(1 for r in control if path_of(r["B"]) == path_of(r["A"]))
        clade_rows.append({
            "clade": "CONTROL_non_arthropod",
            "november_n": len(control),
            "bwa_keeps": cs, "bwa_keeps_pct": pct(cs, len(control)),
            "minimap2_keeps": bs, "minimap2_keeps_pct": pct(bs, len(control)),
            "bwa_lost": len(control) - cs,
            "minimap2_lost": len(control) - bs,
        })
    report["clade_retention"] = clade_rows

    # Depth movement C -> B, arthropod vs rest
    movement = {}
    for lab, sel in (
        ("Arthropoda", lambda r: has_clade(r["A"], "Arthropoda")),
        ("non_Arthropoda", lambda r: not has_clade(r["A"], "Arthropoda")),
    ):
        sub = [r for r in rows if sel(r)]
        if not sub:
            continue
        cd = [depth_of(split_path(path_of(r["C"]))) for r in sub]
        bd = [depth_of(split_path(path_of(r["B"]))) for r in sub]
        movement[lab] = {
            "n": len(sub),
            "deeper": sum(1 for a, b in zip(cd, bd) if b > a),
            "shallower": sum(1 for a, b in zip(cd, bd) if b < a),
            "same": sum(1 for a, b in zip(cd, bd) if b == a),
            "bwa_only_assigned": sum(1 for a, b in zip(cd, bd) if a > 0 and b == 0),
            "minimap2_only_assigned": sum(1 for a, b in zip(cd, bd) if a == 0 and b > 0),
        }
    report["depth_movement_C_to_B"] = movement

    (out / f"{args.label}_aligner_ab.json").write_text(json.dumps(report, indent=2))

    # Per-ASV table for the forensics step and for hand spot-checks
    with (out / f"{args.label}_joined.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["readname", "A_path", "C_path", "B_path",
                    "A_tree", "C_tree", "B_tree", "A_node", "C_node", "B_node",
                    "is_arthropod", "november_assigned",
                    "lost_by_bwa", "lost_by_minimap2"])
        for r in rows:
            a_on = depth_of(split_path(path_of(r["A"]))) > 0
            w.writerow([
                r["readname"], path_of(r["A"]), path_of(r["C"]), path_of(r["B"]),
                (r["A"] or {}).get("tree"), (r["C"] or {}).get("tree"), (r["B"] or {}).get("tree"),
                (r["A"] or {}).get("node"), (r["C"] or {}).get("node"), (r["B"] or {}).get("node"),
                has_clade(r["A"], "Arthropoda"), a_on,
                a_on and depth_of(split_path(path_of(r["C"]))) == 0,
                a_on and depth_of(split_path(path_of(r["B"]))) == 0,
            ])

    print(json.dumps(
        {k: report[k] for k in
         ("coverage", "fidelity_gate_A_vs_C", "comparisons", "arms",
          "clade_retention", "depth_movement_C_to_B")},
        indent=2,
    ))
    if not gate_pass:
        print("\nFIDELITY GATE FAILED -- do not quote C-vs-B as authoritative.",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
