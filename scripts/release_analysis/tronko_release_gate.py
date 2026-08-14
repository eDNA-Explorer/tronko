#!/usr/bin/env python3
"""Evaluate the release gates and emit an explicit PASS/FAIL record.

The point of this script is that a release decision stops depending on someone
reading a report carefully. Each gate is a named, thresholded check over
artifacts the other scripts already produce, and the exit code is the verdict.

Gates:
  G1 determinism      post-vs-post null run byte-identical at cores=1
  G2 provenance       binary SHA256s recorded, post-only guard string present
  G3 reproduces prod  local post-fix vs production assignments.parquet
  G4 BLAST concordance post-agrees-better exceeds pre-agrees-better, per marker
  G5 earned depth     extension_contradicted share in the >=97% identity band
  G6 no clade bias    no well-powered clade inverting the pooled direction
  G7 unassigned       no marker's unassigned rate rising more than 2 pp

G4-G7 thresholds are PROVISIONAL: they are calibrated from the first full run
and then ratified. Until that happens this reports them as `provisional` and
they should not be treated as a release veto on their own. Saying so is part of
the gate's job -- a threshold nobody agreed to is not a standard.

Usage:
    tronko_release_gate.py --config gate_inputs.json --out-dir OUT
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

DEFAULTS = {
    "g3_min_agreement_pct": 99.0,
    "g4_min_post_better_share": 0.50,
    "g5_max_extension_contradicted_share": 0.20,
    "g6_min_clade_n": 30,
    "g7_max_unassigned_increase_pp": 2.0,
}
PROVISIONAL = {"G4", "G5", "G6", "G7"}
GUARD_STRING = "leaf-portion mode) is not supported"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def has_guard(binary: Path) -> bool | None:
    try:
        out = subprocess.run(["strings", str(binary)], capture_output=True,
                             text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return GUARD_STRING in out


def load(path) -> dict | None:
    p = Path(path) if path else None
    if p and p.exists():
        return json.loads(p.read_text())
    return None


def gate(gid: str, title: str, status: str, detail: str, **extra) -> dict:
    return {
        "gate": gid, "title": title, "status": status, "detail": detail,
        "provisional": gid in PROVISIONAL, **extra,
    }


def g1_determinism(cfg) -> dict:
    a, b = cfg.get("null_run_a"), cfg.get("null_run_b")
    if not (a and Path(a).exists() and b and Path(b).exists()):
        return gate("G1", "determinism", "SKIP",
                    "no null-run pair supplied; every delta then needs a noise band")
    same = sha256(Path(a)) == sha256(Path(b))
    return gate("G1", "determinism", "PASS" if same else "FAIL",
                "post-vs-post byte-identical" if same else
                "post-vs-post differ at cores=1; deltas need a noise band")


def g2_provenance(cfg) -> dict:
    bins = cfg.get("binaries", {})
    if not bins:
        return gate("G2", "provenance", "SKIP", "no binaries supplied")
    rec, problems = {}, []
    for name, p in bins.items():
        path = Path(p)
        if not path.exists():
            problems.append(f"{name} missing")
            continue
        g = has_guard(path)
        rec[name] = {"sha256": sha256(path), "guard_string": g}
        if name == "post" and g is False:
            problems.append("post lacks the post-fix guard string")
        if name == "pre" and g is True:
            problems.append("pre carries the post-fix guard string")
    return gate("G2", "provenance", "FAIL" if problems else "PASS",
                "; ".join(problems) if problems else
                "binary hashes recorded and guard strings as expected",
                binaries=rec)


def g3_reproduces_prod(cfg, th) -> dict:
    rows = cfg.get("prod_agreement", [])
    if not rows:
        return gate("G3", "reproduces prod", "SKIP", "no agreement figures supplied")
    bad = [r for r in rows if r.get("agreement_pct", 0) < th["g3_min_agreement_pct"]]
    return gate("G3", "reproduces prod", "FAIL" if bad else "PASS",
                (f"{len(bad)} cell(s) below {th['g3_min_agreement_pct']}%: "
                 + ", ".join(f"{r['cell']} {r['agreement_pct']}%" for r in bad))
                if bad else
                f"all {len(rows)} cell(s) at or above {th['g3_min_agreement_pct']}%",
                cells=rows)


def g4_blast_concordance(cfg, th) -> dict:
    summaries = [(k, load(v)) for k, v in cfg.get("arbiter_summaries", {}).items()]
    summaries = [(k, s) for k, s in summaries if s]
    if not summaries:
        return gate("G4", "BLAST concordance", "SKIP", "no arbiter summaries supplied")
    rows, bad = [], []
    for label, s in summaries:
        h = s.get("headline", {})
        share = h.get("post_better_share")
        rows.append({"cell": label, "post_better_share": share,
                     "decided": (h.get("post_better") or 0) + (h.get("pre_better") or 0)})
        if share is not None and share < th["g4_min_post_better_share"]:
            bad.append(label)
    return gate("G4", "BLAST concordance", "FAIL" if bad else "PASS",
                (f"below {th['g4_min_post_better_share']:.2f} in: " + ", ".join(bad))
                if bad else "post agrees with BLAST better than pre in every cell",
                cells=rows)


def g5_earned_depth(cfg, th) -> dict:
    summaries = [(k, load(v)) for k, v in cfg.get("arbiter_summaries", {}).items()]
    summaries = [(k, s) for k, s in summaries if s]
    if not summaries:
        return gate("G5", "earned depth", "SKIP", "no arbiter summaries supplied")
    rows, bad = [], []
    for label, s in summaries:
        share = s.get("headline", {}).get("extension_contradicted_share_hi_identity")
        rows.append({"cell": label, "extension_contradicted_share_hi_identity": share})
        if share is not None and share > th["g5_max_extension_contradicted_share"]:
            bad.append(f"{label} {share:.3f}")
    return gate("G5", "earned depth", "FAIL" if bad else "PASS",
                ("deeper calls landing off the BLAST-supported lineage above "
                 f"{th['g5_max_extension_contradicted_share']:.2f}: " + ", ".join(bad))
                if bad else
                "deeper calls sit on the lineage the top hit supports",
                cells=rows)


def g6_clade_bias(cfg, th) -> dict:
    bias = load(cfg.get("bias_audit"))
    if not bias:
        return gate("G6", "no clade bias", "SKIP", "no bias audit supplied")
    inversions = [
        r for r in bias.get("direction_inversions", [])
        if (r.get("decided") or 0) >= th["g6_min_clade_n"]
    ]
    return gate("G6", "no clade bias", "FAIL" if inversions else "PASS",
                (f"{len(inversions)} well-powered clade(s) move against the pooled "
                 "direction: " + ", ".join(
                     f"{r['cell']}/{r['clade']} {r['post_better_share']:.2f}"
                     for r in inversions[:5]))
                if inversions else
                "no well-powered clade inverts the pooled direction",
                inversions=inversions)


def g7_unassigned(cfg, th) -> dict:
    anatomies = [(k, load(v)) for k, v in cfg.get("anatomy_summaries", {}).items()]
    anatomies = [(k, a) for k, a in anatomies if a]
    if not anatomies:
        return gate("G7", "unassigned", "SKIP", "no anatomy summaries supplied")
    rows, bad = [], []
    for label, a in anatomies:
        t = a.get("assigned_unassigned_transition", {})
        pre, post = t.get("unassigned_rate_pre_pct"), t.get("unassigned_rate_post_pct")
        if pre is None or post is None:
            continue
        delta = round(post - pre, 3)
        rows.append({"cell": label, "pre_pct": pre, "post_pct": post, "delta_pp": delta})
        if delta > th["g7_max_unassigned_increase_pp"]:
            bad.append(f"{label} +{delta}pp")
    return gate("G7", "unassigned", "FAIL" if bad else "PASS",
                (f"unassigned rate up more than {th['g7_max_unassigned_increase_pp']}pp: "
                 + ", ".join(bad)) if bad else
                "no marker's unassigned rate rose materially",
                cells=rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path,
                    help="JSON naming the artifacts to evaluate")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--fail-on-provisional", action="store_true",
                    help="treat a failing provisional gate as a release veto")
    args = ap.parse_args()

    cfg = json.loads(args.config.read_text())
    th = {**DEFAULTS, **cfg.get("thresholds", {})}
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    gates = [
        g1_determinism(cfg),
        g2_provenance(cfg),
        g3_reproduces_prod(cfg, th),
        g4_blast_concordance(cfg, th),
        g5_earned_depth(cfg, th),
        g6_clade_bias(cfg, th),
        g7_unassigned(cfg, th),
    ]

    hard_fail = [g for g in gates if g["status"] == "FAIL" and not g["provisional"]]
    soft_fail = [g for g in gates if g["status"] == "FAIL" and g["provisional"]]
    skipped = [g for g in gates if g["status"] == "SKIP"]

    verdict = (
        "BLOCKED" if hard_fail
        else "BLOCKED" if (soft_fail and args.fail_on_provisional)
        else "REVIEW" if (soft_fail or skipped)
        else "PASS"
    )

    result = {
        "verdict": verdict,
        "thresholds": th,
        "provisional_gates": sorted(PROVISIONAL),
        "note": (
            "G4-G7 thresholds are provisional until calibrated from a full run "
            "and ratified; a provisional failure is a prompt to look, not a veto."
        ),
        "gates": gates,
    }
    (out / "gate_result.json").write_text(json.dumps(result, indent=2))

    width = max(len(g["title"]) for g in gates)
    print(f"{'gate':<5}{'title':<{width + 2}}{'status':<8}detail")
    print("-" * (width + 60))
    for g in gates:
        mark = g["status"] + ("*" if g["provisional"] else "")
        print(f"{g['gate']:<5}{g['title']:<{width + 2}}{mark:<8}{g['detail']}")
    print("-" * (width + 60))
    print(f"VERDICT: {verdict}")
    if any(g["provisional"] for g in gates):
        print("* provisional threshold, not yet ratified")
    if skipped:
        print(f"{len(skipped)} gate(s) skipped for want of inputs — "
              f"a skipped gate is not a passed gate")

    return 1 if verdict == "BLOCKED" else 0


if __name__ == "__main__":
    sys.exit(main())
