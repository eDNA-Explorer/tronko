#!/usr/bin/env python3
"""Anatomy of what a tronko change did to a real project's assignments.

Answers, per marker/project, without reference to any synthetic data:
  * old vs new assignment -- how many changed, and in which direction
  * assigned vs unassigned -- the full transition matrix
  * placement movement    -- did reads move node, or tree entirely
  * taxon winners/losers  -- which taxa gained or lost reads

Placement comes from Tree_Number / Node_Number, which tronko already emits;
that same signal is what identified the July reference-indexing regression.

Usage:
    tronko_change_anatomy.py --pre PRE --post POST --out-dir OUT [--label NAME]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tronko_blast_arbiter import (  # noqa: E402
    N_RANKS,
    RANK_NAMES,
    classify_change,
    depth_of,
    is_assigned,
    read_assignments,
    split_path,
)


def pct(n: int, d: int) -> float | None:
    return round(100.0 * n / d, 3) if d else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pre", required=True, type=Path)
    ap.add_argument("--post", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="cell")
    ap.add_argument("--top-taxa", type=int, default=40)
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    pre_a = read_assignments(args.pre)
    post_a = read_assignments(args.post)

    common = sorted(set(pre_a) & set(post_a))
    only_pre = len(set(pre_a) - set(post_a))
    only_post = len(set(post_a) - set(pre_a))
    if only_pre or only_post:
        print(
            f"WARNING: readname sets differ (pre-only {only_pre:,}, post-only {only_post:,}); "
            "diffing the intersection",
            file=sys.stderr,
        )

    change_counts: Counter[str] = Counter()
    transition: Counter[tuple[bool, bool]] = Counter()
    placement: Counter[str] = Counter()
    placement_mm: dict[str, list[float]] = defaultdict(list)
    depth_pre: Counter[int] = Counter()
    depth_post: Counter[int] = Counter()
    rank_pre = [0] * N_RANKS
    rank_post = [0] * N_RANKS
    taxa_pre: Counter[str] = Counter()
    taxa_post: Counter[str] = Counter()
    first_div: Counter[str] = Counter()

    for name in common:
        p, q = pre_a[name], post_a[name]
        pr, qr = split_path(p["path"]), split_path(q["path"])
        kind = classify_change(pr, qr)
        change_counts[kind] += 1
        transition[(is_assigned(pr), is_assigned(qr))] += 1

        dp, dq = depth_of(pr), depth_of(qr)
        depth_pre[dp] += 1
        depth_post[dq] += 1
        for i in range(dp):
            rank_pre[i] += 1
        for i in range(dq):
            rank_post[i] += 1

        if is_assigned(pr):
            taxa_pre[";".join(x or "NA" for x in pr[:dp])] += 1
        if is_assigned(qr):
            taxa_post[";".join(x or "NA" for x in qr[:dq])] += 1

        if kind not in ("identical", "both_unassigned"):
            # where did the two paths first part company
            for i in range(N_RANKS):
                if pr[i] != qr[i]:
                    first_div[RANK_NAMES[i]] += 1
                    break

            pt, qt = p.get("tree"), q.get("tree")
            pn, qn = p.get("node"), q.get("node")
            if pt is None or qt is None:
                bucket = "unknown"
            elif pt == qt and pn == qn:
                bucket = "same_tree_same_node"
            elif pt == qt:
                bucket = "same_tree_different_node"
            else:
                bucket = "different_tree"
            placement[bucket] += 1
            pm, qm = p.get("fwd_mm"), q.get("fwd_mm")
            if pm is not None and qm is not None:
                placement_mm[bucket].append(float(qm) - float(pm))

    n = len(common)
    changed = n - change_counts["identical"] - change_counts["both_unassigned"]

    def mm_summary(b: str) -> dict:
        v = placement_mm.get(b) or []
        if not v:
            return {"n": 0}
        v = sorted(v)
        return {
            "n": len(v),
            "mean_fwd_mismatch_delta": round(sum(v) / len(v), 4),
            "median_fwd_mismatch_delta": round(v[len(v) // 2], 4),
            "improved_share": round(sum(1 for x in v if x < 0) / len(v), 4),
        }

    summary = {
        "label": args.label,
        "reads_compared": n,
        "readname_set_mismatch": {"pre_only": only_pre, "post_only": only_post},
        "changed": changed,
        "change_rate_pct": pct(changed, n),
        "change_types": dict(change_counts.most_common()),
        "assigned_unassigned_transition": {
            "assigned_to_assigned": transition[(True, True)],
            "assigned_to_unassigned": transition[(True, False)],
            "unassigned_to_assigned": transition[(False, True)],
            "unassigned_to_unassigned": transition[(False, False)],
            "unassigned_rate_pre_pct": pct(
                transition[(False, True)] + transition[(False, False)], n
            ),
            "unassigned_rate_post_pct": pct(
                transition[(True, False)] + transition[(False, False)], n
            ),
        },
        "placement_movement": {
            b: {"reads": c, "share_of_changed": pct(c, changed), **mm_summary(b)}
            for b, c in placement.most_common()
        },
        "first_divergent_rank": dict(first_div.most_common()),
        "depth": {
            "mean_pre": round(
                sum(k * v for k, v in depth_pre.items()) / max(1, n), 4
            ),
            "mean_post": round(
                sum(k * v for k, v in depth_post.items()) / max(1, n), 4
            ),
            "distribution_pre": dict(sorted(depth_pre.items())),
            "distribution_post": dict(sorted(depth_post.items())),
        },
        "rank_resolution_pct": {
            RANK_NAMES[i]: {
                "pre": pct(rank_pre[i], n),
                "post": pct(rank_post[i], n),
                "delta": round((pct(rank_post[i], n) or 0) - (pct(rank_pre[i], n) or 0), 3),
            }
            for i in range(N_RANKS)
        },
        "taxa": {
            "distinct_pre": len(taxa_pre),
            "distinct_post": len(taxa_post),
            "gained": len(set(taxa_post) - set(taxa_pre)),
            "lost": len(set(taxa_pre) - set(taxa_post)),
        },
    }

    (out / f"{args.label}_anatomy.json").write_text(json.dumps(summary, indent=2))

    # taxon winners and losers, by read count
    deltas = []
    for t in set(taxa_pre) | set(taxa_post):
        a, b = taxa_pre.get(t, 0), taxa_post.get(t, 0)
        deltas.append({"taxon": t, "pre_reads": a, "post_reads": b, "delta": b - a})
    deltas.sort(key=lambda d: d["delta"])
    losers = deltas[: args.top_taxa]
    winners = list(reversed(deltas[-args.top_taxa :]))
    with (out / f"{args.label}_taxon_deltas.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["taxon", "pre_reads", "post_reads", "delta"])
        w.writeheader()
        w.writerows(winners + losers)

    print(json.dumps({k: summary[k] for k in
                      ("label", "reads_compared", "changed", "change_rate_pct",
                       "change_types", "assigned_unassigned_transition",
                       "placement_movement")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
