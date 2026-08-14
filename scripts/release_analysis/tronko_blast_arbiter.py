#!/usr/bin/env python3
"""Adjudicate tronko assignment changes against an independent BLAST call.

The question this answers, on real project data where no ground truth exists:
when an ASV's classification changed between two tronko builds, does the new
call agree better with the ASV's top BLAST hit than the old one did?  And when
the call goes *deeper on the same lineage*, is that extra depth on the lineage
the top hit actually supports, or is it unearned?

The BLAST subject is the tronko database's own reference FASTA, so both the
tronko paths and the BLAST lineages come from one taxonomy and no remapping is
needed.  That also scopes the verdict: this adjudicates *placement*, not
reference coverage.  A taxon missing from the reference is invisible to both
tools alike.

Usage:
    tronko_blast_arbiter.py --pre PRE --post POST \
        --forward F.fasta --reverse R.fasta \
        --marker-fasta ref.fasta --taxonomy ref_taxonomy.txt \
        --out-dir OUT [--label NAME]

PRE/POST accept either tronko TSV output or the parquet prod writes.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

N_RANKS = 7
RANK_NAMES = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
UNASSIGNED_TOKENS = {"unassigned", "", "NA;NA;NA;NA;NA;NA;NA"}

# Percent-identity bands, matching realworld_similarity_bands/compute_bands.py so
# the two analyses stay comparable.
BANDS: tuple[tuple[str, float, float], ...] = (
    (">=99", 99.0, 100.01),
    ("97-99", 97.0, 99.0),
    ("90-97", 90.0, 97.0),
    ("80-90", 80.0, 90.0),
    ("50-80", 50.0, 80.0),
    ("<50/none", -1.0, 50.0),
)


def band_of(pid: float | None) -> str:
    if pid is None:
        return BANDS[-1][0]
    for name, lo, hi in BANDS:
        if lo <= pid < hi:
            return name
    return BANDS[0][0] if pid >= 99.0 else BANDS[-1][0]


# --------------------------------------------------------------------------
# taxonomy / path handling
# --------------------------------------------------------------------------


def split_path(path: str | None) -> list[str | None]:
    """Return a 7-slot lineage, with None for absent or NA ranks.

    Tronko emits variable-length paths and pads missing ranks with the literal
    "NA".  Both mean "not resolved to this rank" and must not be treated as a
    name that can agree or disagree with anything.
    """
    if path is None:
        return [None] * N_RANKS
    text = path.strip()
    if text.lower() in {t.lower() for t in UNASSIGNED_TOKENS}:
        return [None] * N_RANKS
    parts = [p.strip() for p in text.split(";")]
    out: list[str | None] = []
    for i in range(N_RANKS):
        v = parts[i] if i < len(parts) else ""
        out.append(None if v in ("", "NA", "na", "N/A") else v)
    return out


def depth_of(ranks: list[str | None]) -> int:
    """Number of leading resolved ranks (stops at the first gap)."""
    n = 0
    for r in ranks:
        if r is None:
            break
        n += 1
    return n


def is_assigned(ranks: list[str | None]) -> bool:
    return depth_of(ranks) > 0


def agree_depth(a: list[str | None], b: list[str | None]) -> int:
    """Leading ranks where both are resolved and equal."""
    n = 0
    for x, y in zip(a, b):
        if x is None or y is None or x != y:
            break
        n += 1
    return n


def first_conflict(a: list[str | None], b: list[str | None]) -> int | None:
    """Index of the first rank where both resolve but disagree, else None."""
    for i, (x, y) in enumerate(zip(a, b)):
        if x is not None and y is not None and x != y:
            return i
    return None


def load_taxonomy(path: Path) -> dict[str, list[str | None]]:
    tax: dict[str, list[str | None]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            acc, _, lineage = line.partition("\t")
            if not lineage:
                continue
            tax[acc.strip()] = split_path(lineage)
    return tax


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_fasta(path: Path) -> dict[str, str]:
    """Read a (possibly zstd-compressed, possibly wrapped) FASTA."""
    if path.suffix == ".zst":
        data = subprocess.run(
            ["zstd", "-dc", str(path)],
            capture_output=True,
            check=True,
        ).stdout.decode()
        lines = data.splitlines()
    else:
        lines = path.read_text().splitlines()

    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for ln in lines:
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_assignments(path: Path) -> dict[str, dict]:
    """Read tronko output as {readname: {path, tree, node, fwd_mm, rev_mm}}.

    Accepts the TSV tronko-assign writes and the parquet production stores.
    """
    out: dict[str, dict] = {}
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq

        tbl = pq.read_table(path)
        cols = {c.lower(): c for c in tbl.column_names}

        def col(*names):
            for n in names:
                if n in cols:
                    return tbl.column(cols[n]).to_pylist()
            return None

        names = col("readname")
        paths = col("taxonomic_path")
        trees = col("tree_number") or [None] * len(names)
        nodes = col("node_number") or [None] * len(names)
        fmm = col("forward_mismatch") or [None] * len(names)
        rmm = col("reverse_mismatch") or [None] * len(names)

        def dec(v):
            return v.decode() if isinstance(v, (bytes, bytearray)) else v

        for i, nm in enumerate(names):
            out[dec(nm)] = {
                "path": dec(paths[i]),
                "tree": trees[i],
                "node": nodes[i],
                "fwd_mm": fmm[i],
                "rev_mm": rmm[i],
            }
        return out

    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}

        def gi(row, key, cast=None):
            i = idx.get(key)
            if i is None or i >= len(row) or row[i] == "":
                return None
            return cast(row[i]) if cast else row[i]

        for line in fh:
            row = line.rstrip("\n").split("\t")
            nm = row[idx["Readname"]]
            out[nm] = {
                "path": gi(row, "Taxonomic_Path"),
                "tree": gi(row, "Tree_Number", lambda x: int(float(x))),
                "node": gi(row, "Node_Number", lambda x: int(float(x))),
                "fwd_mm": gi(row, "Forward_Mismatch", float),
                "rev_mm": gi(row, "Reverse_Mismatch", float),
            }
    return out


# --------------------------------------------------------------------------
# dereplication
# --------------------------------------------------------------------------


@dataclass
class Asv:
    asv_id: str
    fseq: str
    rseq: str
    reads: list[str] = field(default_factory=list)

    @property
    def abundance(self) -> int:
        return len(self.reads)


def dereplicate(fwd: dict[str, str], rev: dict[str, str]) -> tuple[list[Asv], dict[str, str]]:
    """Collapse read pairs to unique (F, R) sequence pairs.

    Returns the ASV list and a readname -> asv_id index.  Reads whose mate is
    missing are kept with an empty reverse sequence rather than dropped, so the
    abundance totals still reconcile against the input.
    """
    by_seq: dict[tuple[str, str], Asv] = {}
    read_to_asv: dict[str, str] = {}
    for name, fseq in fwd.items():
        rseq = rev.get(name, "")
        key = (fseq, rseq)
        asv = by_seq.get(key)
        if asv is None:
            asv = Asv(asv_id=f"asv_{len(by_seq):07d}", fseq=fseq, rseq=rseq)
            by_seq[key] = asv
        asv.reads.append(name)
        read_to_asv[name] = asv.asv_id
    return list(by_seq.values()), read_to_asv


# --------------------------------------------------------------------------
# BLAST
# --------------------------------------------------------------------------


def ensure_blast_db(marker_fasta: Path, work: Path) -> Path:
    db = work / "blastdb" / marker_fasta.stem
    db.parent.mkdir(parents=True, exist_ok=True)
    if not (db.with_suffix(db.suffix + ".nin").exists() or Path(str(db) + ".nin").exists()):
        subprocess.run(
            ["makeblastdb", "-in", str(marker_fasta), "-dbtype", "nucl", "-out", str(db)],
            check=True,
            capture_output=True,
        )
    return db


def run_blast(
    asvs: list[Asv],
    db: Path,
    work: Path,
    threads: int,
    max_target_seqs: int = 5,
) -> dict[str, dict]:
    """Blast both mates of each ASV; return the best-supported subject per ASV.

    Both mates are queried in one call and their bitscores summed per subject,
    mirroring how tronko scans both mates before deciding.  The reported
    percent identity is the best across the mates that support the winner.
    """
    qpath = work / "blast_query.fasta"
    with qpath.open("w") as fh:
        for a in asvs:
            if a.fseq:
                fh.write(f">{a.asv_id}/1\n{a.fseq}\n")
            if a.rseq:
                fh.write(f">{a.asv_id}/2\n{a.rseq}\n")

    outpath = work / "blast_hits.tsv"
    subprocess.run(
        [
            "blastn",
            "-task", "megablast",
            "-query", str(qpath),
            "-db", str(db),
            "-max_target_seqs", str(max_target_seqs),
            "-evalue", "1e-5",
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid pident length bitscore qcovhsp",
            "-out", str(outpath),
        ],
        check=True,
        capture_output=True,
    )

    # sum bitscore per (asv, subject) across mates, keep best identity
    agg: dict[str, dict[str, dict]] = defaultdict(dict)
    with outpath.open() as fh:
        for line in fh:
            q, s, pid, length, bits, qcov = line.rstrip("\n").split("\t")
            asv_id, _, mate = q.partition("/")
            rec = agg[asv_id].get(s)
            pid_f, bits_f, qcov_f = float(pid), float(bits), float(qcov)
            if rec is None:
                agg[asv_id][s] = {
                    "bits": bits_f,
                    "pident": pid_f,
                    "qcov": qcov_f,
                    "mates": {mate},
                }
            else:
                rec["bits"] += bits_f
                rec["pident"] = max(rec["pident"], pid_f)
                rec["qcov"] = max(rec["qcov"], qcov_f)
                rec["mates"].add(mate)

    best: dict[str, dict] = {}
    for asv_id, subjects in agg.items():
        top_bits = max(r["bits"] for r in subjects.values())
        tied = [s for s, r in subjects.items() if r["bits"] == top_bits]
        s = tied[0]
        r = subjects[s]
        best[asv_id] = {
            "subject": s,
            "pident": r["pident"],
            "qcov": r["qcov"],
            "bits": r["bits"],
            "both_mates": len(r["mates"]) > 1,
            "tied": len(tied) > 1,
            "n_tied": len(tied),
            "tied_subjects": tied if len(tied) > 1 else None,
        }
    return best


# --------------------------------------------------------------------------
# adjudication
# --------------------------------------------------------------------------


def classify_change(pre: list[str | None], post: list[str | None]) -> str:
    pre_on, post_on = is_assigned(pre), is_assigned(post)
    if not pre_on and not post_on:
        return "both_unassigned"
    if pre_on and not post_on:
        return "assigned_to_unassigned"
    if not pre_on and post_on:
        return "unassigned_to_assigned"
    dpre, dpost = depth_of(pre), depth_of(post)
    common = agree_depth(pre, post)
    if common == dpre and common == dpost:
        return "identical"
    if common == dpre and dpost > dpre:
        return "extension"
    if common == dpost and dpre > dpost:
        return "retraction"
    return "divergence"


def adjudicate(
    kind: str,
    pre: list[str | None],
    post: list[str | None],
    blast: list[str | None] | None,
) -> tuple[str, str]:
    """Return (verdict, detail).

    verdict is one of: post_better, pre_better, tie, unverifiable.
    """
    if blast is None or not is_assigned(blast):
        return "unverifiable", "no_blast_lineage"

    a_pre, a_post = agree_depth(pre, blast), agree_depth(post, blast)
    b_depth = depth_of(blast)

    if kind == "extension":
        added_from, added_to = depth_of(pre), depth_of(post)
        if b_depth <= added_from:
            # BLAST cannot speak to the ranks that were added
            return "unverifiable", "blast_too_shallow_for_added_ranks"
        # do the added ranks (as far as BLAST reaches) match?
        checkable = min(added_to, b_depth)
        for i in range(added_from, checkable):
            if post[i] != blast[i]:
                return "pre_better", "extension_contradicted"
        if checkable < added_to:
            return "post_better", "extension_partly_confirmed"
        return "post_better", "extension_confirmed"

    if kind == "retraction":
        # post dropped ranks pre had. Correct if those ranks disagreed with BLAST.
        dropped_from, dropped_to = depth_of(post), depth_of(pre)
        checkable = min(dropped_to, b_depth)
        if checkable <= dropped_from:
            return "unverifiable", "blast_too_shallow_for_dropped_ranks"
        for i in range(dropped_from, checkable):
            if pre[i] != blast[i]:
                return "post_better", "retraction_removed_wrong_ranks"
        return "pre_better", "retraction_lost_correct_ranks"

    if kind == "divergence":
        if a_post > a_pre:
            return "post_better", "divergence_post_closer"
        if a_pre > a_post:
            return "pre_better", "divergence_pre_closer"
        return "tie", "divergence_equal_agreement"

    if kind == "assigned_to_unassigned":
        if a_pre > 0:
            conflict = first_conflict(pre, blast)
            if conflict is None:
                return "pre_better", "lost_a_blast_supported_call"
            return "tie", "lost_a_partly_wrong_call"
        return "post_better", "refused_a_blast_contradicted_call"

    if kind == "unassigned_to_assigned":
        if a_post > 0:
            conflict = first_conflict(post, blast)
            if conflict is None:
                return "post_better", "gained_a_blast_supported_call"
            return "tie", "gained_a_partly_wrong_call"
        return "pre_better", "gained_a_blast_contradicted_call"

    return "tie", kind


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pre", required=True, type=Path)
    ap.add_argument("--post", required=True, type=Path)
    ap.add_argument("--forward", required=True, type=Path)
    ap.add_argument("--reverse", required=True, type=Path)
    ap.add_argument("--marker-fasta", required=True, type=Path)
    ap.add_argument("--taxonomy", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="cell")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument(
        "--control-sample",
        type=int,
        default=2000,
        help="unchanged ASVs to also BLAST, as a concordance baseline",
    )
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    print(f"[{args.label}] loading reads", flush=True)
    fwd, rev = read_fasta(args.forward), read_fasta(args.reverse)
    asvs, read_to_asv = dereplicate(fwd, rev)
    print(f"[{args.label}] {len(fwd):,} read pairs -> {len(asvs):,} ASVs", flush=True)

    print(f"[{args.label}] loading assignments", flush=True)
    pre_a, post_a = read_assignments(args.pre), read_assignments(args.post)

    # Assignment must be a pure function of sequence; verify rather than assume.
    inconsistent = {"pre": 0, "post": 0}
    asv_paths: dict[str, dict[str, list[str | None]]] = {}
    for asv in asvs:
        for tag, table in (("pre", pre_a), ("post", post_a)):
            seen = {table[r]["path"] for r in asv.reads if r in table}
            if len(seen) > 1:
                inconsistent[tag] += 1
            val = next(iter(seen)) if seen else None
            asv_paths.setdefault(asv.asv_id, {})[tag] = split_path(val)

    # placement (tree/node) taken from the first read that has it
    placement: dict[str, dict] = {}
    for asv in asvs:
        rec: dict = {}
        for tag, table in (("pre", pre_a), ("post", post_a)):
            for r in asv.reads:
                if r in table:
                    rec[tag] = table[r]
                    break
        placement[asv.asv_id] = rec

    changed = [
        a for a in asvs
        if classify_change(asv_paths[a.asv_id]["pre"], asv_paths[a.asv_id]["post"])
        not in ("identical", "both_unassigned")
    ]
    print(f"[{args.label}] {len(changed):,} changed ASVs", flush=True)

    import random

    rng = random.Random(args.seed)
    unchanged = [a for a in asvs if a not in set(changed)]
    control = rng.sample(unchanged, min(args.control_sample, len(unchanged)))

    to_blast = changed + control
    print(f"[{args.label}] blasting {len(to_blast):,} ASVs "
          f"({len(changed):,} changed + {len(control):,} control)", flush=True)

    tax = load_taxonomy(args.taxonomy)
    db = ensure_blast_db(args.marker_fasta, out)
    hits = run_blast(to_blast, db, out, args.threads)

    unresolved = 0
    rows = []
    for asv in to_blast:
        pre_p = asv_paths[asv.asv_id]["pre"]
        post_p = asv_paths[asv.asv_id]["post"]
        kind = classify_change(pre_p, post_p)
        hit = hits.get(asv.asv_id)
        blast_lineage = None
        if hit:
            blast_lineage = tax.get(hit["subject"])
            if blast_lineage is None:
                unresolved += 1
        verdict, detail = adjudicate(kind, pre_p, post_p, blast_lineage)
        p = placement[asv.asv_id]
        rows.append(
            {
                "asv_id": asv.asv_id,
                "abundance": asv.abundance,
                "is_control": kind in ("identical", "both_unassigned"),
                "change": kind,
                "verdict": verdict,
                "detail": detail,
                "pre_path": ";".join(x or "NA" for x in pre_p),
                "post_path": ";".join(x or "NA" for x in post_p),
                "pre_depth": depth_of(pre_p),
                "post_depth": depth_of(post_p),
                "blast_subject": hit["subject"] if hit else None,
                "blast_pident": hit["pident"] if hit else None,
                "blast_qcov": hit["qcov"] if hit else None,
                "blast_both_mates": hit["both_mates"] if hit else None,
                "blast_tied": hit["tied"] if hit else None,
                "blast_lineage": (
                    ";".join(x or "NA" for x in blast_lineage) if blast_lineage else None
                ),
                "blast_depth": depth_of(blast_lineage) if blast_lineage else 0,
                "band": band_of(hit["pident"] if hit else None),
                "pre_agree": agree_depth(pre_p, blast_lineage) if blast_lineage else None,
                "post_agree": agree_depth(post_p, blast_lineage) if blast_lineage else None,
                "pre_tree": p.get("pre", {}).get("tree"),
                "post_tree": p.get("post", {}).get("tree"),
                "pre_node": p.get("pre", {}).get("node"),
                "post_node": p.get("post", {}).get("node"),
                "pre_fwd_mm": p.get("pre", {}).get("fwd_mm"),
                "post_fwd_mm": p.get("post", {}).get("fwd_mm"),
            }
        )

    import csv

    detail_path = out / f"{args.label}_asv_verdicts.csv"
    with detail_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = summarise(rows, args.label)
    summary["inputs"] = {
        "read_pairs": len(fwd),
        "asvs": len(asvs),
        "changed_asvs": len(changed),
        "control_asvs": len(control),
        "asv_path_inconsistent": inconsistent,
        "blast_subjects_unresolved_in_taxonomy": unresolved,
    }
    (out / f"{args.label}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["headline"], indent=2))
    print(f"[{args.label}] wrote {detail_path}", flush=True)
    return 0


def summarise(rows: list[dict], label: str) -> dict:
    changed = [r for r in rows if not r["is_control"]]
    control = [r for r in rows if r["is_control"]]

    def tally(subset, key):
        out: dict[str, int] = defaultdict(int)
        for r in subset:
            out[r[key]] += 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    def weighted(subset, key):
        out: dict[str, int] = defaultdict(int)
        for r in subset:
            out[r[key]] += r["abundance"]
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    pb = sum(1 for r in changed if r["verdict"] == "post_better")
    prb = sum(1 for r in changed if r["verdict"] == "pre_better")
    pb_w = sum(r["abundance"] for r in changed if r["verdict"] == "post_better")
    prb_w = sum(r["abundance"] for r in changed if r["verdict"] == "pre_better")

    ext = [r for r in changed if r["change"] == "extension"]
    ext_conf = sum(1 for r in ext if r["detail"].startswith("extension_confirmed")
                   or r["detail"] == "extension_partly_confirmed")
    ext_contra = sum(1 for r in ext if r["detail"] == "extension_contradicted")
    ext_hi = [r for r in ext if r["band"] in (">=99", "97-99")]
    ext_hi_contra = sum(1 for r in ext_hi if r["detail"] == "extension_contradicted")

    by_band: dict[str, dict] = {}
    for name, _, _ in BANDS:
        sub = [r for r in changed if r["band"] == name]
        if not sub:
            continue
        by_band[name] = {
            "changed_asvs": len(sub),
            "post_better": sum(1 for r in sub if r["verdict"] == "post_better"),
            "pre_better": sum(1 for r in sub if r["verdict"] == "pre_better"),
            "tie": sum(1 for r in sub if r["verdict"] == "tie"),
            "unverifiable": sum(1 for r in sub if r["verdict"] == "unverifiable"),
        }

    return {
        "label": label,
        "headline": {
            "changed_asvs": len(changed),
            "post_better": pb,
            "pre_better": prb,
            "post_better_share": round(pb / (pb + prb), 4) if (pb + prb) else None,
            "post_better_share_read_weighted": (
                round(pb_w / (pb_w + prb_w), 4) if (pb_w + prb_w) else None
            ),
            "extensions": len(ext),
            "extension_confirmed": ext_conf,
            "extension_contradicted": ext_contra,
            "extension_contradicted_share_hi_identity": (
                round(ext_hi_contra / len(ext_hi), 4) if ext_hi else None
            ),
        },
        "change_types": tally(changed, "change"),
        "change_types_read_weighted": weighted(changed, "change"),
        "verdicts": tally(changed, "verdict"),
        "verdict_details": tally(changed, "detail"),
        "by_identity_band": by_band,
        "control_agreement": {
            "n": len(control),
            "mean_agree_depth": (
                round(
                    sum(r["post_agree"] for r in control if r["post_agree"] is not None)
                    / max(1, sum(1 for r in control if r["post_agree"] is not None)),
                    3,
                )
            ),
        },
    }


if __name__ == "__main__":
    sys.exit(main())
