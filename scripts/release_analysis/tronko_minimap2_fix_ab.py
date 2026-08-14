"""Marker-driven tronko minimap2 fix A/B: run the matrix and score it.

  ab.py run   --marker 18s_euk_lca --run-root <benchmark_run_dir>
  ab.py score --marker 18s_euk_lca --run-root <benchmark_run_dir>

Matrix is {pre,post} x {fwd,rc} x communities, unmerged paired F+R, at
production parameters. The `rc` arm reverse-complements both mates, which
reproduces the production vert12S condition (query antiparallel to the
reference) and is the only way the strand half of the fix can fire.

Everything marker-specific is derived from MARKER_PRESETS.
"""

import argparse
import hashlib
import json
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path("/Users/ryanmartin/edna-explorer-data-pipelines")
sys.path.insert(0, str(REPO / "projects/assignment_benchmarks/src"))
sys.path.insert(0, str(REPO / "libraries/taxonomy-classifier-tuning/src"))

import polars as pl  # noqa: E402

from assignment_benchmarks.infrastructure.tronko_taxonomy_remap import (  # noqa: E402
    build_runtime_tronko_remapper,
    translate_tronko_output_path,
)
from taxonomy_classifier_tuning.domain.metrics import (  # noqa: E402
    compute_per_rank_metrics,
    compute_taxonomic_distance,
    compute_truth_depth_asv_metrics,
)
from taxonomy_classifier_tuning.domain.models import TAXONOMIC_RANKS  # noqa: E402
from taxonomy_classifier_tuning.pipeline.markers import MARKER_PRESETS  # noqa: E402

BASE = Path("/Users/ryanmartin/tronko-ab")
BIN = BASE / "bin"
CACHE = BASE / "cache"
TAXDUMP = Path("/Users/ryanmartin/rcrux-py/data/taxdump")
UNASSIGNED = "Unassigned"
NA_PATH = ";".join(["NA"] * len(TAXONOMIC_RANKS))

# Production parameters (also the binary's compiled-in defaults at both
# commits, tronko-assign.c:939-956). cores=1 because tronko is ~3.2%
# non-deterministic above -C 1, which would sit under the signal.
PROD_ARGS = [
    "-6", "--number-of-cores", "1", "--Cinterval", "0.02", "--aligner", "minimap2",
    "-u", "0.0001", "--max-leaf-matches", "10",
    "--best-leaf-threshold", "-0.1", "--best-leaf-max-votes", "10",
]

COMP = str.maketrans(
    "ACGTUacgtuRYSWKMBDHVNryswkmbdhvn", "TGCAAtgcaaYRSWMKVHDBNyrswmkvhdbn"
)

# tronko's option struct uses fixed 200-byte buffers for -f/-a/-o
# (global.h:227). Longer paths silently overrun the NUL terminator into the
# adjacent field; that is what made the earlier sweep fail with exit 255.
MAX_TRONKO_PATH = 199


def preset(marker: str):
    if marker not in MARKER_PRESETS:
        raise SystemExit(f"unknown marker {marker!r}; choose from {sorted(MARKER_PRESETS)}")
    return MARKER_PRESETS[marker]


def workdir(marker: str) -> Path:
    """Short-path symlink farm, so no tronko argument exceeds 200 bytes."""
    return BASE / "w" / marker


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, header, chunks = [], None, []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: Path, records, rc: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            out = seq.translate(COMP)[::-1] if rc else seq
            fh.write(f">{header}\n{out}\n")


def discover(marker: str, run_root: Path) -> tuple[Path, list[str]]:
    """Link DB + QC dirs into the short-path farm; return (farm, communities)."""
    farm = workdir(marker)
    farm.mkdir(parents=True, exist_ok=True)
    shared = run_root / "data" / "shared"
    comms = sorted(p.name for p in shared.glob("unified_*") if p.is_dir())
    kept = []
    for c in comms:
        idx = c.split("_")[-1]
        db = run_root / "data" / "unified_holdout" / c / "tool_dbs" / "tronko_db"
        qc = shared / c / "qc_edna_damage"
        if not db.exists() or not qc.exists():
            print(f"  skip {c}: db={db.exists()} qc={qc.exists()}")
            continue
        for link, target in ((f"db{idx}", db), (f"q{idx}", qc)):
            lp = farm / link
            if lp.is_symlink() or lp.exists():
                lp.unlink()
            lp.symlink_to(target)
        kept.append(idx)
    return farm, kept


def _mate_paths(farm: Path, marker_name: str, idx: str) -> tuple[Path, Path]:
    d = farm / f"q{idx}" / "fasta_gen" / "tronko"
    return d / f"{marker_name}_paired_F.fasta", d / f"{marker_name}_paired_R.fasta"


def cmd_run(marker: str, run_root: Path, out: Path) -> None:
    p = preset(marker)
    farm, comms = discover(marker, run_root)
    if not comms:
        raise SystemExit("no communities with both a tronko_db and QC output")
    out.mkdir(parents=True, exist_ok=True)

    print(f"=== {marker}: preparing query FASTAs (unwrapped; rc = revcomp both mates) ===")
    for idx in comms:
        f_src, r_src = _mate_paths(farm, p.marker.name, idx)
        for arm, do_rc in (("fwd", False), ("rc", True)):
            for mate, src in (("F", f_src), ("R", r_src)):
                write_fasta(out / "queries" / arm / f"c{idx}_{mate}.fasta",
                            read_fasta(src), rc=do_rc)
        print(f"  community {idx}: {len(read_fasta(f_src)):,} pairs")

    prov = {
        "marker": marker,
        "commits": {"pre": "f3bdfac", "post": "530caf9"},
        "binaries": {b: {"sha256": hashlib.sha256((BIN / f"tronko-assign.{b}").read_bytes()).hexdigest()}
                     for b in ("pre", "post")},
        "prod_args": PROD_ARGS,
        "communities": comms,
        "cells": [],
    }

    print(f"\n=== running matrix: {{pre,post}} x {{fwd,rc}} x {len(comms)} communities ===")
    for binary in ("pre", "post"):
        for arm in ("fwd", "rc"):
            for idx in comms:
                exp = out / f"{binary}__{arm}__c{idx}"
                exp.mkdir(parents=True, exist_ok=True)
                tsv = exp / "tronko_output.tsv"
                db = farm / f"db{idx}"
                cmd = [str(BIN / f"tronko-assign.{binary}"), "-r",
                       "-f", str(db / "reference_tree.trkb"),
                       "-a", str(db / "marker.fasta"), "-p", "-z", "-w",
                       "-1", str(out / "queries" / arm / f"c{idx}_F.fasta"),
                       "-2", str(out / "queries" / arm / f"c{idx}_R.fasta"),
                       *PROD_ARGS, "-o", str(tsv)]
                for i, tok in enumerate(cmd):
                    if cmd[i - 1] in ("-f", "-a", "-o") and len(tok) > MAX_TRONKO_PATH:
                        raise SystemExit(f"path too long for tronko ({len(tok)}): {tok}")

                start = time.time()
                proc = subprocess.run(cmd, capture_output=True, text=True)
                elapsed = time.time() - start
                if proc.returncode != 0:
                    # tronko printf()s fatal errors to stdout, not stderr.
                    raise SystemExit(f"[{binary}/{arm}/c{idx}] exit {proc.returncode}\n"
                                     f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
                rows = sum(1 for _ in tsv.open()) - 1
                digest = hashlib.sha256(tsv.read_bytes()).hexdigest()
                print(f"  {binary:4} {arm:3} c{idx}  {elapsed:6.1f}s  {rows:,} rows  {digest[:12]}",
                      flush=True)
                prov["cells"].append({"binary": binary, "arm": arm, "community": idx,
                                      "rows": rows, "seconds": round(elapsed, 1),
                                      "sha256": digest})
    (out / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"\nwrote {out / 'provenance.json'}")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def split_path(path):
    if not path or path.strip().lower().startswith("unassigned"):
        return [None] * len(TAXONOMIC_RANKS)
    parts = [p.strip() for p in path.split(";")]
    parts = [p if p and not p.upper().startswith("NA") else None for p in parts]
    parts += [None] * (len(TAXONOMIC_RANKS) - len(parts))
    return parts[: len(TAXONOMIC_RANKS)]


def load_cell(out: Path, binary: str, arm: str, idx: str, remapper) -> pl.DataFrame:
    df = pl.read_csv(out / f"{binary}__{arm}__c{idx}" / "tronko_output.tsv",
                     separator="\t", infer_schema_length=0)
    memo, translated = {}, []
    for raw in df["Taxonomic_Path"].to_list():
        key = raw or ""
        if key not in memo:
            memo[key] = (UNASSIGNED if not key or key.lower().startswith("unassigned")
                         else translate_tronko_output_path(
                             key, remapper, ambiguity_policy="rank-na").translated_path)
        translated.append(memo[key])
    return pl.DataFrame({
        "read_id": df["Readname"].to_list(),
        "pred_path": translated,
        "fwd_mm": [_num(x) for x in df["Forward_Mismatch"].to_list()],
        "rev_mm": [_num(x) for x in df["Reverse_Mismatch"].to_list()],
    })


def paired_frame(cell: pl.DataFrame, truth: dict) -> pl.DataFrame:
    rows = cell.filter(pl.col("read_id").is_in(list(truth.keys())))
    ids = rows["read_id"].to_list()
    t = [split_path(truth[r]) for r in ids]
    p = [split_path(x) for x in rows["pred_path"].to_list()]
    data = {"read_id": ids}
    for i, rank in enumerate(TAXONOMIC_RANKS):
        data[f"true_{rank}"] = [v[i] for v in t]
        data[f"pred_{rank}"] = [v[i] for v in p]
    return pl.DataFrame(data, schema_overrides={
        **{f"true_{r}": pl.Utf8 for r in TAXONOMIC_RANKS},
        **{f"pred_{r}": pl.Utf8 for r in TAXONOMIC_RANKS}})


def scoreable(p) -> bool:
    """A path counts only if it survived the taxonomy remap. Shallow legacy
    paths are genuinely ambiguous and come back all-NA; treating those as real
    answers would score 'pre unresolvable -> post resolvable' as improvement."""
    return bool(p) and p != UNASSIGNED and p != NA_PATH


def churn(pre: pl.DataFrame, post: pl.DataFrame, truth: dict) -> dict:
    a = dict(zip(pre["read_id"], pre["pred_path"]))
    b = dict(zip(post["read_id"], post["pred_path"]))
    bk = {"improved": 0, "regressed": 0, "lateral": 0, "newly_assigned": 0,
          "newly_unassigned": 0, "unresolvable_both": 0}
    n_both = n_changed = 0
    for k in truth:
        if k not in a or k not in b:
            continue
        pa, pb = a[k], b[k]
        if scoreable(pa) and not scoreable(pb):
            bk["newly_unassigned"] += 1; continue
        if not scoreable(pa) and scoreable(pb):
            bk["newly_assigned"] += 1; continue
        if not scoreable(pa) and not scoreable(pb):
            bk["unresolvable_both"] += 1; continue
        n_both += 1
        if pa == pb:
            continue
        n_changed += 1
        t = split_path(truth[k])
        da = compute_taxonomic_distance(t, split_path(pa))
        db = compute_taxonomic_distance(t, split_path(pb))
        bk["improved" if db < da else "regressed" if db > da else "lateral"] += 1
    moved = bk["improved"] + bk["regressed"]
    return {**bk, "scoreable_in_both": n_both, "path_changed": n_changed,
            "churn_pct": round(100 * n_changed / n_both, 2) if n_both else None,
            "toward_truth_ratio": round(bk["improved"] / moved, 4) if moved else None}


def cmd_score(marker: str, run_root: Path, out: Path) -> None:
    p = preset(marker)
    farm, comms = discover(marker, run_root)
    CACHE.mkdir(parents=True, exist_ok=True)
    print("building taxonomy remapper (shared by every cell, so it cancels)...")
    remapper = build_runtime_tronko_remapper(
        marker=p.marker.name, main_rcrux_db_dir=p.marker.cruxv2_dir,
        taxdump_dir=TAXDUMP, cache_dir=CACHE,
        reference_stem=p.marker.reference_file_stem)
    print(f"  cache_key={remapper.cache_key[:16]}")

    report = {"marker": marker, "remapper_cache_key": remapper.cache_key, "communities": {}}
    for idx in comms:
        gt = json.loads((farm / f"q{idx}" / "asv_ground_truth.json").read_text())
        truth = gt["paired_asv_ground_truth"]
        split = json.loads(((farm / f"q{idx}").resolve().parent / "reads"
                            / "ground_truth.json").read_text()).get("unified_split", {})
        ov = {s.lower(): "genus" for s in split.get("species_holdout_species", [])}
        ov.update({s.lower(): "family" for s in split.get("genus_holdout_species", [])})

        lengths: dict[str, int] = {}
        for mate in ("F", "R"):
            h = None
            for line in (out / "queries" / "fwd" / f"c{idx}_{mate}.fasta").open():
                line = line.rstrip()
                if line.startswith(">"):
                    h = line[1:].replace(f"_paired_{mate}_", "_paired_F_")
                elif h:
                    lengths[h] = lengths.get(h, 0) + len(line)

        entry = {"n_truth": len(truth), "cells": {}}
        cells = {}
        for binary in ("pre", "post"):
            for arm in ("fwd", "rc"):
                cell = load_cell(out, binary, arm, idx, remapper)
                cells[(binary, arm)] = cell
                rates = [(f + r) / lengths[i] for i, f, r in
                         zip(cell["read_id"], cell["fwd_mm"], cell["rev_mm"])
                         if lengths.get(i) and f is not None and r is not None]
                acc = {"truth_depth_f1": round(compute_truth_depth_asv_metrics(
                    paired_frame(cell, truth), eval_rank_override=ov).f1, 4)}
                for rank in ("phylum", "class", "order", "family", "genus", "species"):
                    m = compute_per_rank_metrics(paired_frame(cell, truth), rank)
                    acc[rank] = {"f1": round(m.f1, 4), "precision": round(m.precision, 4),
                                 "recall": round(m.recall, 4)}
                entry["cells"][f"{binary}_{arm}"] = {
                    "strand": {"median_divergence": round(statistics.median(rates), 4),
                               "frac_above_25pct": round(sum(1 for x in rates if x > 0.25) / len(rates), 4),
                               "n": len(rates)} if rates else {},
                    "accuracy": acc}
        entry["churn_fwd"] = churn(cells[("pre", "fwd")], cells[("post", "fwd")], truth)
        entry["churn_rc"] = churn(cells[("pre", "rc")], cells[("post", "rc")], truth)
        report["communities"][idx] = entry
        print(f"  community {idx} scored")

    (out / "tronko_minimap2_ab.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out / 'tronko_minimap2_ab.json'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--marker", required=True)
    ap.add_argument("--run-root", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    out = args.out or (BASE / "runs" / args.marker)
    (cmd_run if args.action == "run" else cmd_score)(args.marker, args.run_root, out)


if __name__ == "__main__":
    main()
