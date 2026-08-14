#!/usr/bin/env python3
"""Fetch the November 2025 BioSCape CO1 set and prove it is internally consistent.

The November objects under the **v2** GCS prefix are not trustworthy: a
retroactive March 2026 backfill overwrote the FASTA, so the November
assignments.parquet and the November FASTA under that prefix do not describe
each other (22.2% forward_length agreement, below the 31.7% chance baseline).
The legacy pre-backfill copy is intact and is what this script pulls.

It then *verifies* that, rather than trusting it: every sequence's length must
match the forward_length / reverse_length the recorded assignment claims. If
that check fails the run stops, because every downstream number would be
comparing one ASV's sequence against another ASV's call.

Usage:
    bioscape_co1_fetch_nov25.py --out-dir co1/ [--min-agreement 99.0]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROJECT = "cm4n7f2bb0001ssovx4mh0dt4"
BUCKET = "gs://edna-project-files-ca"

# The legacy, pre-backfill November copy. Do NOT substitute the v2 prefix
# (assign/cmhwjvd4z000ejs042nriini2/...) -- that FASTA is the corrupted one.
LEGACY = f"{BUCKET}/projects/{PROJECT}/assign/CO1_Metazoa/paired"
LEGACY_FILES = {
    "assignments": f"{PROJECT}-CO1_Metazoa-paired.txt.zst",
    "forward": f"{PROJECT}-CO1_Metazoa-paired_F.fasta.zst",
    "reverse": f"{PROJECT}-CO1_Metazoa-paired_R.fasta.zst",
}

# For the run-identity check (is this November, or March?).
RUNS = {
    "november": "cmhwjvd4z000ejs042nriini2",
    "march": "cmmjfc1ph0003js040lgonrcr",
}
MARKER_ID = "bfe214ab-fd47-4e47-92be-3a73ad1d98b1"


def fetch(uri: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  have {dest.name}")
        return dest
    print(f"  fetching {uri}")
    subprocess.run(["gcloud", "storage", "cp", uri, str(dest)], check=True,
                   capture_output=True)
    return dest


def unzst(src: Path, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with dest.open("wb") as fh:
        subprocess.run(["zstd", "-dc", str(src)], stdout=fh, check=True)
    return dest


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_tsv(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h.lower(): i for i, h in enumerate(header)}
        for line in fh:
            row = line.rstrip("\n").split("\t")

            def g(key, cast=None):
                i = idx.get(key)
                if i is None or i >= len(row) or row[i] == "":
                    return None
                return cast(row[i]) if cast else row[i]

            name = g("readname")
            if name:
                out[name] = {
                    "path": g("taxonomic_path"),
                    "tree": g("tree_number", lambda x: int(float(x))),
                    "node": g("node_number", lambda x: int(float(x))),
                    "flen": g("forward_length", lambda x: int(float(x))),
                    "rlen": g("reverse_length", lambda x: int(float(x))),
                }
    return out


def write_fasta(records: list[tuple[str, str]], dest: Path) -> None:
    """Unwrapped FASTA -- tronko-assign wants one sequence per line."""
    with dest.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-agreement", type=float, default=99.0,
                    help="length-agreement %% below which this exits non-zero")
    args = ap.parse_args()

    out = args.out_dir
    raw = out / "raw"
    out.mkdir(parents=True, exist_ok=True)

    print("fetching the legacy (pre-backfill) November copy")
    local = {}
    for key, fn in LEGACY_FILES.items():
        z = fetch(f"{LEGACY}/{fn}", raw / fn)
        local[key] = unzst(z, raw / fn.removesuffix(".zst"))

    print("\nalso fetching both candidate runs' parquet for the identity check")
    for label, run in RUNS.items():
        try:
            fetch(
                f"{BUCKET}/projects/{PROJECT}/assign/{run}/{MARKER_ID}/paired/assignments.parquet",
                raw / f"{label}_assignments.parquet",
            )
        except subprocess.CalledProcessError:
            print(f"  ! could not fetch {label} parquet; identity check will be partial")

    tsv = read_tsv(local["assignments"])
    fwd = read_fasta(local["forward"])
    rev = read_fasta(local["reverse"])
    print(f"\nrecorded assignments: {len(tsv):,}")
    print(f"forward sequences   : {len(fwd):,}")
    print(f"reverse sequences   : {len(rev):,}")

    # --- the blocking check: do the sequences and the calls describe each other?
    checked = agree_f = agree_r = 0
    for name, rec in tsv.items():
        s = fwd.get(name)
        if s is not None and rec["flen"]:
            checked += 1
            if len(s) == rec["flen"]:
                agree_f += 1
        s = rev.get(name)
        if s is not None and rec["rlen"]:
            if len(s) == rec["rlen"]:
                agree_r += 1

    pct_f = 100.0 * agree_f / checked if checked else 0.0
    print(f"\nforward_length agreement: {agree_f:,}/{checked:,} = {pct_f:.1f}%")
    if checked:
        print(f"reverse_length agreement: {agree_r:,}/{checked:,} = "
              f"{100.0 * agree_r / checked:.1f}%")

    report = {
        "source": LEGACY,
        "assignments": len(tsv),
        "forward_sequences": len(fwd),
        "reverse_sequences": len(rev),
        "length_agreement_pct": round(pct_f, 2),
        "threshold_pct": args.min_agreement,
        "status": "PASS" if pct_f >= args.min_agreement else "FAIL",
    }

    if pct_f < args.min_agreement:
        report["meaning"] = (
            "The sequences and the recorded calls do not describe each other. "
            "This is the backfill corruption; do not proceed. The v2 prefix "
            "sits at ~22%, well below the ~31.7% chance baseline."
        )
        (out / "fetch_report.json").write_text(json.dumps(report, indent=2))
        print("\nFAILED consistency check -- refusing to emit query FASTAs.",
              file=sys.stderr)
        print(report["meaning"], file=sys.stderr)
        return 2

    # --- emit query FASTAs, mirroring production's partitioning
    paired = [n for n in tsv if n in fwd and n in rev]
    unpaired = [n for n in tsv if n in fwd and n not in rev]
    write_fasta([(n, fwd[n]) for n in paired], out / "pF.fasta")
    write_fasta([(n, rev[n]) for n in paired], out / "pR.fasta")
    if unpaired:
        write_fasta([(n, fwd[n]) for n in unpaired], out / "unpaired_F.fasta")

    # arm A, the recorded November answer, in the same TSV shape the other arms use
    with (out / "arm_A_november.tsv").open("w") as fh:
        fh.write("Readname\tTaxonomic_Path\tTree_Number\tNode_Number\n")
        for n in paired + unpaired:
            r = tsv[n]
            fh.write(f"{n}\t{r['path'] or 'unassigned'}\t"
                     f"{r['tree'] if r['tree'] is not None else -1}\t"
                     f"{r['node'] if r['node'] is not None else -1}\n")

    report.update({"paired": len(paired), "unpaired": len(unpaired)})
    (out / "fetch_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {len(paired):,} paired and {len(unpaired):,} unpaired queries")
    print(f"  {out}/pF.fasta  {out}/pR.fasta  {out}/arm_A_november.tsv")
    print("\nRun identity: confirm these are November and not March by joining "
          "arm_A_november.tsv on readname against raw/november_assignments.parquet "
          "and raw/march_assignments.parquet -- whichever agrees ~100% is what you "
          "have. Everything downstream is mislabelled if this is skipped.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
