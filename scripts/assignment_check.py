#!/usr/bin/env python3
"""
Check how a finished build assigns reads, with synthetic holdout communities.

The design, reads, holdout references and scoring are community-lib's
(edna-explorer-data-pipelines, libraries/community-lib), the corpus design the
classifier benchmarks use:

- species-disjoint communities drawn from the build's own reference, each split
  into four tiers: normal, haplotype (the query's own record removed), species
  holdout and genus holdout;
- per community, a Tronko reference ablated from this build: the build less
  that community's three holdout tiers, rebuilding only the partitions that
  lose records;
- per community and error profile (illumina_typical, edna_damage), read pairs
  from every query, assigned with tronko-assign the way production does
  (paired, -z -w -6, LCA cut-off 10) against that community's reference;
- each read scored at the rank its truth can be resolved to, per tier.

The build is never modified; everything is written under --work.

Usage:
    python assignment_check.py --build-dir OUT --marker MARKER --inputs inputs.json \
        --threads N --work DIR [--communities 6] [--richness 150] \
        [--reads-per-query 20] --out result.json
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from community_lib import check
from taxonomy_classifier_harness.tronko.database import TronkoBuildConfig

# production's tronko-assign settings (qc-tools-lib LCA_CUTOFF)
LCA_CUTOFF = 10


def assign(reference: Path, marker: str, reads: check.MemberReads, threads: int, out: Path) -> dict:
    """tronko-assign's call for every read, keyed by the read's truth id."""
    command = [
        "tronko-assign", "-r", "-f", str(reference / "reference_tree.trkb"),
        "-a", str(reference / f"{marker}.fasta"), "-q",
        "-p", "-z", "-w", "-1", str(reads.r1), "-2", str(reads.r2),
        "-6", "--number-of-cores", str(threads), "--Cinterval", str(LCA_CUTOFF),
        "-o", str(out),
    ]  # fmt: skip
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"tronko-assign failed: {completed.stderr.strip()[-2000:]}")
    ids = {row["read_id"] for row in check.load_read_truth(reads.truth)}
    predictions = {}
    with open(out) as handle:
        next(handle)  # header
        for line in handle:
            name, path = line.rstrip("\n").split("\t")[:2]
            # tronko-assign keeps the whole FASTQ header, words joined with "_";
            # the truth has its first word, which itself holds "_".
            while name not in ids and "_" in name:
                name = name.rsplit("_", 1)[0]
            predictions[name] = path
    return predictions


def run(
    build_dir: Path,
    marker: str,
    inputs: dict,
    threads: int,
    work: Path,
    communities: int,
    richness: int,
    reads_per_query: int,
) -> dict:
    primers = {"forward_primer": inputs["forward_primer"], "reverse_primer": inputs["reverse_primer"]}
    design = check.design_communities(
        build_dir / f"{marker}.fasta",
        build_dir / f"{marker}_taxonomy.txt",
        work / "design",
        marker=marker,
        richness=richness,
        n_communities=communities,
        **primers,
    )
    config = TronkoBuildConfig(primer_name=marker, threads=threads, famsa_threads=threads)
    members, scores = [], {profile: {} for profile in design.profiles}
    for member in design.members:
        holdout = check.build_holdout_tronko_reference(
            member, build_dir, work / "reference" / member.name, marker=marker, build_config=config
        )
        members.append(
            {
                "name": member.name,
                "seed": member.seed,
                "species_per_tier": {tier: len(species) for tier, species in member.tiers.items()},
                "removed_records": holdout["removed_records"],
                # the ablated reference is not kept; this is what remains of its check
                "audit": holdout["audit"],
            }
        )
        for profile in design.profiles:
            reads = check.write_member_reads(
                member, profile, work / "reads" / member.name / profile,
                reads_per_query=reads_per_query, **primers,
            )  # fmt: skip
            predictions = assign(
                Path(holdout["dir"]), marker, reads, threads,
                work / "reads" / member.name / profile / "assigned.txt",
            )  # fmt: skip
            scores[profile][member.name] = check.score_reads(reads.truth, predictions)

    summary = {}
    for profile, by_member in scores.items():
        tiers = next(iter(by_member.values()))["by_tier"]
        summary[profile] = {
            tier: {
                key: sum(s["by_tier"][tier][key] for s in by_member.values()) / len(by_member)
                for key in ("precision", "recall", "f1")
            }
            | {
                "reads": sum(s["by_tier"][tier]["reads"] for s in by_member.values()),
                "fp": sum(s["by_tier"][tier]["fp"] for s in by_member.values()),
            }
            for tier in tiers
        }
    return {
        "settings": {
            "communities": len(design.members),
            "richness": richness,
            "reads_per_query": reads_per_query,
            "profiles": list(design.profiles),
            "lca_cutoff": LCA_CUTOFF,
            **primers,
        },
        "communities": members,
        "summary": summary,
        "scores": scores,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--inputs", type=Path, required=True, help="fetch's inputs.json (the primers)")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--communities", type=int, default=6)
    parser.add_argument("--richness", type=int, default=150)
    parser.add_argument("--reads-per-query", type=int, default=20)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    inputs = json.loads(args.inputs.read_text())
    if not (inputs.get("forward_primer") and inputs.get("reverse_primer")):
        raise SystemExit(f"{args.inputs} records no primers; fetch with --forward-primer and --reverse-primer")
    result = run(
        args.build_dir, args.marker, inputs, args.threads, args.work, args.communities,
        args.richness, args.reads_per_query,
    )  # fmt: skip
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
