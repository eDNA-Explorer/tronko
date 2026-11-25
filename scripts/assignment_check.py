#!/usr/bin/env python3
"""
Check that a finished build assigns reads correctly, with a synthetic community.

The community, its reads and the scoring are synthetic-fastq-toolkit's, the
code the assignment benchmark uses (edna-explorer-data-pipelines,
libraries/synthetic-fastq-toolkit): species drawn from the build's own
reference with equal abundance, read pairs as QC hands them to assignment, and
the benchmark's per-rank metrics. Equal abundance gives each species the same
weight in the scores, as the classifier benchmarks do; a skewed community lets
a few species decide them. tronko-assign assigns them the way
production does (paired, -z -w -6, LCA cut-off 10).

The reads come from sequences in the database, so this shows the database and
tronko-assign work together, not how well it does on species it has not seen.

Usage:
    python assignment_check.py --build-dir OUT --marker MARKER --threads N \
        [--species 100] [--reads 5000] [--read-length 150] [--seed 1] --out result.json
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from synthetic_fastq.benchmark_reads import load_read_truth, simulate_community
from synthetic_fastq.community import RankAbundanceModel
from synthetic_fastq.scoring.reads import score_reads

# production's tronko-assign settings (qc-tools-lib LCA_CUTOFF)
LCA_CUTOFF = 10


def run(
    build_dir: Path,
    marker: str,
    threads: int,
    species: int,
    reads: int,
    read_length: int,
    seed: int,
) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        community = simulate_community(
            reference_dir=build_dir,
            marker=marker,
            num_species=species,
            num_reads=reads,
            read_length=read_length,
            output_dir=work,
            seed=seed,
            abundance_model=RankAbundanceModel.UNIFORM,
        )
        output = work / "assigned.txt"
        command = [
            "tronko-assign", "-r", "-f", str(build_dir / "reference_tree.trkb"),
            "-a", str(build_dir / f"{marker}.fasta"), "-q",
            "-p", "-z", "-w", "-1", str(community.reads.r1_path), "-2", str(community.reads.r2_path),
            "-6", "--number-of-cores", str(threads), "--Cinterval", str(LCA_CUTOFF),
            "-o", str(output),
        ]  # fmt: skip
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"tronko-assign failed: {completed.stderr.strip()[-2000:]}")
        truth = load_read_truth(community.reads.read_truth_path)
        predictions = {}
        with open(output) as handle:
            next(handle)  # header
            for line in handle:
                name, path = line.rstrip("\n").split("\t")[:2]
                # tronko-assign keeps the whole FASTQ header, joined with "_"
                # ("<id>_1:N:0:ATCG"); the truth has the id alone.
                predictions[name if name in truth else name.rsplit("_", 1)[0]] = path

    scores = score_reads(truth, predictions)
    return {
        "settings": {
            "species": len(community.abundances),
            "reads": len(truth),
            "read_length": read_length,
            "seed": seed,
            "lca_cutoff": LCA_CUTOFF,
        },
        "by_rank": {
            rank: {
                "correct": score.correct,
                "unresolved": score.unresolved,
                "wrong": score.wrong,
                "precision": score.metrics.precision,
                "recall": score.metrics.recall,
                "f1": score.metrics.f1,
                "mcc": score.metrics.mcc,
            }
            for rank, score in scores.items()
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--species", type=int, default=100)
    parser.add_argument("--reads", type=int, default=5000)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run(
        args.build_dir, args.marker, args.threads, args.species, args.reads,
        args.read_length, args.seed,
    )  # fmt: skip
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["by_rank"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
