#!/usr/bin/env python3
"""
Check that a finished build assigns reads correctly, with a synthetic community.

Built from synthetic-fastq-toolkit (edna-explorer-data-pipelines,
libraries/synthetic-fastq-toolkit), with the read profile of the velocitaxa
classifier benchmark's read arm: species drawn from the build's own reference,
up to two query sequences per species (the first distinct ones by accession),
and the same number of read pairs from every query, MiSeq, 300 bp, with
sequencing errors. Every query weighs the same, so no few species decide the
score. The reads are what QC hands to assignment: no primers, no adapters.
tronko-assign assigns them the way production does (paired, -z -w -6, LCA
cut-off 10), and they are scored with the benchmark's per-rank metrics.

The reads come from sequences in the database, so this shows the database and
tronko-assign work together, not how well it does on species it has not seen.

Usage:
    python assignment_check.py --build-dir OUT --marker MARKER --threads N \
        [--species 100] [--queries-per-species 2] [--reads-per-query 20] \
        [--read-length 300] [--seed 1] --out result.json
"""

import argparse
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from synthetic_fastq.reference import CruxV2Database, load_cruxv2_database, subsample_species
from synthetic_fastq.scoring.reads import score_reads
from synthetic_fastq.simulate_reads import simulate_amplicon_paired_reads_with_sources

# production's tronko-assign settings (qc-tools-lib LCA_CUTOFF)
LCA_CUTOFF = 10
INSTRUMENT = "MiSeq"


def queries(db: CruxV2Database, species: int, per_species: int, seed: int) -> list[str]:
    """Accessions to read from: up to per_species distinct sequences of each drawn species."""
    sequence = dict(db.sequences)
    named = {a: r for a, r in db.taxonomy.items() if a in sequence and r.species not in ("", "NA")}
    drawable = CruxV2Database(
        marker=db.marker,
        sequences=[(a, s) for a, s in db.sequences if a in named],
        taxonomy=named,
        source_dir=db.source_dir,
    )
    community = subsample_species(drawable, species, random.Random(seed))
    chosen = []
    for _, accessions in sorted(community.species_to_accessions.items()):
        seen = set()
        for accession in sorted(accessions):
            if sequence[accession] not in seen and len(seen) < per_species:
                seen.add(sequence[accession])
                chosen.append(accession)
    return chosen


def run(
    build_dir: Path,
    marker: str,
    threads: int,
    species: int,
    queries_per_species: int,
    reads_per_query: int,
    read_length: int,
    seed: int,
) -> dict:
    db = load_cruxv2_database(build_dir, marker)
    sequence = dict(db.sequences)
    chosen = queries(db, species, queries_per_species, seed)
    rng = random.Random(seed)
    truth, r1_all, r2_all = {}, [], []
    for accession in chosen:
        r1, r2, _ = simulate_amplicon_paired_reads_with_sources(
            [(accession, sequence[accession])],
            read_length=read_length,
            num_reads=reads_per_query,
            instrument=INSTRUMENT,
            clean=False,
            rng=rng,
            r1_adapter_seq="",
            r2_adapter_seq="",
            pad_to_read_length=False,
        )
        for (_, s1, q1), (_, s2, q2) in zip(r1, r2, strict=True):
            # A one-word name: tronko-assign keeps it as is.
            name = f"q{len(truth)}"
            truth[name] = db.taxonomy[accession].taxonomic_path
            r1_all.append(f"@{name}\n{s1}\n+\n{q1}\n")
            r2_all.append(f"@{name}\n{s2}\n+\n{q2}\n")

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / "R1.fastq").write_text("".join(r1_all))
        (work / "R2.fastq").write_text("".join(r2_all))
        output = work / "assigned.txt"
        command = [
            "tronko-assign", "-r", "-f", str(build_dir / "reference_tree.trkb"),
            "-a", str(build_dir / f"{marker}.fasta"), "-q",
            "-p", "-z", "-w", "-1", str(work / "R1.fastq"), "-2", str(work / "R2.fastq"),
            "-6", "--number-of-cores", str(threads), "--Cinterval", str(LCA_CUTOFF),
            "-o", str(output),
        ]  # fmt: skip
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"tronko-assign failed: {completed.stderr.strip()[-2000:]}")
        with open(output) as handle:
            next(handle)  # header
            predictions = dict(line.rstrip("\n").split("\t")[:2] for line in handle)

    scores = score_reads(truth, predictions)
    return {
        "settings": {
            "species": len({db.taxonomy[a].species for a in chosen}),
            "queries": len(chosen),
            "reads_per_query": reads_per_query,
            "reads": len(truth),
            "read_length": read_length,
            "instrument": INSTRUMENT,
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
    parser.add_argument("--queries-per-species", type=int, default=2)
    parser.add_argument("--reads-per-query", type=int, default=20)
    parser.add_argument("--read-length", type=int, default=300)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run(
        args.build_dir, args.marker, args.threads, args.species, args.queries_per_species,
        args.reads_per_query, args.read_length, args.seed,
    )  # fmt: skip
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["by_rank"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
