#!/usr/bin/env python3
"""
Tronko-build parameter sweep: build databases with different configs
and compute pre-ablation quality metrics.

Usage:
    python sweep_tronko_build.py --input-dir <merged_clusters> --output-dir <sweep_output>

Each config produces a separate database directory. Quality metrics are
computed from the reference_tree.txt + taxonomy files WITHOUT running
the full ablation pipeline — giving fast signal on tree quality.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional
import math
import gzip


# ── Sweep configurations ─────────────────────────────────────────────

@dataclass
class BuildConfig:
    """A single tronko-build parameter configuration."""
    name: str
    # Partitioning strategy: "sp" or "minleaf"
    partition_mode: str = "sp"
    # SP-score threshold (only used when partition_mode="sp")
    sp_threshold: float = 0.1
    # Min-leaf count (only used when partition_mode="minleaf")
    min_leaf: int = 50
    # Whether to use FastTree (-a flag)
    use_fasttree: bool = True
    # Whether to adjust for missing data (True=default, False=-g flag)
    missing_data: bool = True

    def tronko_flags(self) -> list[str]:
        """Return the tronko-build-specific flags for this config."""
        flags = []
        if self.partition_mode == "sp":
            flags.extend(["-s", "-u", str(self.sp_threshold)])
        elif self.partition_mode == "minleaf":
            flags.extend(["-v", "-f", str(self.min_leaf)])
        elif self.partition_mode == "none":
            flags.extend(["-v", "-f", "999999"])
        if self.use_fasttree:
            flags.append("-a")
        if not self.missing_data:
            flags.append("-g")
        return flags

    def short_desc(self) -> str:
        if self.partition_mode == "sp":
            return f"SP={self.sp_threshold}"
        elif self.partition_mode == "minleaf":
            return f"MinLeaf={self.min_leaf}"
        else:
            return "NoPartition"


# Phase 1 configs: vary partitioning strategy (cheap to compare)
SWEEP_CONFIGS = [
    # SP-score threshold sweep
    BuildConfig(name="sp_0.05", partition_mode="sp", sp_threshold=0.05),
    BuildConfig(name="sp_0.10", partition_mode="sp", sp_threshold=0.10),  # baseline
    BuildConfig(name="sp_0.20", partition_mode="sp", sp_threshold=0.20),
    BuildConfig(name="sp_0.50", partition_mode="sp", sp_threshold=0.50),
    # Min-leaf threshold sweep (alternative strategy)
    BuildConfig(name="minleaf_50",  partition_mode="minleaf", min_leaf=50),
    BuildConfig(name="minleaf_200", partition_mode="minleaf", min_leaf=200),
    # Missing data toggle (with baseline SP)
    BuildConfig(name="sp_0.10_no_missdata", partition_mode="sp", sp_threshold=0.10, missing_data=False),
]


# ── Quality metrics ──────────────────────────────────────────────────

@dataclass
class TreeStats:
    """Statistics for a single tree/partition."""
    tree_idx: int
    num_leaves: int
    num_nodes: int
    num_bases: int
    # Taxonomic metrics (populated later)
    genus_purity: float = 0.0    # fraction of leaves in the dominant genus
    family_purity: float = 0.0   # fraction of leaves in the dominant family
    num_genera: int = 0
    num_families: int = 0
    num_species: int = 0


@dataclass
class QualityReport:
    """Aggregate quality metrics for one database build."""
    config_name: str
    build_time_seconds: float = 0.0
    reference_tree_size_mb: float = 0.0

    # Partition statistics
    num_trees: int = 0
    total_leaves: int = 0
    total_nodes: int = 0
    leaves_per_tree_min: int = 0
    leaves_per_tree_max: int = 0
    leaves_per_tree_mean: float = 0.0
    leaves_per_tree_median: float = 0.0
    leaves_per_tree_p25: float = 0.0
    leaves_per_tree_p75: float = 0.0
    singleton_trees: int = 0         # trees with only 1 leaf
    small_trees: int = 0             # trees with <= 3 leaves

    # Taxonomic purity (averaged across trees, weighted by leaf count)
    mean_genus_purity: float = 0.0   # how genus-pure are trees?
    mean_family_purity: float = 0.0  # how family-pure are trees?
    mean_genera_per_tree: float = 0.0
    mean_species_per_tree: float = 0.0

    # Monophyly (global)
    genus_monophyly_rate: float = 0.0   # fraction of genera that are in a single tree
    family_monophyly_rate: float = 0.0  # fraction of families that are in a single tree

    # Taxonomic coverage
    total_genera: int = 0
    total_families: int = 0
    total_species: int = 0

    # Taxonomic entropy (lower = more uniform, higher = more mixed)
    mean_genus_entropy: float = 0.0


def percentile(sorted_vals: list, p: float) -> float:
    """Compute the p-th percentile (0-100) of a sorted list."""
    if not sorted_vals:
        return 0.0
    k = (len(sorted_vals) - 1) * p / 100.0
    f = int(k)
    c = f + 1 if f + 1 < len(sorted_vals) else f
    d = k - f
    return sorted_vals[f] + d * (sorted_vals[c] - sorted_vals[f])


def shannon_entropy(counts: dict) -> float:
    """Shannon entropy of a frequency distribution."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    ent = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            ent -= p * math.log2(p)
    return ent


def parse_reference_tree_header(ref_tree_path: Path) -> dict:
    """Parse the header of reference_tree.txt to get basic stats.

    Format:
        line 1: numberOfTrees
        line 2: max_nodename
        line 3: max_tax_name
        line 4: max_lineTaxonomy
        then per-tree: numbase\troot\tnumspec
        then taxonomy blocks
        then node+posterior data
    """
    opener = gzip.open if ref_tree_path.suffix == ".gz" else open
    with opener(ref_tree_path, "rt") as f:
        num_trees = int(f.readline().strip())
        max_nodename = int(f.readline().strip())
        max_tax_name = int(f.readline().strip())
        max_line_taxonomy = int(f.readline().strip())

        trees = []
        for i in range(num_trees):
            parts = f.readline().strip().split("\t")
            numbase = int(parts[0])
            root = int(parts[1])
            numspec = int(parts[2])
            trees.append({
                "idx": i,
                "numbase": numbase,
                "root": root,
                "numspec": numspec,
            })

        # Parse taxonomy blocks: for each tree, numspec lines of 7 semicolon-separated ranks
        taxonomies = []  # list of list of 7-tuples
        for i in range(num_trees):
            tree_tax = []
            for j in range(trees[i]["numspec"]):
                line = f.readline().strip()
                ranks = line.replace("\t", ";").split(";")
                # Pad to 7 if needed
                while len(ranks) < 7:
                    ranks.append("")
                tree_tax.append(ranks[:7])
            taxonomies.append(tree_tax)

    return {
        "num_trees": num_trees,
        "trees": trees,
        "taxonomies": taxonomies,
    }


def compute_quality_metrics(
    ref_tree_path: Path,
    config_name: str,
    build_time: float = 0.0,
) -> QualityReport:
    """Compute quality metrics from a reference_tree.txt file."""
    report = QualityReport(config_name=config_name, build_time_seconds=build_time)

    # File size
    report.reference_tree_size_mb = ref_tree_path.stat().st_size / (1024 * 1024)

    print(f"  Parsing reference tree header...")
    data = parse_reference_tree_header(ref_tree_path)
    report.num_trees = data["num_trees"]

    trees = data["trees"]
    taxonomies = data["taxonomies"]

    # ── Partition size statistics ─────────────────────────────────
    leaf_counts = sorted([t["numspec"] for t in trees])
    report.total_leaves = sum(leaf_counts)
    report.total_nodes = sum(2 * t["numspec"] - 1 for t in trees)
    report.leaves_per_tree_min = leaf_counts[0] if leaf_counts else 0
    report.leaves_per_tree_max = leaf_counts[-1] if leaf_counts else 0
    report.leaves_per_tree_mean = sum(leaf_counts) / len(leaf_counts) if leaf_counts else 0
    report.leaves_per_tree_median = percentile(leaf_counts, 50)
    report.leaves_per_tree_p25 = percentile(leaf_counts, 25)
    report.leaves_per_tree_p75 = percentile(leaf_counts, 75)
    report.singleton_trees = sum(1 for c in leaf_counts if c == 1)
    report.small_trees = sum(1 for c in leaf_counts if c <= 3)

    # ── Taxonomic purity per tree ─────────────────────────────────
    # Track which trees each genus/family appears in (for monophyly)
    genus_to_trees: dict[str, set[int]] = defaultdict(set)
    family_to_trees: dict[str, set[int]] = defaultdict(set)
    all_genera = set()
    all_families = set()
    all_species = set()

    genus_purity_sum = 0.0
    family_purity_sum = 0.0
    genera_per_tree_sum = 0.0
    species_per_tree_sum = 0.0
    genus_entropy_sum = 0.0
    weight_sum = 0

    for i, tree_tax in enumerate(taxonomies):
        numspec = trees[i]["numspec"]
        if numspec == 0:
            continue

        # reference_tree.txt stores ranks REVERSED: species(0);genus(1);family(2);order(3);class(4);phylum(5);domain(6)
        genera = Counter()
        families = Counter()
        species_set = set()

        for ranks in tree_tax:
            genus = ranks[1] if len(ranks) > 1 else ""
            family = ranks[2] if len(ranks) > 2 else ""
            species = ranks[0] if len(ranks) > 0 else ""

            if genus:
                genera[genus] += 1
                genus_to_trees[genus].add(i)
                all_genera.add(genus)
            if family:
                families[family] += 1
                family_to_trees[family].add(i)
                all_families.add(family)
            if species:
                species_set.add(species)
                all_species.add(species)

        # Purity = fraction of leaves belonging to the most common taxon
        if genera:
            genus_purity = max(genera.values()) / numspec
            genus_purity_sum += genus_purity * numspec
        if families:
            family_purity = max(families.values()) / numspec
            family_purity_sum += family_purity * numspec

        genera_per_tree_sum += len(genera)
        species_per_tree_sum += len(species_set)
        genus_entropy_sum += shannon_entropy(genera)
        weight_sum += numspec

    if weight_sum > 0:
        report.mean_genus_purity = genus_purity_sum / weight_sum
        report.mean_family_purity = family_purity_sum / weight_sum
    report.mean_genera_per_tree = genera_per_tree_sum / max(report.num_trees, 1)
    report.mean_species_per_tree = species_per_tree_sum / max(report.num_trees, 1)
    report.mean_genus_entropy = genus_entropy_sum / max(report.num_trees, 1)

    # ── Monophyly ─────────────────────────────────────────────────
    # A genus is "monophyletic" if all its members appear in a single tree
    if all_genera:
        mono_genera = sum(1 for g, ts in genus_to_trees.items() if len(ts) == 1)
        report.genus_monophyly_rate = mono_genera / len(all_genera)
    if all_families:
        mono_families = sum(1 for f, ts in family_to_trees.items() if len(ts) == 1)
        report.family_monophyly_rate = mono_families / len(all_families)

    report.total_genera = len(all_genera)
    report.total_families = len(all_families)
    report.total_species = len(all_species)

    return report


# ── Build runner ─────────────────────────────────────────────────────

def count_clusters(cluster_dir: Path) -> int:
    """Count the number of *_MSA.fasta files in a directory."""
    return sum(1 for f in cluster_dir.iterdir() if f.name.endswith("_MSA.fasta"))


def run_tronko_build(
    config: BuildConfig,
    cluster_dir: Path,
    num_clusters: int,
    output_dir: Path,
    tronko_binary: str,
    threads: int,
) -> float:
    """Run tronko-build with the given config. Returns build time in seconds."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        tronko_binary,
        "-y",
        "-e", str(cluster_dir),
        "-n", str(num_clusters),
        "-d", str(output_dir),
    ]
    cmd.extend(config.tronko_flags())
    cmd.extend(["-c", str(threads)])

    print(f"\n{'='*70}")
    print(f"Config: {config.name} ({config.short_desc()})")
    print(f"Command: {' '.join(cmd)}")
    print(f"Output: {output_dir}")
    print(f"{'='*70}")

    stdout_log = output_dir / "_build_stdout.log"
    stderr_log = output_dir / "_build_stderr.log"

    t0 = time.monotonic()
    with open(stdout_log, "w") as f_out, open(stderr_log, "w") as f_err:
        proc = subprocess.Popen(cmd, stdout=f_out, stderr=f_err)

        # Poll for progress
        progress_file = output_dir / "_progress.txt"
        last_msg = ""
        while proc.poll() is None:
            time.sleep(10)
            if progress_file.exists():
                try:
                    cur = progress_file.read_text().strip()
                    if cur != last_msg:
                        print(f"  [{config.name}] {cur}")
                        last_msg = cur
                except OSError:
                    pass

    elapsed = time.monotonic() - t0

    if proc.returncode != 0:
        stderr_text = stderr_log.read_text() if stderr_log.exists() else ""
        print(f"  FAILED (exit code {proc.returncode})")
        if stderr_text:
            print(f"  stderr: {stderr_text[:500]}")
        return -1.0

    print(f"  Completed in {elapsed:.1f}s")

    # Print last few lines of stdout
    if stdout_log.exists():
        lines = stdout_log.read_text().strip().splitlines()
        for line in lines[-5:]:
            print(f"  > {line}")

    # Clean up progress file
    progress_file.unlink(missing_ok=True)

    # Convert reference_tree.txt to .trkb (required by tronko-assign)
    ref_txt = output_dir / "reference_tree.txt"
    ref_trkb = output_dir / "reference_tree.trkb"
    tronko_convert = shutil.which("tronko-convert")
    if ref_txt.exists() and tronko_convert:
        print(f"  Converting reference tree to .trkb...")
        conv = subprocess.run(
            [tronko_convert, "-i", str(ref_txt), "-o", str(ref_trkb)],
            capture_output=True, text=True,
        )
        if conv.returncode == 0:
            print(f"  .trkb created ({ref_trkb.stat().st_size / 1024 / 1024:.1f} MB)")
        else:
            print(f"  WARNING: tronko-convert failed: {conv.stderr[:200]}")

    return elapsed


def find_reference_tree(output_dir: Path) -> Optional[Path]:
    """Find the reference_tree.txt file in the output directory."""
    for name in ["reference_tree.txt", "reference_tree.txt.gz", "reference_tree.trkb"]:
        p = output_dir / name
        if p.exists():
            return p
    return None


# ── Report formatting ────────────────────────────────────────────────

def print_comparison_table(reports: list[QualityReport]):
    """Print a comparison table of all configs."""
    if not reports:
        return

    print(f"\n{'='*120}")
    print("PARAMETER SWEEP RESULTS")
    print(f"{'='*120}")

    # Header
    header = (
        f"{'Config':<25s} "
        f"{'Trees':>6s} "
        f"{'Leaves':>7s} "
        f"{'Med':>5s} "
        f"{'P25':>5s} "
        f"{'P75':>5s} "
        f"{'Max':>5s} "
        f"{'Small':>5s} "
        f"{'GenPur':>6s} "
        f"{'FamPur':>6s} "
        f"{'GenMono':>7s} "
        f"{'Entropy':>7s} "
        f"{'Size MB':>7s} "
        f"{'Time':>6s}"
    )
    print(header)
    print("-" * 120)

    for r in reports:
        line = (
            f"{r.config_name:<25s} "
            f"{r.num_trees:>6d} "
            f"{r.total_leaves:>7d} "
            f"{r.leaves_per_tree_median:>5.0f} "
            f"{r.leaves_per_tree_p25:>5.0f} "
            f"{r.leaves_per_tree_p75:>5.0f} "
            f"{r.leaves_per_tree_max:>5d} "
            f"{r.small_trees:>5d} "
            f"{r.mean_genus_purity:>6.3f} "
            f"{r.mean_family_purity:>6.3f} "
            f"{r.genus_monophyly_rate:>7.3f} "
            f"{r.mean_genus_entropy:>7.3f} "
            f"{r.reference_tree_size_mb:>7.1f} "
            f"{r.build_time_seconds:>5.0f}s"
        )
        print(line)

    print(f"\n{'='*120}")
    print("METRIC EXPLANATIONS:")
    print("  Trees    = Number of final partitions/trees in the database")
    print("  Leaves   = Total leaf nodes across all trees")
    print("  Med/P25/P75/Max = Leaf count distribution per tree")
    print("  Small    = Trees with <= 3 leaves (limited phylogenetic signal)")
    print("  GenPur   = Mean genus purity (1.0 = each tree is one genus)")
    print("  FamPur   = Mean family purity (1.0 = each tree is one family)")
    print("  GenMono  = Genus monophyly rate (fraction of genera in a single tree)")
    print("  Entropy  = Mean genus Shannon entropy per tree (lower = more pure)")
    print("  Size MB  = reference_tree.txt file size")
    print("  Time     = Build time in seconds")
    print()
    print("INTERPRETATION:")
    print("  Higher GenPur/FamPur  → trees are taxonomically coherent (good for LCA)")
    print("  Higher GenMono        → genera aren't split across trees (good for recall)")
    print("  Lower Entropy         → trees are more genus-uniform")
    print("  Fewer Small trees     → more phylogenetic signal per partition")
    print("  Trade-off: more trees → higher purity but lower monophyly (genera get split)")
    print(f"{'='*120}\n")


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Sweep tronko-build parameters and compute quality metrics"
    )
    parser.add_argument(
        "--input-dir", "-i", required=True,
        help="Directory containing pre-built clusters (*_MSA.fasta, *_taxonomy.txt, RAxML_bestTree.*.reroot)"
    )
    parser.add_argument(
        "--output-dir", "-o", required=True,
        help="Base output directory (each config gets a subdirectory)"
    )
    parser.add_argument(
        "--tronko-build", "-b", default=None,
        help="Path to tronko-build binary (default: auto-detect)"
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=0,
        help="Number of threads (0 = auto-detect, default: 0)"
    )
    parser.add_argument(
        "--configs", "-c", nargs="*", default=None,
        help="Run only these configs (by name). Default: all"
    )
    parser.add_argument(
        "--metrics-only", "-m", action="store_true",
        help="Skip builds, only compute metrics on existing outputs"
    )
    parser.add_argument(
        "--json", "-j", default=None,
        help="Write results as JSON to this file"
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_base = Path(args.output_dir)
    output_base.mkdir(parents=True, exist_ok=True)

    # Find tronko-build binary
    tronko_bin = args.tronko_build
    if tronko_bin is None:
        for candidate in [
            os.path.expanduser("~/.local/bin/tronko-build"),
            shutil.which("tronko-build"),
        ]:
            if candidate and os.path.isfile(candidate):
                tronko_bin = candidate
                break
    if tronko_bin is None and not args.metrics_only:
        print("ERROR: tronko-build binary not found. Use --tronko-build to specify.")
        sys.exit(1)

    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)

    # Count input clusters
    num_clusters = count_clusters(input_dir)
    print(f"Input directory: {input_dir}")
    print(f"Number of input clusters: {num_clusters}")
    print(f"tronko-build: {tronko_bin}")
    print(f"Threads: {threads}")

    # Filter configs if requested
    configs = SWEEP_CONFIGS
    if args.configs:
        configs = [c for c in configs if c.name in args.configs]
        if not configs:
            print(f"ERROR: No matching configs found. Available: {[c.name for c in SWEEP_CONFIGS]}")
            sys.exit(1)

    print(f"\nConfigs to run ({len(configs)}):")
    for c in configs:
        print(f"  {c.name}: {c.short_desc()}, fasttree={c.use_fasttree}, missing_data={c.missing_data}")
        print(f"    flags: {' '.join(c.tronko_flags())}")

    # ── Run builds ────────────────────────────────────────────────
    reports = []
    for config in configs:
        config_dir = output_base / config.name
        ref_tree = find_reference_tree(config_dir)

        if args.metrics_only and ref_tree:
            print(f"\n[metrics-only] {config.name}: using existing {ref_tree.name}")
            build_time = 0.0
        elif args.metrics_only:
            print(f"\n[metrics-only] {config.name}: no reference tree found, skipping")
            continue
        else:
            build_time = run_tronko_build(
                config, input_dir, num_clusters, config_dir, tronko_bin, threads,
            )
            if build_time < 0:
                print(f"  Skipping metrics for {config.name} (build failed)")
                continue
            ref_tree = find_reference_tree(config_dir)

        if ref_tree is None:
            print(f"  WARNING: No reference_tree.txt found in {config_dir}")
            continue

        # Compute quality metrics
        print(f"  Computing quality metrics for {config.name}...")
        report = compute_quality_metrics(ref_tree, config.name, build_time)
        reports.append(report)
        print(f"  {report.num_trees} trees, {report.total_leaves} leaves, "
              f"genus_purity={report.mean_genus_purity:.3f}, "
              f"genus_mono={report.genus_monophyly_rate:.3f}")

    # ── Output ────────────────────────────────────────────────────
    print_comparison_table(reports)

    if args.json:
        json_path = Path(args.json)
        with open(json_path, "w") as f:
            json.dump([asdict(r) for r in reports], f, indent=2)
        print(f"Results written to {json_path}")

    # Save summary alongside outputs
    summary_path = output_base / "sweep_summary.json"
    with open(summary_path, "w") as f:
        json.dump([asdict(r) for r in reports], f, indent=2)
    print(f"Summary saved to {summary_path}")


if __name__ == "__main__":
    main()
