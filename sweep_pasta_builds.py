#!/usr/bin/env python3
"""
Sweep PASTA database build parameters with caching at each pipeline stage.

Pipeline stages and their cache keys:
  Stage 1: PASTA alignment+tree  → key: (input_fasta, iter_limit, max_subproblem_size, num_cpus)
  Stage 2: Tree decomposition    → key: (pasta_tree, max_size, max_diam, min_size, strategy)
  Stage 3: Per-partition pipeline → key: (decomposition, fasttree_gamma)
  Stage 4: tronko-build          → key: (partitions, sp_threshold)

Changing a downstream parameter reuses all upstream cached results.

Usage:
    python3 sweep_pasta_builds.py --config sweep_config.yaml --outdir databases/pasta_sweep
    python3 sweep_pasta_builds.py --list  # show what would run without executing
"""

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from itertools import product
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Defaults / config
# ---------------------------------------------------------------------------

DEFAULT_CONFIG = {
    # Input files
    "input_fasta": None,       # required
    "input_taxonomy": None,    # required
    "tronko_build": None,      # path to tronko-build binary

    # PASTA parameters to sweep (Stage 1)
    "pasta_iter_limits": [3],
    "pasta_max_subproblem_sizes": [200],
    "pasta_num_cpus": 4,
    "pasta_treeshrink": [False],

    # Decomposition parameters to sweep (Stage 2)
    "decomp_max_sizes": [1000],
    "decomp_max_diams": [0],        # 0 = no limit
    "decomp_min_size": 3,
    "decomp_strategy": "centroid",

    # Per-partition parameters (Stage 3)
    "fasttree_gamma": [False, True],
    "famsa_threads": 4,
    "partition_workers": 3,

    # tronko-build parameters (Stage 4)
    "sp_thresholds": [0],           # 0 = no SP partitioning
    "tronko_threads": 12,

    # General
    "threads": 12,
}


def load_config(path):
    """Load YAML config, falling back to JSON."""
    if path.endswith(".yaml") or path.endswith(".yml"):
        try:
            import yaml
            with open(path) as f:
                return yaml.safe_load(f)
        except ImportError:
            log.error("PyYAML not installed. Use JSON config or: pip install pyyaml")
            sys.exit(1)
    else:
        with open(path) as f:
            return json.load(f)


def cache_key(*args):
    """Create a short hash from arguments for cache directory naming."""
    h = hashlib.md5(json.dumps(args, sort_keys=True).encode()).hexdigest()[:8]
    return h


def run_cmd(cmd, label="", **kwargs):
    """Run a command, log it, return CompletedProcess."""
    cmd_str = " ".join(str(c) for c in cmd)
    log.info(f"  [{label}] {cmd_str}")
    t0 = time.time()
    result = subprocess.run(cmd, **kwargs)
    elapsed = time.time() - t0
    if result.returncode != 0:
        log.error(f"  [{label}] FAILED (exit {result.returncode}, {elapsed:.1f}s)")
        if hasattr(result, 'stderr') and result.stderr:
            log.error(f"  stderr: {result.stderr[:500]}")
    else:
        log.info(f"  [{label}] done ({elapsed:.1f}s)")
    return result


# ---------------------------------------------------------------------------
# Stage 1: PASTA
# ---------------------------------------------------------------------------

def stage1_pasta(cfg, pasta_dir, iter_limit, max_subproblem_size, treeshrink):
    """Run PASTA. Returns path to final tree."""
    done_marker = os.path.join(pasta_dir, ".stage1_done")
    tree_path = os.path.join(pasta_dir, "pasta_final.tre")

    if os.path.exists(done_marker):
        log.info(f"  Stage 1 CACHED: {pasta_dir}")
        return tree_path

    os.makedirs(pasta_dir, exist_ok=True)

    # Build PASTA command
    cmd = [
        "python3", "-m", "pasta",
        "--input", cfg["input_fasta"],
        "--datatype", "dna",
        "--iter-limit", str(iter_limit),
        "--max-subproblem-size", str(max_subproblem_size),
        "--num-cpus", str(cfg.get("pasta_num_cpus", 4)),
        "--job", "pasta_sweep",
        "--output-directory", pasta_dir,
        "-o", pasta_dir,
    ]
    if treeshrink:
        cmd.append("--treeshrink")

    result = run_cmd(cmd, label="PASTA")
    if result.returncode != 0:
        return None

    # Find the output tree (PASTA names it based on --job)
    for f in os.listdir(pasta_dir):
        if f.endswith(".tre") and "temp" not in f and "topo" not in f:
            # Fix: strip quotes, copy to canonical name
            src = os.path.join(pasta_dir, f)
            with open(src) as fh:
                content = fh.read().replace("'", "")
            with open(tree_path, "w") as fh:
                fh.write(content)
            break
    else:
        log.error("PASTA did not produce a tree file")
        return None

    with open(done_marker, "w") as f:
        f.write(json.dumps({
            "iter_limit": iter_limit,
            "max_subproblem_size": max_subproblem_size,
            "treeshrink": treeshrink,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }, indent=2))

    return tree_path


# ---------------------------------------------------------------------------
# Stage 2: Decomposition
# ---------------------------------------------------------------------------

def stage2_decompose(cfg, decomp_dir, tree_path, max_size, max_diam, min_size, strategy):
    """Decompose PASTA tree into partitions. Returns partition dir."""
    done_marker = os.path.join(decomp_dir, ".stage2_done")

    if os.path.exists(done_marker):
        log.info(f"  Stage 2 CACHED: {decomp_dir}")
        return decomp_dir

    os.makedirs(decomp_dir, exist_ok=True)

    # Use partition_and_build.py with --skip-tronko
    script = os.path.join(os.path.dirname(__file__),
                          "tronko-build", "partition_and_build.py")

    cmd = [
        "python3", script,
        "--tree", tree_path,
        "--fasta", cfg["input_fasta"],
        "--taxonomy", cfg["input_taxonomy"],
        "--outdir", decomp_dir,
        "--min-size", str(min_size),
        "--strategy", strategy,
        "--threads", str(cfg.get("threads", 12)),
        "--famsa-threads", str(cfg.get("famsa_threads", 4)),
        "--workers", str(cfg.get("partition_workers", 3)),
        "--skip-tronko",
    ]
    if max_size > 0:
        cmd.extend(["--max-size", str(max_size)])
    if max_diam > 0:
        cmd.extend(["--max-diam", str(max_diam)])

    result = run_cmd(cmd, label="decompose")
    if result.returncode != 0:
        return None

    with open(done_marker, "w") as f:
        f.write(json.dumps({
            "max_size": max_size,
            "max_diam": max_diam,
            "min_size": min_size,
            "strategy": strategy,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }, indent=2))

    return decomp_dir


# ---------------------------------------------------------------------------
# Stage 3: Per-partition FAMSA + FastTree (already done by partition_and_build.py)
# For gamma reoptimization, we post-process the trees
# ---------------------------------------------------------------------------

def stage3_gamma_reopt(cfg, partitions_dir, gamma_dir, use_gamma):
    """Optionally re-estimate branch lengths with -gamma."""
    done_marker = os.path.join(gamma_dir, ".stage3_done")

    if os.path.exists(done_marker):
        log.info(f"  Stage 3 CACHED: {gamma_dir}")
        return gamma_dir

    os.makedirs(gamma_dir, exist_ok=True)

    if not use_gamma:
        # Just symlink/copy the partition files as-is
        for f in os.listdir(partitions_dir):
            if f.startswith("."):
                continue
            src = os.path.join(partitions_dir, f)
            dst = os.path.join(gamma_dir, f)
            if os.path.isfile(src) and not os.path.exists(dst):
                os.symlink(os.path.abspath(src), dst)
    else:
        # Copy MSA, taxonomy, and reroot files; re-run FastTree -gamma on trees
        for f in os.listdir(partitions_dir):
            if f.startswith("."):
                continue
            src = os.path.join(partitions_dir, f)
            dst = os.path.join(gamma_dir, f)

            if f.endswith("_MSA.fasta") or f.endswith("_taxonomy.txt"):
                if not os.path.exists(dst):
                    os.symlink(os.path.abspath(src), dst)
            elif f.startswith("RAxML_bestTree.") and f.endswith(".reroot"):
                # Re-run FastTree with -gamma to reoptimize branch lengths
                # FastTree -gamma -nt -gtr -nome -mllen -intree <tree> <alignment>
                partition_name = f.replace("RAxML_bestTree.", "").replace(".reroot", "")
                msa = os.path.join(partitions_dir, f"{partition_name}_MSA.fasta")
                if os.path.exists(msa):
                    gamma_cmd = [
                        "FastTree", "-nt", "-gtr", "-gamma",
                        "-nome", "-mllen",
                        "-intree", src,
                        msa,
                    ]
                    result = subprocess.run(
                        gamma_cmd,
                        capture_output=True, text=True
                    )
                    if result.returncode == 0 and result.stdout.strip():
                        with open(dst, "w") as fh:
                            fh.write(result.stdout)
                    else:
                        # Fallback: use original tree
                        os.symlink(os.path.abspath(src), dst)
                else:
                    os.symlink(os.path.abspath(src), dst)

    with open(done_marker, "w") as f:
        f.write(json.dumps({
            "use_gamma": use_gamma,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }, indent=2))

    return gamma_dir


# ---------------------------------------------------------------------------
# Stage 4: tronko-build
# ---------------------------------------------------------------------------

def stage4_tronko_build(cfg, partitions_dir, output_dir, sp_threshold):
    """Run tronko-build on partitions."""
    done_marker = os.path.join(output_dir, ".stage4_done")
    ref_tree = os.path.join(output_dir, "reference_tree.txt")

    if os.path.exists(done_marker) and os.path.exists(ref_tree):
        log.info(f"  Stage 4 CACHED: {output_dir}")
        return output_dir

    os.makedirs(output_dir, exist_ok=True)

    # Count partitions
    n_partitions = len([f for f in os.listdir(partitions_dir)
                        if f.startswith("RAxML_bestTree.") and f.endswith(".reroot")])

    tronko_bin = cfg.get("tronko_build") or "tronko-build"

    cmd = [
        tronko_bin,
        "-y",
        "-e", partitions_dir,
        "-n", str(n_partitions),
        "-d", output_dir,
        "-E",  # export subtrees for ablation
        "-c", str(cfg.get("tronko_threads", 12)),
    ]
    if sp_threshold > 0:
        cmd.extend(["-s", "-u", str(sp_threshold), "-a"])

    result = run_cmd(cmd, label="tronko-build")
    if result.returncode != 0:
        return None

    if not os.path.exists(ref_tree):
        log.error(f"tronko-build did not produce {ref_tree}")
        return None

    # Build marker.fasta and BWA index
    _build_marker_and_index(cfg, partitions_dir, output_dir, n_partitions, sp_threshold)

    with open(done_marker, "w") as f:
        f.write(json.dumps({
            "sp_threshold": sp_threshold,
            "n_partitions": n_partitions,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }, indent=2))

    return output_dir


def _build_marker_and_index(cfg, partitions_dir, output_dir, n_partitions, sp_threshold):
    """Build marker.fasta from partition FASTAs and create BWA index."""
    marker_path = os.path.join(output_dir, "marker.fasta")
    tax_path = os.path.join(output_dir, "marker_taxonomy.txt")

    # If SP partitioning created final_partitions.txt, use that
    final_parts = os.path.join(output_dir, "final_partitions.txt")
    if sp_threshold > 0 and os.path.exists(final_parts):
        # Use final_partitions.txt to pick which partitions to include
        with open(final_parts) as f:
            part_nums = [line.strip() for line in f if line.strip()]
        source_dir = output_dir  # SP-partitioned files are in output_dir
    else:
        # No SP partitioning — use all original partitions
        part_nums = list(range(n_partitions))
        source_dir = partitions_dir

    # Concatenate raw (gap-free) FASTAs
    with open(marker_path, "w") as out:
        for pn in part_nums:
            # Try unaligned first, then strip gaps from MSA
            for pattern in [f"partition{pn}.fasta", f"partition{pn}_unaligned.fasta"]:
                fpath = os.path.join(source_dir, pattern)
                if os.path.exists(fpath):
                    with open(fpath) as fin:
                        for line in fin:
                            if line.startswith(">"):
                                out.write(line)
                            else:
                                out.write(line.replace("-", ""))
                    break
            else:
                # Fall back to MSA with gap stripping
                msa = os.path.join(source_dir, f"partition{pn}_MSA.fasta")
                if os.path.exists(msa):
                    with open(msa) as fin:
                        for line in fin:
                            if line.startswith(">"):
                                out.write(line)
                            else:
                                out.write(line.replace("-", ""))

    # Build taxonomy
    with open(tax_path, "w") as out:
        for pn in part_nums:
            tpath = os.path.join(source_dir, f"partition{pn}_taxonomy.txt")
            if os.path.exists(tpath):
                with open(tpath) as fin:
                    out.write(fin.read())

    # BWA index
    run_cmd(["bwa", "index", marker_path], label="bwa-index")

    # Copy input files for provenance
    for src_name, dst_name in [("input_fasta", "input.fasta"), ("input_taxonomy", "input_taxonomy.txt")]:
        src = cfg.get(src_name)
        dst = os.path.join(output_dir, dst_name)
        if src and os.path.exists(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)


# ---------------------------------------------------------------------------
# Sweep orchestrator
# ---------------------------------------------------------------------------

def generate_sweep_configs(cfg):
    """Generate all parameter combinations."""
    combos = []

    for iter_limit in cfg.get("pasta_iter_limits", [3]):
        for max_sub in cfg.get("pasta_max_subproblem_sizes", [200]):
            for treeshrink in cfg.get("pasta_treeshrink", [False]):
                for max_size in cfg.get("decomp_max_sizes", [1000]):
                    for max_diam in cfg.get("decomp_max_diams", [0]):
                        for gamma in cfg.get("fasttree_gamma", [False]):
                            for sp in cfg.get("sp_thresholds", [0]):
                                combos.append({
                                    "iter_limit": iter_limit,
                                    "max_subproblem_size": max_sub,
                                    "treeshrink": treeshrink,
                                    "max_size": max_size,
                                    "max_diam": max_diam,
                                    "gamma": gamma,
                                    "sp_threshold": sp,
                                })

    return combos


def combo_name(combo):
    """Human-readable name for a parameter combination."""
    parts = [f"pasta_iter{combo['iter_limit']}_sub{combo['max_subproblem_size']}"]
    if combo["treeshrink"]:
        parts.append("treeshrink")
    parts.append(f"maxsize{combo['max_size']}")
    if combo["max_diam"] > 0:
        parts.append(f"maxdiam{combo['max_diam']}")
    if combo["gamma"]:
        parts.append("gamma")
    if combo["sp_threshold"] > 0:
        parts.append(f"sp{combo['sp_threshold']}")
    else:
        parts.append("nofilter")
    return "_".join(parts)


def run_sweep(cfg, outdir, dry_run=False):
    """Run the full parameter sweep with caching."""
    combos = generate_sweep_configs(cfg)
    log.info(f"Sweep: {len(combos)} configurations")

    results = []
    for i, combo in enumerate(combos):
        name = combo_name(combo)
        log.info(f"\n{'='*60}")
        log.info(f"Config {i+1}/{len(combos)}: {name}")
        log.info(f"{'='*60}")

        if dry_run:
            log.info(f"  [DRY RUN] Would build: {name}")
            results.append({"name": name, "status": "dry_run", **combo})
            continue

        # Stage 1: PASTA (cached by iter_limit + max_subproblem_size + treeshrink)
        pasta_key = cache_key(
            cfg["input_fasta"],
            combo["iter_limit"],
            combo["max_subproblem_size"],
            combo["treeshrink"],
        )
        pasta_dir = os.path.join(outdir, "cache", f"pasta_{pasta_key}")
        tree_path = stage1_pasta(
            cfg, pasta_dir,
            combo["iter_limit"],
            combo["max_subproblem_size"],
            combo["treeshrink"],
        )
        if tree_path is None:
            results.append({"name": name, "status": "failed_stage1", **combo})
            continue

        # Stage 2: Decomposition (cached by tree + decomp params)
        decomp_key = cache_key(
            pasta_key,
            combo["max_size"],
            combo["max_diam"],
            cfg.get("decomp_min_size", 3),
            cfg.get("decomp_strategy", "centroid"),
        )
        decomp_dir = os.path.join(outdir, "cache", f"decomp_{decomp_key}")
        partitions_dir = stage2_decompose(
            cfg, decomp_dir, tree_path,
            combo["max_size"],
            combo["max_diam"],
            cfg.get("decomp_min_size", 3),
            cfg.get("decomp_strategy", "centroid"),
        )
        if partitions_dir is None:
            results.append({"name": name, "status": "failed_stage2", **combo})
            continue

        # Stage 3: Gamma reoptimization (cached by decomp + gamma flag)
        gamma_key = cache_key(decomp_key, combo["gamma"])
        gamma_dir = os.path.join(outdir, "cache", f"gamma_{gamma_key}")
        final_partitions = stage3_gamma_reopt(
            cfg, partitions_dir, gamma_dir, combo["gamma"],
        )
        if final_partitions is None:
            results.append({"name": name, "status": "failed_stage3", **combo})
            continue

        # Stage 4: tronko-build
        db_dir = os.path.join(outdir, "dbs", name)
        result_dir = stage4_tronko_build(
            cfg, final_partitions, db_dir, combo["sp_threshold"],
        )
        if result_dir is None:
            results.append({"name": name, "status": "failed_stage4", **combo})
            continue

        # Record result
        ref_tree = os.path.join(db_dir, "reference_tree.txt")
        size_mb = os.path.getsize(ref_tree) / 1e6 if os.path.exists(ref_tree) else 0
        n_trees = 0
        if os.path.exists(ref_tree):
            with open(ref_tree) as f:
                n_trees = int(f.readline().strip())

        results.append({
            "name": name,
            "status": "success",
            "db_dir": db_dir,
            "size_mb": round(size_mb, 1),
            "n_trees": n_trees,
            **combo,
        })
        log.info(f"  SUCCESS: {name} → {n_trees} trees, {size_mb:.1f} MB")

    # Write results summary
    os.makedirs(outdir, exist_ok=True)
    summary_path = os.path.join(outdir, "sweep_results.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    log.info(f"\nSweep complete. Results: {summary_path}")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Sweep PASTA database build parameters")
    parser.add_argument("--config", help="Config file (YAML or JSON)")
    parser.add_argument("--outdir", default="databases/pasta_sweep", help="Output directory")
    parser.add_argument("--list", action="store_true", help="List configs without running")

    # Quick CLI overrides (avoid needing a config file for simple sweeps)
    parser.add_argument("--fasta", help="Input FASTA")
    parser.add_argument("--taxonomy", help="Input taxonomy")
    parser.add_argument("--tronko-build", help="Path to tronko-build binary")
    parser.add_argument("--pasta-iters", help="Comma-separated iter limits (e.g., 3,5,10)")
    parser.add_argument("--max-sizes", help="Comma-separated max sizes (e.g., 500,1000,2000)")
    parser.add_argument("--threads", type=int, default=12)

    args = parser.parse_args()

    # Load config
    if args.config:
        cfg = {**DEFAULT_CONFIG, **load_config(args.config)}
    else:
        cfg = dict(DEFAULT_CONFIG)

    # Apply CLI overrides
    if args.fasta:
        cfg["input_fasta"] = args.fasta
    if args.taxonomy:
        cfg["input_taxonomy"] = args.taxonomy
    if args.tronko_build:
        cfg["tronko_build"] = args.tronko_build
    if args.pasta_iters:
        cfg["pasta_iter_limits"] = [int(x) for x in args.pasta_iters.split(",")]
    if args.max_sizes:
        cfg["decomp_max_sizes"] = [int(x) for x in args.max_sizes.split(",")]
    if args.threads:
        cfg["threads"] = args.threads
        cfg["pasta_num_cpus"] = args.threads

    # Validate
    if not cfg.get("input_fasta"):
        parser.error("--fasta or config with input_fasta required")
    if not cfg.get("input_taxonomy"):
        parser.error("--taxonomy or config with input_taxonomy required")

    run_sweep(cfg, args.outdir, dry_run=args.list)


if __name__ == "__main__":
    main()
