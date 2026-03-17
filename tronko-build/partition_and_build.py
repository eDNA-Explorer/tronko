#!/usr/bin/env python3
"""
Partition a PASTA tree into tronko-build-compatible partitions,
then run FAMSA + FastTree per partition and invoke tronko-build.
"""

import argparse
import os
import sys
import subprocess
import time
import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_fasta(path):
    """Load FASTA into dict {name: sequence}."""
    d = {}
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    d[name] = ''.join(parts)
                name = line[1:].rstrip()
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        d[name] = ''.join(parts)
    return d


def load_taxonomy(path):
    """Load taxonomy into dict {name: taxonomy_string}."""
    d = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if '\t' in line:
                name, tax = line.split('\t', 1)
                d[name] = tax
    return d


def get_group_leaves(group_root):
    """Collect leaf taxon labels belonging to this group.

    Traverse descendants of group_root, stopping at children
    that are marked (they belong to a different group).
    """
    leaves = []
    stack = [group_root]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            leaves.append(node.taxon.label)
        else:
            for child in node.child_nodes():
                if not child.marked:
                    stack.append(child)
    return leaves


def unwrap_fasta(path):
    """Convert wrapped FASTA to single-line-per-sequence in-place."""
    with open(path) as f:
        lines = f.readlines()
    with open(path, 'w') as f:
        seq_parts = []
        for line in lines:
            if line.startswith('>'):
                if seq_parts:
                    f.write(''.join(seq_parts) + '\n')
                    seq_parts = []
                f.write(line)
            else:
                seq_parts.append(line.rstrip('\n'))
        if seq_parts:
            f.write(''.join(seq_parts) + '\n')


def strip_quotes_from_newick(path):
    """Remove single quotes from Newick leaf names."""
    with open(path) as f:
        content = f.read()
    if "'" in content:
        content = content.replace("'", "")
        with open(path, 'w') as f:
            f.write(content)


def count_newick_leaves(path):
    """Count leaf nodes in a Newick file."""
    with open(path) as f:
        tree = f.read()
    # Leaves are tokens between (/, and : that aren't empty
    return len(re.findall(r'[\(,]([^:,\(\)]+?):', tree))


# ---------------------------------------------------------------------------
# Per-partition pipeline (runs in subprocess pool)
# ---------------------------------------------------------------------------

def process_partition(outdir, partition_idx, expected_count, famsa_threads):
    """Run FAMSA, unwrap, FastTree, nw_reroot for one partition.

    Returns (partition_idx, success, elapsed_time, error_msg).
    """
    name = f"partition{partition_idx}"
    unaligned = os.path.join(outdir, f"{name}_unaligned.fasta")
    msa_path = os.path.join(outdir, f"{name}_MSA.fasta")
    tree_raw = os.path.join(outdir, f"{name}_fasttree_raw.tre")
    tree_final = os.path.join(outdir, f"RAxML_bestTree.{name}.reroot")

    t0 = time.time()

    # Step 1: FAMSA alignment
    result = subprocess.run(
        ["famsa", "-t", str(famsa_threads), unaligned, msa_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return (partition_idx, False, 0, f"FAMSA failed: {result.stderr[:500]}")

    # Step 2: Unwrap FASTA to single-line-per-sequence
    unwrap_fasta(msa_path)

    # Step 3: FastTree (flags match sweep: -nt -gtr only)
    with open(tree_raw, 'w') as tree_out:
        result = subprocess.run(
            ["FastTree", "-nt", "-gtr", msa_path],
            stdout=tree_out, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        return (partition_idx, False, 0, f"FastTree failed: {result.stderr[:500]}")

    # Step 4: Strip quotes and support values from FastTree output
    with open(tree_raw) as f:
        tree_data = f.read()
    tree_data = tree_data.replace("'", "")
    tree_data = re.sub(r'\)([0-9][0-9.eE+-]*):', '):', tree_data)
    with open(tree_raw, 'w') as f:
        f.write(tree_data)

    # Step 5: nw_reroot (midpoint)
    with open(tree_final, 'w') as reroot_out:
        result = subprocess.run(
            ["nw_reroot", tree_raw],
            stdout=reroot_out, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        return (partition_idx, False, 0, f"nw_reroot failed: {result.stderr[:500]}")

    # Validate rerooted tree has correct leaf count
    leaf_count = count_newick_leaves(tree_final)
    if leaf_count != expected_count:
        return (partition_idx, False, 0,
                f"nw_reroot leaf count mismatch: expected {expected_count}, got {leaf_count}")

    # Clean up intermediates
    os.unlink(unaligned)
    os.unlink(tree_raw)

    elapsed = time.time() - t0
    return (partition_idx, True, elapsed, None)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Partition a PASTA tree and build tronko-build database")
    parser.add_argument("--tree", required=True, help="PASTA Newick tree")
    parser.add_argument("--fasta", required=True, help="Unaligned FASTA")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy file")
    parser.add_argument("--outdir", required=True, help="Partition output directory")
    parser.add_argument("--max-size", type=int, default=0,
                        help="Max leaves per partition (default: 0 = no limit)")
    parser.add_argument("--max-diam", type=float, default=0,
                        help="Max tree diameter per partition (default: 0 = no limit)")
    parser.add_argument("--min-size", type=int, default=3,
                        help="Min leaves per partition (default: 3)")
    parser.add_argument("--strategy", default="centroid",
                        choices=["centroid", "midpoint"],
                        help="Decomposition strategy (default: centroid)")
    parser.add_argument("--sp-threshold", type=float, default=0,
                        help="SP score threshold for tronko-build partitioning (default: 0 = no SP partitioning)")
    parser.add_argument("--threads", type=int, default=12,
                        help="Total threads (default: 12)")
    parser.add_argument("--famsa-threads", type=int, default=4,
                        help="Threads per FAMSA call (default: 4)")
    parser.add_argument("--workers", type=int, default=0,
                        help="Parallel partitions (default: threads // famsa_threads)")
    parser.add_argument("--tronko-build", default=None,
                        help="Path to tronko-build binary")
    parser.add_argument("--tronko-outdir", default=None,
                        help="Output dir for tronko-build -d")
    parser.add_argument("--skip-tronko", action="store_true",
                        help="Only partition, don't run tronko-build")
    args = parser.parse_args()

    if args.workers <= 0:
        args.workers = max(1, args.threads // args.famsa_threads)
    if args.tronko_outdir is None:
        args.tronko_outdir = os.path.join(args.outdir, "tronko_db")

    os.makedirs(args.outdir, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(os.path.join(args.outdir, "partition_and_build.log"))
        ]
    )

    # ------------------------------------------------------------------
    # 1. Load inputs
    # ------------------------------------------------------------------
    logging.info("Loading tree...")
    import dendropy
    tree = dendropy.Tree.get(path=args.tree, schema="newick",
                             preserve_underscores=True)
    total_leaves = len(tree.leaf_nodes())
    logging.info(f"Tree loaded: {total_leaves} leaves")

    logging.info("Loading FASTA...")
    fasta_dict = load_fasta(args.fasta)
    logging.info(f"FASTA loaded: {len(fasta_dict)} sequences")

    logging.info("Loading taxonomy...")
    taxonomy_dict = load_taxonomy(args.taxonomy)
    logging.info(f"Taxonomy loaded: {len(taxonomy_dict)} entries")

    # ------------------------------------------------------------------
    # 2. Decompose tree
    # ------------------------------------------------------------------
    decomp_params = f"strategy={args.strategy}, min_size={args.min_size}"
    if args.max_diam > 0:
        decomp_params += f", max_diam={args.max_diam}"
    if args.max_size > 0:
        decomp_params += f", max_size={args.max_size}"
    logging.info(f"Decomposing tree ({decomp_params})...")

    # Import PASTA's decompose function
    pasta_dir = os.environ.get("PASTA_DIR", os.path.expanduser("~/pasta"))
    sys.path.insert(0, pasta_dir)
    from pasta.decompose_lib import decompose_by_diameter

    decomp_kwargs = dict(strategy=args.strategy, min_size=args.min_size)
    if args.max_size > 0:
        decomp_kwargs['max_size'] = args.max_size
    if args.max_diam > 0:
        decomp_kwargs['max_diam'] = args.max_diam
    tree_map = decompose_by_diameter(tree, **decomp_kwargs)

    # Extract leaves per group
    # decompose_by_diameter returns a dict {name: node} normally, but returns
    # a list of (node, name) tuples when no splitting occurs
    partitions = []
    total_assigned = 0
    if isinstance(tree_map, dict):
        items = [(name, node) for name, node in sorted(tree_map.items())]
    else:
        # List of (node, name) tuples
        items = [(name, node) for node, name in tree_map]
    for group_name, group_root in items:
        leaves = get_group_leaves(group_root)
        partitions.append(leaves)
        total_assigned += len(leaves)

    assert total_assigned == total_leaves, \
        f"Partition sum ({total_assigned}) != total leaves ({total_leaves})"

    sizes = [len(p) for p in partitions]
    logging.info(f"Decomposed into {len(partitions)} partitions")
    logging.info(f"  sizes: min={min(sizes)}, max={max(sizes)}, "
                 f"mean={sum(sizes)/len(sizes):.0f}, "
                 f"median={sorted(sizes)[len(sizes)//2]}")

    # ------------------------------------------------------------------
    # 3. Write partition files
    # ------------------------------------------------------------------
    logging.info("Writing partition files...")
    partition_info = []  # (index, leaf_count)
    missing_fasta = 0
    missing_tax = 0

    for i, leaves in enumerate(partitions):
        name = f"partition{i}"

        # Unaligned FASTA (temporary, input to FAMSA)
        fasta_path = os.path.join(args.outdir, f"{name}_unaligned.fasta")
        with open(fasta_path, 'w') as f:
            for leaf in leaves:
                seq = fasta_dict.get(leaf)
                if seq is None:
                    missing_fasta += 1
                    continue
                f.write(f">{leaf}\n{seq}\n")

        # Taxonomy
        tax_path = os.path.join(args.outdir, f"{name}_taxonomy.txt")
        with open(tax_path, 'w') as f:
            for leaf in leaves:
                tax = taxonomy_dict.get(leaf)
                if tax is None:
                    missing_tax += 1
                    continue
                f.write(f"{leaf}\t{tax}\n")

        partition_info.append((i, len(leaves)))

    if missing_fasta:
        logging.warning(f"{missing_fasta} leaves missing from FASTA")
    if missing_tax:
        logging.warning(f"{missing_tax} leaves missing from taxonomy")

    logging.info(f"Written {len(partition_info)} partition file sets")

    # ------------------------------------------------------------------
    # 4. Run FAMSA + FastTree + nw_reroot in parallel
    # ------------------------------------------------------------------
    logging.info(f"Processing partitions ({args.workers} workers, "
                 f"{args.famsa_threads} FAMSA threads each)...")

    failed = []
    completed = 0
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_partition, args.outdir, idx, count,
                            args.famsa_threads): idx
            for idx, count in partition_info
        }
        for future in as_completed(futures):
            idx, success, elapsed, err = future.result()
            completed += 1
            if success:
                logging.info(f"  Partition {idx} done ({elapsed:.1f}s) "
                             f"[{completed}/{len(partition_info)}]")
            else:
                logging.error(f"  Partition {idx} FAILED: {err}")
                failed.append(idx)

    if failed:
        logging.error(f"{len(failed)} partitions failed: {failed}")
        sys.exit(1)

    logging.info(f"All {len(partition_info)} partitions processed successfully")

    # ------------------------------------------------------------------
    # 5. Run tronko-build
    # ------------------------------------------------------------------
    if args.skip_tronko:
        logging.info("Skipping tronko-build (--skip-tronko)")
        return

    os.makedirs(args.tronko_outdir, exist_ok=True)

    tronko_bin = args.tronko_build
    if tronko_bin is None:
        # Try to find it relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        candidate = os.path.join(script_dir, "..", "..", "tronko-build")
        if os.path.isfile(candidate):
            tronko_bin = candidate
        else:
            logging.error("tronko-build binary not found. Use --tronko-build.")
            sys.exit(1)

    cmd = [
        tronko_bin,
        "-y",
        "-e", args.outdir,
        "-n", str(len(partition_info)),
        "-d", args.tronko_outdir
    ]
    if args.sp_threshold > 0:
        cmd.extend(["-s", "-u", str(args.sp_threshold), "-a"])
    cmd.append("-E")
    logging.info(f"Running tronko-build: {' '.join(cmd)}")
    t0 = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - t0

    if result.returncode != 0:
        logging.error(f"tronko-build failed with exit code {result.returncode}")
        sys.exit(1)

    ref_tree = os.path.join(args.tronko_outdir, "reference_tree.txt")
    if os.path.isfile(ref_tree):
        size_mb = os.path.getsize(ref_tree) / 1e6
        logging.info(f"tronko-build completed in {elapsed:.0f}s. "
                     f"reference_tree.txt: {size_mb:.1f} MB")
    else:
        logging.error("tronko-build completed but reference_tree.txt not found")
        sys.exit(1)


if __name__ == "__main__":
    main()
