#!/usr/bin/env python3
"""Remap a finished tronko database onto a different taxonomy, without rebuilding.

Everything expensive in a tronko build is taxonomy-independent:

  * AncestralClust receives no taxonomy at all (build-tronko-db.sh passes only
    the FASTA) -- clustering is purely sequence-based.
  * FAMSA alignment and tree inference are sequence-only.
  * tronko-build's SP partitioning reads only the alignment matrix
    (calculateSPArr in tronko-build.c touches m->msa / numspec / numbase and
    nothing else).
  * The per-node posteriors are nucleotide likelihoods -- sequence-derived.

Taxonomy enters only as *labels*: a contiguous block of leaf lineages near the
top of reference_tree.txt, plus the per-partition and marker taxonomy files.
So an LCA and a species database built from the same FASTA differ only in those
labels, and one can be produced from the other by rewriting them -- minutes
instead of the hours (or days) that re-running Step 4 costs.

reference_tree.txt layout (writer: printtree.c printTreeFile,
reader: readreference.c) -- verified empirically against a real database:

    line 1            numberOfTrees
    line 2            max_nodename        longest leaf/accession name
    line 3            max_tax_name        longest single rank field
    line 4            max_lineTaxonomy    longest full lineage line
    next N lines      numbase \t root \t numspec   (one per tree)
    taxonomy block    sum(numspec) lines, each 7 ranks joined by ';'
    remainder         posteriors -- the bulk of the file, left untouched

Ordering rule for the taxonomy block (this is the crux, and is verified by
--verify): trees appear in final_partitions.txt order, and within a tree the
leaves appear in partition{N}_MSA.fasta record order -- NOT the order of
partition{N}_taxonomy.txt. Lineages are written reversed (species first,
domain last), the opposite of the rCRUX input files.

Note on the species variant: rCRUX species_taxonomy.txt carries multiple rows
per accession (one per taxonomic path). taxonomyArr[i][j] holds exactly one
lineage per leaf, so the format cannot represent multiple paths -- one must be
chosen. This script takes the first row per accession and reports how many
accessions had alternatives, rather than silently discarding them.

Usage:
    remap-tronko-taxonomy.py <db_dir> <new_taxonomy.tsv> [--out DIR] [--verify]

    --verify  rebuild the taxonomy block from the database's OWN current
              taxonomy and diff it against what is already in the file. A clean
              run proves the ordering rule holds for this database before any
              remap is trusted. Writes nothing.
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys


def log(msg):
    print(msg, flush=True)


def load_taxonomy(path):
    """accession -> lineage. First row wins; report accessions with alternates."""
    tax = {}
    dupes = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc, lin = parts[0], parts[1]
            if acc in tax:
                if tax[acc] != lin:
                    dupes += 1
                continue
            tax[acc] = lin
    return tax, dupes


def msa_accessions(db, part):
    """Leaf order for a partition == its MSA record order."""
    accs = []
    with open(os.path.join(db, f"partition{part}_MSA.fasta")) as f:
        for line in f:
            if line.startswith(">"):
                accs.append(line[1:].strip().split()[0])
    return accs


def read_partitions(db):
    with open(os.path.join(db, "final_partitions.txt")) as f:
        return [l.strip() for l in f if l.strip()]


def open_tree(db):
    """reference_tree.txt, plain or gzipped."""
    plain = os.path.join(db, "reference_tree.txt")
    gz = plain + ".gz"
    if os.path.exists(plain):
        return open(plain, "rt"), plain, False
    if os.path.exists(gz):
        return gzip.open(gz, "rt"), gz, True
    sys.exit(f"ERROR: no reference_tree.txt(.gz) in {db}")


def build_block(db, parts, tax, missing_limit=20):
    """The new taxonomy block, in tree order then MSA order, reversed."""
    lines = []
    missing = []
    for p in parts:
        for acc in msa_accessions(db, p):
            lin = tax.get(acc)
            if lin is None:
                missing.append(acc)
                if len(missing) <= missing_limit:
                    pass
                lin = ";".join(["NA"] * 7)
            lines.append(";".join(reversed(lin.split(";"))))
    return lines, missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    ap.add_argument("new_taxonomy", nargs="?")
    ap.add_argument("--out", help="write remapped DB here (default: in place)")
    ap.add_argument("--verify", action="store_true",
                    help="check the ordering rule against the existing file; write nothing")
    args = ap.parse_args()

    db = os.path.abspath(args.db)
    parts = read_partitions(db)
    log(f"partitions (trees): {len(parts)}")

    # ---- verify mode: reconstruct from the DB's own taxonomy ----------------
    if args.verify:
        tax = {}
        for p in parts:
            with open(os.path.join(db, f"partition{p}_taxonomy.txt")) as f:
                for line in f:
                    a, l = line.rstrip("\n").split("\t", 1)
                    tax.setdefault(a, l)
        expected, missing = build_block(db, parts, tax)
        fh, path, _ = open_tree(db)
        with fh:
            n_trees = int(fh.readline())
            fh.readline(); fh.readline(); fh.readline()
            for _ in range(n_trees):
                fh.readline()
            bad = 0
            for i, want in enumerate(expected):
                got = fh.readline().rstrip("\n")
                if got != want:
                    bad += 1
                    if bad <= 5:
                        log(f"  line {i}: expected {want!r} got {got!r}")
        log(f"leaves checked : {len(expected)}")
        log(f"missing lookups: {len(missing)}")
        log("VERIFY: PASS — ordering rule holds" if bad == 0
            else f"VERIFY: FAIL — {bad} mismatched leaves")
        return 0 if bad == 0 else 1

    if not args.new_taxonomy:
        sys.exit("ERROR: new_taxonomy is required unless --verify")

    tax, dupes = load_taxonomy(args.new_taxonomy)
    log(f"taxonomy entries  : {len(tax)}")
    if dupes:
        log(f"NOTE: {dupes} extra rows for accessions that already had one "
            f"(multi-path species taxonomy); first row kept per accession.")

    out = os.path.abspath(args.out) if args.out else db
    if out != db:
        log(f"copying database -> {out}")
        shutil.copytree(db, out, symlinks=True, dirs_exist_ok=True)

    block, missing = build_block(out, parts, tax)
    log(f"leaves            : {len(block)}")
    if missing:
        log(f"WARNING: {len(missing)} leaves absent from the new taxonomy "
            f"(written as NA), e.g. {missing[:3]}")

    max_tax_name = max((len(f) for l in block for f in l.split(";")), default=0)
    max_line = max((len(l) for l in block), default=0)
    log(f"max_tax_name      : {max_tax_name}")
    log(f"max_lineTaxonomy  : {max_line}")

    # ---- rewrite reference_tree: header + triplets + new block + posteriors --
    fh, path, was_gz = open_tree(out)
    tmp = os.path.join(out, "reference_tree.txt.remap.tmp")
    written = 0
    with fh, open(tmp, "wt") as w:
        n_trees = int(fh.readline())
        max_nodename = fh.readline().rstrip("\n")
        fh.readline()  # old max_tax_name
        fh.readline()  # old max_lineTaxonomy
        w.write(f"{n_trees}\n{max_nodename}\n{max_tax_name}\n{max_line}\n")
        n_leaves = 0
        for _ in range(n_trees):
            line = fh.readline()
            w.write(line)
            n_leaves += int(line.split("\t")[2])
        if n_leaves != len(block):
            sys.exit(f"ERROR: tree header says {n_leaves} leaves, "
                     f"MSA files give {len(block)} — refusing to write")
        for l in block:
            w.write(l + "\n")
        for _ in range(n_leaves):      # discard the old block
            fh.readline()
        shutil.copyfileobj(fh, w, 1024 * 1024)   # posteriors, untouched
        written = w.tell()
    log(f"reference_tree written: {written/1e9:.2f} GB")

    if was_gz:
        log("recompressing reference_tree.txt.gz ...")
        with open(tmp, "rb") as src, gzip.open(path, "wb", compresslevel=6) as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        os.remove(tmp)
    else:
        os.replace(tmp, path)

    # ---- per-partition taxonomy files --------------------------------------
    # Rewrite lineages in place, preserving each file's existing accession
    # order. tronko-build writes these in its own order, which is neither
    # sorted nor MSA order; regenerating them in MSA order would leave a
    # remapped database gratuitously different from a natively built one
    # (same content, different line order) and make verification harder.
    # Only the reference_tree taxonomy block follows MSA order.
    log("rewriting partition taxonomy files (preserving line order) ...")
    NA7 = ";".join(["NA"] * 7)
    for p in parts:
        path = os.path.join(out, f"partition{p}_taxonomy.txt")
        with open(path) as f:
            accs = [l.rstrip("\n").split("\t", 1)[0] for l in f if l.strip()]
        with open(path, "w") as f:
            for acc in accs:
                f.write(f"{acc}\t{tax.get(acc, NA7)}\n")

    # ---- marker taxonomy (concat in final_partitions order) ----------------
    primer = None
    for f in os.listdir(out):
        if f.endswith("_taxonomy.txt") and not f.startswith(("partition", "input", "marker")):
            primer = f[: -len("_taxonomy.txt")]
            break
    if primer:
        log(f"rewriting {primer}_taxonomy.txt ...")
        with open(os.path.join(out, f"{primer}_taxonomy.txt"), "w") as w:
            for p in parts:
                with open(os.path.join(out, f"partition{p}_taxonomy.txt")) as r:
                    shutil.copyfileobj(r, w)

    # ---- input_taxonomy.txt (harness "is this built" marker) ---------------
    shutil.copyfile(args.new_taxonomy, os.path.join(out, "input_taxonomy.txt"))

    # ---- .trkb ------------------------------------------------------------
    trkb = os.path.join(out, "reference_tree.trkb")
    plain = os.path.join(out, "reference_tree.txt")
    cleanup = False
    if not os.path.exists(plain):
        log("decompressing for tronko-convert ...")
        with gzip.open(plain + ".gz", "rb") as src, open(plain, "wb") as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        cleanup = True
    log("running tronko-convert ...")
    rc = subprocess.call(["tronko-convert", "-i", plain, "-o", trkb])
    if cleanup:
        os.remove(plain)
    if rc != 0:
        sys.exit(f"ERROR: tronko-convert failed rc={rc}")

    log("")
    log("DONE. Sequences, alignments, trees, partitions, FASTA and BWA index")
    log("are unchanged — only taxonomy labels were rewritten.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
