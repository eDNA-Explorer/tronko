#!/usr/bin/env python3
"""Patch reference_tree.txt to fix ;;;;;;  taxonomy entries caused by
colon-to-underscore mismatch between Newick leaf names and taxonomy accessions.

Instead of a full tronko-build rebuild, this:
1. Extracts leaf node names from the node data section
2. Looks up correct taxonomy from master taxonomy file (with colon normalization)
3. Rewrites the taxonomy section in-place
"""
import sys
import os
import re
import subprocess

def main():
    if len(sys.argv) < 3:
        print("Usage: patch_taxonomy.py <reference_tree.txt> <master_taxonomy.txt>")
        sys.exit(1)

    db_path = sys.argv[1]
    tax_path = sys.argv[2]

    # Step 1: Read master taxonomy file
    # Format: accession<TAB>domain;phylum;class;order;family;genus;species
    print("Reading master taxonomy file...")
    master_tax = {}  # normalized_name -> taxonomy_string
    with open(tax_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            parts = line.split('\t', 1)
            if len(parts) != 2:
                continue
            accession, taxonomy = parts
            # Store with original name
            master_tax[accession] = taxonomy
            # Also store with colons replaced by underscores
            if ':' in accession:
                normalized = accession.replace(':', '_')
                master_tax[normalized] = taxonomy

    print(f"  Loaded {len(master_tax)} taxonomy entries (including normalized)")

    # Step 2: Parse database header
    print("Parsing database header...")
    with open(db_path, 'r') as f:
        ntrees = int(f.readline().strip())
        max_nodename = int(f.readline().strip())
        max_tax_name = int(f.readline().strip())
        max_linetax = int(f.readline().strip())
        header_lines = 4

        numbase = []
        root = []
        numspec = []
        for i in range(ntrees):
            parts = f.readline().strip().split('\t')
            numbase.append(int(parts[0]))
            root.append(int(parts[1]))
            numspec.append(int(parts[2]))

    meta_lines = ntrees
    tax_start = header_lines + meta_lines  # 0-indexed line number where taxonomy begins
    total_tax = sum(numspec)
    tax_end = tax_start + total_tax  # 0-indexed line number where node data begins

    print(f"  {ntrees} trees, {total_tax} taxonomy entries")
    print(f"  Taxonomy section: lines {tax_start+1} to {tax_end} (1-indexed)")

    # Step 3: Extract leaf node lines using grep (fast)
    print("Extracting leaf node names from node data section...")
    result = subprocess.run(
        ['grep', '-E', r'^[0-9]+\t[0-9]+\t-1\t-1\t', db_path],
        capture_output=True, text=True
    )
    leaf_lines = result.stdout.strip().split('\n')
    print(f"  Found {len(leaf_lines)} leaf nodes")

    # Build mapping: (treeNum, taxIdx0) -> leafName
    leaf_map = {}
    for line in leaf_lines:
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        tree_num = int(fields[0])
        tax_idx0 = int(fields[6])
        leaf_name = fields[8]
        leaf_map[(tree_num, tax_idx0)] = leaf_name

    print(f"  Mapped {len(leaf_map)} (tree, speciesIdx) -> leafName entries")

    # Step 4: Read the taxonomy section, identify ;;;;;;  lines, look up corrections
    print("Scanning taxonomy section for ;;;;;;  entries...")
    corrections = {}  # line_number (0-indexed) -> corrected_taxonomy_string
    empty_tax = ';;;;;;'
    missed = 0
    fixed = 0

    line_idx = tax_start
    for tree_i in range(ntrees):
        for spec_j in range(numspec[tree_i]):
            # This taxonomy line corresponds to (tree_i, spec_j)
            leaf_name = leaf_map.get((tree_i, spec_j))
            if leaf_name is not None:
                tax = master_tax.get(leaf_name)
                if tax is not None:
                    corrections[line_idx] = tax
                    fixed += 1
                else:
                    missed += 1
            line_idx += 1

    print(f"  Can fix {fixed} entries, {missed} still unresolvable")

    # Step 5: Rewrite the file with corrections
    print("Rewriting database with corrected taxonomy...")
    tmp_path = db_path + '.patched'
    lines_written = 0
    tax_fixed = 0

    with open(db_path, 'r') as fin, open(tmp_path, 'w') as fout:
        for line_num, line in enumerate(fin):
            if line_num in corrections and line.rstrip('\n') == empty_tax:
                fout.write(corrections[line_num] + '\n')
                tax_fixed += 1
            else:
                fout.write(line)
            lines_written += 1
            if lines_written % 5000000 == 0:
                print(f"  Processed {lines_written} lines...")

    print(f"  Fixed {tax_fixed} taxonomy entries")
    print(f"  Total lines written: {lines_written}")

    # Step 6: Replace original with patched
    backup_path = db_path + '.bak'
    os.rename(db_path, backup_path)
    os.rename(tmp_path, db_path)
    print(f"  Original backed up to {backup_path}")
    print("Done!")

if __name__ == '__main__':
    main()
