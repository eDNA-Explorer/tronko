#!/bin/bash
# Detect "orphaned top-level clusters" in a finished tronko database.
#
# BUG: tronko-build's createNewRoots() aborts a split when any of the 3 candidate
# sub-partitions would hold < 4 sequences (tronko-build.c:1069). It marks the node
# a leaf and returns BEFORE printPartitionsToFileArr() (:1096), so no partition
# files are written. For a non-root node that is harmless — its parent already
# wrote them. For a TOP-LEVEL cluster there is no parent, so its data is never
# materialised into the DB dir and stays only in the tronko-build input dir.
#
# Downstream: tree_list.txt lists every leaf tree (so those entries point back at
# the input dir), but final_partitions.txt only lists *partition-numbered* leaves.
# build-tronko-db.sh step 5 builds the DB FASTA by looping final_partitions.txt,
# so the orphaned clusters' sequences never reach the FASTA or the BWA index.
#
# This is detectable purely from build outputs — no logs, no cache required.
#
# Usage: check-orphaned-partitions.sh <db_dir> [<db_dir> ...]

set -uo pipefail

check_one() {
    local db="$1"
    local tl="$db/tree_list.txt"
    local fp="$db/final_partitions.txt"

    echo "=== $db"

    if [[ ! -f "$tl" ]]; then
        echo "  SKIP: no tree_list.txt (not a partitioned tronko build)"
        return
    fi

    local n_trees n_parts n_orphan
    n_trees=$(wc -l < "$tl")
    if [[ -f "$fp" ]]; then n_parts=$(wc -l < "$fp"); else n_parts=0; fi

    # An orphan is any tree_list line whose filename is NOT partition-named.
    # That is exactly the set build-tronko-db.sh step 5 cannot see.
    mapfile -t orphans < <(grep -vE 'RAxML_bestTree\.partition[0-9]+\.reroot' "$tl" || true)
    n_orphan=${#orphans[@]}

    printf "  trees=%s  partitions=%s  orphaned=%s\n" "$n_trees" "$n_parts" "$n_orphan"

    if (( n_orphan == 0 )); then
        echo "  OK: no orphaned clusters — DB FASTA covers every leaf tree."
        return
    fi

    # Volume matters more than count: the <4 guard can trip on a LARGE cluster
    # whose min-variance split is lopsided (e.g. [10000, 2, 5000]), which would
    # orphan a big chunk of the reference set, not just a handful of sequences.
    local total=0 unknown=0 line id dir src n
    echo "  affected clusters:"
    for line in "${orphans[@]}"; do
        [[ -z "$line" ]] && continue
        id=$(basename "$line"); id=${id#RAxML_bestTree.}; id=${id%.reroot}
        dir=$(dirname "$line")
        src="$dir/$id.fasta"
        if [[ -f "$src" ]]; then
            n=$(grep -c '^>' "$src")
            total=$((total + n))
            printf "    cluster %-10s %8s seqs   (%s)\n" "$id" "$n" "$dir"
        else
            unknown=$((unknown + 1))
            printf "    cluster %-10s %8s          (source dir missing: %s)\n" "$id" "?" "$dir"
        fi
    done

    echo "  --------"
    if (( unknown == 0 )); then
        local dbfa dbn
        dbfa=$(find "$db" -maxdepth 1 -name '*.fasta' ! -name 'partition*' ! -name 'marker.fasta' | head -1)
        if [[ -n "$dbfa" ]]; then
            dbn=$(grep -c '^>' "$dbfa")
            printf "  MISSING from DB FASTA: %s seqs  (DB has %s; should have %s)\n" \
                   "$total" "$dbn" "$((dbn + total))"
            awk -v t="$total" -v d="$dbn" 'BEGIN{printf "  impact: %.4f%% of the reference set\n", 100*t/(d+t)}'
        else
            printf "  MISSING from DB FASTA: %s seqs\n" "$total"
        fi
    else
        echo "  $unknown cluster(s) have no source dir — volume unknown from here."
        echo "  Recoverable anyway: the leaf names live in reference_tree.txt, and the"
        echo "  sequences can be pulled from the original rCRUX input FASTA by name."
    fi
    echo "  FIXABLE POST-HOC: yes. reference_tree.txt/.trkb already contain all"
    echo "  $n_trees trees at correct indices, so no rebuild is needed — only the"
    echo "  FASTA/BWA index and the two list files need repair."
    return 1
}

if [[ $# -eq 0 ]]; then
    echo "usage: $0 <db_dir> [<db_dir> ...]" >&2
    exit 2
fi
# Non-zero exit when any database has orphans, so callers can branch on it.
# Without this the script is informational only and a build wrapper cannot tell
# a clean database from one that is silently missing sequences.
rc=0
for d in "$@"; do check_one "$d" || rc=1; done
exit $rc
