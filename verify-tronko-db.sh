#!/bin/bash
# Verify a finished tronko database. Standalone: works on any database that is
# already built, without rebuilding and without the build cache.
#
# Exists because two real defects shipped in databases built this session and
# were invisible to every check that existed at the time:
#
#   1. The benchmark harness resolves the reference FASTA with
#      sorted(glob("*.fasta")) whenever marker= is not passed -- which is every
#      production caller. It skips BWA suffixes and partition* but NOT
#      input.fasta, so any primer name sorting after "input" silently hands
#      tronko-assign the raw unaligned rCRUX input, with a BWA index that does
#      not match it and leaf coordinates unrelated to reference_tree.trkb.
#      (vert12S: "input" < "vert12S", so it lost. 16Smamm and MiFish escaped
#      only because "1" < "i".)
#
#   2. Reference coverage was checked by record COUNT. AncestralClust can drop
#      some sequences while duplicating others across clusters, and the two
#      nearly cancel: vert12S reported 63,183 records against 63,298 inputs --
#      a plausible-looking 0.2% gap -- while actually holding only 56,260
#      distinct accessions and missing 7,038 (11.1%). Coverage must be
#      compared as SETS.
#
# Usage: verify-tronko-db.sh <db_dir> [<db_dir> ...] [--max-missing-pct N]
#                            [--no-smoke] [--json]
#
# Exit 0 only if every database passes every abort-level check.

set -uo pipefail

LIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tronko-db.sh"
[ -f "$LIB" ] || { echo "ERROR: missing $LIB" >&2; exit 1; }
# shellcheck source=lib/tronko-db.sh
source "$LIB"

MAX_MISSING_PCT=1.0
RUN_SMOKE=1
WRITE_JSON=1
DBS=()

while [ $# -gt 0 ]; do
    case "$1" in
        --max-missing-pct) MAX_MISSING_PCT="$2"; shift 2 ;;
        --no-smoke)        RUN_SMOKE=0; shift ;;
        --no-json)         WRITE_JSON=0; shift ;;
        -h|--help)         sed -n '2,26p' "$0"; exit 0 ;;
        -*)                die "unknown option: $1" ;;
        *)                 DBS+=("$1"); shift ;;
    esac
done
[ ${#DBS[@]} -gt 0 ] || die "usage: $0 <db_dir> [<db_dir> ...]"

export PATH="$TRONKO_REPO_DIR/bin:$TRONKO_REPO_DIR/tronko-build:$TRONKO_REPO_DIR/tronko-convert:$TRONKO_REPO_DIR/tronko-assign:$PATH"

OVERALL=0

# ---------------------------------------------------------------------------
verify_one() {
    local db fails=0 warns=0
    db=$(cd "$1" 2>/dev/null && pwd) || { echo "ERROR: no such directory: $1" >&2; return 1; }
    echo "=== $db"

    local -a FAIL_MSGS=() WARN_MSGS=()
    fail() { FAIL_MSGS+=("$1"); printf '  FAIL  %s\n' "$1"; fails=$((fails+1)); }
    pass() { printf '  ok    %s\n' "$1"; }
    note() { WARN_MSGS+=("$1"); printf '  warn  %s\n' "$1"; warns=$((warns+1)); }

    # -- primer -------------------------------------------------------------
    local primer
    if ! primer=$(db_primer "$db"); then
        fail "cannot determine primer name (no non-alias .fasta found)"
        echo "  RESULT: FAIL ($fails)"; return 1
    fi
    pass "primer: $primer"

    # -- artifacts ----------------------------------------------------------
    # The assign-time set plus the ablation set. A missing input_taxonomy.txt
    # makes the harness treat the database as "not built" and silently drop
    # tronko from --tools, so it is an abort, not a warning.
    local f
    for f in reference_tree.trkb reference_tree.txt.gz tree_list.txt final_partitions.txt \
             "${primer}.fasta" "${primer}_taxonomy.txt" \
             marker.fasta marker_taxonomy.txt input.fasta input_taxonomy.txt \
             build_metadata.json; do
        [ -e "$db/$f" ] || fail "missing artifact: $f"
    done
    for f in "${BWA_EXTS[@]}"; do
        [ -e "$db/${primer}.fasta.$f" ] || fail "missing BWA index: ${primer}.fasta.$f"
        [ -e "$db/marker.fasta.$f" ]    || fail "missing alias: marker.fasta.$f"
    done
    local g
    for g in 'partition*.fasta' 'partition*_MSA.fasta' 'partition*_taxonomy.txt' \
             'RAxML_bestTree.partition*.reroot'; do
        compgen -G "$db/$g" > /dev/null || fail "no files matching $g"
    done
    [ $fails -eq 0 ] && pass "artifact checklist"

    # -- harness resolution simulation (defect 1) ---------------------------
    local res hfa htx
    res=$(harness_resolution "$db")
    hfa=${res%% *}; htx=${res##* }
    if [ "$hfa" != "${primer}.fasta" ]; then
        fail "harness would load '$hfa' as the reference FASTA, not ${primer}.fasta"
    else
        pass "harness resolves FASTA -> $hfa"
    fi
    if [ "$htx" != "${primer}_taxonomy.txt" ]; then
        fail "harness would load '$htx' as the taxonomy, not ${primer}_taxonomy.txt"
    else
        pass "harness resolves taxonomy -> $htx"
    fi

    # -- alias integrity (dereferenced content, not readlink) ---------------
    # dvc add materialises symlinks into real files; after an in-place repair a
    # previously-created alias becomes a stale COPY holding pre-repair content.
    local a b
    if [ -e "$db/marker.fasta" ] && [ -e "$db/${primer}.fasta" ]; then
        a=$(md5sum "$db/marker.fasta" | cut -d' ' -f1)
        b=$(md5sum "$db/${primer}.fasta" | cut -d' ' -f1)
        [ "$a" = "$b" ] && pass "marker.fasta content matches ${primer}.fasta" \
                        || fail "marker.fasta is STALE (content differs from ${primer}.fasta)"
    fi
    if [ -e "$db/marker_taxonomy.txt" ] && [ -e "$db/${primer}_taxonomy.txt" ]; then
        a=$(md5sum "$db/marker_taxonomy.txt" | cut -d' ' -f1)
        b=$(md5sum "$db/${primer}_taxonomy.txt" | cut -d' ' -f1)
        [ "$a" = "$b" ] && pass "marker_taxonomy.txt content matches" \
                        || fail "marker_taxonomy.txt is STALE"
    fi
    if [ -e "$db/${primer}.fasta" ] && [ -e "$db/${primer}.fasta.bwt" ]; then
        [ "$db/${primer}.fasta.bwt" -nt "$db/${primer}.fasta" ] \
            && pass "BWA index newer than FASTA" \
            || fail "BWA index is OLDER than ${primer}.fasta -- reindex required"
    fi

    # -- FASTA / taxonomy agreement ----------------------------------------
    local tmp; tmp=$(mktemp -d)
    local nrec ndist ntax
    nrec=$(grep -c '^>' "$db/${primer}.fasta" 2>/dev/null || echo 0)
    grep '^>' "$db/${primer}.fasta" 2>/dev/null | sed 's/^>//' | awk '{print $1}' | sort > "$tmp/db_all"
    sort -u "$tmp/db_all" > "$tmp/db_uniq"
    ndist=$(wc -l < "$tmp/db_uniq")
    ntax=$(wc -l < "$db/${primer}_taxonomy.txt" 2>/dev/null || echo 0)
    [ "$nrec" = "$ntax" ] && pass "FASTA records == taxonomy rows ($nrec)" \
                          || fail "FASTA records ($nrec) != taxonomy rows ($ntax)"

    # -- coverage, as SETS (defect 2) --------------------------------------
    local ninput nmissing dupes pct
    if [ -f "$db/input.fasta" ]; then
        grep '^>' "$db/input.fasta" | sed 's/^>//' | awk '{print $1}' | sort -u > "$tmp/in_uniq"
        ninput=$(wc -l < "$tmp/in_uniq")
        nmissing=$(comm -23 "$tmp/in_uniq" "$tmp/db_uniq" | wc -l)
        dupes=$((nrec - ndist))
        pct=$(awk -v m="$nmissing" -v n="$ninput" 'BEGIN{printf "%.3f", n?100*m/n:0}')
        printf '  ---   coverage: input_distinct=%s db_records=%s db_distinct=%s duplicated=%s missing=%s (%s%%)\n' \
               "$ninput" "$nrec" "$ndist" "$dupes" "$nmissing" "$pct"
        if awk -v p="$pct" -v m="$MAX_MISSING_PCT" 'BEGIN{exit !(p>m)}'; then
            fail "coverage: $nmissing of $ninput input accessions absent (${pct}% > ${MAX_MISSING_PCT}%)"
        else
            pass "coverage within ${MAX_MISSING_PCT}%"
        fi
        [ "$dupes" -gt 0 ] && note "$dupes duplicated records in the reference FASTA (AncestralClust emits duplicates across clusters)"
    else
        note "no input.fasta -- coverage cannot be checked"
        ninput=0; nmissing=0; dupes=0; pct=0
    fi

    # -- tree_list / final_partitions --------------------------------------
    local ntree npart nlisted nmiss_path
    ntree=$(wc -l < "$db/tree_list.txt" 2>/dev/null || echo 0)
    npart=$(wc -l < "$db/final_partitions.txt" 2>/dev/null || echo 0)
    nlisted=$(grep -cE 'RAxML_bestTree\.partition[0-9]+\.reroot' "$db/tree_list.txt" 2>/dev/null || echo 0)
    [ "$npart" = "$nlisted" ] && pass "final_partitions ($npart) == partition-named tree_list lines" \
                              || fail "final_partitions ($npart) != partition-named tree_list lines ($nlisted)"
    [ "$ntree" = "$nlisted" ] && pass "no orphaned clusters ($ntree trees)" \
                              || fail "$((ntree - nlisted)) orphaned cluster(s): trees present but absent from the reference FASTA"
    # Every referenced tree file must exist. Relative entries are resolved
    # against $db and then the repo root: repair-orphaned-partitions.sh wrote
    # its rewritten lines using whatever $DB it was invoked with, so a repair
    # run from the repo root leaves a handful of repo-root-relative lines among
    # otherwise absolute ones. That is a consistency problem, not a missing
    # file, and the two are worth distinguishing.
    local nrel
    # grep -c prints 0 AND exits 1 when nothing matches, so `|| echo 0` would
    # emit a second zero. Count with awk instead.
    nrel=$(awk '$0 !~ /^\//{c++} END{print c+0}' "$db/tree_list.txt" 2>/dev/null)
    nmiss_path=$(awk -v d="$db" -v r="$TRONKO_REPO_DIR" '
        { p=$0
          if (p ~ /^\//) { if (system("[ -e \""p"\" ]")) c++ }
          else if (system("[ -e \""d"/"p"\" ]") && system("[ -e \""r"/"p"\" ]")) c++ }
        END{print c+0}' "$db/tree_list.txt" 2>/dev/null)
    [ "$nmiss_path" = "0" ] && pass "all tree_list.txt entries exist on disk" \
                            || fail "$nmiss_path tree_list.txt entry/entries do not exist on disk"
    [ "$nrel" = "0" ] && pass "tree_list.txt paths are uniformly absolute" \
                      || note "$nrel tree_list.txt entries are relative while the rest are absolute (from an in-place repair)"

    # -- self-assignment smoke test ----------------------------------------
    local sm_total=0 sm_dom=0 sm_fam=0 sm_exact=0
    if [ "$RUN_SMOKE" = "1" ] && [ -f "$db/reference_tree.trkb" ] && [ "$fails" -eq 0 ]; then
        if command -v tronko-assign >/dev/null 2>&1; then
            python3 - "$db/${primer}.fasta" "$tmp/smoke.fasta" 200 <<'PY'
import sys
src, out, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs=[]; hdr=None; seq=[]
for line in open(src):
    if line.startswith('>'):
        if hdr: recs.append((hdr, ''.join(seq)))
        hdr=line.strip(); seq=[]
    else: seq.append(line.strip())
if hdr: recs.append((hdr, ''.join(seq)))
step=max(1, len(recs)//n)
with open(out,'w') as o:
    for h,s in recs[::step][:n]:
        o.write(h+'\n'+s+'\n')
PY
            ( cd "$tmp" && timeout 1800 tronko-assign -r \
                -f "$db/reference_tree.trkb" -a "$db/${primer}.fasta" \
                -o "$tmp/smoke_results.txt" -s -g "$tmp/smoke.fasta" \
                -6 -C 1 -c 5 >/dev/null 2>&1 )
            if [ -s "$tmp/smoke_results.txt" ]; then
                read -r sm_total sm_exact sm_dom sm_fam < <(
                  python3 - "$db/${primer}_taxonomy.txt" "$tmp/smoke_results.txt" <<'PY'
import sys
truth={}
for line in open(sys.argv[1]):
    p=line.rstrip('\n').split('\t',1)
    if len(p)==2: truth[p[0]]=p[1]
tot=ex=dom=fam=0
f=open(sys.argv[2]); next(f,None)
for line in f:
    c=line.rstrip('\n').split('\t')
    if len(c)<2: continue
    t=truth.get(c[0])
    if t is None: continue
    tot+=1
    tr=t.split(';'); ar=c[1].split(';')
    # reference_tree stores lineages reversed (species-first); assignment
    # output is domain-first. Compare by matching orientation empirically.
    if tr and ar and tr[0]!=ar[0] and tr[-1]==ar[0]:
        tr=list(reversed(tr))
    if t==c[1]: ex+=1
    if tr and ar and tr[0]==ar[0]: dom+=1
    if len(tr)>=5 and len(ar)>=5 and tr[4]==ar[4]: fam+=1
print(tot, ex, dom, fam)
PY
                )
                if [ "${sm_total:-0}" -gt 0 ]; then
                    printf '  ---   smoke: n=%s domain=%s%% family=%s%% exact=%s%%\n' \
                        "$sm_total" \
                        "$(awk -v a="$sm_dom" -v b="$sm_total" 'BEGIN{printf "%.1f",100*a/b}')" \
                        "$(awk -v a="$sm_fam" -v b="$sm_total" 'BEGIN{printf "%.1f",100*a/b}')" \
                        "$(awk -v a="$sm_exact" -v b="$sm_total" 'BEGIN{printf "%.1f",100*a/b}')"
                    awk -v a="$sm_dom" -v b="$sm_total" 'BEGIN{exit !(100*a/b < 90)}' \
                        && note "domain-level self-assignment below 90%"
                else
                    note "smoke test produced no comparable rows"
                fi
            else
                note "smoke test produced no output"
            fi
        else
            note "tronko-assign not on PATH -- smoke test skipped"
        fi
    fi

    # -- report -------------------------------------------------------------
    if [ "$WRITE_JSON" = "1" ]; then
        python3 - "$db" "$primer" "$ninput" "$nrec" "$ndist" "$dupes" "$nmissing" \
                 "$ntree" "$npart" "$fails" "$warns" "$sm_total" "$sm_dom" "$hfa" <<'PY' 2>/dev/null || true
import json,sys,datetime
d,primer,ninput,nrec,ndist,dup,miss,ntree,npart,fails,warns,smn,smd,hfa = sys.argv[1:15]
i=lambda x:int(x or 0)
json.dump({
 "verified_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
 "primer": primer,
 "harness_resolves_fasta": hfa,
 "input_distinct": i(ninput), "db_records": i(nrec), "db_distinct": i(ndist),
 "duplicated_records": i(dup), "missing_accessions": i(miss),
 "trees": i(ntree), "final_partitions": i(npart),
 "smoke_n": i(smn), "smoke_domain": i(smd),
 "failures": i(fails), "warnings": i(warns),
 "result": "PASS" if i(fails)==0 else "FAIL",
}, open(f"{d}/verify_report.json","w"), indent=2)
PY
    fi
    rm -rf "$tmp"

    if [ "$fails" -eq 0 ]; then
        echo "  RESULT: PASS ($warns warning(s))"
        return 0
    fi
    echo "  RESULT: FAIL ($fails failure(s), $warns warning(s))"
    return 1
}

for d in "${DBS[@]}"; do
    verify_one "$d" || OVERALL=1
    echo
done
exit $OVERALL
