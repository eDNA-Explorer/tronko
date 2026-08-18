#!/usr/bin/env bash
# Append a status line every 2 minutes so progress is reconstructable
# after any disconnect. Runs until the build process exits.

BASE="$HOME/tronko-build-branch/databases/CO1_mlCOIintF_Fol-degen-rev/lca/ac/default"
OUT="$BASE/sp0.10"
MON="$OUT/monitor.log"

printf '%-20s %6s %7s %9s %7s %7s %6s  %s\n' \
    TIMESTAMP STEP CLUSTERS TREES RSS_GB DISK_GB CPU% NOTE >> "$MON"

while true; do
    ts=$(date -Is)

    read -r cpu rss <<< "$(ps -eo pcpu,rss,comm --no-headers \
        | awk '$3=="ancestralclust"||$3=="tronko-build"{c+=$1; if($2>r)r=$2} END{print c+0, r+0}')"
    rss_gb=$(awk -v r="$rss" 'BEGIN{printf "%.1f", r/1048576}')

    clusters=$(ls "$BASE/cache/ancestralclust/"*.fasta 2>/dev/null | wc -l)
    trees=$(ls "$BASE/cache/newick/"*.reroot 2>/dev/null | wc -l)
    parts=$(ls "$OUT"/partition*.fasta 2>/dev/null | grep -vc '_MSA\|_taxonomy')
    disk=$(df --output=used -BG / | tail -1 | tr -dc '0-9')

    if   [[ -f "$OUT/reference_tree.txt" || -f "$OUT/reference_tree.txt.gz" ]]; then step=5
    elif [[ "$parts" -gt 0 ]];                                                  then step=4
    elif [[ "$trees" -gt 0 ]];                                                  then step=2
    else                                                                             step=1
    fi

    note=""
    grep -q 'slot limit reached' "$OUT/build.log" 2>/dev/null && note="CAP-HIT"
    n_att=$(grep -c 'ancestralclust attempt' "$OUT/build.log" 2>/dev/null || echo 0)
    [[ "$n_att" -gt 1 ]] && note="$note AC-RETRY:$n_att"

    printf '%-20s %6s %7s %9s %7s %7s %6s  %s\n' \
        "$ts" "$step" "$clusters" "$trees/$parts" "$rss_gb" "$disk" "$cpu" "$note" >> "$MON"

    # stop once the build is gone and status.txt has an EXIT marker
    if ! pgrep -x ancestralclust >/dev/null && ! pgrep -x tronko-build >/dev/null \
       && grep -q '^EXIT' "$OUT/status.txt" 2>/dev/null; then
        echo "$(date -Is)  watcher exiting: $(tail -1 "$OUT/status.txt")" >> "$MON"
        break
    fi
    sleep 120
done
