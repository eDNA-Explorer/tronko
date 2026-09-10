#!/bin/bash
# Compare memory usage across multiple tronko-assign runs
# Usage: ./compare_memlogs.sh log1.tsv [log2.tsv ...]

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 log1.tsv [log2.tsv ...]"
    exit 1
fi

echo "Run,Peak_RSS_MB,Time_to_Peak_s,Post_Batch_RSS_MB,Final_RSS_MB,Total_Time_s,CPU_User_s,CPU_Sys_s"
for log in "$@"; do
    if [ ! -f "$log" ]; then
        echo "Warning: File not found: $log" >&2
        continue
    fi

    run_name=$(basename "$log" .tsv)

    # Skip header lines (starting with # or containing 'wall_time')
    # Some platforms do not expose ru_maxrss. Use the larger of the reported
    # process peak and the highest phase-sampled current RSS.
    peak_rss=$(awk -F'\t' 'NR>1 && !/^#/ && !/wall_time/ {
        value = ($5 > $3 ? $5 : $3); if (value > peak) peak = value
    } END {printf "%.1f", peak}' "$log")
    peak_time=$(awk -F'\t' -v peak="$peak_rss" 'NR>1 && !/^#/ && !/wall_time/ {
        value = ($5 > $3 ? $5 : $3); if (sprintf("%.1f", value) == peak) {print $1; exit}
    }' "$log")
    post_batch_rss=$(awk -F'\t' '$2=="BATCH_FREED" {value=$3} END {
        if (value == "") print "N/A"; else print value
    }' "$log")
    final_line=$(awk -F'\t' '/FINAL/ {print}' "$log" | tail -1)

    if [ -n "$final_line" ]; then
        final_rss=$(echo "$final_line" | cut -f3)
        total_time=$(echo "$final_line" | cut -f1)
        cpu_user=$(echo "$final_line" | cut -f6)
        cpu_sys=$(echo "$final_line" | cut -f7)
    else
        final_rss="N/A"
        total_time="N/A"
        cpu_user="N/A"
        cpu_sys="N/A"
    fi

    echo "$run_name,$peak_rss,$peak_time,$post_batch_rss,$final_rss,$total_time,$cpu_user,$cpu_sys"
done
