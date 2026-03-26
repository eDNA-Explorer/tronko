#!/bin/bash
# Monitor pasta-3builds.sh progress
LOG="/home/ubuntu/tronko-fork/pasta-3builds.log"
MONITOR_LOG="/home/ubuntu/tronko-fork/monitor_pasta.log"

echo "=== $(date) ===" >> "$MONITOR_LOG"

# Check if process is running
if pgrep -f "pasta-3builds" > /dev/null || pgrep -f "run_pasta.py" > /dev/null; then
    echo "STATUS: RUNNING" >> "$MONITOR_LOG"
else
    echo "STATUS: NOT RUNNING" >> "$MONITOR_LOG"
fi

# Memory usage
free -h | grep Mem | awk '{print "MEMORY: " $3 " used / " $2 " total (" $7 " available)"}' >> "$MONITOR_LOG"

# Current marker (last marker header seen in log)
CURRENT_MARKER=$(grep -E "^#+ " "$LOG" 2>/dev/null | tail -1)
echo "CURRENT MARKER: ${CURRENT_MARKER:-unknown}" >> "$MONITOR_LOG"

# PASTA step (last "=== " line)
CURRENT_STEP=$(grep "^=== " "$LOG" 2>/dev/null | tail -1)
echo "CURRENT STEP: ${CURRENT_STEP:-unknown}" >> "$MONITOR_LOG"

# VeryFastTree / FastTree progress
VFT_PID=$(pgrep -f "VeryFastTree" 2>/dev/null | head -1)
FT_PID=$(pgrep -f "fasttreeMP" 2>/dev/null | head -1)
if [ -n "$VFT_PID" ]; then
    VFT_RSS=$(ps -o rss= -p "$VFT_PID" 2>/dev/null | awk '{printf "%.1f GB", $1/1048576}')
    VFT_CPU=$(ps -o %cpu= -p "$VFT_PID" 2>/dev/null)
    echo "TREE: VeryFastTree (PID $VFT_PID) CPU=${VFT_CPU}% MEM=${VFT_RSS}" >> "$MONITOR_LOG"
elif [ -n "$FT_PID" ]; then
    FT_RSS=$(ps -o rss= -p "$FT_PID" 2>/dev/null | awk '{printf "%.1f GB", $1/1048576}')
    FT_CPU=$(ps -o %cpu= -p "$FT_PID" 2>/dev/null)
    echo "TREE: FastTreeMP (PID $FT_PID) CPU=${FT_CPU}% MEM=${FT_RSS}" >> "$MONITOR_LOG"
else
    echo "TREE: no tree estimator running" >> "$MONITOR_LOG"
fi

# MAFFT processes (alignment step)
MAFFT_COUNT=$(pgrep -c mafft 2>/dev/null || echo 0)
if [ "$MAFFT_COUNT" -gt 0 ]; then
    echo "MAFFT PROCESSES: $MAFFT_COUNT" >> "$MONITOR_LOG"
fi

# tronko-build progress
TB_PID=$(pgrep -f "tronko-build" 2>/dev/null | head -1)
if [ -n "$TB_PID" ]; then
    echo "TRONKO-BUILD: running (PID $TB_PID)" >> "$MONITOR_LOG"
fi

# Last log lines
LAST_LINE=$(tail -1 "$LOG" 2>/dev/null)
echo "LAST LOG: ${LAST_LINE:-(empty)}" >> "$MONITOR_LOG"

# Check for errors
LAST_ERROR=$(grep -i "error\|traceback\|killed\|MemoryError" "$LOG" 2>/dev/null | tail -1)
if [ -n "$LAST_ERROR" ]; then
    echo "LAST ERROR: $LAST_ERROR" >> "$MONITOR_LOG"
fi

echo "" >> "$MONITOR_LOG"
