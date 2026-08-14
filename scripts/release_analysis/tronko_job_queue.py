#!/usr/bin/env python3
"""Serial tronko job runner with a memory gate.

Exists because running two CO1-scale tronko jobs concurrently OOM'd a 48 GB
machine: tronko-assign + minimap2 against CO1_Metazoa peaks at ~24.6 GB.

The gate is the point. It waits for memory to be *actually released* before
starting the next job, rather than for an output file to appear -- a
finished-looking output file does not mean the writing process has exited and
freed its pages, and that assumption is precisely what caused the OOM.

Jobs are described in a small JSON file so a run is reproducible and resumable:

    [
      {"name": "jalama_vert12S/post",
       "binary": "/path/bin/tronko-assign.post",
       "db": "/path/db/vert12S",
       "marker_fasta": "vert12S.fasta",
       "cell": "/path/work/jalama_vert12S",
       "args": ["--aligner", "minimap2"]}
    ]

Usage:
    tronko_job_queue.py --jobs jobs.json --out-dir runs/ [--max-rss-gb 8]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Production parameters; also the compiled-in defaults at both commits
# (tronko-assign.c:939-956). cores=1 because tronko is ~3.2% non-deterministic
# above -C 1, which would sit on top of the signal being measured.
PROD_ARGS = [
    "-6", "--number-of-cores", "1", "--Cinterval", "0.02", "--aligner", "minimap2",
    "-u", "0.0001", "--max-leaf-matches", "10",
    "--best-leaf-threshold", "-0.1", "--best-leaf-max-votes", "10",
]

# tronko's option struct uses fixed 200-byte buffers for -f/-a/-o (global.h:227)
# and sscanf("%s") overflows rather than truncating.
MAX_TRONKO_PATH = 199


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def tronko_rss_bytes() -> int:
    """Total RSS held by any live tronko-assign process."""
    try:
        out = subprocess.run(
            ["ps", "-Ao", "rss=,comm="], capture_output=True, text=True, check=True
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 0
    total = 0
    for line in out.splitlines():
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        rss_kb, comm = parts
        if "tronko-assign" in comm:
            try:
                total += int(rss_kb) * 1024
            except ValueError:
                pass
    return total


def wait_for_memory(limit_bytes: int, poll_s: float, timeout_s: float, log) -> bool:
    """Block until live tronko processes hold less than limit_bytes."""
    waited = 0.0
    while True:
        held = tronko_rss_bytes()
        if held < limit_bytes:
            if waited:
                log(f"    memory gate cleared after {waited:.0f}s "
                    f"({held / 2**30:.2f} GB still held)")
            return True
        if waited == 0:
            log(f"    memory gate: {held / 2**30:.2f} GB held by live tronko "
                f"processes, waiting for < {limit_bytes / 2**30:.2f} GB")
        if waited >= timeout_s:
            log(f"    memory gate TIMED OUT after {waited:.0f}s "
                f"({held / 2**30:.2f} GB still held)")
            return False
        time.sleep(poll_s)
        waited += poll_s


def count_queries(cell: Path) -> int | None:
    f = cell / "pF.fasta"
    if not f.exists():
        return None
    n = 0
    with f.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def peak_rss_bytes(time_output: str) -> int | None:
    """Parse the peak RSS out of `/usr/bin/time -l` (macOS) or `-v` (GNU)."""
    for line in time_output.splitlines():
        s = line.strip()
        if "maximum resident set size" in s.lower():
            tok = s.split()[0]
            if tok.isdigit():
                return int(tok)  # macOS reports bytes
        if "Maximum resident set size" in s and "kbytes" in s:
            tok = s.rsplit(":", 1)[-1].strip()
            if tok.isdigit():
                return int(tok) * 1024  # GNU reports kbytes
    return None


def run_job(job: dict, out_dir: Path, args, log) -> dict:
    name = job["name"]
    cell = Path(job["cell"])
    safe = name.replace("/", "__")
    tsv = out_dir / f"{safe}.tsv"
    logf = out_dir / f"{safe}.log"

    want = count_queries(cell)
    record: dict = {"name": name, "started": now(), "output": str(tsv)}

    # A non-empty output is not proof of completion; a killed run leaves a
    # partial TSV. Only skip when the row count matches exactly (+1 header).
    if tsv.exists() and tsv.stat().st_size > 0 and want is not None:
        have = sum(1 for _ in tsv.open())
        if have == want + 1:
            log(f"  {name}: already complete ({have} lines), skipping")
            return {**record, "status": "skipped_complete", "rows": have}
        log(f"  {name}: partial ({have}/{want + 1} lines), redoing")
        tsv.rename(tsv.with_suffix(f".partial.{int(time.time())}.tsv"))

    db = Path(job["db"])
    cmd_core = [
        job["binary"], "-r",
        "-f", str(db / "reference_tree.trkb"),
        "-a", str(db / job["marker_fasta"]),
        "-p", "-z", "-w",
        "-1", str(cell / "pF.fasta"),
        "-2", str(cell / "pR.fasta"),
        *job.get("args", PROD_ARGS),
        "-o", str(tsv),
    ]
    for flag, val in (("-f", cmd_core[3]), ("-a", cmd_core[5]), ("-o", str(tsv))):
        if len(val) > MAX_TRONKO_PATH:
            raise SystemExit(
                f"{name}: {flag} path is {len(val)} bytes, over tronko's "
                f"{MAX_TRONKO_PATH}-byte limit; use a shorter working directory"
            )

    timer = shutil.which("gtime") or "/usr/bin/time"
    cmd = [timer, "-l", *cmd_core] if Path(timer).exists() else cmd_core

    if not wait_for_memory(
        int(args.max_rss_gb * 2**30), args.poll_seconds, args.gate_timeout_seconds, log
    ):
        return {**record, "status": "gate_timeout",
                "note": "another tronko process never released its memory; "
                        "this is a queue bug, not something to retry blindly"}

    log(f"  {name}: start {now()}")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    logf.write_text((proc.stdout or "") + "\n--- stderr ---\n" + (proc.stderr or ""))

    rows = sum(1 for _ in tsv.open()) if tsv.exists() else 0
    peak = peak_rss_bytes(proc.stderr or "")
    complete = want is None or rows == want + 1

    # tronko printf()s fatal errors to stdout, not stderr.
    detail = (proc.stderr or "").strip() or (proc.stdout or "").strip()
    status = (
        "ok" if proc.returncode == 0 and complete
        else "incomplete" if proc.returncode == 0
        else "failed"
    )
    if status != "ok":
        log(f"  {name}: {status} rc={proc.returncode} rows={rows}/{(want or 0) + 1}")
        if detail:
            log(f"    {detail.splitlines()[-1][:200]}")
    else:
        log(f"  {name}: ok {rows} rows in {elapsed / 60:.1f} min"
            + (f", peak RSS {peak / 2**30:.2f} GB" if peak else ""))

    # index-integrity guard: f0a608a prints this when the reference silently truncates
    if "exceeds one index batch" in (proc.stderr or "") + (proc.stdout or ""):
        log(f"  {name}: !! INDEX BATCH OVERFLOW in stderr -- this is the bug that "
            f"broke the July run; the output is not trustworthy")
        status = "index_overflow"

    return {
        **record, "finished": now(), "status": status,
        "returncode": proc.returncode, "rows": rows,
        "expected_rows": (want + 1) if want is not None else None,
        "elapsed_seconds": round(elapsed, 1),
        "peak_rss_bytes": peak,
        "peak_rss_gb": round(peak / 2**30, 2) if peak else None,
        "command": cmd_core,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--max-rss-gb", type=float, default=8.0,
                    help="start the next job only when live tronko processes "
                         "hold less than this")
    ap.add_argument("--poll-seconds", type=float, default=10.0)
    ap.add_argument("--gate-timeout-seconds", type=float, default=3600.0)
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    jobs = json.loads(args.jobs.read_text())

    lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        lines.append(msg)

    log(f"queue: {len(jobs)} jobs, serial, memory gate at {args.max_rss_gb} GB")
    results = []
    for i, job in enumerate(jobs, 1):
        log(f"[{i}/{len(jobs)}] {job['name']}")
        results.append(run_job(job, out, args, log))
        (out / "queue_results.json").write_text(json.dumps(results, indent=2))

    bad = [r for r in results if r["status"] not in ("ok", "skipped_complete")]
    peaks = [r["peak_rss_gb"] for r in results if r.get("peak_rss_gb")]
    log("")
    log(f"done: {len(results) - len(bad)}/{len(results)} ok"
        + (f", peak RSS observed {max(peaks):.2f} GB" if peaks else ""))
    for r in bad:
        log(f"  FAILED {r['name']}: {r['status']}")
    (out / "queue_log.txt").write_text("\n".join(lines) + "\n")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
