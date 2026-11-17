"""Tests for scripts/release_tronko.sh: fetch → build → publish in one command.

The build tools (famsa, VeryFastTree, nw_reroot, ancestralclust, bwa,
tronko-build, tronko-convert) are small stand-ins that write what the real
tools write, so the script's own work is tested: clustering choice, per-
cluster files, flags, final partitions, params.json, build.log, and the real
fetch and publish against a local git + DVC repo.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
MARKER = "CO1_SCRIPT"
RCRUX_TAG = f"{MARKER}_v1.0.0"
N_SEQS = 6

STUBS = {
    # Wraps sequences at 3 bases, as FAMSA wraps long lines.
    "famsa": r"""#!/bin/sh
echo "$*" > "$CALLS_DIR/famsa.args"
for a; do in="$out"; out="$a"; done
awk '/^>/ { print; next } { s = $0; while (length(s) > 3) { print substr(s, 1, 3); s = substr(s, 4) } print s }' "$in" > "$out"
""",
    "VeryFastTree": r"""#!/bin/sh
for a; do msa="$a"; done
printf '(%s);\n' "$(grep '^>' "$msa" | cut -c2- | paste -sd, -)"
""",
    "nw_reroot": '#!/bin/sh\ncat "$1"\n',
    "bwa": '#!/bin/sh\nfor e in amb ann bwt pac sa; do echo index > "$2.$e"; done\n',
    # Splits the input into two clusters, alternating records.
    "ancestralclust": r"""#!/bin/sh
while [ $# -gt 0 ]; do case "$1" in -i) in="$2"; shift 2;; -d) d="$2"; shift 2;; *) shift;; esac; done
awk -v d="$d" '/^>/ { n++ } { print > (d (n % 2) ".fasta") }' "$in"
echo "$*" > "$d/args"
""",
    # One final partition per cluster, gap-free, as tronko-build -y writes them.
    "tronko-build": r"""#!/bin/sh
echo "$*" > "$CALLS_DIR/tronko-build.args"
while [ $# -gt 0 ]; do case "$1" in -e) e="$2"; shift 2;; -n) n="$2"; shift 2;; -d) d="$2"; shift 2;; *) shift;; esac; done
i=0; : > "$d/final_partitions.txt"
while [ "$i" -lt "$n" ]; do
  tr -d '-' < "$e/${i}_MSA.fasta" > "$d/partition$i.fasta"
  cp "$e/${i}_taxonomy.txt" "$d/partition${i}_taxonomy.txt"
  echo "$i" >> "$d/final_partitions.txt"; i=$((i + 1))
done
echo "tree of $n clusters" > "$d/reference_tree.txt"
""",
    "tronko-convert": r"""#!/bin/sh
while [ $# -gt 0 ]; do case "$1" in -i) in="$2"; shift 2;; -o) out="$2"; shift 2;; *) shift;; esac; done
{ echo TRKB; cat "$in"; } > "$out"
""",
}


def git(cwd: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.fixture
def env(tmp_path, monkeypatch):
    for key in ("AUTHOR", "COMMITTER"):
        monkeypatch.setenv(f"GIT_{key}_NAME", "test")
        monkeypatch.setenv(f"GIT_{key}_EMAIL", "test@example.com")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    monkeypatch.setenv(
        "PATH",
        os.pathsep.join([str(bin_dir), str(Path(sys.executable).parent), os.environ["PATH"]]),
    )
    monkeypatch.setenv("PYTHON", sys.executable)
    monkeypatch.setenv("TRONKO_BUILD", str(bin_dir / "tronko-build"))
    monkeypatch.setenv("TRONKO_CONVERT", str(bin_dir / "tronko-convert"))
    monkeypatch.setenv("CALLS_DIR", str(tmp_path))

    # A tronko checkout holding the scripts, committed, so its HEAD is the sha.
    tronko = tmp_path / "tronko"
    (tronko / "scripts").mkdir(parents=True)
    for name in ("release_tronko.sh", "tronko_release.py", "tronko_build.conf"):
        shutil.copy2(SCRIPTS / name, tronko / "scripts" / name)
    git(tmp_path, "init", "-q", str(tronko))
    git(tronko, "add", "-A")
    git(tronko, "commit", "-q", "-m", "tronko")

    # rcrux-py on a bare 'GitHub', with a tagged curated release.
    origin, dev, remote = tmp_path / "github.git", tmp_path / "rcrux-dev", tmp_path / "remote"
    remote.mkdir()
    git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    git(tmp_path, "clone", "-q", str(origin), str(dev))
    git(dev, "checkout", "-q", "-b", "main")
    subprocess.run(["dvc", "init", "-q"], cwd=dev, check=True)
    subprocess.run(["dvc", "remote", "add", "-d", "rcrux", str(remote)], cwd=dev, check=True)
    curated = dev / "databases" / MARKER / "curated"
    curated.mkdir(parents=True)
    ids = [f"seq_{i:03d}" for i in range(N_SEQS)]
    (curated / "lca.fasta").write_text(
        "".join(f">{s}\nACGTACGT{'A' * i}\n" for i, s in enumerate(ids))
    )
    (curated / "lca_taxonomy.txt").write_text(
        "".join(f"{s}\tEukaryota;Arthropoda;Insecta;Diptera;F;G;S{i}\n" for i, s in enumerate(ids))
    )
    subprocess.run(["dvc", "add", "-q", f"databases/{MARKER}/curated"], cwd=dev, check=True)
    git(dev, "add", "-A")
    git(dev, "commit", "-q", "-m", "rcrux release")
    subprocess.run(["dvc", "push", "-q", "-r", "rcrux"], cwd=dev, check=True)
    git(dev, "tag", RCRUX_TAG)
    git(dev, "push", "-q", "-u", "origin", "main", "--tags")

    checkout = tmp_path / "rcrux-checkout"
    git(tmp_path, "clone", "-q", str(origin), str(checkout))
    return {"tmp": tmp_path, "tronko": tronko, "origin": origin, "checkout": checkout}


def run(env, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "bash", str(env["tronko"] / "scripts" / "release_tronko.sh"),
            "--marker", MARKER, "--rcrux-tag", RCRUX_TAG,
            "--rcrux-checkout", str(env["checkout"]), "--rcrux-repo", str(env["origin"]),
            "--work", str(env["tmp"] / "work"), "--threads", "2", *extra,
        ],
        capture_output=True, text=True, check=False,
    )  # fmt: skip


def published(env) -> dict:
    return json.loads(
        (env["checkout"] / "databases" / MARKER / "tronko" / "metadata.json").read_text()
    )


def test_one_cluster_build_is_published_with_its_settings(env):
    result = run(env, "--push")
    assert result.returncode == 0, result.stdout + result.stderr
    checkout, work = env["checkout"], env["tmp"] / "work"

    branch = git(checkout, "rev-parse", "--abbrev-ref", "HEAD")
    assert branch == f"publish/{MARKER}-tronko-v1.0.0"
    assert git(checkout, "ls-remote", "--heads", "origin", branch)
    meta = published(env)
    assert meta["tronko_sha"] == git(env["tronko"], "rev-parse", "HEAD")
    assert meta["rcrux_tag"] == RCRUX_TAG

    params = meta["build_params"]
    assert (params["input_sequences"], params["clusters"], params["clustered"]) == (
        N_SEQS,
        1,
        False,
    )
    assert params["tree_builder"] == "VeryFastTree"
    assert params["config_files"] == ["scripts/tronko_build.conf"]
    assert params["settings"]["PARTITION"] == "sum_of_pairs"
    assert params["settings"]["SP_THRESHOLD"] == "0.1"
    assert params["settings"]["CLUSTER_CUTOFF"] == "25000"
    args = (env["tmp"] / "tronko-build.args").read_text().split()
    assert args[:2] == ["-y", "-e"] and args[3:5] == ["-n", "1"]
    assert args[7:] == ["-s", "-u", "0.1", "-a", "-c", "2"]

    # FAMSA's wrapped lines were joined; the final FASTA is gap-free, every record.
    msa = (work / "build" / "merged" / "0_MSA.fasta").read_text().splitlines()
    assert len(msa) == 2 * N_SEQS
    out = (work / "out" / f"{MARKER}.fasta").read_text()
    assert out.count(">") == N_SEQS and "-" not in out
    assert (work / "out" / f"{MARKER}_taxonomy.txt").read_text().count("\n") == N_SEQS
    log = (work / "build.log").read_text()
    assert "6 sequences, at most 25000: one cluster" in log and "tronko-build" in log


def test_large_input_is_clustered_and_min_leaf_partitioning_is_passed(env):
    result = run(env, "--cluster-cutoff", "3", "--min-leaf", "500")
    assert result.returncode == 0, result.stdout + result.stderr
    params = published(env)["build_params"]
    assert params["clustered"] is True and params["clusters"] == 2
    assert (params["settings"]["PARTITION"], params["settings"]["MIN_LEAF"]) == ("min_leaf", "500")
    assert params["overrides"] == [
        "set:THREADS=2",
        "set:CLUSTER_CUTOFF=3",
        "set:PARTITION=min_leaf",
        "set:MIN_LEAF=500",
    ]
    args = (env["tmp"] / "tronko-build.args").read_text().split()
    assert args[3:5] == ["-n", "2"] and "-v -f 500" in " ".join(args)
    merged = env["tmp"] / "work" / "build" / "merged"
    taxonomy_rows = sum(len(p.read_text().splitlines()) for p in merged.glob("*_taxonomy.txt"))
    assert taxonomy_rows == N_SEQS, "every sequence keeps its taxonomy row in its cluster"
    out = (env["tmp"] / "work" / "out" / f"{MARKER}.fasta").read_text()
    assert out.count(">") == N_SEQS


def test_uncommitted_tronko_changes_stop_it(env):
    (env["tronko"] / "scripts" / "tronko_release.py").write_text("# edited\n")
    result = run(env)
    assert result.returncode != 0
    assert "uncommitted changes" in result.stderr


def test_missing_required_arguments(env):
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--marker", MARKER],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode != 0 and "needs --rcrux-tag" in result.stderr


def commit(repo: Path, path: str, text: str) -> None:
    (repo / path).parent.mkdir(parents=True, exist_ok=True)
    (repo / path).write_text(text)
    git(repo, "add", path)
    git(repo, "commit", "-q", "-m", f"edit {path}")


def test_a_marker_config_is_picked_up_and_reaches_the_tools(env):
    commit(
        env["tronko"], f"scripts/tronko_build.d/{MARKER}.conf",
        "# this marker partitions by leaf count\nPARTITION=min_leaf\nMIN_LEAF = 321\n"
        'FAMSA_ARGS="--refine-mode on"\nTRONKO_BUILD_ARGS=\n',
    )  # fmt: skip
    result = run(env)
    assert result.returncode == 0, result.stdout + result.stderr
    params = published(env)["build_params"]
    assert params["config_files"] == [
        "scripts/tronko_build.conf",
        f"scripts/tronko_build.d/{MARKER}.conf",
    ]
    assert (params["settings"]["MIN_LEAF"], params["settings"]["FAMSA_ARGS"]) == (
        "321",
        "--refine-mode on",
    )
    args = (env["tmp"] / "tronko-build.args").read_text().split()
    assert "-v -f 321" in " ".join(args) and "-a" not in args, "TRONKO_BUILD_ARGS emptied"
    assert (env["tmp"] / "famsa.args").read_text().split()[:4] == ["-t", "2", "--refine-mode", "on"]


@pytest.mark.parametrize(
    ("order", "expected"),
    [(["--config", "{conf}", "--set", "SP_THRESHOLD=0.3"], "0.3"),
     (["--set", "SP_THRESHOLD=0.3", "--config", "{conf}"], "0.25")],
)  # fmt: skip
def test_later_overrides_win(env, order, expected):
    conf = env["tmp"] / "run.conf"
    conf.write_text("SP_THRESHOLD='0.25'\n")
    args = [a.format(conf=conf) for a in order]
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--print-config", *args],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert f"SP_THRESHOLD={expected}" in result.stdout.splitlines()


@pytest.mark.parametrize(
    ("line", "message"),
    [("SP_THRESHHOLD=0.2", "unknown setting SP_THRESHHOLD"),
     ("PARTITION=leaves", "PARTITION must be"),
     ("CLUSTER_CUTOFF=lots", "CLUSTER_CUTOFF must be a whole number"),
     ("just words", "expected KEY=VALUE")],
)  # fmt: skip
def test_bad_settings_are_refused(env, line, message):
    conf = env["tmp"] / "bad.conf"
    conf.write_text(line + "\n")
    result = run(env, "--config", str(conf))
    assert result.returncode != 0 and message in result.stderr


def test_config_values_are_never_run_as_shell(env):
    conf = env["tmp"] / "run.conf"
    marker_file = env["tmp"] / "ran"
    conf.write_text(f"FAMSA_ARGS=$(touch {marker_file})\n")
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--print-config",
         "--config", str(conf)],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert f"FAMSA_ARGS=$(touch {marker_file})" in result.stdout
    assert not marker_file.exists()


def test_steps_rebuild_after_an_edit_then_publish_separately(env):
    first = run(env, "--steps", "fetch,build")
    assert first.returncode == 0, first.stdout + first.stderr
    assert "Publish" not in first.stdout
    assert git(env["checkout"], "rev-parse", "--abbrev-ref", "HEAD") == "main"

    # Edit a setting, rebuild without fetching again, then publish.
    commit(env["tronko"], "scripts/tronko_build.conf",
           (env["tronko"] / "scripts" / "tronko_build.conf").read_text()
           .replace("SP_THRESHOLD=0.1", "SP_THRESHOLD=0.2"))  # fmt: skip
    rebuilt = run(env, "--steps", "build")
    assert rebuilt.returncode == 0, rebuilt.stdout + rebuilt.stderr
    assert "Fetch" not in rebuilt.stdout
    published_ = run(env, "--steps", "publish")
    assert published_.returncode == 0, published_.stdout + published_.stderr
    meta = published(env)
    assert meta["build_params"]["settings"]["SP_THRESHOLD"] == "0.2"
    assert meta["tronko_sha"] == git(env["tronko"], "rev-parse", "HEAD")


def test_publish_alone_needs_a_finished_build(env):
    result = run(env, "--steps", "publish")
    assert result.returncode != 0 and "no finished build" in result.stderr
