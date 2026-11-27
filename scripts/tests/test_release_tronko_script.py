"""Tests for scripts/release_tronko.sh: fetch → build → publish in one command.

build-tronko-db.sh is a stand-in that records its arguments and writes the
files the real one writes, so the script's own work is tested: settings to
arguments, the taxonomy columns, params.json, build.log, and the real fetch
and publish against a local git + DVC repo. One test at the end runs the real
build-tronko-db.sh on a small input when its tools are installed.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import copy_repos, lineage_nodes, write_taxdump

SCRIPTS = Path(__file__).resolve().parents[1]
REPO = SCRIPTS.parent
MARKER = "CO1_SCRIPT"
RCRUX_TAG = f"{MARKER}_v1.0.0"
N_SEQS = 6

# Copies the input through, as the real build's outputs; DROP_LAST=1 leaves the last
# sequence out. A taxonomy with taxids gets the taxid files, as the real build writes them.
BUILD_STUB = r"""#!/bin/bash
printf '%s\n' "$@" > "$CALLS_DIR/build.args"
while [ $# -gt 0 ]; do case "$1" in -f) f="$2"; shift 2;; -t) t="$2"; shift 2;; -o) o="$2"; shift 2;; -p) p="$2"; shift 2;; --taxdump) d="$2"; shift 2;; *) shift;; esac; done
mkdir -p "$o"
echo "stand-in build of $p"
if [ -n "${DROP_LAST:-}" ]; then sed '$d' "$f" | sed '$d' > "$o/$p.fasta"; else cp "$f" "$o/$p.fasta"; fi
if [ -n "${d:-}" ]; then
    cut -f1,3 "$t" > "$o/${p}_taxonomy.txt"
    cut -f1,2 "$t" > "$o/${p}_taxids.txt"
    python3 "$LINEAGE_TAXIDS" --marker "$p" --build-dir "$o" --taxdump "$d"
else
    cp "$t" "$o/${p}_taxonomy.txt"
fi
for e in amb ann bwt pac sa bwt.2bit.64; do echo index > "$o/$p.fasta.$e"; done
echo TRKB > "$o/reference_tree.trkb"
echo '{"pipeline": "stand-in", "tree_tool": "veryfasttree"}' > "$o/build_metadata.json"
echo "$o/RAxML_bestTree.partition0.reroot" > "$o/tree_list.txt"
echo tree | gzip > "$o/reference_tree.txt.gz"
cp "$o/${p}_taxonomy.txt" "$o/partition0_taxonomy.txt"
"""


def git(cwd: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.fixture(scope="session")
def repos_template(tmp_path_factory):
    """tronko and rcrux-py as git + DVC repositories with their 'GitHub's and remotes,
    built once; each test works on its own copy (the env fixture)."""
    tmp_path = tmp_path_factory.mktemp("repos").resolve()
    with pytest.MonkeyPatch.context() as patch:
        for key in ("AUTHOR", "COMMITTER"):
            patch.setenv(f"GIT_{key}_NAME", "test")
            patch.setenv(f"GIT_{key}_EMAIL", "test@example.com")
        patch.setenv("PATH", os.pathsep.join([str(Path(sys.executable).parent), os.environ["PATH"]]))
        _build_repos(tmp_path)
    return tmp_path


@pytest.fixture
def env(tmp_path, monkeypatch, repos_template):
    for key in ("AUTHOR", "COMMITTER"):
        monkeypatch.setenv(f"GIT_{key}_NAME", "test")
        monkeypatch.setenv(f"GIT_{key}_EMAIL", "test@example.com")
    copy_repos(repos_template, tmp_path)
    monkeypatch.setenv(
        "PATH", os.pathsep.join([str(Path(sys.executable).parent), os.environ["PATH"]])
    )
    monkeypatch.setenv("PYTHON", sys.executable)
    monkeypatch.setenv("BUILD_TRONKO_DB", str(tmp_path / "build-stub.sh"))
    monkeypatch.setenv("CALLS_DIR", str(tmp_path))
    monkeypatch.setenv("LINEAGE_TAXIDS", str(SCRIPTS / "lineage_taxids.py"))
    return {
        "tmp": tmp_path,
        "tronko": tmp_path / "tronko",
        # rcrux-py, read in place: a clone of it could not follow the relative remote
        "origin": tmp_path / "rcrux-dev",
        "checkout": tmp_path / "tronko-release",
        "tronko_remote": tmp_path / "tronko-remote",
        "rcrux_remote": tmp_path / "remote",
    }


def _build_repos(tmp_path: Path) -> None:
    # The DVC remotes are named relative to each repository, so a copy of these
    # directories reads and writes its own.
    stub = tmp_path / "build-stub.sh"
    stub.write_text(BUILD_STUB)

    # A tronko checkout holding the scripts, committed, so its HEAD is the sha.
    tronko = tmp_path / "tronko"
    (tronko / "scripts").mkdir(parents=True)
    for name in ("release_tronko.sh", "tronko_release.py", "lineage_taxids.py", "rcrux_primers.py", "tronko_build.conf"):
        shutil.copy2(SCRIPTS / name, tronko / "scripts" / name)
    (tronko / "lib").mkdir()
    shutil.copy2(REPO / "lib" / "tronko-db.sh", tronko / "lib" / "tronko-db.sh")
    git(tmp_path, "init", "-q", "-b", "main", str(tronko))
    subprocess.run(["dvc", "init", "-q"], cwd=tronko, check=True)
    tronko_remote = tmp_path / "tronko-remote"
    tronko_remote.mkdir()
    subprocess.run(
        ["dvc", "remote", "add", "-d", "tronko", "../tronko-remote"], cwd=tronko, check=True
    )
    git(tronko, "add", "-A")
    git(tronko, "commit", "-q", "-m", "tronko")
    # tronko's 'GitHub', and the second checkout a release is committed on.
    tronko_origin = tmp_path / "tronko-github.git"
    git(tmp_path, "clone", "-q", "--bare", str(tronko), str(tronko_origin))
    checkout = tmp_path / "tronko-release"
    git(tmp_path, "clone", "-q", str(tronko_origin), str(checkout))

    # rcrux-py on a bare 'GitHub', with a tagged curated release.
    origin, dev, remote = tmp_path / "github.git", tmp_path / "rcrux-dev", tmp_path / "remote"
    remote.mkdir()
    git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    git(tmp_path, "clone", "-q", str(origin), str(dev))
    git(dev, "checkout", "-q", "-b", "main")
    subprocess.run(["dvc", "init", "-q"], cwd=dev, check=True)
    subprocess.run(["dvc", "remote", "add", "-d", "rcrux", "../remote"], cwd=dev, check=True)
    curated = dev / "databases" / MARKER / "curated"
    curated.mkdir(parents=True)
    ids = [f"seq_{i:03d}" for i in range(N_SEQS)]
    (curated / "lca.fasta").write_text(
        "".join(f">{s}\nACGTACGT{'A' * i}\n" for i, s in enumerate(ids))
    )
    (curated / "lca_taxonomy.txt").write_text(
        "".join(
            f"{s}\t{7000 + i}\tEukaryota;Arthropoda;Insecta;Diptera;F;G;S{i}\n"
            for i, s in enumerate(ids)
        )
    )
    subprocess.run(["dvc", "add", "-q", f"databases/{MARKER}/curated"], cwd=dev, check=True)
    # The taxdump: S0..S4 are taxids 7000..7004 under genus G (6995); it lacks S5 (7005).
    nodes = lineage_nodes(6990, "Eukaryota;Arthropoda;Insecta;Diptera;F;G;S0")[:-1]
    write_taxdump(dev / "data" / "taxdump", nodes + [(7000 + i, 6995, "species", f"S{i}") for i in range(N_SEQS - 1)])
    subprocess.run(["dvc", "add", "-q", "data"], cwd=dev, check=True)
    git(dev, "add", "-A")
    git(dev, "commit", "-q", "-m", "rcrux release")
    subprocess.run(["dvc", "push", "-q", "-r", "rcrux"], cwd=dev, check=True)
    git(dev, "tag", RCRUX_TAG)
    git(dev, "push", "-q", "-u", "origin", "main", "--tags")



def run(env, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "bash", str(env["tronko"] / "scripts" / "release_tronko.sh"),
            "--marker", MARKER, "--rcrux-tag", RCRUX_TAG,
            "--tronko-checkout", str(env["checkout"]), "--rcrux-repo", str(env["origin"]),
            "--work", str(env["tmp"] / "work"), "--threads", "2", *extra,
            "--set", "ASSIGNMENT_CHECK_READS_PER_QUERY=0",
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
    # The release is stored in tronko's DVC remote; rcrux's holds only its own tier.
    rcrux_objects = sum(1 for p in env["rcrux_remote"].rglob("*") if p.is_file())
    assert any(env["tronko_remote"].rglob("*")) and rcrux_objects == 7  # curated/ and data/taxdump
    meta = published(env)
    assert meta["tronko_sha"] == git(env["tronko"], "rev-parse", "HEAD")
    assert meta["rcrux_tag"] == RCRUX_TAG

    params = meta["build_params"]
    assert params["build"] == {"pipeline": "stand-in", "tree_tool": "veryfasttree"}
    assert params["config_files"] == ["scripts/tronko_build.conf"]
    assert params["settings"]["PARTITION"] == "sum_of_pairs"
    assert params["settings"]["SP_THRESHOLD"] == "0.1"
    assert params["settings"]["CLUSTER_CUTOFF"] == "25000"
    args = build_args(env)
    assert args[args.index("-T") + 1] == "2"
    assert args[args.index("-s") + 1] == "0.1" and "-L" in args and "--min-leaf" not in args
    assert args[args.index("--tree-tool") + 1] == "veryfasttree"
    assert "-J" not in args  # no cap: Step 2's thread shares decide

    # The build gets rCRUX's taxonomy as it is, with the taxdump, and publishes the
    # taxonomy as accession and lineage with the taxids beside it.
    assert args[args.index("-t") + 1].endswith("lca_taxonomy.txt")
    assert args[args.index("--taxdump") + 1] == str(work / "input" / "taxdump")
    released = checkout / "databases" / MARKER / "tronko"
    assert (released / f"{MARKER}_taxonomy.txt").read_text().splitlines()[0] == (
        "seq_000\tEukaryota;Arthropoda;Insecta;Diptera;F;G;S0"
    )
    taxids = (released / f"{MARKER}_taxids.txt").read_text().splitlines()
    assert taxids[0] == "seq_000\t7000" and len(taxids) == N_SEQS
    assert f"{MARKER}_taxids.txt" in meta["files"]
    # and the taxid of every lineage path tronko-assign can print
    table = (released / f"{MARKER}_lineage_taxids.txt").read_text().splitlines()
    assert "Eukaryota;Arthropoda;Insecta;Diptera;F;G\tgenus\t6995" in table
    assert "Eukaryota;Arthropoda;Insecta;Diptera;F;G;S3\tspecies\t7003" in table
    assert len(table) == 6 + N_SEQS and f"{MARKER}_lineage_taxids.txt" in meta["files"]
    # A path without a taxid is logged, so the published build.log and the QC report show it.
    missing = f"WARNING: no NCBI taxid for lineage Eukaryota;Arthropoda;Insecta;Diptera;F;G;S{N_SEQS - 1}"
    assert missing in (released / "build.log").read_text()
    assert missing in git(checkout, "log", "-1", "--format=%b")
    assert (work / "out" / f"{MARKER}.fasta").read_text().count(">") == N_SEQS
    assert "stand-in build" in (work / "build.log").read_text()


def build_args(env) -> list[str]:
    return (env["tmp"] / "build.args").read_text().split("\n")[:-1]


def test_settings_become_build_arguments(env):
    result = run(
        env, "--cluster-cutoff", "3", "--min-leaf", "500",
        "--set", "SP_FORMULA=corrected", "--set", "TREE_ARGS=-gamma -fastexp 2",
        "--set", "CLUSTER_JOBS=4",
    )  # fmt: skip
    assert result.returncode == 0, result.stdout + result.stderr
    params = published(env)["build_params"]
    assert (params["settings"]["PARTITION"], params["settings"]["MIN_LEAF"]) == ("min_leaf", "500")
    assert params["settings"]["TREE_ARGS"] == "-gamma -fastexp 2"
    assert params["overrides"][:4] == [
        "set:THREADS=2",
        "set:CLUSTER_CUTOFF=3",
        "set:PARTITION=min_leaf",
        "set:MIN_LEAF=500",
    ]
    args = build_args(env)
    assert args[args.index("-C") + 1] == "3"
    assert args[args.index("--min-leaf") + 1] == "500" and "-s" not in args and "-L" not in args
    assert args[args.index("--tree-args") + 1] == "-gamma -fastexp 2"
    assert args[args.index("-J") + 1] == "4"


def test_a_build_that_loses_sequences_is_not_published(env, monkeypatch):
    monkeypatch.setenv("DROP_LAST", "1")
    result = run(env)
    assert result.returncode != 0
    assert f"only {N_SEQS - 1} of {N_SEQS} sequences" in result.stderr


def test_a_bad_taxonomy_stops_before_the_build(env):
    work = env["tmp"] / "work"
    assert run(env, "--steps", "fetch").returncode == 0
    taxonomy = work / "input" / "lca_taxonomy.txt"
    taxonomy.write_text(taxonomy.read_text().replace(";S0\n", "\n", 1))
    result = run(env, "--steps", "build")
    assert result.returncode != 0 and "every lineage needs 7 ranks" in result.stderr
    assert not (env["tmp"] / "build.args").exists(), "the build never started"


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


def test_a_marker_config_is_picked_up_and_reaches_the_build(env):
    commit(
        env["tronko"], f"scripts/tronko_build.d/{MARKER}.conf",
        "# this marker partitions by leaf count\nPARTITION=min_leaf\nMIN_LEAF = 321\n"
        'FAMSA_ARGS="--refine-mode on"\n',
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
    args = build_args(env)
    assert args[args.index("--min-leaf") + 1] == "321"
    assert args[args.index("--famsa-args") + 1] == "--refine-mode on"


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
     ("FASTTREE_ARGS=-nt", "unknown setting FASTTREE_ARGS"),
     ("CLUSTER_CUTOFF=lots", "CLUSTER_CUTOFF must be a whole number"),
     ("SP_THRESHOLD=.1", "SP_THRESHOLD must be a number like 0.1"),
     ("CLUSTER_JOBS=0", "CLUSTER_JOBS must be at least 1"),
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


REAL_TOOLS = (
    "famsa",
    "nw_reroot",
    "bwa",
    "bwa-mem3",
    "raxmlHPC-PTHREADS",
    "fasta2phyml.pl",
    "make",
    "gcc",
)


@pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in REAL_TOOLS),
    reason="needs " + ", ".join(REAL_TOOLS) + " on PATH",
)
def test_the_real_build_makes_a_database_from_a_small_input(tmp_path, monkeypatch):
    """release_tronko.sh's build step with the real build-tronko-db.sh and
    tronko-build, compiled from a clone of this checkout's HEAD."""
    example = REPO / "tronko-build" / "example_datasets" / "single_tree"
    records = (example / "Charadriiformes.fasta").read_text().split(">")[1:41]
    ids = [record.split("\n", 1)[0] for record in records]
    lineage = dict(
        line.rstrip("\n").split("\t", 1)
        for line in (example / "Charadriiformes_taxonomy.txt").read_text().splitlines()
    )
    work = tmp_path / "work"
    (work / "input").mkdir(parents=True)
    (work / "input" / "lca.fasta").write_text("".join(">" + record for record in records))
    (work / "input" / "lca_taxonomy.txt").write_text(
        "".join(f"{i}\t1\t{lineage[i]}\n" for i in ids)
    )
    write_taxdump(work / "input" / "taxdump", [])  # taxid 1 is the root, which names no rank
    (work / "input" / "inputs.json").write_text(
        json.dumps(
            {"rcrux_tag": "EXAMPLE_v1.0.0", "fasta": "lca.fasta", "taxonomy": "lca_taxonomy.txt"},
            indent=2,
        )
    )
    tronko = tmp_path / "tronko"
    git(tmp_path, "clone", "-q", str(REPO), str(tronko))
    monkeypatch.setenv("PYTHON", sys.executable)
    monkeypatch.delenv("BUILD_TRONKO_DB", raising=False)

    result = subprocess.run(
        ["bash", str(tronko / "scripts" / "release_tronko.sh"), "--marker", "EXAMPLE",
         "--work", str(work), "--steps", "build", "--threads", "2",
         "--set", "TREE_TOOL=raxml", "--set", "TRONKO_BUILD_ARGS=--one-thread-below 0",
         "--set", "ASSIGNMENT_CHECK_READS_PER_QUERY=0"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    out = work / "out"
    assert (out / "EXAMPLE.fasta").read_text().count(">") == len(ids)
    assert (out / "reference_tree.trkb").stat().st_size > 0
    assert (out / "EXAMPLE.fasta.bwt").stat().st_size > 0
    assert (out / "EXAMPLE.fasta.bwt.2bit.64").stat().st_size > 0
    taxonomy = (out / "EXAMPLE_taxonomy.txt").read_text().splitlines()
    assert len(taxonomy) == len(ids)
    assert all(line.split("\t")[1:] == [lineage[line.split("\t")[0]]] for line in taxonomy)
    # one taxid row per taxonomy row, in the same order, and a row per lineage path
    assert (out / "EXAMPLE_taxids.txt").read_text().splitlines() == [
        line.split("\t")[0] + "\t1" for line in taxonomy
    ]
    paths = {";".join(line.split("\t")[1].split(";")[: i + 1]) for line in taxonomy for i in range(7)}
    table = (out / "EXAMPLE_lineage_taxids.txt").read_text().splitlines()
    assert {row.split("\t")[0] for row in table} == paths
    assert all(row.endswith("\t\t") for row in table)  # taxid 1 is the root, which names no rank
    assert "lineage taxids:" in (work / "build.log").read_text()
    build = json.loads((work / "params.json").read_text())["build"]
    assert build["tree_tool"] == "raxml" and build["reference_sequences"] == len(ids)
    assert build["tronko_build_args"] == "--one-thread-below 0"
    assert build["rcrux_tag"] == "EXAMPLE_v1.0.0"


@pytest.mark.parametrize(
    ("settings", "message"),
    [(["ASSIGNMENT_CHECK_SPECIES=0"], "ASSIGNMENT_CHECK_SPECIES must be at least 1"),
     (["ASSIGNMENT_CHECK_READS_PER_QUERY=x"], "ASSIGNMENT_CHECK_READS_PER_QUERY must be a whole number")],
)  # fmt: skip
def test_bad_assignment_check_settings_are_refused(env, settings, message):
    # Not through run(), which turns the check off.
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--print-config",
         "--marker", MARKER, *(a for s in settings for a in ("--set", s))],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode != 0 and message in result.stderr


@pytest.mark.parametrize("reads", ["0", "00"])
def test_zero_reads_skip_the_check_whatever_the_species(env, reads):
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--print-config",
         "--marker", MARKER, "--set", "ASSIGNMENT_CHECK_SPECIES=0",
         "--set", f"ASSIGNMENT_CHECK_READS_PER_QUERY={reads}"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode == 0, result.stderr


def test_the_build_stops_at_preflight_when_the_check_cannot_score(env, tmp_path, monkeypatch):
    """A Python without the toolkit (or with a crashing polars) stops before the build."""
    fake = tmp_path / "python"
    fake.write_text('#!/bin/bash\ncase "$*" in *polars*) exit 139 ;; esac\nexec "' + sys.executable + '" "$@"\n')
    fake.chmod(0o755)
    monkeypatch.setenv("PYTHON", str(fake))
    result = subprocess.run(
        ["bash", str(env["tronko"] / "scripts" / "release_tronko.sh"), "--marker", MARKER,
         "--rcrux-tag", RCRUX_TAG, "--tronko-checkout", str(env["checkout"]),
         "--rcrux-repo", str(env["origin"]), "--work", str(env["tmp"] / "work"), "--threads", "2"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode != 0
    assert "the assignment check needs community-lib and a working polars" in result.stderr
    assert "==> Build" not in result.stdout


@pytest.mark.parametrize(
    ("taxonomy", "taxdump", "message"),
    [("a\t5\tE;P;C;O;F;G;S\n", False, "needs --taxdump DIR"),
     ("a\tE;P;C;O;F;G;S\n", True, "--taxdump needs a taxonomy with taxids"),
     ("a\tE\nb\t5\tE\n", False, "must have 2 or 3 tab-separated columns")],
)  # fmt: skip
def test_build_script_checks_its_taxonomy_and_taxdump_first(tmp_path, taxonomy, taxdump, message):
    """build-tronko-db.sh stops on these before it looks for its tools."""
    (tmp_path / "in.fasta").write_text(">a\nACGT\n")
    (tmp_path / "tax.txt").write_text(taxonomy)
    extra = ["--taxdump", str(write_taxdump(tmp_path / "taxdump", []))] if taxdump else []
    result = subprocess.run(
        ["bash", str(REPO / "build-tronko-db.sh"), "-f", str(tmp_path / "in.fasta"),
         "-t", str(tmp_path / "tax.txt"), "-o", str(tmp_path / "out"), *extra],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert result.returncode == 1 and message in result.stderr
