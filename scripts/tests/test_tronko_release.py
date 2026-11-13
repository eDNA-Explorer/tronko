import datetime
import json
import subprocess
from pathlib import Path

import pytest

import tronko_release as tr

MARKER = "CO1_Test"


def git(cwd: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


def dvc(cwd: Path, *args: str) -> None:
    subprocess.run(["dvc", *args], cwd=cwd, check=True, capture_output=True, text=True)


@pytest.fixture(autouse=True)
def git_identity(monkeypatch):
    for key in ("AUTHOR", "COMMITTER"):
        monkeypatch.setenv(f"GIT_{key}_NAME", "test")
        monkeypatch.setenv(f"GIT_{key}_EMAIL", "test@example.com")
    monkeypatch.setenv("TRONKO_REF", "abc1234")


def write_build(build_dir: Path, marker: str = MARKER) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    for name in tr.required_files(marker):
        (build_dir / name).write_text(f"{name}\n")
    (build_dir / f"{marker}.fasta").write_text(">a\nACGT\n>b\nACGA\n")
    (build_dir / f"{marker}_taxonomy.txt").write_text("a\tE;P;C;O;F;G;S\nb\tE;P;C;O;F;G;T\n")
    return build_dir


def write_params(tmp_path: Path) -> Path:
    path = tmp_path / "params.json"
    path.write_text(json.dumps({"sp_threshold": 0.1, "partition_mode": "-y"}))
    return path


def test_required_files_are_the_eight_assign_reads():
    assert tr.required_files("M") == [
        "M.fasta",
        "M.fasta.amb",
        "M.fasta.ann",
        "M.fasta.bwt",
        "M.fasta.pac",
        "M.fasta.sa",
        "reference_tree.trkb",
        "M_taxonomy.txt",
    ]


def test_dvc_out_md5_reads_the_single_output():
    text = "outs:\n- md5: d2710fa6bb72e164c0ebb623c7ea906a.dir\n  path: curated\n"
    assert tr.dvc_out_md5(text) == "d2710fa6bb72e164c0ebb623c7ea906a.dir"
    with pytest.raises(ValueError):
        tr.dvc_out_md5("outs: []\n")


@pytest.mark.parametrize(
    ("listing", "expected"),
    [
        (
            ["lca.fasta", "lca_taxonomy.txt", "species.fasta", "species_taxonomy.txt"],
            ("lca.fasta", "lca_taxonomy.txt"),
        ),
        (
            [f"{MARKER}_lca.fasta", f"{MARKER}_lca_taxonomy.txt"],
            (f"{MARKER}_lca.fasta", f"{MARKER}_lca_taxonomy.txt"),
        ),
    ],
)
def test_curated_input_names_handles_both_layouts(listing, expected):
    assert tr.curated_input_names(listing, MARKER) == expected


@pytest.mark.parametrize(
    "listing",
    [
        ["lca.fasta"],  # taxonomy missing
        ["species.fasta", "species_taxonomy.txt"],  # only the species tier
        [
            "lca.fasta",
            "lca_taxonomy.txt",
            f"{MARKER}_lca.fasta",
            f"{MARKER}_lca_taxonomy.txt",
        ],  # ambiguous
    ],
)
def test_curated_input_names_rejects_missing_or_ambiguous(listing):
    with pytest.raises(ValueError):
        tr.curated_input_names(listing, MARKER)


def test_validate_build_names_missing_and_empty_files(tmp_path):
    build = write_build(tmp_path / "build")
    (build / "reference_tree.trkb").unlink()
    (build / f"{MARKER}.fasta.sa").write_text("")
    with pytest.raises(ValueError, match="fasta.sa empty; reference_tree.trkb missing"):
        tr.validate_build(build, MARKER)


@pytest.mark.parametrize("version", ["1.0.0", "v1.0", "v1.0.0-rc1", "latest"])
def test_versions_must_be_semver_with_v(version):
    with pytest.raises(ValueError):
        tr.tronko_tag(MARKER, version)


def test_build_metadata_records_inputs_outputs_and_tool(tmp_path):
    build = write_build(tmp_path / "build")
    inputs = {
        "marker": MARKER,
        "rcrux_tag": f"{MARKER}_v1.0.0",
        "rcrux_commit": "c0ffee",
        "input_tier": "lca",
        "input_md5": "abc.dir",
        "input_files": {"curated/lca.fasta": "m1"},
    }
    now = datetime.datetime(2026, 11, 13, tzinfo=datetime.timezone.utc)
    meta = tr.build_metadata(MARKER, "v1.0.0", build, inputs, "abc1234", {"-y": True}, now)
    assert meta["tronko_tag"] == f"{MARKER}_tronko_v1.0.0"
    assert meta["tronko_sha"] == "abc1234"
    assert meta["rcrux_tag"] == f"{MARKER}_v1.0.0"
    assert meta["input_md5"] == "abc.dir"
    assert set(meta["files"]) == set(tr.required_files(MARKER))
    assert meta["row_counts"] == {"fasta_records": 2, "taxonomy_rows": 2}
    assert meta["build_params"] == {"-y": True}
    assert meta["published_at"] == now.isoformat()

    with pytest.raises(ValueError, match="not 'Other'"):
        tr.build_metadata(
            "Other", "v1.0.0", write_build(tmp_path / "o", "Other"), inputs, "s", {"-y": 1}, now
        )


@pytest.mark.parametrize("params", [None, {}, [], "-y"])
def test_build_metadata_requires_build_settings(tmp_path, params):
    inputs = {"marker": MARKER}
    now = datetime.datetime(2026, 11, 13, tzinfo=datetime.timezone.utc)
    with pytest.raises(ValueError, match="Build settings are required"):
        tr.build_metadata(MARKER, "v1.0.0", write_build(tmp_path / "b"), inputs, "s", params, now)


@pytest.fixture
def rcrux_repo(tmp_path):
    """An rcrux-py-shaped repo: curated/ tracked by DVC on remote 'rcrux', tagged."""
    origin = tmp_path / "rcrux-origin"
    remote = tmp_path / "remote-rcrux"
    decoy = tmp_path / "remote-gcs"
    for path in (origin, remote, decoy):
        path.mkdir()
    git(origin, "init", "-q", "-b", "main")
    dvc(origin, "init", "-q")
    dvc(origin, "remote", "add", "-d", "rcrux", str(remote))
    dvc(origin, "remote", "add", "gcs", str(decoy))
    curated = origin / "databases" / MARKER / "curated"
    curated.mkdir(parents=True)
    (curated / "species.fasta").write_text(">s1\nACGTACGT\n")
    (curated / "species_taxonomy.txt").write_text("s1\tE;P;C;O;F;G;S\n")
    (curated / "lca.fasta").write_text(">l1\nACGT\n")
    (curated / "lca_taxonomy.txt").write_text("l1\tE;P;C;O;F;G;NA\n")
    dvc(origin, "add", f"databases/{MARKER}/curated")
    git(origin, "add", "-A")
    git(origin, "commit", "-q", "-m", "rcrux build")
    dvc(origin, "push", "-r", "rcrux")
    git(origin, "tag", f"{MARKER}_v1.0.0")
    # Drop the local cache so fetch must read the remote.
    subprocess.run(["rm", "-rf", str(origin / ".dvc" / "cache")], check=True)
    return {"origin": origin, "remote": remote, "decoy": decoy}


def test_fetch_publish_tag_round_trip_uses_only_the_rcrux_remote(tmp_path, rcrux_repo):
    origin, remote, decoy = rcrux_repo["origin"], rcrux_repo["remote"], rcrux_repo["decoy"]

    work = tmp_path / "work"
    inputs = tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(origin))
    assert inputs["rcrux_commit"] == git(origin, "rev-parse", f"{MARKER}_v1.0.0")
    assert inputs["input_md5"].endswith(".dir")
    assert inputs["input_tier"] == "lca"
    assert set(inputs["input_files"]) == {"curated/lca.fasta", "curated/lca_taxonomy.txt"}
    assert (work / "lca.fasta").read_text() == ">l1\nACGT\n"
    assert json.loads((work / "inputs.json").read_text()) == inputs

    checkout = tmp_path / "rcrux-checkout"
    git(tmp_path, "clone", "-q", str(origin), str(checkout))
    build = write_build(tmp_path / "build")
    log = tmp_path / "build.log"
    log.write_text("tronko-build -y -u 0.1\n")
    params = write_params(tmp_path)

    missing = tmp_path / "missing.json"
    with pytest.raises(ValueError, match="not found"):
        tr.publish(MARKER, build, work / "inputs.json", log, missing, "v1.0.0", checkout)

    meta = tr.publish(MARKER, build, work / "inputs.json", log, params, "v1.0.0", checkout)
    assert meta["build_params"] == {"sp_threshold": 0.1, "partition_mode": "-y"}
    published = checkout / "databases" / MARKER / "tronko"
    assert json.loads((published / "build_params.json").read_text()) == meta["build_params"]
    assert git(checkout, "rev-parse", "--abbrev-ref", "HEAD") == f"publish/{MARKER}-tronko-v1.0.0"
    assert git(checkout, "status", "--porcelain") == ""
    assert f"databases/{MARKER}/tronko.dvc" in git(checkout, "show", "--name-only", "--format=")
    assert meta["rcrux_tag"] == f"{MARKER}_v1.0.0"
    assert any(remote.rglob("*")), "build was not pushed to the rcrux remote"
    assert not any(decoy.rglob("*")), "versioning must never touch the gcs remote"

    # A second publish of the same version is refused once the tag exists.
    git(checkout, "checkout", "-q", "main")
    git(checkout, "merge", "-q", "--no-ff", "-m", "merge", f"publish/{MARKER}-tronko-v1.0.0")
    with pytest.raises(ValueError, match="carries .* v1.0.0, not v2.0.0"):
        tr.tag(MARKER, "v2.0.0", checkout, ref="main", push=False)
    assert tr.tag(MARKER, "v1.0.0", checkout, ref="main", push=False) == f"{MARKER}_tronko_v1.0.0"
    assert git(checkout, "rev-parse", f"{MARKER}_tronko_v1.0.0") == git(
        checkout, "rev-parse", "main"
    )
    with pytest.raises(ValueError, match="already exists"):
        tr.tag(MARKER, "v1.0.0", checkout, ref="main", push=False)


def test_publish_refuses_existing_tronko_without_replace(tmp_path, rcrux_repo):
    work = tmp_path / "work"
    tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(rcrux_repo["origin"]))
    checkout = tmp_path / "checkout"
    git(tmp_path, "clone", "-q", str(rcrux_repo["origin"]), str(checkout))
    (checkout / "databases" / MARKER / "tronko").mkdir()
    log = tmp_path / "build.log"
    log.write_text("log\n")
    with pytest.raises(ValueError, match="--replace"):
        tr.publish(
            MARKER,
            write_build(tmp_path / "build"),
            work / "inputs.json",
            log,
            write_params(tmp_path),
            "v1.0.0",
            checkout,
        )


def test_fetch_rejects_unknown_tag(tmp_path, rcrux_repo):
    with pytest.raises(ValueError, match="not found"):
        tr.fetch(MARKER, "nope_v9.9.9", tmp_path / "w", rcrux_repo=str(rcrux_repo["origin"]))


@pytest.mark.parametrize(
    ("existing", "bump", "expected"),
    [
        ([], "minor", "v1.0.0"),
        (["v1.0.0"], "minor", "v1.1.0"),
        (["v1.0.0", "v1.10.0", "v1.2.0"], "minor", "v1.11.0"),
        (["v1.2.3"], "patch", "v1.2.4"),
        (["v1.2.3"], "major", "v2.0.0"),
    ],
)
def test_next_version(existing, bump, expected):
    assert tr.next_version(existing, bump) == expected


def test_tronko_versions_match_their_marker_exactly():
    tags = [
        f"{MARKER}_tronko_v1.0.0",
        f"{MARKER}_v9.0.0",  # rcrux release, not tronko
        f"{MARKER}_v2_tronko_v5.0.0",  # another marker
        f"{MARKER}_tronko_v1.1",  # not a version
    ]
    assert tr.tronko_versions(tags, MARKER) == ["v1.0.0"]


def _publish_ready(tmp_path, rcrux_repo):
    work = tmp_path / "work"
    tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(rcrux_repo["origin"]))
    checkout = tmp_path / "checkout"
    git(tmp_path, "clone", "-q", str(rcrux_repo["origin"]), str(checkout))
    log = tmp_path / "build.log"
    log.write_text("log\n")
    return work / "inputs.json", checkout, log, write_build(tmp_path / "build")


def test_publish_picks_the_version_and_stamps_the_pointer(tmp_path, rcrux_repo):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    meta = tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert meta["tronko_version"] == "v1.0.0"
    pointer = (checkout / "databases" / MARKER / "tronko.dvc").read_text()
    assert tr.read_pointer_version(pointer) == "v1.0.0"
    assert f"rcrux_tag: {MARKER}_v1.0.0" in pointer
    assert pointer in git(checkout, "show", f"HEAD:databases/{MARKER}/tronko.dvc") + "\n"


def test_publish_bumps_past_existing_tronko_tags(tmp_path, rcrux_repo):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    git(checkout, "tag", f"{MARKER}_tronko_v1.0.0")
    git(checkout, "tag", f"{MARKER}_tronko_v1.3.0")
    params = write_params(tmp_path)
    meta = tr.publish(MARKER, build, inputs, log, params, None, checkout, bump="patch")
    assert meta["tronko_version"] == "v1.3.1"
    assert meta["tronko_tag"] == f"{MARKER}_tronko_v1.3.1"


def test_publish_refuses_a_version_already_tagged(tmp_path, rcrux_repo):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    git(checkout, "tag", f"{MARKER}_tronko_v1.0.0")
    with pytest.raises(ValueError, match="already exists"):
        tr.publish(MARKER, build, inputs, log, write_params(tmp_path), "v1.0.0", checkout)
