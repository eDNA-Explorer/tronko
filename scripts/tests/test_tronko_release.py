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


def test_required_files_include_both_aligner_indexes():
    assert tr.required_files("M") == [
        "M.fasta",
        "M.fasta.amb",
        "M.fasta.ann",
        "M.fasta.bwt",
        "M.fasta.pac",
        "M.fasta.sa",
        "M.fasta.bwt.2bit.64",
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
    # The tronko repository, where a build is published: remote 'tronko'.
    tronko = tmp_path / "tronko-origin"
    tronko_remote = tmp_path / "remote-tronko"
    for path in (tronko, tronko_remote):
        path.mkdir()
    git(tronko, "init", "-q", "-b", "main")
    dvc(tronko, "init", "-q")
    dvc(tronko, "remote", "add", "-d", "tronko", str(tronko_remote))
    dvc(tronko, "remote", "add", "gcs", str(decoy))
    git(tronko, "add", "-A")
    git(tronko, "commit", "-q", "-m", "tronko")
    return {
        "origin": origin,
        "remote": remote,
        "decoy": decoy,
        "tronko": tronko,
        "tronko_remote": tronko_remote,
    }


def test_fetch_publish_tag_round_trip_reads_rcrux_and_publishes_to_tronko(tmp_path, rcrux_repo):
    origin, remote, decoy = rcrux_repo["origin"], rcrux_repo["remote"], rcrux_repo["decoy"]
    rcrux_objects = sorted(p for p in remote.rglob("*") if p.is_file())

    work = tmp_path / "work"
    inputs = tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(origin))
    assert inputs["rcrux_commit"] == git(origin, "rev-parse", f"{MARKER}_v1.0.0")
    assert inputs["input_md5"].endswith(".dir")
    assert inputs["input_tier"] == "lca"
    assert set(inputs["input_files"]) == {"curated/lca.fasta", "curated/lca_taxonomy.txt"}
    assert (work / "lca.fasta").read_text() == ">l1\nACGT\n"
    assert json.loads((work / "inputs.json").read_text()) == inputs

    checkout = tmp_path / "tronko-checkout"
    git(tmp_path, "clone", "-q", str(rcrux_repo["tronko"]), str(checkout))
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
    assert any(rcrux_repo["tronko_remote"].rglob("*")), "build was not pushed to tronko's remote"
    assert sorted(p for p in remote.rglob("*") if p.is_file()) == rcrux_objects, (
        "publish must not write to the rcrux remote"
    )
    assert not any(decoy.rglob("*")), "versioning must never touch the gcs remote"

    # A second publish of the same version is refused once the tag exists.
    git(checkout, "checkout", "-q", "main")
    git(checkout, "merge", "-q", "--no-ff", "-m", "merge", f"publish/{MARKER}-tronko-v1.0.0")
    with pytest.raises(ValueError, match="carries .* v1.0.0, not v2.0.0"):
        tr.tag(MARKER, "v2.0.0", checkout, ref="main", push=False)
    assert tr.tag(MARKER, None, checkout, ref="main", push=False) == f"{MARKER}_tronko_v1.0.0"
    assert git(checkout, "rev-parse", f"{MARKER}_tronko_v1.0.0") == git(
        checkout, "rev-parse", "main"
    )
    with pytest.raises(ValueError, match="already exists"):
        tr.tag(MARKER, "v1.0.0", checkout, ref="main", push=False)


def _merge_into_origin(origin: Path, checkout: Path, branch: str) -> None:
    git(checkout, "push", "-q", "origin", branch)
    git(origin, "merge", "-q", "--no-ff", "-m", f"merge {branch}", branch)


def test_publish_refuses_existing_tronko_without_replace(tmp_path, rcrux_repo):
    origin = rcrux_repo["tronko"]
    work = tmp_path / "work"
    tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(rcrux_repo["origin"]))
    checkout = tmp_path / "checkout"
    git(tmp_path, "clone", "-q", str(origin), str(checkout))
    log = tmp_path / "build.log"
    log.write_text("log\n")
    args = (MARKER, write_build(tmp_path / "build"), work / "inputs.json", log)
    tr.publish(*args, write_params(tmp_path), None, checkout)
    _merge_into_origin(origin, checkout, f"publish/{MARKER}-tronko-v1.0.0")
    git(checkout, "checkout", "-q", "main")
    with pytest.raises(ValueError, match="--replace"):
        tr.publish(*args, write_params(tmp_path), None, checkout)
    # With --replace, the next version; main's stamped v1.0.0 counts though untagged.
    meta = tr.publish(*args, write_params(tmp_path), None, checkout, replace=True)
    assert meta["tronko_version"] == "v1.1.0"


def test_publish_branches_from_the_latest_origin_main(tmp_path, rcrux_repo):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    origin = rcrux_repo["tronko"]
    (origin / "NEWS.md").write_text("merged after the clone\n")
    git(origin, "add", "NEWS.md")
    git(origin, "commit", "-q", "-m", "later work on main")
    git(checkout, "checkout", "-q", "-b", "someone-elses-branch")
    tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert (checkout / "NEWS.md").exists()
    assert git(checkout, "rev-parse", "HEAD~1") == git(origin, "rev-parse", "main")


def test_a_failed_publish_rolls_back_and_retries_at_the_same_version(
    tmp_path, rcrux_repo, monkeypatch
):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    real_run = tr.run

    def failing_push(command, cwd=None):
        if command[:2] == ["dvc", "push"]:
            raise RuntimeError("dvc push failed: no credentials")
        return real_run(command, cwd)

    monkeypatch.setattr(tr, "run", failing_push)
    with pytest.raises(RuntimeError, match="no credentials"):
        tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    monkeypatch.setattr(tr, "run", real_run)
    assert git(checkout, "rev-parse", "--abbrev-ref", "HEAD") == "main"
    leftovers = git(checkout, "status", "--porcelain", "--ignored").splitlines()
    assert [line for line in leftovers if not line.startswith("!! .dvc/")] == []
    assert git(checkout, "branch", "--list", f"publish/{MARKER}-tronko-*") == ""
    meta = tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert meta["tronko_version"] == "v1.0.0"


@pytest.mark.parametrize(("env", "arg"), [(None, None), (None, "not-a-sha"), ("main", None)])
def test_publish_requires_the_tronko_commit(tmp_path, rcrux_repo, monkeypatch, env, arg):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    monkeypatch.delenv("TRONKO_REF", raising=False)
    if env:
        monkeypatch.setenv("TRONKO_REF", env)
    with pytest.raises(ValueError, match="tronko"):
        tr.publish(
            MARKER, build, inputs, log, write_params(tmp_path), None, checkout, tronko_sha=arg
        )


def test_publish_records_an_explicit_tronko_commit(tmp_path, rcrux_repo, monkeypatch):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    monkeypatch.delenv("TRONKO_REF", raising=False)
    meta = tr.publish(
        MARKER, build, inputs, log, write_params(tmp_path), None, checkout, tronko_sha="cbe20c3"
    )
    assert meta["tronko_sha"] == "cbe20c3"


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
    git(tmp_path, "clone", "-q", str(rcrux_repo["tronko"]), str(checkout))
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


def test_tag_changed_tags_a_merged_release_once(tmp_path, rcrux_repo):
    """What the Tag tronko releases workflow runs on a push to main."""
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    before = git(checkout, "rev-parse", "main")
    tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    git(checkout, "checkout", "-q", "main")
    git(checkout, "merge", "-q", "--no-ff", "-m", "merge", f"publish/{MARKER}-tronko-v1.0.0")
    assert tr.tag_changed(checkout, before, "main", push=False) == [f"{MARKER}_tronko_v1.0.0"]
    assert git(checkout, "rev-parse", f"{MARKER}_tronko_v1.0.0") == git(
        checkout, "rev-parse", "main"
    )
    assert tr.tag_changed(checkout, before, "main", push=False) == [], "a rerun tags nothing"


def _report_metadata(tmp_path, taxonomy_rows=2):
    return {
        "marker": MARKER,
        "tronko_version": "v1.1.0",
        "tronko_sha": "a" * 40,
        "rcrux_tag": f"{MARKER}_v1.0.0",
        "rcrux_commit": "b" * 40,
        "input_tier": "lca",
        "files": {f"{MARKER}.fasta": "c" * 64},
        "row_counts": {"fasta_records": 2, "taxonomy_rows": taxonomy_rows},
        "build_params": {
            "overrides": ["set:CLUSTER_JOBS=8"],
            "settings": {"PARTITION": "sum_of_pairs", "SP_THRESHOLD": "0.1"},
            "build": {"input_sequences": 2, "ac_clusters": 1, "n_trees": 3},
        },
    }


def test_release_report_passes_a_consistent_build_and_describes_it(tmp_path):
    build = write_build(tmp_path / "build")
    (build / "final_partitions.txt").write_text("0\n")
    (build / "partition0_taxonomy.txt").write_text("a\tx\nb\ty\n")
    log = tmp_path / "build.log"
    # As tronko-build writes it: score lines cut by its other output, a 0.05 parent
    # split into three, one of them under the threshold but too small to split, and
    # an input cluster too small to split at all.
    log.write_text(
        "  Launched cluster 0 (PID 1, 2 seqs)\n"
        "SPscore: 0.050000\nSPscore: 0.200000\nSPscore: 0.300000\n"
        "SPscore:   writeNewick start (root=0, numspec=11)\n080000\n"
        "SPscore: -  createNewRoots: partition too small, returning early\n"
        "  createNewRoots: partition too small, returning early\n"
        "materialising unsplit cluster 1 as partition9\n"
    )
    report = tr.release_report(_report_metadata(tmp_path), build, log, None)
    assert report.startswith("All checks passed")
    assert "settings changed from the recipe: CLUSTER_JOBS=8" in report
    assert "sequences per cluster: 2 / 2 / 2; 1 clusters" in report
    assert "sequences per tree: 2 / 2 / 2; 1 trees" in report
    assert "trees under the 0.1 sum-of-pairs threshold, too small to split further: 1 of 3" in report
    assert "2 species" in report and "sequences with no species name: 0 of 2" in report


def test_release_report_leads_with_failed_checks(tmp_path):
    build = write_build(tmp_path / "build")
    log = tmp_path / "build.log"
    log.write_text("ERROR: 1 cluster(s) failed during Step 2\n")
    previous = _report_metadata(tmp_path) | {"tronko_version": "v1.0.0"}
    previous["row_counts"] = {"fasta_records": 10, "taxonomy_rows": 10}
    report = tr.release_report(_report_metadata(tmp_path, taxonomy_rows=3), build, log, previous)
    assert report.startswith("**Checks failed:**")
    assert "2 FASTA records but 3 taxonomy rows" in report
    assert "the build log has 1 error lines" in report
    assert "more than 10% fewer sequences than v1.0.0" in report
    assert "| sequences | 2 | 10 | -8 |" in report


def test_publish_puts_the_report_in_the_commit_the_pr_is_made_from(tmp_path, rcrux_repo):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert git(checkout, "log", "-1", "--format=%s").startswith(f"Publish {MARKER} tronko v1.0.0")
    assert "## Taxonomy" in git(checkout, "log", "-1", "--format=%b")
