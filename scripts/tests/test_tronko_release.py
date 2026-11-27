import datetime
import json
import subprocess
from pathlib import Path

import pytest

import lineage_taxids as lt
import tronko_release as tr
from conftest import lineage_nodes, write_taxdump

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
    # rcrux-py tracks the NCBI taxdump its lineages came from in data/.
    write_taxdump(origin / "data" / "taxdump", lineage_nodes(100, "E;P;C;O;F;G;S"))
    dvc(origin, "add", "data")
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
    assert set(inputs["input_files"]) == {
        "curated/lca.fasta", "curated/lca_taxonomy.txt",
        *(f"data/taxdump/{name}" for name in lt.TAXDUMP_FILES),
    }  # fmt: skip
    assert (work / "taxdump" / "nodes.dmp").is_file()
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
    monkeypatch.setattr(tr, "PUSH_RETRY_SECONDS", 0)
    with pytest.raises(RuntimeError, match="no credentials"):
        tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    monkeypatch.setattr(tr, "run", real_run)
    assert git(checkout, "rev-parse", "--abbrev-ref", "HEAD") == "main"
    leftovers = git(checkout, "status", "--porcelain", "--ignored").splitlines()
    assert [line for line in leftovers if not line.startswith("!! .dvc/")] == []
    assert git(checkout, "branch", "--list", f"publish/{MARKER}-tronko-*") == ""
    meta = tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert meta["tronko_version"] == "v1.0.0"


def test_publish_pushes_again_after_a_failed_upload(tmp_path, rcrux_repo, monkeypatch):
    inputs, checkout, log, build = _publish_ready(tmp_path, rcrux_repo)
    real_run = tr.run
    pushes = []

    def flaky_push(command, cwd=None):
        if command[:2] == ["dvc", "push"]:
            pushes.append(command)
            if len(pushes) == 1:
                raise RuntimeError("An error occurred () when calling the CompleteMultipartUpload operation")
        return real_run(command, cwd)

    monkeypatch.setattr(tr, "run", flaky_push)
    monkeypatch.setattr(tr, "PUSH_RETRY_SECONDS", 0)
    meta = tr.publish(MARKER, build, inputs, log, write_params(tmp_path), None, checkout)
    assert len(pushes) == 2
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
    # named for tronko: the rcrux release beside it has its own version
    assert "tronko_release_version: v1.0.0" in pointer
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
    (build / "partition0_taxonomy.txt").write_text("a\tE;P;C;O;F;G;S\nb\tE;P;C;O;F;G;T\n")
    (build / "RAxML_bestTree.partition0.reroot").write_text("(a:0.1,b:0.2);\n")
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
    assert "1 trees; sequences per tree: 2 / 2 / 2" in report
    assert "| 1 to 9 | 1 | 2 |" in report
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


def test_release_report_shows_the_assignment_check_and_flags_false_calls(tmp_path):
    build = write_build(tmp_path / "build")
    log = tmp_path / "build.log"
    log.write_text("")

    def tier(reads, fp, f1):
        return {"reads": reads, "fp": fp, "precision": f1, "recall": f1, "f1": f1}

    check = {
        "settings": {
            "communities": 6, "richness": 150, "reads_per_query": 20, "lca_cutoff": 10,
            "data_pipelines_commit": "1026bc049" + "0" * 31,
        },  # fmt: skip
        "communities": [
            {"removed_records": n, "audit": {"unplaced_orphans": o}} for n, o in ((30, 0), (45, 2))
        ],
        "summary": {
            "illumina_typical": {"normal": tier(1000, 20, 0.9), "species_holdout": tier(500, 100, 0.4)},
            "edna_damage": {"normal": tier(1000, 50, 0.7), "species_holdout": tier(500, 120, 0.3)},
        },
    }
    (build / tr.ASSIGNMENT_CHECK_FILENAME).write_text(json.dumps(check))
    assert tr.ASSIGNMENT_CHECK_FILENAME in tr.published_files(build, MARKER)
    report = tr.release_report(_report_metadata(tmp_path), build, log, None)
    assert report.startswith("**Checks failed:**")
    # only the normal tier of illumina_typical decides: held-out tiers expect false calls
    assert "20 false calls on 1,000 reads of species the database holds" in report
    assert report.count("false calls on") == 1
    assert "| normal | 1,000 | 20 | 0.900 | 0.900 | 0.900 |" in report
    assert "| species held out | 500 | 100 | 0.400 | 0.400 | 0.400 |" in report
    assert "**edna_damage**" in report
    assert "community-lib from edna-explorer-data-pipelines `1026bc049000`." in report
    assert "Held-out records removed per community: 30 to 45." in report
    assert "Kept records not placed back in a tree per community: 0 to 2." in report

def test_fetch_from_a_pointer_reads_the_remote_by_hash_alone(tmp_path, rcrux_repo):
    """What the release VM does: the tag's commit, curated.dvc and remote settings
    from files, the data from the remote, no access to the rcrux repository."""
    origin, tag = rcrux_repo["origin"], f"{MARKER}_v1.0.0"
    pointer = tmp_path / "pointer"
    pointer.mkdir()
    (pointer / "commit").write_text(git(origin, "rev-parse", f"{tag}^{{commit}}") + "\n")
    (pointer / "curated.dvc").write_text(git(origin, "show", f"{tag}:databases/{MARKER}/curated.dvc"))
    (pointer / "data.dvc").write_text(git(origin, "show", f"{tag}:data.dvc"))
    (pointer / "dvc_config").write_text(git(origin, "show", f"{tag}:.dvc/config"))

    by_git = tr.fetch(MARKER, tag, tmp_path / "git", rcrux_repo=str(origin))
    by_pointer = tr.fetch(MARKER, tag, tmp_path / "pointer-out", rcrux_repo=str(origin), pointer_dir=pointer)
    assert by_pointer == by_git
    for name in ("lca.fasta", "taxdump/nodes.dmp"):
        assert (tmp_path / "pointer-out" / name).read_bytes() == (tmp_path / "git" / name).read_bytes()

    # A file that is not the one the listing names is refused.
    listings = [json.loads(p.read_text()) for p in rcrux_repo["remote"].rglob("*.dir")]
    listing = next(entries for entries in listings if any(e["relpath"] == "lca.fasta" for e in entries))
    lca = next(e for e in listing if e["relpath"] == "lca.fasta")
    stored = rcrux_repo["remote"] / "files" / "md5" / lca["md5"][:2] / lca["md5"][2:]
    stored.chmod(0o644)  # DVC stores objects read-only
    stored.write_text(">x\nAAAA\n")
    with pytest.raises(ValueError, match="is not md5"):
        tr.fetch(MARKER, tag, tmp_path / "bad", rcrux_repo=str(origin), pointer_dir=pointer)


def test_tree_report_counts_clades_against_the_taxonomy(tmp_path):
    # Tree 0: S is one clade, T is split by U, and genus G is the whole tree.
    # Tree 1: S again, alone. Family F is in both trees.
    trees = {
        "0": ("(((s1:0.1,s2:0.1)90:0.2,(t1:0.1,u1:0.3):0.1):0.05,t2:0.4);\n",
              {"s1": "S", "s2": "S", "t1": "T", "u1": "U", "t2": "T"}),
        "1": ("(s3:0.1,s4:0.5);\n", {"s3": "S", "s4": "S"}),
    }
    for number, (newick, species) in trees.items():
        (tmp_path / f"RAxML_bestTree.partition{number}.reroot").write_text(newick)
        (tmp_path / f"partition{number}_taxonomy.txt").write_text(
            "".join(f"{leaf}\tE;P;C;O;F;G;{name}\n" for leaf, name in species.items())
        )
    (tmp_path / "final_partitions.txt").write_text("0\n1\n")
    report = "\n".join(tr.tree_report(tmp_path))
    assert "2 trees; sequences per tree: 2 / 3.5 / 5" in report
    assert "| species | 2 of 3 (67%) | 1 of 3 (33%) | 1 of 2 (50%) |" in report
    assert "| genus | 2 of 2 (100%) | 1 of 1 (100%) | 2 of 2 (100%) |" in report
    assert "- T: 2 sequences in tree 0, whose smallest clade holding them has 5" in report


def test_clades_reads_inner_labels_and_trailing_newlines():
    leaves, spans, lengths = tr._clades("((a:1,b:2)90:0.5,c:1);\n")
    assert leaves == ["a", "b", "c"]
    assert spans == {(0, 0), (1, 1), (2, 2), (0, 1), (0, 2)}
    assert lengths == [1, 2, 0.5, 1]


def _taxid_tree(tmp_path: Path) -> Path:
    """Eukaryota..Gen with species a and b, a subspecies of a, subgenus Gen holding
    species c, a species of order O with no family or genus, and a homonym genus
    Gen under Bacteria; taxid 99 was merged into b."""
    nodes = [
        *lineage_nodes(10, "Eukaryota;P;C;O;F;Gen;Gen a"),
        (17, 15, "species", "Gen b"),
        (18, 16, "subspecies", "Gen a subsp. x"),
        (19, 15, "subgenus", "Gen"),
        (20, 19, "species", "Gen c"),
        (21, 13, "species", "O sp. 1"),
        *lineage_nodes(30, "Bacteria;Q;D;R;H;Gen;Gen z"),
    ]
    return write_taxdump(tmp_path / "taxdump", nodes, merged=[(99, 17)])


def test_lineage_taxids_names_every_path_tronko_assign_can_print(tmp_path):
    taxdump = _taxid_tree(tmp_path)
    records = [
        ("r1", "18", "Eukaryota;P;C;O;F;Gen;Gen a"),  # a subspecies taxid: the species is a
        ("r2", "99", "Eukaryota;P;C;O;F;Gen;Gen b"),  # merged into b
        ("r3", "20", "Eukaryota;P;C;O;F;Gen;Gen c"),  # under subgenus Gen: the genus is 15
        ("r4", "16,17,21", "Eukaryota;P;C;O;F;Gen;NA"),  # rCRUX's LCA: 21 names no F or Gen
        ("r5", "36", "Bacteria;Q;D;R;H;Gen;Gen z"),  # the homonym genus
        ("r6", "16", "NA;NA;NA;NA;NA;NA;NA"),
        ("r7", "16", "Eukaryota;P;C;O;F;Gen;Not here"),
    ]
    (tmp_path / "tax.txt").write_text("".join(f"{a}\t{lineage}\n" for a, _, lineage in records))
    (tmp_path / "ids.txt").write_text("".join(f"{a}\t{ids}\n" for a, ids, _ in records))

    rows, conflicts = lt.lineage_taxids(tmp_path / "tax.txt", tmp_path / "ids.txt", taxdump)
    table = {path: (rank, taxid) for path, rank, taxid in rows}
    assert conflicts == 0
    assert table["Eukaryota"] == ("domain", "10")
    assert table["Eukaryota;P;C;O;F;Gen"] == ("genus", "15")
    assert table["Eukaryota;P;C;O;F;Gen;Gen a"] == ("species", "16")
    assert table["Eukaryota;P;C;O;F;Gen;Gen b"] == ("species", "17")
    assert table["Eukaryota;P;C;O;F;Gen;Gen c"] == ("species", "20")
    assert table["Eukaryota;P;C;O;F;Gen;NA"] == ("genus", "15")
    assert table["Bacteria;Q;D;R;H;Gen"] == ("genus", "35")
    assert table["NA;NA;NA;NA;NA;NA;NA"] == ("", "")
    assert table["Eukaryota;P;C;O;F;Gen;Not here"] == ("", "")
    # every prefix of every lineage, once
    assert len(rows) == len({";".join(l.split(";")[: i + 1]) for _, _, l in records for i in range(7)})


def test_write_lineage_taxids_warns_of_paths_without_a_taxid(tmp_path):
    taxdump = _taxid_tree(tmp_path)
    build = tmp_path / "build"
    build.mkdir()
    (build / f"{MARKER}_taxonomy.txt").write_text("r1\tEukaryota;P;C;O;F;Gen;Not here\nr2\tNA;NA;NA;NA;NA;NA;NA\n")
    (build / f"{MARKER}_taxids.txt").write_text("r1\t16\nr2\t16\n")
    lines = lt.write_lineage_taxids(build, MARKER, taxdump)
    assert lines[0] == "lineage taxids: 14 paths, 1 without a taxid, 0 with records disagreeing"
    assert lines[1:] == ["WARNING: no NCBI taxid for lineage Eukaryota;P;C;O;F;Gen;Not here"]
    assert "Eukaryota;P;C;O;F;Gen\tgenus\t15\n" in (build / f"{MARKER}_lineage_taxids.txt").read_text()


def test_lineage_taxids_refuses_taxids_out_of_step_with_the_taxonomy(tmp_path):
    (tmp_path / "tax.txt").write_text("a\tE\nb\tE\n")
    (tmp_path / "ids.txt").write_text("b\t1\na\t1\n")
    with pytest.raises(ValueError, match="disagree at a"):
        lt.lineage_taxids(tmp_path / "tax.txt", tmp_path / "ids.txt", _taxid_tree(tmp_path))


def test_add_lineage_taxids_adds_the_table_to_a_published_release(tmp_path, rcrux_repo):
    work = tmp_path / "work"
    tr.fetch(MARKER, f"{MARKER}_v1.0.0", work, rcrux_repo=str(rcrux_repo["origin"]))
    checkout = tmp_path / "tronko-checkout"
    git(tmp_path, "clone", "-q", str(rcrux_repo["tronko"]), str(checkout))
    build = write_build(tmp_path / "build")
    (build / f"{MARKER}_taxids.txt").write_text("a\t106\nb\t105\n")
    log = tmp_path / "build.log"
    log.write_text("tronko-build -y -u 0.1\n")
    tr.publish(MARKER, build, work / "inputs.json", log, write_params(tmp_path), "v1.0.0", checkout)
    pointer = checkout / "databases" / MARKER / "tronko.dvc"
    before = pointer.read_text()
    subprocess.run(["rm", "-rf", str(checkout / "databases" / MARKER / "tronko"), str(checkout / ".dvc" / "cache")], check=True)

    lines = tr.add_lineage_taxids(MARKER, checkout, tmp_path / "add", rcrux_repo=str(rcrux_repo["origin"]))
    assert lines == ["lineage taxids: 8 paths, 1 without a taxid, 0 with records disagreeing",
                     "WARNING: no NCBI taxid for lineage E;P;C;O;F;G;T"]  # fmt: skip
    after = pointer.read_text()
    assert after != before
    assert tr.read_pointer_version(after) == "v1.0.0" and f"rcrux_tag: {MARKER}_v1.0.0" in after
    assert git(checkout, "status", "--porcelain") == ""
    assert git(checkout, "log", "-1", "--format=%s") == f"Release: add {MARKER}_lineage_taxids.txt"
    released = checkout / "databases" / MARKER / "tronko"
    meta = json.loads((released / "metadata.json").read_text())
    assert f"{MARKER}_lineage_taxids.txt" in meta["files"]
    assert "data/taxdump/nodes.dmp" in meta["input_files"]
    assert "E;P;C;O;F;G;S\tspecies\t106\n" in (released / f"{MARKER}_lineage_taxids.txt").read_text()
    # what was pushed is what a fresh pull gets
    subprocess.run(["rm", "-rf", str(released), str(checkout / ".dvc" / "cache")], check=True)
    dvc(checkout, "pull", "-r", "tronko", f"databases/{MARKER}/tronko.dvc")
    assert (released / f"{MARKER}_lineage_taxids.txt").is_file()


def test_release_report_counts_lineage_paths_with_a_taxid(tmp_path):
    build = write_build(tmp_path / "build")
    (build / f"{MARKER}_lineage_taxids.txt").write_text("E\tdomain\t1\nE;P\tphylum\t2\nNA\t\t\n")
    log = tmp_path / "build.log"
    log.write_text("x\n")
    meta = {"marker": MARKER, "rcrux_tag": "t", "input_tier": "lca", "rcrux_commit": "c" * 12,
            "tronko_sha": "a" * 12, "row_counts": {"fasta_records": 2, "taxonomy_rows": 2},
            "build_params": {}, "files": {}}  # fmt: skip
    assert "- lineage paths with an NCBI taxid: 2 of 3" in tr.release_report(meta, build, log, None)
