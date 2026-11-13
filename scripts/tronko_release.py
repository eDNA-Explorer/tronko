#!/usr/bin/env python3
"""
Version a tronko reference build against the rcrux-py release it was built from.

The rcrux build is versioned by a git tag ({marker}_v1.0.0) whose .dvc pointers
name the exact curated sequences. This script gives the tronko build the same
treatment, with DVC as the only storage on both sides:

1. fetch    Read the curated FASTA + taxonomy from the rcrux DVC remote at the
            rcrux tag (DVCFileSystem, remote "rcrux"), and record which bytes
            were read in inputs.json.
2. (build)  Run tronko-build however the marker needs; this script does not
            build trees.
3. publish  Check the eight files tronko-assign needs, write metadata.json
            beside them with the build log, then in an rcrux-py checkout:
            dvc add databases/{marker}/tronko, dvc push -r rcrux, and commit
            tronko.dvc on a publish branch.
4. tag      After that branch merges, tag {marker}_tronko_v{version} on the
            commit that carries tronko.dvc and push the tag.

Usage:
    python tronko_release.py fetch --marker MARKER --rcrux-tag TAG --input-tier species --out DIR
    python tronko_release.py publish --marker MARKER --build-dir DIR --inputs DIR/inputs.json \
        --build-log build.log --tronko-version v1.0.0 --rcrux-checkout PATH
    python tronko_release.py tag --marker MARKER --tronko-version v1.0.0 --rcrux-checkout PATH
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

SCHEMA_VERSION = "1.0"
DEFAULT_RCRUX_REPO = "https://github.com/eDNA-Explorer/rcrux-py"
# The only DVC remote this script reads from or pushes to. rcrux-py also
# configures a "gcs" remote; versioning never touches it.
DVC_REMOTE = "rcrux"
INPUT_TIERS = ("species", "lca")
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")
VERSION_PATTERN = re.compile(r"v\d+\.\d+\.\d+")
METADATA_FILENAME = "metadata.json"
INPUTS_FILENAME = "inputs.json"
BUILD_LOG_FILENAME = "build.log"


def required_files(marker: str) -> list[str]:
    """The eight files tronko-assign reads, the same set RefDB ingestion checks."""
    fasta = f"{marker}.fasta"
    return [
        fasta,
        *(fasta + ext for ext in BWA_INDEX_EXTENSIONS),
        "reference_tree.trkb",
        f"{marker}_taxonomy.txt",
    ]


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_fasta_records(path: Path) -> int:
    with open(path, "rb") as handle:
        return sum(1 for line in handle if line.startswith(b">"))


def count_lines(path: Path) -> int:
    with open(path, "rb") as handle:
        return sum(1 for line in handle if line.strip())


def run(command: list[str], cwd: Path | None = None) -> str:
    """Run a command, returning stdout; raise with its stderr on failure."""
    result = subprocess.run(command, check=False, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{' '.join(command)} failed ({result.returncode}): {result.stderr.strip()}"
        )
    return result.stdout.strip()


def check_version(version: str) -> str:
    if not VERSION_PATTERN.fullmatch(version):
        raise ValueError(f"Version must look like v1.0.0, got {version!r}")
    return version


def tronko_tag(marker: str, version: str) -> str:
    return f"{marker}_tronko_{check_version(version)}"


# ---------------------------------------------------------------------------
# fetch
# ---------------------------------------------------------------------------


def dvc_out_md5(dvc_file_text: str) -> str:
    """The md5 a .dvc pointer records for its single output."""
    outs = (yaml.safe_load(dvc_file_text) or {}).get("outs") or []
    if len(outs) != 1 or not outs[0].get("md5"):
        raise ValueError("Expected a .dvc file with exactly one output carrying an md5")
    return outs[0]["md5"]


def curated_input_names(listing: list[str], marker: str, tier: str) -> tuple[str, str]:
    """Pick the tier's FASTA and taxonomy from a curated/ listing.

    rcrux writes unprefixed names (species.fasta); its legacy layout prefixes
    them with the marker (MiFish_species.fasta). Exactly one layout must match.
    """
    names = set(listing)
    candidates = [
        (f"{tier}.fasta", f"{tier}_taxonomy.txt"),
        (f"{marker}_{tier}.fasta", f"{marker}_{tier}_taxonomy.txt"),
    ]
    found = [pair for pair in candidates if pair[0] in names and pair[1] in names]
    if len(found) != 1:
        raise ValueError(
            f"Expected one {tier} FASTA + taxonomy pair in curated/, found "
            f"{len(found)} among: {', '.join(sorted(names))}"
        )
    return found[0]


def resolve_tag_commit(repo: str, tag: str) -> str:
    """The commit a tag points to, for annotated and lightweight tags alike."""
    out = run(["git", "ls-remote", repo, f"refs/tags/{tag}", f"refs/tags/{tag}^{{}}"])
    refs = dict(line.split("\t")[::-1] for line in out.splitlines() if line)
    commit = refs.get(f"refs/tags/{tag}^{{}}") or refs.get(f"refs/tags/{tag}")
    if not commit:
        raise ValueError(f"Tag {tag} not found in {repo}")
    return commit


def fetch(
    marker: str,
    rcrux_tag: str,
    input_tier: str,
    out_dir: Path,
    rcrux_repo: str = DEFAULT_RCRUX_REPO,
) -> dict:
    """Download the curated inputs at an rcrux tag and write inputs.json."""
    from dvc.api import DVCFileSystem

    if input_tier not in INPUT_TIERS:
        raise ValueError(f"input tier must be one of {INPUT_TIERS}, got {input_tier!r}")
    commit = resolve_tag_commit(rcrux_repo, rcrux_tag)
    # Pin the commit, not the tag name, so a moved tag cannot change what is read.
    fs = DVCFileSystem(rcrux_repo, rev=commit, remote=DVC_REMOTE)
    marker_dir = f"databases/{marker}"
    with fs.open(f"{marker_dir}/curated.dvc", "r") as handle:
        curated_md5 = dvc_out_md5(handle.read())

    listing = [Path(p).name for p in fs.ls(f"{marker_dir}/curated", detail=False)]
    fasta_name, taxonomy_name = curated_input_names(listing, marker, input_tier)

    out_dir.mkdir(parents=True, exist_ok=True)
    input_files = {}
    for name in (fasta_name, taxonomy_name):
        local = out_dir / name
        fs.get_file(f"{marker_dir}/curated/{name}", str(local))
        if local.stat().st_size == 0:
            raise ValueError(f"curated/{name} at {rcrux_tag} is empty")
        input_files[f"curated/{name}"] = file_digest(local, "md5")

    inputs = {
        "schema_version": SCHEMA_VERSION,
        "marker": marker,
        "rcrux_repo": rcrux_repo,
        "rcrux_tag": rcrux_tag,
        "rcrux_commit": commit,
        "dvc_remote": DVC_REMOTE,
        "input_tier": input_tier,
        "input_md5": curated_md5,
        "input_files": input_files,
        "fasta": fasta_name,
        "taxonomy": taxonomy_name,
    }
    (out_dir / INPUTS_FILENAME).write_text(json.dumps(inputs, indent=2) + "\n")
    return inputs


# ---------------------------------------------------------------------------
# publish
# ---------------------------------------------------------------------------


def validate_build(build_dir: Path, marker: str) -> None:
    problems = []
    for name in required_files(marker):
        path = build_dir / name
        if not path.is_file():
            problems.append(f"{name} missing")
        elif path.stat().st_size == 0:
            problems.append(f"{name} empty")
    if problems:
        raise ValueError(f"Build in {build_dir} is incomplete: " + "; ".join(problems))


def tronko_sha() -> tuple[str, bool]:
    """The tronko commit that built this, and whether it was pinned explicitly.

    TRONKO_REF is the Dockerfile ARG the build image is pinned with; outside
    that image, fall back to this checkout's HEAD and refuse local changes.
    """
    pinned = os.environ.get("TRONKO_REF")
    if pinned:
        return pinned, True
    repo = Path(__file__).resolve().parent.parent
    if run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=repo):
        raise RuntimeError(
            "tronko checkout has uncommitted changes; set TRONKO_REF to the commit "
            "that built the database, or commit first"
        )
    return run(["git", "rev-parse", "HEAD"], cwd=repo), False


def build_metadata(
    marker: str,
    version: str,
    build_dir: Path,
    inputs: dict,
    sha: str,
    build_params: dict | None,
    now: datetime.datetime,
) -> dict:
    if inputs.get("marker") != marker:
        raise ValueError(f"inputs.json is for {inputs.get('marker')!r}, not {marker!r}")
    files = {name: file_digest(build_dir / name, "sha256") for name in required_files(marker)}
    tree = build_dir / "reference_tree.trkb"
    return {
        "schema_version": SCHEMA_VERSION,
        "marker": marker,
        "tronko_version": check_version(version),
        "tronko_tag": tronko_tag(marker, version),
        "tronko_sha": sha,
        "rcrux_tag": inputs["rcrux_tag"],
        "rcrux_commit": inputs["rcrux_commit"],
        "input_tier": inputs["input_tier"],
        "input_md5": inputs["input_md5"],
        "input_files": inputs["input_files"],
        "files": files,
        "row_counts": {
            "fasta_records": count_fasta_records(build_dir / f"{marker}.fasta"),
            "taxonomy_rows": count_lines(build_dir / f"{marker}_taxonomy.txt"),
        },
        "build_params": build_params or {},
        "built_at": datetime.datetime.fromtimestamp(
            tree.stat().st_mtime, tz=datetime.timezone.utc
        ).isoformat(),
        "published_at": now.isoformat(),
    }


def publish(
    marker: str,
    build_dir: Path,
    inputs_path: Path,
    build_log: Path,
    version: str,
    rcrux_checkout: Path,
    build_params_path: Path | None = None,
    replace: bool = False,
    push_branch: bool = False,
    dry_run: bool = False,
) -> dict:
    """Track a finished build in rcrux-py's DVC and commit it on a branch."""
    check_version(version)
    validate_build(build_dir, marker)
    if not build_log.is_file() or build_log.stat().st_size == 0:
        raise ValueError(f"Build log {build_log} is missing or empty")
    inputs = json.loads(inputs_path.read_text())
    params = json.loads(build_params_path.read_text()) if build_params_path else None
    sha, _ = tronko_sha()
    metadata = build_metadata(
        marker,
        version,
        build_dir,
        inputs,
        sha,
        params,
        datetime.datetime.now(datetime.timezone.utc),
    )

    if run(["git", "status", "--porcelain"], cwd=rcrux_checkout):
        raise RuntimeError(f"{rcrux_checkout} has uncommitted changes")
    marker_dir = rcrux_checkout / "databases" / marker
    if not marker_dir.is_dir():
        raise ValueError(f"{marker_dir} does not exist; publish the rcrux build first")
    dest = marker_dir / "tronko"
    if dest.exists() and not replace:
        raise ValueError(f"{dest} already exists; pass --replace for a new version")
    tag = metadata["tronko_tag"]
    if run(["git", "tag", "--list", tag], cwd=rcrux_checkout):
        raise ValueError(f"Tag {tag} already exists in {rcrux_checkout}")
    branch = f"publish/{marker}-tronko-{version}"
    if dry_run:
        print(json.dumps(metadata, indent=2))
        print(f"would track {dest} with DVC, push to {DVC_REMOTE}, commit on {branch}")
        return metadata

    run(["git", "checkout", "-b", branch], cwd=rcrux_checkout)
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir()
    for name in required_files(marker):
        shutil.copy2(build_dir / name, dest / name)
    shutil.copy2(build_log, dest / BUILD_LOG_FILENAME)
    (dest / METADATA_FILENAME).write_text(json.dumps(metadata, indent=2) + "\n")

    target = f"databases/{marker}/tronko"
    run(["dvc", "add", target], cwd=rcrux_checkout)
    run(["dvc", "push", "-r", DVC_REMOTE, f"{target}.dvc"], cwd=rcrux_checkout)
    run(
        ["git", "add", f"{target}.dvc", f"databases/{marker}/.gitignore"],
        cwd=rcrux_checkout,
    )
    run(
        [
            "git",
            "commit",
            "-m",
            f"Publish {marker} tronko {version} built from {inputs['rcrux_tag']}",
        ],
        cwd=rcrux_checkout,
    )
    if push_branch:
        run(["git", "push", "-u", "origin", branch], cwd=rcrux_checkout)
    return metadata


# ---------------------------------------------------------------------------
# tag
# ---------------------------------------------------------------------------


def tag(
    marker: str,
    version: str,
    rcrux_checkout: Path,
    ref: str = "origin/main",
    push: bool = True,
) -> str:
    """Tag the merged commit that carries this build's tronko.dvc."""
    name = tronko_tag(marker, version)
    if push:
        run(["git", "fetch", "origin", "--tags"], cwd=rcrux_checkout)
    if run(["git", "tag", "--list", name], cwd=rcrux_checkout):
        raise ValueError(f"Tag {name} already exists")
    commit = run(["git", "rev-parse", f"{ref}^{{commit}}"], cwd=rcrux_checkout)
    pointer = f"databases/{marker}/tronko.dvc"
    if subprocess.run(
        ["git", "cat-file", "-e", f"{commit}:{pointer}"], check=False, cwd=rcrux_checkout
    ).returncode:
        raise ValueError(f"{ref} has no {pointer}; merge the publish branch first")
    # The tag must point at this version's tree, not a later publish of the marker.
    from dvc.api import DVCFileSystem

    fs = DVCFileSystem(str(rcrux_checkout), rev=commit, remote=DVC_REMOTE)
    with fs.open(f"databases/{marker}/tronko/{METADATA_FILENAME}", "r") as handle:
        published = json.load(handle).get("tronko_version")
    if published != version:
        raise ValueError(f"{ref} carries {marker} tronko {published}, not {version}")
    run(["git", "tag", name, commit], cwd=rcrux_checkout)
    if push:
        run(["git", "push", "origin", name], cwd=rcrux_checkout)
    return name


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser("fetch", help="Download curated inputs at an rcrux tag")
    p_fetch.add_argument("--marker", required=True)
    p_fetch.add_argument("--rcrux-tag", required=True, help="e.g. CO1_BF3_BR2_v1.0.0")
    p_fetch.add_argument("--input-tier", required=True, choices=INPUT_TIERS)
    p_fetch.add_argument("--out", required=True, type=Path)
    p_fetch.add_argument("--rcrux-repo", default=DEFAULT_RCRUX_REPO)

    p_pub = sub.add_parser("publish", help="Track a finished build in rcrux-py DVC")
    p_pub.add_argument("--marker", required=True)
    p_pub.add_argument("--build-dir", required=True, type=Path)
    p_pub.add_argument("--inputs", required=True, type=Path, help="inputs.json from fetch")
    p_pub.add_argument("--build-log", required=True, type=Path)
    p_pub.add_argument("--build-params", type=Path, help="JSON of the build settings")
    p_pub.add_argument("--tronko-version", required=True, help="e.g. v1.0.0")
    p_pub.add_argument("--rcrux-checkout", required=True, type=Path)
    p_pub.add_argument("--replace", action="store_true", help="Replace an earlier tronko/")
    p_pub.add_argument("--push-branch", action="store_true")
    p_pub.add_argument("--dry-run", action="store_true")

    p_tag = sub.add_parser("tag", help="Tag the merged commit carrying tronko.dvc")
    p_tag.add_argument("--marker", required=True)
    p_tag.add_argument("--tronko-version", required=True)
    p_tag.add_argument("--rcrux-checkout", required=True, type=Path)
    p_tag.add_argument("--ref", default="origin/main")

    args = parser.parse_args(argv)
    try:
        if args.command == "fetch":
            inputs = fetch(args.marker, args.rcrux_tag, args.input_tier, args.out, args.rcrux_repo)
            print(f"Wrote {args.out / INPUTS_FILENAME} (input_md5 {inputs['input_md5']})")
        elif args.command == "publish":
            metadata = publish(
                args.marker,
                args.build_dir,
                args.inputs,
                args.build_log,
                args.tronko_version,
                args.rcrux_checkout,
                args.build_params,
                args.replace,
                args.push_branch,
                args.dry_run,
            )
            if not args.dry_run:
                print(
                    f"Committed {args.marker} tronko {metadata['tronko_version']}; "
                    f"open a PR, then run: tronko_release.py tag ..."
                )
        else:
            print(f"Tagged {tag(args.marker, args.tronko_version, args.rcrux_checkout, args.ref)}")
    except (ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
