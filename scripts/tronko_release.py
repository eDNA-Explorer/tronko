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
            beside them with the build log and the build settings (always
            required), then in an rcrux-py checkout:
            dvc add databases/{marker}/tronko, dvc push -r rcrux, and commit
            tronko.dvc on a publish branch.
            The version is the next after this marker's tronko tags unless
            given, and is stamped into tronko.dvc as meta.release_version.
4. tag      Automatic: when the branch merges, rcrux-py's Tag releases
            workflow tags {marker}_tronko_v{version}. `tag` does the same by
            hand, for a merge the workflow did not tag.

Usage:
    python tronko_release.py fetch --marker MARKER --rcrux-tag TAG --out DIR
    python tronko_release.py publish --marker MARKER --build-dir DIR --inputs DIR/inputs.json \
        --build-log build.log --build-params params.json --rcrux-checkout PATH \
        --tronko-sha <commit that built the tree, or set TRONKO_REF> \
        [--tronko-version v2.0.0 | --bump major|minor|patch]
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
# tronko builds from rcrux's LCA tier: curated/lca.fasta + lca_taxonomy.txt.
INPUT_TIER = "lca"
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")
VERSION_PATTERN = re.compile(r"v\d+\.\d+\.\d+")
METADATA_FILENAME = "metadata.json"
INPUTS_FILENAME = "inputs.json"
BUILD_LOG_FILENAME = "build.log"
BUILD_PARAMS_FILENAME = "build_params.json"


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


# Versions are chosen the way rcrux-py's `rcrux release` chooses them, and the
# chosen one is stamped into tronko.dvc as meta.release_version so rcrux-py's
# Tag releases workflow can tag the merge from git alone.
BUMPS = ("major", "minor", "patch")
FIRST_VERSION = "v1.0.0"
META_KEY = "release_version"


def tronko_versions(tags: list[str], marker: str) -> list[str]:
    """This marker's tronko release versions among ``tags``, exact names only."""
    pattern = re.compile(rf"{re.escape(marker)}_tronko_(v\d+\.\d+\.\d+)")
    return [m.group(1) for t in tags if (m := pattern.fullmatch(t))]


def next_version(existing: list[str], bump: str = "minor") -> str:
    """The version after the highest existing one; v1.0.0 for a first release."""
    if bump not in BUMPS:
        raise ValueError(f"bump must be one of {BUMPS}, got {bump!r}")
    if not existing:
        return FIRST_VERSION
    major, minor, patch = max(tuple(int(n) for n in v[1:].split(".")) for v in existing)
    if bump == "major":
        return f"v{major + 1}.0.0"
    if bump == "minor":
        return f"v{major}.{minor + 1}.0"
    return f"v{major}.{minor}.{patch + 1}"


def known_tags(repo: Path) -> list[str]:
    """Local tags plus origin's, so a release cut elsewhere is not reused."""
    tags = set(run(["git", "tag", "--list"], cwd=repo).splitlines())
    try:
        remote = run(["git", "ls-remote", "--tags", "origin"], cwd=repo)
    except RuntimeError:
        remote = ""
    for line in remote.splitlines():
        tags.add(line.split("\t")[-1].removeprefix("refs/tags/").removesuffix("^{}"))
    return sorted(t for t in tags if t)


def read_pointer_version(text: str) -> str | None:
    meta = (yaml.safe_load(text) or {}).get("meta")
    version = meta.get(META_KEY) if isinstance(meta, dict) else None
    return str(version) if version is not None else None


def stamp_pointer(pointer: Path, version: str, **extra: str) -> None:
    """Record the release version (and any extra keys) in a .dvc pointer's meta."""
    data = yaml.safe_load(pointer.read_text()) or {}
    current = data.get("meta")
    meta = dict(current) if isinstance(current, dict) else {}
    meta.update({META_KEY: check_version(version), **extra})
    data["meta"] = meta
    pointer.write_text(yaml.safe_dump(data, sort_keys=False))


def choose_version(
    rcrux_checkout: Path, marker: str, requested: str | None, bump: str, base: str = "HEAD"
) -> str:
    """The requested version if free, else the next after this marker's tronko tags.

    A version stamped in ``base``'s tronko.dvc but not yet tagged (merged, tag
    pending) counts too, so versions only increase.
    """
    tagged = tronko_versions(known_tags(rcrux_checkout), marker)
    if requested is not None:
        check_version(requested)
        if requested in tagged:
            raise ValueError(f"Tag {tronko_tag(marker, requested)} already exists")
        return requested
    on_base = subprocess.run(
        ["git", "show", f"{base}:databases/{marker}/tronko.dvc"],
        cwd=rcrux_checkout,
        capture_output=True,
        text=True,
        check=False,
    )
    pending = read_pointer_version(on_base.stdout) if on_base.returncode == 0 else None
    return next_version(tagged + ([pending] if pending else []), bump)


# ---------------------------------------------------------------------------
# fetch
# ---------------------------------------------------------------------------


def dvc_out_md5(dvc_file_text: str) -> str:
    """The md5 a .dvc pointer records for its single output."""
    outs = (yaml.safe_load(dvc_file_text) or {}).get("outs") or []
    if len(outs) != 1 or not outs[0].get("md5"):
        raise ValueError("Expected a .dvc file with exactly one output carrying an md5")
    return outs[0]["md5"]


def curated_input_names(listing: list[str], marker: str) -> tuple[str, str]:
    """Pick the LCA FASTA and taxonomy from a curated/ listing.

    rcrux writes unprefixed names (lca.fasta); its legacy layout prefixes
    them with the marker (MiFish_lca.fasta). Exactly one layout must match.
    """
    names = set(listing)
    tier = INPUT_TIER
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
    out_dir: Path,
    rcrux_repo: str = DEFAULT_RCRUX_REPO,
) -> dict:
    """Download the curated inputs at an rcrux tag and write inputs.json."""
    from dvc.api import DVCFileSystem

    commit = resolve_tag_commit(rcrux_repo, rcrux_tag)
    # Pin the commit, not the tag name, so a moved tag cannot change what is read.
    fs = DVCFileSystem(rcrux_repo, rev=commit, remote=DVC_REMOTE)
    marker_dir = f"databases/{marker}"
    with fs.open(f"{marker_dir}/curated.dvc", "r") as handle:
        curated_md5 = dvc_out_md5(handle.read())

    listing = [Path(p).name for p in fs.ls(f"{marker_dir}/curated", detail=False)]
    fasta_name, taxonomy_name = curated_input_names(listing, marker)

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
        "input_tier": INPUT_TIER,
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


SHA_PATTERN = re.compile(r"[0-9a-f]{7,40}")


def resolve_tronko_sha(explicit: str | None) -> str:
    """The tronko commit that built the tree, which must be stated.

    It is never taken from the checkout running publish: the tree may have
    been built from an image or another checkout. TRONKO_REF is the Dockerfile
    ARG the build image is pinned with.
    """
    sha = explicit or os.environ.get("TRONKO_REF")
    if not sha:
        raise ValueError(
            "Say which tronko commit built the tree: --tronko-sha <commit> or TRONKO_REF"
        )
    if not SHA_PATTERN.fullmatch(sha):
        raise ValueError(f"tronko sha must be a git commit (7-40 hex characters), got {sha!r}")
    return sha


def build_metadata(
    marker: str,
    version: str,
    build_dir: Path,
    inputs: dict,
    sha: str,
    build_params: dict,
    now: datetime.datetime,
) -> dict:
    # A release without its build settings cannot be reproduced; never publish one.
    if not isinstance(build_params, dict) or not build_params:
        raise ValueError("Build settings are required: a non-empty JSON object")
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
        "build_params": build_params,
        "built_at": datetime.datetime.fromtimestamp(
            tree.stat().st_mtime, tz=datetime.timezone.utc
        ).isoformat(),
        "published_at": now.isoformat(),
    }


def _branch_exists(repo: Path, branch: str) -> bool:
    return (
        subprocess.run(
            ["git", "rev-parse", "--verify", "--quiet", f"refs/heads/{branch}"],
            cwd=repo,
            capture_output=True,
            check=False,
        ).returncode
        == 0
    )


def _roll_back(repo: Path, start: str, branch: str, marker: str) -> None:
    """Undo a publish that failed after creating its branch: back to ``start``.

    The checkout was clean before publish began, so discarding is safe; a
    tronko/ the start commit tracks is restored from the DVC cache.
    """
    run(["git", "checkout", "-f", start], cwd=repo)
    if _branch_exists(repo, branch):
        run(["git", "branch", "-D", branch], cwd=repo)
    pointer = repo / "databases" / marker / "tronko.dvc"
    if pointer.is_file() and not run(["git", "ls-files", str(pointer.relative_to(repo))], cwd=repo):
        pointer.unlink()
    if pointer.is_file():
        subprocess.run(
            ["dvc", "checkout", str(pointer.relative_to(repo))],
            cwd=repo,
            capture_output=True,
            check=False,
        )
    elif (repo / "databases" / marker / "tronko").exists():
        shutil.rmtree(repo / "databases" / marker / "tronko")


def publish(
    marker: str,
    build_dir: Path,
    inputs_path: Path,
    build_log: Path,
    build_params_path: Path,
    version: str | None,
    rcrux_checkout: Path,
    replace: bool = False,
    push_branch: bool = False,
    dry_run: bool = False,
    bump: str = "minor",
    tronko_sha: str | None = None,
    base_branch: str = "main",
) -> dict:
    """Track a finished build in rcrux-py's DVC and commit it on a branch.

    The branch starts from the latest ``origin/<base_branch>``, whatever the
    checkout had checked out. With no ``version``, the next one after this
    marker's tronko tags is used. If any step after the branch is created
    fails, the checkout is put back as it was, so publish can simply be rerun.
    """
    if version is not None:
        check_version(version)
    validate_build(build_dir, marker)
    if not build_log.is_file() or build_log.stat().st_size == 0:
        raise ValueError(f"Build log {build_log} is missing or empty")
    inputs = json.loads(inputs_path.read_text())
    if not build_params_path.is_file():
        raise ValueError(f"Build settings {build_params_path} not found")
    params = json.loads(build_params_path.read_text())
    sha = resolve_tronko_sha(tronko_sha)

    repo = rcrux_checkout
    # Untracked files are fine: publish commits only the files it names, and a
    # tronko/ left over from an earlier publish is untracked once main is back.
    if run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=repo):
        raise RuntimeError(f"{repo} has uncommitted changes to tracked files")
    run(["git", "fetch", "origin", base_branch], cwd=repo)
    base = f"origin/{base_branch}"
    marker_path = f"databases/{marker}"
    if not run(["git", "ls-tree", "--name-only", base, f"{marker_path}/curated.dvc"], cwd=repo):
        raise ValueError(f"{base} has no {marker_path}/curated.dvc; publish the rcrux build first")
    already_released = run(
        ["git", "ls-tree", "--name-only", base, f"{marker_path}/tronko.dvc"], cwd=repo
    )
    if already_released and not replace:
        raise ValueError(
            f"{base} already has a tronko release of {marker}; pass --replace to publish "
            "the next version"
        )
    version = choose_version(repo, marker, version, bump, base)
    metadata = build_metadata(
        marker,
        version,
        build_dir,
        inputs,
        sha,
        params,
        datetime.datetime.now(datetime.timezone.utc),
    )

    branch = f"publish/{marker}-tronko-{version}"
    if _branch_exists(repo, branch):
        raise ValueError(f"Branch {branch} already exists in {repo}; delete it or publish again")
    dest = repo / marker_path / "tronko"
    if dry_run:
        print(json.dumps(metadata, indent=2))
        print(f"would track {dest} with DVC, push to {DVC_REMOTE}, commit on {branch} from {base}")
        return metadata

    current = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo)
    start = current if current != "HEAD" else run(["git", "rev-parse", "HEAD"], cwd=repo)
    run(["git", "checkout", "-b", branch, base], cwd=repo)
    try:
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        for name in required_files(marker):
            shutil.copy2(build_dir / name, dest / name)
        shutil.copy2(build_log, dest / BUILD_LOG_FILENAME)
        (dest / BUILD_PARAMS_FILENAME).write_text(json.dumps(params, indent=2) + "\n")
        (dest / METADATA_FILENAME).write_text(json.dumps(metadata, indent=2) + "\n")

        target = f"{marker_path}/tronko"
        run(["dvc", "add", target], cwd=repo)
        stamp_pointer(repo / f"{target}.dvc", version, rcrux_tag=inputs["rcrux_tag"])
        run(["dvc", "push", "-r", DVC_REMOTE, f"{target}.dvc"], cwd=repo)
        run(["git", "add", f"{target}.dvc", f"{marker_path}/.gitignore"], cwd=repo)
        run(
            [
                "git",
                "commit",
                "-m",
                f"Publish {marker} tronko {version} built from {inputs['rcrux_tag']}",
            ],
            cwd=repo,
        )
    except Exception:
        _roll_back(repo, start, branch, marker)
        raise
    if push_branch:
        run(["git", "push", "-u", "origin", branch], cwd=repo)
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
    p_fetch.add_argument("--out", required=True, type=Path)
    p_fetch.add_argument("--rcrux-repo", default=DEFAULT_RCRUX_REPO)

    p_pub = sub.add_parser("publish", help="Track a finished build in rcrux-py DVC")
    p_pub.add_argument("--marker", required=True)
    p_pub.add_argument("--build-dir", required=True, type=Path)
    p_pub.add_argument("--inputs", required=True, type=Path, help="inputs.json from fetch")
    p_pub.add_argument("--build-log", required=True, type=Path)
    p_pub.add_argument(
        "--build-params", required=True, type=Path, help="JSON of every build setting used"
    )
    p_pub.add_argument(
        "--tronko-version",
        help="e.g. v2.0.0. Default: the next after this marker's tronko tags (see --bump)",
    )
    p_pub.add_argument("--bump", choices=BUMPS, default="minor")
    p_pub.add_argument(
        "--tronko-sha",
        help="Commit of the tronko that built the tree (default: $TRONKO_REF; one is required)",
    )
    p_pub.add_argument("--base-branch", default="main", help="rcrux-py branch to start from")
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
            inputs = fetch(args.marker, args.rcrux_tag, args.out, args.rcrux_repo)
            print(f"Wrote {args.out / INPUTS_FILENAME} (input_md5 {inputs['input_md5']})")
        elif args.command == "publish":
            metadata = publish(
                args.marker,
                args.build_dir,
                args.inputs,
                args.build_log,
                args.build_params,
                args.tronko_version,
                args.rcrux_checkout,
                args.replace,
                args.push_branch,
                args.dry_run,
                args.bump,
                args.tronko_sha,
                args.base_branch,
            )
            if not args.dry_run:
                print(
                    f"Committed {args.marker} tronko {metadata['tronko_version']}; open a PR. "
                    f"Merging it tags {metadata['tronko_tag']} (rcrux-py Tag releases workflow)."
                )
        else:
            print(f"Tagged {tag(args.marker, args.tronko_version, args.rcrux_checkout, args.ref)}")
    except (ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
