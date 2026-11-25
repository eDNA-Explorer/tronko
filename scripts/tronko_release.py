#!/usr/bin/env python3
"""
Version a tronko reference build against the rcrux-py release it was built from.

The rcrux build is versioned by a git tag in rcrux-py ({marker}_v1.0.0) whose
.dvc pointers name the exact curated sequences. This script gives the tronko
build the same treatment in the tronko repository, with DVC as the only
storage on both sides:

1. fetch    Read the curated FASTA + taxonomy from the rcrux DVC remote at the
            rcrux tag (DVCFileSystem, remote "rcrux"), and record which bytes
            were read in inputs.json.
2. (build)  Run tronko-build however the marker needs; this script does not
            build trees.
3. publish  Check the nine files tronko-assign needs, write metadata.json
            beside them with the build log and the build settings (always
            required), then in a tronko checkout:
            dvc add databases/{marker}/tronko, dvc push -r tronko, and commit
            tronko.dvc on a publish branch.
            The version is the next after this marker's tronko tags unless
            given, and is stamped into tronko.dvc as meta.release_version,
            beside the rcrux tag the build came from.
4. tag      Automatic: when the branch merges, the tronko repository's "Tag
            tronko releases" workflow tags {marker}_tronko_v{version} there.
            `tag` does the same by hand, for a merge the workflow did not tag.

Usage:
    python tronko_release.py fetch --marker MARKER --rcrux-tag TAG --out DIR
    python tronko_release.py publish --marker MARKER --build-dir DIR --inputs DIR/inputs.json \
        --build-log build.log --build-params params.json --tronko-checkout PATH \
        --tronko-sha <commit that built the tree, or set TRONKO_REF> \
        [--tronko-version v2.0.0 | --bump major|minor|patch]
    python tronko_release.py tag --marker MARKER [--tronko-version v1.0.0] [--tronko-checkout PATH]
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

import yaml

SCHEMA_VERSION = "1.0"
DEFAULT_RCRUX_REPO = "https://github.com/eDNA-Explorer/rcrux-py"
# Inputs are read from rcrux-py's "rcrux" DVC remote; the build is pushed to
# the tronko repository's "tronko" remote. Both repositories also configure a
# "gcs" remote; versioning never touches it.
RCRUX_REMOTE = "rcrux"
TRONKO_REMOTE = "tronko"
# tronko builds from rcrux's LCA tier: curated/lca.fasta + lca_taxonomy.txt.
INPUT_TIER = "lca"
# bwa's index, then the one file BWA-MEM3 adds (it shares .amb .ann .pac);
# tronko-assign aligns with BWA-MEM3.
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa", ".bwt.2bit.64")
VERSION_PATTERN = re.compile(r"v\d+\.\d+\.\d+")
METADATA_FILENAME = "metadata.json"
INPUTS_FILENAME = "inputs.json"
BUILD_LOG_FILENAME = "build.log"
BUILD_PARAMS_FILENAME = "build_params.json"
ASSIGNMENT_CHECK_FILENAME = "assignment_check.json"
# Reads drawn from the database itself should not be assigned a wrong name.
ASSIGNMENT_WRONG_WARNING = 0.01


def required_files(marker: str) -> list[str]:
    """The nine files of a database: RefDB ingestion's eight, plus the BWA-MEM3 index."""
    fasta = f"{marker}.fasta"
    return [
        fasta,
        *(fasta + ext for ext in BWA_INDEX_EXTENSIONS),
        "reference_tree.trkb",
        f"{marker}_taxonomy.txt",
    ]


def published_files(build_dir: Path, marker: str) -> list[str]:
    """The required files, plus <marker>_taxids.txt when the build wrote one
    (rCRUX's taxid per accession; a 2-column input has none)."""
    optional = (f"{marker}_taxids.txt", ASSIGNMENT_CHECK_FILENAME)
    return required_files(marker) + [name for name in optional if (build_dir / name).is_file()]


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
# chosen one is stamped into tronko.dvc as meta.release_version so the merge
# can be tagged from git alone.
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
    checkout: Path, marker: str, requested: str | None, bump: str, base: str = "HEAD"
) -> str:
    """The requested version if free, else the next after this marker's tronko tags.

    A version stamped in ``base``'s tronko.dvc but not yet tagged (merged, tag
    pending) counts too, so versions only increase.
    """
    tagged = tronko_versions(known_tags(checkout), marker)
    if requested is not None:
        check_version(requested)
        if requested in tagged:
            raise ValueError(f"Tag {tronko_tag(marker, requested)} already exists")
        return requested
    on_base = subprocess.run(
        ["git", "show", f"{base}:databases/{marker}/tronko.dvc"],
        cwd=checkout,
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
    fs = DVCFileSystem(rcrux_repo, rev=commit, remote=RCRUX_REMOTE)
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
        "dvc_remote": RCRUX_REMOTE,
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
    been built from an image or another checkout.
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
    files = {
        name: file_digest(build_dir / name, "sha256") for name in published_files(build_dir, marker)
    }
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


TAXONOMY_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
# A release this much smaller than the one before it is flagged for the reviewer.
SHRINK_WARNING = 0.10


def _spread(values: list[float]) -> str:
    if not values:
        return "none"
    return f"{min(values):,g} / {statistics.median(values):,g} / {max(values):,g}"


def assignment_section(check: dict) -> list[str]:
    settings, reads = check["settings"], check["settings"]["reads"]
    return [
        "## Assignment check",
        f"{reads:,} read pairs of up to {settings['read_length']} bp, {settings['reads_per_query']} from each "
        f"of {settings['queries']} sequences of {settings['species']} species drawn from this database, "
        f"assigned with tronko-assign "
        f"(LCA cut-off {settings['lca_cutoff']}) and scored as the assignment benchmark scores:",
        "",
        "| | correct | unresolved | wrong | precision | recall | F1 |",
        "|---|---|---|---|---|---|---|",
        *(
            f"| {rank} | "
            + " | ".join(f"{s[k]:,} ({s[k] / reads:.1%})" for k in ("correct", "unresolved", "wrong"))
            + " | "
            + " | ".join(f"{s[k]:.3f}" for k in ("precision", "recall", "f1"))
            + " |"
            for rank, s in check["by_rank"].items()
        ),
        "",
        "Unresolved: assigned no further than a higher rank. The reads come from sequences "
        "in the database, so this does not measure accuracy on unseen species.",
        "",
    ]


def release_report(metadata: dict, build_dir: Path, build_log: Path, previous: dict | None) -> str:
    """The release PR's description: what a reviewer checks before merging.

    Everything comes from files the build already wrote. Checks that fail are
    listed first; they do not stop the publish.
    """
    marker = metadata["marker"]
    build = metadata["build_params"].get("build", {})
    settings = metadata["build_params"].get("settings", {})
    counts = metadata["row_counts"]
    log = build_log.read_text(errors="replace").splitlines()

    cluster_sizes = [
        int(m.group(1)) for line in log if (m := re.search(r"Launched cluster \d+ \(.* (\d+) seqs\)", line))
    ]
    # A tree is under the sum-of-pairs threshold when tronko-build wanted to split it
    # and could not: every refused split, less the input clusters refused before any
    # scoring. (Its SPscore lines also cover partitions it split and dropped, and
    # arrive cut up by its other output, so they are not read.)
    weak_trees = sum("partition too small" in line for line in log) - sum(
        "materialising unsplit cluster" in line for line in log
    )
    # tronko-build's list of the partitions that became trees; absent from a build made elsewhere.
    final = build_dir / "final_partitions.txt"
    partition_sizes = [
        count_lines(build_dir / f"partition{number}_taxonomy.txt")
        for number in (final.read_text().split() if final.is_file() else [])
    ]
    problems = [line.strip() for line in log if re.search(r"\bERROR\b|Traceback", line)]
    notices = [line.strip() for line in log if re.search(r"\bWARNING\b", line, re.IGNORECASE)]
    with open(build_dir / f"{marker}_taxonomy.txt") as handle:
        lineages = [line.rstrip("\n").split("\t")[1].split(";") for line in handle if line.strip()]
    named = [
        {lineage[i] for lineage in lineages if i < len(lineage) and lineage[i] not in ("", "NA")}
        for i in range(len(TAXONOMY_RANKS))
    ]
    no_species = sum(1 for lineage in lineages if len(lineage) < 7 or lineage[6] in ("", "NA"))

    warnings = []
    n_in = build.get("input_sequences")
    if n_in is not None and n_in != counts["fasta_records"]:
        warnings.append(f"{n_in:,} sequences went in but {counts['fasta_records']:,} came out")
    if counts["fasta_records"] != counts["taxonomy_rows"]:
        warnings.append(
            f"{counts['fasta_records']:,} FASTA records but {counts['taxonomy_rows']:,} taxonomy rows"
        )
    if problems:
        warnings.append(f"the build log has {len(problems)} error lines, the first: `{problems[0][:200]}`")
    check_path = build_dir / ASSIGNMENT_CHECK_FILENAME
    check = json.loads(check_path.read_text()) if check_path.is_file() else None
    for rank, result in (check or {}).get("by_rank", {}).items():
        if result["wrong"] > check["settings"]["reads"] * ASSIGNMENT_WRONG_WARNING:
            warnings.append(
                f"the assignment check put {result['wrong']:,} of {check['settings']['reads']:,} "
                f"reads in the wrong {rank}"
            )
    if previous and counts["fasta_records"] < previous["row_counts"]["fasta_records"] * (1 - SHRINK_WARNING):
        warnings.append(
            f"more than {SHRINK_WARNING:.0%} fewer sequences than {previous['tronko_version']}"
        )

    def change(now, before):
        return "" if before is None else f" | {before:,} | {now - before:+,}"

    old = previous or {}
    old_build = old.get("build_params", {}).get("build", {})
    was = f" | {old['tronko_version']} | change" if previous else ""
    lines = [
        *(
            ["**Checks failed:**", *(f"- {w}" for w in warnings)]
            if warnings
            else ["All checks passed: sequence counts agree and the build log has no errors."]
        ),
        "",
        "## Build",
        f"- rcrux release: `{metadata['rcrux_tag']}` ({metadata['input_tier']} tier, "
        f"commit `{metadata['rcrux_commit'][:12]}`)",
        f"- tronko commit: `{metadata['tronko_sha'][:12]}`",
        f"- settings changed from the recipe: "
        f"{', '.join(o.removeprefix('set:') for o in metadata['build_params'].get('overrides', [])) or 'none'}",
        f"- {build.get('threads', '?')} threads, {build.get('total_build_time_seconds', '?')} s "
        f"({build.get('tronko_build_time_seconds', '?')} s in tronko-build)",
        "",
        "## Size",
        f"| | this release{was} |",
        "|---|---|" + ("---|---|" if previous else ""),
        f"| sequences | {counts['fasta_records']:,}"
        f"{change(counts['fasta_records'], old.get('row_counts', {}).get('fasta_records'))} |",
        *(
            f"| {label} | {build[key]:,}{change(build[key], old_build.get(key) if previous else None)} |"
            for key, label in (("ac_clusters", "clusters"), ("n_trees", "trees"))
            if isinstance(build.get(key), int)
        ),
        "",
        "## Shape (smallest / median / largest)",
        f"- sequences per cluster: {_spread(cluster_sizes)}; {len(cluster_sizes)} clusters",
        f"- sequences per tree: {_spread(partition_sizes)}; {len(partition_sizes):,} trees",
        *(
            [
                f"- trees under the {settings.get('SP_THRESHOLD')} sum-of-pairs threshold, "
                f"too small to split further: {weak_trees:,} of {build.get('n_trees', 0):,}"
            ]
            if settings.get("PARTITION") == "sum_of_pairs"
            else []
        ),
        "",
        *(assignment_section(check) if check else []),
        "## Taxonomy",
        "- distinct names: " + ", ".join(f"{len(names):,} {rank}" for rank, names in zip(TAXONOMY_RANKS, named)),
        f"- sequences with no species name: {no_species:,} of {len(lineages):,}",
        "",
        f"## Build log warnings: {len(notices)}",
        *(f"- `{line[:200]}`" for line in sorted(set(notices))[:10]),
        "",
        "<details><summary>Every setting</summary>",
        "",
        *(f"- `{key}={value}`" for key, value in settings.items()),
        "",
        "</details>",
        "",
        "<details><summary>Files (sha256)</summary>",
        "",
        *(
            f"- `{name}` {(build_dir / name).stat().st_size:,} bytes `{digest}`"
            for name, digest in metadata["files"].items()
        ),
        "",
        "</details>",
    ]
    return "\n".join(lines)


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
    # A marker's first publish creates its directory and the .gitignore dvc add writes.
    ignore = repo / "databases" / marker / ".gitignore"
    if ignore.is_file() and not run(["git", "ls-files", str(ignore.relative_to(repo))], cwd=repo):
        ignore.unlink()
    for directory in (ignore.parent, ignore.parent.parent):
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()


def publish(
    marker: str,
    build_dir: Path,
    inputs_path: Path,
    build_log: Path,
    build_params_path: Path,
    version: str | None,
    tronko_checkout: Path,
    replace: bool = False,
    push_branch: bool = False,
    dry_run: bool = False,
    bump: str = "minor",
    tronko_sha: str | None = None,
    base_branch: str = "main",
) -> dict:
    """Track a finished build in the tronko repository's DVC and commit it on a branch.

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

    repo = tronko_checkout
    # Untracked files are fine: publish commits only the files it names, and a
    # tronko/ left over from an earlier publish is untracked once main is back.
    if run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=repo):
        raise RuntimeError(f"{repo} has uncommitted changes to tracked files")
    run(["git", "fetch", "origin", base_branch], cwd=repo)
    base = f"origin/{base_branch}"
    marker_path = f"databases/{marker}"
    already_released = run(
        ["git", "ls-tree", "--name-only", base, f"{marker_path}/tronko.dvc"], cwd=repo
    )
    if already_released and not replace:
        raise ValueError(
            f"{base} already has a tronko release of {marker}; pass --replace to publish "
            "the next version"
        )
    previous = None
    if already_released:
        from dvc.api import DVCFileSystem

        with DVCFileSystem(str(repo), rev=base).open(f"{marker_path}/tronko/{METADATA_FILENAME}") as handle:
            previous = json.load(handle)
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
    report = release_report(metadata, build_dir, build_log, previous)
    if dry_run:
        print(json.dumps(metadata, indent=2))
        print(report)
        print(
            f"would track {dest} with DVC, push to {TRONKO_REMOTE}, commit on {branch} from {base}"
        )
        return metadata

    current = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo)
    start = current if current != "HEAD" else run(["git", "rev-parse", "HEAD"], cwd=repo)
    run(["git", "checkout", "-b", branch, base], cwd=repo)
    try:
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        for name in published_files(build_dir, marker):
            shutil.copy2(build_dir / name, dest / name)
        shutil.copy2(build_log, dest / BUILD_LOG_FILENAME)
        (dest / BUILD_PARAMS_FILENAME).write_text(json.dumps(params, indent=2) + "\n")
        (dest / METADATA_FILENAME).write_text(json.dumps(metadata, indent=2) + "\n")

        target = f"{marker_path}/tronko"
        run(["dvc", "add", target], cwd=repo)
        stamp_pointer(repo / f"{target}.dvc", version, rcrux_tag=inputs["rcrux_tag"])
        run(["dvc", "push", "-r", TRONKO_REMOTE, f"{target}.dvc"], cwd=repo)
        run(["git", "add", f"{target}.dvc", f"{marker_path}/.gitignore"], cwd=repo)
        run(
            [
                "git",
                "commit",
                "-m",
                f"Publish {marker} tronko {version} built from {inputs['rcrux_tag']}",
                "-m",
                report,
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
    version: str | None,
    tronko_checkout: Path,
    ref: str = "origin/main",
    push: bool = True,
) -> str:
    """Tag the merged commit that carries this build's tronko.dvc.

    The version is the one publish stamped into the pointer; a ``version``
    given here must agree with it, so a tag never names a later publish.
    """
    if push:
        run(["git", "fetch", "origin", "--tags"], cwd=tronko_checkout)
    commit = run(["git", "rev-parse", f"{ref}^{{commit}}"], cwd=tronko_checkout)
    pointer = f"databases/{marker}/tronko.dvc"
    shown = subprocess.run(
        ["git", "show", f"{commit}:{pointer}"],
        check=False,
        cwd=tronko_checkout,
        capture_output=True,
        text=True,
    )
    if shown.returncode:
        raise ValueError(f"{ref} has no {pointer}; merge the publish branch first")
    published = read_pointer_version(shown.stdout)
    if published is None:
        raise ValueError(f"{pointer} at {ref} has no meta.{META_KEY}; publish stamps it")
    if version is not None and published != version:
        raise ValueError(f"{ref} carries {marker} tronko {published}, not {version}")
    name = tronko_tag(marker, published)
    if run(["git", "tag", "--list", name], cwd=tronko_checkout):
        raise ValueError(f"Tag {name} already exists")
    run(["git", "tag", name, commit], cwd=tronko_checkout)
    if push:
        run(["git", "push", "origin", name], cwd=tronko_checkout)
    return name


def tag_changed(tronko_checkout: Path, before: str, ref: str, push: bool = True) -> list[str]:
    """Tag every tronko.dvc that changed between ``before`` and ``ref``.

    This is what the "Tag tronko releases" workflow runs on a push to main. A
    pointer whose tag already exists is left alone, so a rerun is harmless.
    """
    changed = run(
        ["git", "diff", "--name-only", before, ref, "--", "databases/*/tronko.dvc"],
        cwd=tronko_checkout,
    ).splitlines()
    names = []
    for path in changed:
        try:
            names.append(tag(Path(path).parent.name, None, tronko_checkout, ref, push))
        except ValueError as error:
            if "already exists" not in str(error) and "has no" not in str(error):
                raise
            print(f"skipped {path}: {error}")
    return names


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser("fetch", help="Download curated inputs at an rcrux tag")
    p_fetch.add_argument("--marker", required=True)
    p_fetch.add_argument("--rcrux-tag", required=True, help="e.g. CO1_BF3_BR2_v1.0.0")
    p_fetch.add_argument("--out", required=True, type=Path)
    p_fetch.add_argument("--rcrux-repo", default=DEFAULT_RCRUX_REPO)

    p_pub = sub.add_parser("publish", help="Track a finished build in the tronko repo's DVC")
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
    p_pub.add_argument("--base-branch", default="main", help="tronko branch to start from")
    p_pub.add_argument(
        "--tronko-checkout",
        required=True,
        type=Path,
        help="A clean tronko checkout to commit the release on (not the one running this)",
    )
    p_pub.add_argument("--replace", action="store_true", help="Replace an earlier tronko/")
    p_pub.add_argument("--push-branch", action="store_true")
    p_pub.add_argument("--dry-run", action="store_true")

    p_tag = sub.add_parser("tag", help="Tag the merged commit carrying tronko.dvc")
    p_tag.add_argument("--marker", help="Tag this marker's release")
    p_tag.add_argument("--tronko-version", help="Must match the version publish stamped")
    p_tag.add_argument(
        "--changed-since",
        metavar="COMMIT",
        help="Instead of --marker: tag every tronko.dvc changed between COMMIT and --ref",
    )
    p_tag.add_argument("--tronko-checkout", type=Path, default=Path("."))
    p_tag.add_argument("--ref", default="origin/main")
    p_tag.add_argument("--no-push", action="store_true", help="Create the tag locally only")

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
                args.tronko_checkout,
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
                    f"Merging it tags {metadata['tronko_tag']} (Tag tronko releases workflow)."
                )
        elif args.changed_since:
            names = tag_changed(
                args.tronko_checkout, args.changed_since, args.ref, not args.no_push
            )
            print(f"Tagged {', '.join(names) or 'nothing'}")
        elif args.marker:
            name = tag(
                args.marker, args.tronko_version, args.tronko_checkout, args.ref, not args.no_push
            )
            print(f"Tagged {name}")
        else:
            parser.error("tag needs --marker or --changed-since")
    except (ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
