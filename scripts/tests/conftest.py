import shutil
from pathlib import Path


def write_taxdump(directory: Path, nodes: list[tuple[int, int, str, str]], merged: tuple = ()) -> Path:
    """An NCBI taxdump of (taxid, parent, rank, scientific name) nodes under root 1."""
    directory.mkdir(parents=True, exist_ok=True)
    rows = [(1, 1, "no rank", "root"), *nodes]
    (directory / "nodes.dmp").write_text(
        "".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\t0\t|\n" for t, p, r, _ in rows)
    )
    (directory / "names.dmp").write_text(
        "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in rows)
    )
    (directory / "merged.dmp").write_text("".join(f"{old}\t|\t{new}\t|\n" for old, new in merged))
    return directory


def lineage_nodes(first: int, lineage: str) -> list[tuple[int, int, str, str]]:
    """Nodes for one seven-rank lineage, taxids from first down, its species last."""
    ranks = ("domain", "phylum", "class", "order", "family", "genus", "species")
    return [
        (first + i, first + i - 1 if i else 1, rank, name)
        for i, (rank, name) in enumerate(zip(ranks, lineage.split(";")))
    ]


def copy_repos(template: Path, dest: Path) -> None:
    """A test's own copy of repositories built once for the session: git and DVC
    setup is seconds per repository, a copy is not. Their configs name each other
    by absolute path, so those are pointed at the copy."""
    shutil.copytree(template, dest, dirs_exist_ok=True, symlinks=True)
    for config in [*dest.glob("*/.git/config"), *dest.glob("*/.dvc/config"), *dest.glob("*.git/config")]:
        config.write_text(config.read_text().replace(str(template), str(dest)))
