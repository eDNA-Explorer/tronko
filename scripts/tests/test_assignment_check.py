import json
from pathlib import Path

import pytest

pytest.importorskip("synthetic_fastq")
import assignment_check as ac  # noqa: E402

# a2 repeats a1's sequence, so S1 has one distinct query; S2 has three, of which two are read.
SEQUENCES = {
    "a1": "ACGT" * 30, "a2": "ACGT" * 30, "b1": "TTGA" * 30, "b2": "TTGC" * 30,
    "b3": "TAGC" * 30, "c1": "GGCA" * 30, "n1": "CCCC" * 30,
}  # fmt: skip
LINEAGES = {
    "a1": "E;P;C;O;F1;G1;S1", "a2": "E;P;C;O;F1;G1;S1", "b1": "E;P;C;O;F1;G1;S2",
    "b2": "E;P;C;O;F1;G1;S2", "b3": "E;P;C;O;F1;G1;S2", "c1": "E;P;C;O;F2;G2;S3",
    "n1": "E;P;C;O;F2;G2;NA",
}  # fmt: skip


def write_reference(directory: Path) -> None:
    (directory / "M.fasta").write_text("".join(f">{a}\n{s}\n" for a, s in SEQUENCES.items()))
    (directory / "M_taxonomy.txt").write_text("".join(f"{a}\t{p}\n" for a, p in LINEAGES.items()))


def test_queries_take_up_to_two_distinct_sequences_of_named_species(tmp_path):
    write_reference(tmp_path)
    db = ac.load_cruxv2_database(tmp_path, "M")
    assert ac.queries(db, species=3, per_species=2, seed=1) == ["a1", "b1", "b2", "c1"]


def test_run_reads_every_query_equally_and_scores_by_rank(tmp_path, monkeypatch):
    write_reference(tmp_path)
    calls = []

    def assign(command, **_):
        # A stand-in tronko-assign: right to species for S1, genus only for S2,
        # the wrong family for S3. Sources are told apart by R1's first bases.
        calls.append(command)
        answer = {"ACGT": "E;P;C;O;F1;G1;S1", "TTGA": "E;P;C;O;F1;G1", "TTGC": "E;P;C;O;F1;G1",
                  "GGCA": "E;P;C;O;F1;G1;S1"}  # fmt: skip
        lines = Path(command[command.index("-1") + 1]).read_text().splitlines()
        with open(command[command.index("-o") + 1], "w") as handle:
            handle.write("Readname\tTaxonomic_Path\tScore\n")
            for header, read in zip(lines[0::4], lines[1::4], strict=True):
                handle.write(f"{header[1:]}\t{answer[read[:4]]}\t0\n")
        return type("Done", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr(ac.subprocess, "run", assign)
    result = ac.run(tmp_path, "M", threads=2, species=3, queries_per_species=2,
                    reads_per_query=10, read_length=300, seed=1)  # fmt: skip

    command = calls[0]
    for flag in ("-r", "-q", "-p", "-z", "-w", "-6"):
        assert flag in command
    assert command[command.index("--Cinterval") + 1] == "10"
    assert (result["settings"]["queries"], result["settings"]["reads"]) == (4, 40)
    species = result["by_rank"]["species"]
    assert (species["correct"], species["unresolved"], species["wrong"]) == (10, 20, 10)
    assert result["by_rank"]["family"]["wrong"] == 10
    json.dumps(result)
