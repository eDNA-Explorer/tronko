import gzip
import json
from pathlib import Path

import pytest

pytest.importorskip("synthetic_fastq")
import assignment_check as ac  # noqa: E402

LINEAGES = {
    "a1": "E;P;C;O;F1;G1;S1",
    "a2": "E;P;C;O;F1;G1;S2",
    "a3": "E;P;C;O;F2;G2;S3",
}


def test_run_assigns_the_community_and_scores_it_by_rank(tmp_path, monkeypatch):
    sequences = {"a1": "ACGT" * 30, "a2": "TTGA" * 30, "a3": "GGCA" * 30}
    (tmp_path / "M.fasta").write_text("".join(f">{a}\n{s}\n" for a, s in sequences.items()))
    (tmp_path / "M_taxonomy.txt").write_text("".join(f"{a}\t{p}\n" for a, p in LINEAGES.items()))
    calls = []

    def assign(command, **_):
        # A stand-in tronko-assign: right to species for S1, genus only for S2,
        # the wrong family for S3.
        calls.append(command)
        r1 = Path(command[command.index("-1") + 1])
        with gzip.open(r1.parent / "read_truth.tsv.gz", "rt") as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle][1:]
        out = {"a1": LINEAGES["a1"], "a2": "E;P;C;O;F1;G1", "a3": "E;P;C;O;F1;G1;S1"}
        with open(command[command.index("-o") + 1], "w") as handle:
            handle.write("Readname\tTaxonomic_Path\tScore\n")
            # named as tronko-assign names them: the whole header, joined with "_"
            handle.writelines(
                f"{read}_1:N:0:ATCG\t{out[accession]}\t0\n" for read, accession, _ in rows
            )
        return type("Done", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr(ac.subprocess, "run", assign)
    result = ac.run(tmp_path, "M", threads=2, species=3, reads=300, read_length=150, seed=1)

    command = calls[0]
    for flag in ("-r", "-q", "-p", "-z", "-w", "-6"):
        assert flag in command
    assert command[command.index("--Cinterval") + 1] == "10"
    assert result["settings"]["species"] == 3 and result["settings"]["reads"] == 300
    species, family = result["by_rank"]["species"], result["by_rank"]["family"]
    assert species["correct"] + species["unresolved"] + species["wrong"] == 300
    assert species["correct"] > 0 and species["unresolved"] > 0 and species["wrong"] > 0
    assert family["unresolved"] == 0 and family["wrong"] == species["wrong"]
    assert 0 < species["precision"] < 1
    json.dumps(result)
