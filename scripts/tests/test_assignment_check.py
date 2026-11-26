import json
from pathlib import Path

import pytest

pytest.importorskip("community_lib")
import assignment_check as ac  # noqa: E402

TRUTH_IDS = ["asv_aaa_0_PE_0:1:SYNTHFC:1:1101:0:0", "asv_bbb_0_PE_3:1:SYNTHFC:1:1101:0:3"]


def test_assign_matches_tronko_assign_names_to_truth_ids(tmp_path, monkeypatch):
    """tronko-assign joins the whole FASTQ header with "_"; ids themselves hold "_"."""
    reads = ac.check.MemberReads(tmp_path / "r1", tmp_path / "r2", tmp_path / "truth", 2, 0)
    monkeypatch.setattr(ac.check, "load_read_truth", lambda _: [{"read_id": i} for i in TRUTH_IDS])
    calls = []

    def fake(command, **_):
        calls.append(command)
        Path(command[command.index("-o") + 1]).write_text(
            "Readname\tTaxonomic_Path\tScore\n"
            f"{TRUTH_IDS[0]}_1:N:0:ATCG\tE;P;C;O;F;G;S\t0\n"  # with the header's comment
            f"{TRUTH_IDS[1]}\tE;P;C\t0\n"  # without one
        )
        return type("Done", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr(ac.subprocess, "run", fake)
    predictions = ac.assign(tmp_path / "ref", "M", reads, 4, tmp_path / "out.txt")
    assert predictions == {TRUTH_IDS[0]: "E;P;C;O;F;G;S", TRUTH_IDS[1]: "E;P;C"}
    command = calls[0]
    assert command[command.index("-f") + 1] == str(tmp_path / "ref" / "reference_tree.trkb")
    assert command[command.index("-a") + 1] == str(tmp_path / "ref" / "M.fasta")
    for flag in ("-r", "-q", "-p", "-z", "-w", "-6"):
        assert flag in command
    assert command[command.index("--Cinterval") + 1] == "10"


def test_run_scores_every_community_against_its_own_holdout_reference(tmp_path, monkeypatch):
    members = [
        ac.check.CommunityMember(f"unified_{i}", 42 + i, {"normal": ("a",), "species_holdout": ("b", "c")},
                                 tmp_path / "f.fasta", tmp_path / "f.txt", {})
        for i in range(2)
    ]  # fmt: skip
    design = ac.check.CommunitySet(tmp_path, "M", ("illumina_typical", "edna_damage"), tuple(members))
    seen = {"design": [], "reference": [], "reads": [], "assign": []}

    def design_communities(fasta, taxonomy, out, **kwargs):
        seen["design"].append((fasta.name, taxonomy.name, kwargs))
        return design

    def build_reference(member, master, out, **kwargs):
        seen["reference"].append((member.name, master, kwargs["build_config"].threads))
        return {"dir": str(out), "removed_records": 10 + member.seed, "audit": {"unplaced_orphans": member.seed}}

    def write_reads(member, profile, out, **kwargs):
        seen["reads"].append((member.name, profile, kwargs))
        return ac.check.MemberReads(out / "r1", out / "r2", out / "truth", 1, 0)

    def score(truth, predictions):
        f1 = 0.9 if "illumina" in str(truth) else 0.5
        tier = {"reads": 10, "fp": 1, "precision": f1, "recall": f1, "f1": f1}
        return {"by_tier": {"normal": tier, "species_holdout": tier}}

    monkeypatch.setattr(ac.check, "design_communities", design_communities)
    monkeypatch.setattr(ac.check, "build_holdout_tronko_reference", build_reference)
    monkeypatch.setattr(ac.check, "write_member_reads", write_reads)
    monkeypatch.setattr(ac.check, "score_reads", score)
    monkeypatch.setattr(ac, "assign", lambda ref, *a: seen["assign"].append(ref.name) or {})

    inputs = {"forward_primer": "ACGT", "reverse_primer": "TTGG"}
    result = ac.run(tmp_path / "out", "M", inputs, 8, tmp_path / "work", 2, 150, 20)

    assert seen["design"] == [("M.fasta", "M_taxonomy.txt", {
        "marker": "M", "richness": 150, "n_communities": 2,
        "forward_primer": "ACGT", "reverse_primer": "TTGG"})]  # fmt: skip
    assert [r[0] for r in seen["reference"]] == ["unified_0", "unified_1"]
    assert all(r[1] == tmp_path / "out" and r[2] == 8 for r in seen["reference"])
    assert [(r[0], r[1]) for r in seen["reads"]] == [
        (m, p) for m in ("unified_0", "unified_1") for p in ("illumina_typical", "edna_damage")
    ]
    assert seen["reads"][0][2]["reads_per_query"] == 20
    assert seen["assign"] == ["unified_0", "unified_0", "unified_1", "unified_1"]
    assert [c["removed_records"] for c in result["communities"]] == [52, 53]
    assert [c["audit"] for c in result["communities"]] == [{"unplaced_orphans": 42}, {"unplaced_orphans": 43}]
    assert result["summary"]["illumina_typical"]["normal"] == {
        "precision": 0.9, "recall": 0.9, "f1": 0.9, "reads": 20, "fp": 2}  # fmt: skip
    json.dumps(result)
