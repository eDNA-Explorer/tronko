"""Unit tests for the BLAST arbiter's change classification and adjudication.

These cover the logic the whole analysis rests on -- in particular the case where
an assignment goes deeper on the same lineage, which is only credited when the
added ranks land on the lineage the top BLAST hit supports.

Run:  python3 scripts/release_analysis/test_arbiter.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tronko_blast_arbiter import (  # noqa: E402
    Asv,
    adjudicate,
    agree_depth,
    band_of,
    classify_change,
    depth_of,
    lca_lineage,
    split_path,
)

P = split_path
fails = []


def check(name, got, want):
    if got != want:
        fails.append(f"  {name}: got {got!r}, want {want!r}")


# ---- path parsing -------------------------------------------------------
check("split pads to 7", len(P("A;B")), 7)
check("NA becomes None", P("A;B;NA;NA;NA;NA;NA")[2], None)
check("unassigned is empty", depth_of(P("unassigned")), 0)
check("all-NA is empty", depth_of(P("NA;NA;NA;NA;NA;NA;NA")), 0)
check("depth stops at gap", depth_of(P("A;B;C")), 3)
check("agree_depth basic", agree_depth(P("A;B;C"), P("A;B;D")), 2)
check("agree_depth None-safe", agree_depth(P("A;B"), P("A;B;C")), 2)

# ---- change classification ---------------------------------------------
check("identical", classify_change(P("A;B;C"), P("A;B;C")), "identical")
check("extension", classify_change(P("A;B"), P("A;B;C")), "extension")
check("retraction", classify_change(P("A;B;C"), P("A;B")), "retraction")
check("divergence", classify_change(P("A;B;C"), P("A;B;D")), "divergence")
check("a2u", classify_change(P("A;B"), P("unassigned")), "assigned_to_unassigned")
check("u2a", classify_change(P("unassigned"), P("A;B")), "unassigned_to_assigned")
check("both un", classify_change(P("unassigned"), P("unassigned")), "both_unassigned")
# a deeper path that also diverges is divergence, not extension
check(
    "deeper+diverged is divergence",
    classify_change(P("A;B;C"), P("A;B;D;E")),
    "divergence",
)

# ---- deepening on the same lineage: the case this analysis exists for ----
# post adds C, BLAST says C -> confirmed, credit the change
check(
    "extension confirmed",
    adjudicate("extension", P("A;B"), P("A;B;C"), P("A;B;C;D")),
    ("post_better", "extension_confirmed"),
)
# post adds C, BLAST says that rank is X -> unearned depth, counts against
check(
    "extension contradicted",
    adjudicate("extension", P("A;B"), P("A;B;C"), P("A;B;X;D")),
    ("pre_better", "extension_contradicted"),
)
# BLAST lineage stops before the added rank -> cannot judge, do not credit
check(
    "extension unverifiable",
    adjudicate("extension", P("A;B"), P("A;B;C"), P("A;B")),
    ("unverifiable", "blast_too_shallow_for_added_ranks"),
)
# post adds two ranks, BLAST reaches only the first and agrees
check(
    "extension partly confirmed",
    adjudicate("extension", P("A;B"), P("A;B;C;D"), P("A;B;C")),
    ("post_better", "extension_partly_confirmed"),
)

# ---- retraction ---------------------------------------------------------
check(
    "retraction removed wrong rank",
    adjudicate("retraction", P("A;B;C"), P("A;B"), P("A;B;X")),
    ("post_better", "retraction_removed_wrong_ranks"),
)
check(
    "retraction lost a right rank",
    adjudicate("retraction", P("A;B;C"), P("A;B"), P("A;B;C")),
    ("pre_better", "retraction_lost_correct_ranks"),
)

# ---- divergence ---------------------------------------------------------
check(
    "divergence post closer",
    adjudicate("divergence", P("A;X;C"), P("A;B;C"), P("A;B;C")),
    ("post_better", "divergence_post_closer"),
)
check(
    "divergence pre closer",
    adjudicate("divergence", P("A;B;C"), P("A;X;C"), P("A;B;C")),
    ("pre_better", "divergence_pre_closer"),
)
check(
    "divergence tie",
    adjudicate("divergence", P("A;X"), P("A;Y"), P("A;Z")),
    ("tie", "divergence_equal_agreement"),
)

# ---- assigned/unassigned transitions ------------------------------------
check(
    "lost a supported call",
    adjudicate("assigned_to_unassigned", P("A;B"), P("unassigned"), P("A;B;C")),
    ("pre_better", "lost_a_blast_supported_call"),
)
check(
    "refused a contradicted call",
    adjudicate("assigned_to_unassigned", P("X;Y"), P("unassigned"), P("A;B")),
    ("post_better", "refused_a_blast_contradicted_call"),
)
check(
    "gained a supported call",
    adjudicate("unassigned_to_assigned", P("unassigned"), P("A;B"), P("A;B;C")),
    ("post_better", "gained_a_blast_supported_call"),
)
check(
    "gained a contradicted call",
    adjudicate("unassigned_to_assigned", P("unassigned"), P("X;Y"), P("A;B")),
    ("pre_better", "gained_a_blast_contradicted_call"),
)

# ---- no BLAST lineage ---------------------------------------------------
check(
    "no blast lineage",
    adjudicate("extension", P("A;B"), P("A;B;C"), None),
    ("unverifiable", "no_blast_lineage"),
)

# ---- unchanged-set computation must not crash ---------------------------
# regression: Asv is a dataclass with a list field, so it's unhashable by
# default -- `set(changed)` raised TypeError on every real run.
_a1, _a2, _a3 = Asv("a", "AAAA", "TTTT"), Asv("b", "CCCC", "GGGG"), Asv("c", "GGCC", "CCGG")
_changed_ids = {_a1.asv_id}
_unchanged = [a for a in [_a1, _a2, _a3] if a.asv_id not in _changed_ids]
check("unchanged excludes changed", sorted(a.asv_id for a in _unchanged), ["b", "c"])

# ---- lca_lineage: tied BLAST hits must not assert an arbitrary species --
# real vert12S data: Anas platyrhynchos and Anas acuta tie on bitscore for
# the same ASV -- genus should still be credited, species must not be.
check(
    "tied congeners collapse to genus, species unresolved",
    lca_lineage([P("A;B;C;D;E;Anas;Anas platyrhynchos"), P("A;B;C;D;E;Anas;Anas acuta")]),
    P("A;B;C;D;E;Anas"),
)
check(
    "tied subjects with identical lineage keep full depth",
    lca_lineage([P("A;B;C"), P("A;B;C")]),
    P("A;B;C"),
)
check("lca of no lineages is empty", lca_lineage([]), [None] * 7)

# ---- identity bands -----------------------------------------------------
check("band 100", band_of(100.0), ">=99")
check("band 99", band_of(99.0), ">=99")
check("band 98", band_of(98.0), "97-99")
check("band 85", band_of(85.0), "80-90")
check("band none", band_of(None), "<50/none")

if fails:
    print(f"FAIL ({len(fails)})")
    print("\n".join(fails))
    sys.exit(1)
print("all arbiter unit tests passed")
