from rcrux_primers import primers

CONFIG = '''
VERT12S_PRIMER = PrimerConfig(
    name="12S_12SV5-F_12SV5-R",
    forward="TAGAACAGGCTCCTCTAG",
    reverse="TTAGATACCCCACTATGC",
)
'''


def test_primers_from_the_build_logs_echoed_command():
    log = "      --forward-primer     CGGTTGGGGTGACCTCGGA\n      --reverse-primer     GCTGTTATCCCTAGGGTAACT\n"
    assert primers("16S_x", log, "") == ("CGGTTGGGGTGACCTCGGA", "GCTGTTATCCCTAGGGTAACT")


def test_primers_from_the_build_logs_summary_line():
    log = "[INFO] Primers: fwd=CCHGAYATRGCHTTYCCHCG, rev=TCDGGRTGNCCRAARAAYCA, amplicon=[300-600bp]\n"
    assert primers("CO1_x", log, "") == ("CCHGAYATRGCHTTYCCHCG", "TCDGGRTGNCCRAARAAYCA")


def test_primers_from_the_config_preset_when_the_log_has_none():
    assert primers("12S_12SV5-F_12SV5-R", "no primers here", CONFIG) == (
        "TAGAACAGGCTCCTCTAG", "TTAGATACCCCACTATGC")  # fmt: skip


def test_no_primers_when_neither_records_them():
    assert primers("MarVer3", "", CONFIG) is None
