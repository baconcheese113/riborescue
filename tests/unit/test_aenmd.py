import pandas as pd

from riborescue.variants.aenmd import (
    aenmd_verdicts,
    mane_ensembl_by_gene,
    model_agreement,
    read_aenmd_rules,
)


def _rules() -> pd.DataFrame:
    return pd.DataFrame(
        [
            # V1: annotated on its MANE transcript, escapes by the last-exon rule.
            _rule("1:000000100|C|T", "ENST0001", is_last=True),
            # a row that keeps ENST0002 in the annotated set, but for a different variant than V2.
            _rule("2:000000999|G|A", "ENST0002"),
            # V5: annotated on its MANE transcript, no rule fires — decay.
            _rule("5:000000500|C|T", "ENST0005"),
        ]
    )


def _rule(key: str, transcript: str, **fired: bool) -> dict:
    row = {
        "key": key,
        "transcript": transcript,
        "is_ptc": True,
        "is_last": False,
        "is_penultimate": False,
        "is_css_proximal": False,
        "is_single": False,
        "is_407plus": False,
    }
    return row | fired


def _nonsense() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"variant_id": "V1", "gene_id": 1, "chrom": "1", "pos": 100, "ref": "C", "alt": "T"},
            {"variant_id": "V2", "gene_id": 2, "chrom": "2", "pos": 200, "ref": "G", "alt": "A"},
            {"variant_id": "V3", "gene_id": 3, "chrom": "3", "pos": 300, "ref": "C", "alt": "T"},
            {"variant_id": "V4", "gene_id": 4, "chrom": "4", "pos": 400, "ref": "C", "alt": "T"},
            {"variant_id": "V5", "gene_id": 5, "chrom": "5", "pos": 500, "ref": "C", "alt": "T"},
        ]
    )


# gene -> MANE Ensembl transcript, with gene 4 deliberately absent.
_MANE = {1: "ENST0001", 2: "ENST0002", 3: "ENST0003", 5: "ENST0005"}


def test_a_variant_annotated_on_its_mane_transcript_takes_aenmds_escape_call():
    out = aenmd_verdicts(_rules(), _nonsense(), _MANE).set_index("variant_id")
    assert out.loc["V1", "aenmd_available"]
    assert out.loc["V1", "aenmd_escape"]
    assert out.loc["V1", "aenmd_is_last"]


def test_a_variant_with_no_rule_firing_is_available_but_decays():
    out = aenmd_verdicts(_rules(), _nonsense(), _MANE).set_index("variant_id")
    assert out.loc["V5", "aenmd_available"]
    assert not out.loc["V5", "aenmd_escape"]


def test_a_transcript_aenmd_never_annotated_is_a_build_gap_not_a_variant_claim():
    out = aenmd_verdicts(_rules(), _nonsense(), _MANE).set_index("variant_id")
    assert not out.loc["V3", "aenmd_available"]
    assert out.loc["V3", "reason"] == "transcript_absent"


def test_a_transcript_aenmd_carried_but_did_not_annotate_here_is_variant_filtered():
    out = aenmd_verdicts(_rules(), _nonsense(), _MANE).set_index("variant_id")
    assert not out.loc["V2", "aenmd_available"]
    assert out.loc["V2", "reason"] == "variant_filtered"


def test_a_gene_with_no_mane_ensembl_pairing_is_marked_as_such():
    out = aenmd_verdicts(_rules(), _nonsense(), _MANE).set_index("variant_id")
    assert not out.loc["V4", "aenmd_available"]
    assert out.loc["V4", "reason"] == "no_mane_ensembl"


def test_agreement_counts_only_variants_aenmd_scored():
    verdicts = aenmd_verdicts(_rules(), _nonsense(), _MANE)
    predictors = pd.DataFrame(
        [
            {"variant_id": "V1", "escape_guideline": True, "escape_full_rules": True},
            {"variant_id": "V3", "escape_guideline": False, "escape_full_rules": True},
            {"variant_id": "V5", "escape_guideline": False, "escape_full_rules": False},
        ]
    )
    agree = model_agreement(predictors, verdicts)
    # V1 (both escape) and V5 (both decay) are the only variants aenmd scored and the rules cover.
    assert agree["both_available"] == 2
    assert agree["full_rules_vs_aenmd_agree"] == 2
    assert agree["full_rules_vs_aenmd_agree_fraction"] == 1.0


def test_read_rules_strips_the_transcript_version(tmp_path):
    path = tmp_path / "aenmd.tsv"
    path.write_text(
        "key\ttranscript\tis_ptc\tis_last\tis_penultimate\tis_css_proximal\tis_single\tis_407plus\n"
        "1:000000100|C|T\tENST0001.7\tTrue\tTrue\tFalse\tFalse\tFalse\tFalse\n"
    )
    assert read_aenmd_rules(str(path))["transcript"].iloc[0] == "ENST0001"


def test_mane_map_keeps_only_select_and_strips_versions(tmp_path):
    path = tmp_path / "summary.txt"
    path.write_text(
        "#NCBI_GeneID\tEnsembl_nuc\tMANE_status\n"
        "GeneID:1\tENST0001.7\tMANE Select\n"
        "GeneID:9\tENST0009.2\tMANE Plus Clinical\n"
    )
    mapping = mane_ensembl_by_gene(str(path))
    assert mapping == {1: "ENST0001"}
