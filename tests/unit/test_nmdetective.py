import pandas as pd

from riborescue.variants.nmdetective import nmdetective_summary, read_nmdetective


def _verdicts() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"variant_id": "E1", "nmd_available": True, "reason": "", "nmd_efficiency": -0.4},
            {"variant_id": "E2", "nmd_available": True, "reason": "", "nmd_efficiency": -0.2},
            {"variant_id": "D1", "nmd_available": True, "reason": "", "nmd_efficiency": 0.6},
            {"variant_id": "D2", "nmd_available": True, "reason": "", "nmd_efficiency": 0.4},
            {"variant_id": "U1", "nmd_available": False, "reason": "encode_error:ValueError"},
        ]
    )


def _predictors() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"variant_id": "E1", "escape_guideline": True, "escape_full_rules": True},
            {"variant_id": "E2", "escape_guideline": False, "escape_full_rules": True},
            {"variant_id": "D1", "escape_guideline": False, "escape_full_rules": False},
            {"variant_id": "D2", "escape_guideline": False, "escape_full_rules": False},
            {"variant_id": "U1", "escape_guideline": False, "escape_full_rules": False},
        ]
    )


def test_the_summary_counts_only_the_variants_the_model_scored():
    summary = nmdetective_summary(_predictors(), _verdicts())
    assert summary["available"] == 4
    assert summary["both_available"] == 4
    assert summary["available_fraction"] == 0.8


def test_escaping_stops_carry_lower_efficiency_giving_a_positive_separation():
    summary = nmdetective_summary(_predictors(), _verdicts())
    full = summary["full_rules"]
    # full-rules escape = {E1:-0.4, E2:-0.2} mean -0.3; decay = {D1:0.6, D2:0.4} mean 0.5.
    assert full["escape_mean_efficiency"] == -0.3
    assert full["decay_mean_efficiency"] == 0.5
    assert full["separation"] == 0.8  # decay minus escape, positive = concordant with the rules


def test_an_all_unavailable_table_yields_no_separation():
    verdicts = _verdicts().assign(nmd_available=False)
    summary = nmdetective_summary(_predictors(), verdicts)
    assert summary["both_available"] == 0
    assert "full_rules" not in summary


def test_read_tolerates_a_table_without_the_efficiency_column(tmp_path):
    path = tmp_path / "nmdet.tsv"
    path.write_text("variant_id\tnmd_available\treason\nV1\tFalse\tno_mane_ensembl\n")
    table = read_nmdetective(str(path))
    assert "nmd_efficiency" in table.columns
