import pandas as pd

from riborescue.variants.nmd import disagreement_atlas, nmd_predictors


def _contexts(rows: list[dict]) -> pd.DataFrame:
    base = {
        "variant_id": "V",
        "gene_symbol": "G",
        "scoreable": True,
        "in_last_exon": False,
        "nt_to_last_junction": 500.0,
        "nt_from_start": 500,
        "ptc_exon_length": 100,
    }
    return pd.DataFrame([base | row for row in rows])


def test_a_last_exon_stop_escapes_by_both_predictors_and_they_agree():
    row = nmd_predictors(
        _contexts([{"in_last_exon": True, "nt_to_last_junction": float("nan")}])
    ).iloc[0]
    assert row["escape_guideline"]
    assert row["escape_full_rules"]
    assert not row["predictors_disagree"]


def test_within_the_last_junction_window_escapes_by_the_guideline():
    row = nmd_predictors(_contexts([{"nt_to_last_junction": 40.0}])).iloc[0]
    assert row["rule_within_last_junction"]
    assert row["escape_guideline"]
    assert not row["predictors_disagree"]


def test_a_start_proximal_stop_splits_the_predictors():
    # Far from the last junction (guideline: decay), but close to the start (full rules: escape).
    row = nmd_predictors(_contexts([{"nt_from_start": 90, "nt_to_last_junction": 800.0}])).iloc[0]
    assert not row["escape_guideline"]
    assert row["escape_full_rules"]
    assert row["predictors_disagree"]
    assert row["rule_start_proximal"]


def test_a_long_exon_stop_splits_the_predictors():
    row = nmd_predictors(_contexts([{"ptc_exon_length": 900, "nt_to_last_junction": 800.0}])).iloc[
        0
    ]
    assert not row["escape_guideline"]
    assert row["escape_full_rules"]
    assert row["predictors_disagree"]
    assert row["rule_long_exon"]


def test_a_stop_no_rule_catches_is_predicted_to_decay_by_both():
    row = nmd_predictors(
        _contexts([{"nt_from_start": 500, "nt_to_last_junction": 800.0, "ptc_exon_length": 100}])
    ).iloc[0]
    assert not row["escape_guideline"]
    assert not row["escape_full_rules"]
    assert not row["predictors_disagree"]


def test_unscoreable_variants_get_no_verdict():
    contexts = _contexts([{"variant_id": "V1"}])
    contexts.loc[0, "scoreable"] = False
    assert nmd_predictors(contexts).empty


def test_the_atlas_attributes_each_disagreement_to_the_driving_rule():
    predictors = nmd_predictors(
        _contexts(
            [
                {"variant_id": "A", "in_last_exon": True, "nt_to_last_junction": float("nan")},
                {"variant_id": "B", "nt_from_start": 90, "nt_to_last_junction": 800.0},
                {"variant_id": "C", "ptc_exon_length": 900, "nt_to_last_junction": 800.0},
                {
                    "variant_id": "D",
                    "nt_from_start": 90,
                    "ptc_exon_length": 900,
                    "nt_to_last_junction": 800.0,
                },
                {"variant_id": "E", "nt_to_last_junction": 800.0},
            ]
        )
    )
    atlas = disagreement_atlas(predictors)
    assert atlas["scoreable"] == 5
    assert atlas["disagree"] == 3  # B, C, D
    assert atlas["driven_by_start_proximal"] == 1  # B
    assert atlas["driven_by_long_exon"] == 1  # C
    assert atlas["driven_by_both"] == 1  # D
    assert atlas["escape_guideline"] == 1  # only A
    assert atlas["escape_full_rules"] == 4  # A, B, C, D
