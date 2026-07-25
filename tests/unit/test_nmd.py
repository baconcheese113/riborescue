import pandas as pd

from riborescue.variants.nmd import (
    LAST_JUNCTION_NT,
    LONG_EXON_NT,
    START_PROXIMAL_NT,
    disagreement_atlas,
    nmd_predictors,
)


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


def _one(row: dict) -> pd.Series:
    """The single predictor row for one hand-built context, over the shared defaults."""

    return nmd_predictors(_contexts([row])).iloc[0]


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


# --- Adversarial boundary tests: each threshold is one comparison from flipping a case, so pin the
# exactly-at-threshold behaviour against the module constants rather than against magic numbers. ---


def test_last_junction_window_is_inclusive_at_the_threshold_and_open_just_past_it():
    inside = _one({"nt_to_last_junction": float(LAST_JUNCTION_NT)})
    assert inside["rule_within_last_junction"]  # exactly at the threshold is within (inclusive)
    assert inside["escape_guideline"]
    past = _one({"nt_to_last_junction": float(LAST_JUNCTION_NT + 1)})
    assert not past["rule_within_last_junction"]  # one nt further is outside
    assert not past["escape_guideline"]


def test_a_last_exon_stop_is_not_reported_as_within_the_junction_window():
    # A stop in the last exon sits 3' of the last junction: its distance is negative, and the rule
    # must not fire on it — the last-exon flag already carries the escape, and "within the window"
    # means 5' of the junction, not past it.
    row = _one({"in_last_exon": True, "nt_to_last_junction": -3.0})
    assert not row["rule_within_last_junction"]
    assert row["escape_guideline"]  # still escapes, by the last-exon rule


def test_a_stop_exactly_at_the_junction_is_last_exon_not_within_window():
    # Distance 0 is the first base of the last exon: last-exon, not "within 50 nt upstream".
    row = _one({"in_last_exon": True, "nt_to_last_junction": 0.0})
    assert not row["rule_within_last_junction"]
    assert row["escape_guideline"]


def test_start_proximal_is_open_at_the_threshold():
    # "Within the first 150 nt" is codons 1-50 (offsets 0,3,...,147); offset 150 is codon 51, out.
    assert _one({"nt_from_start": START_PROXIMAL_NT - 3, "nt_to_last_junction": 800.0})[
        "rule_start_proximal"
    ]
    assert not _one({"nt_from_start": START_PROXIMAL_NT, "nt_to_last_junction": 800.0})[
        "rule_start_proximal"
    ]


def test_long_exon_is_open_at_the_threshold_and_fires_one_nt_past_it():
    assert not _one({"ptc_exon_length": LONG_EXON_NT, "nt_to_last_junction": 800.0})[
        "rule_long_exon"
    ]  # exactly 407 nt is not "longer than 407 nt"
    assert _one({"ptc_exon_length": LONG_EXON_NT + 1, "nt_to_last_junction": 800.0})[
        "rule_long_exon"
    ]


def test_predictors_stay_a_superset_at_every_boundary():
    # full_rules is a strict superset of guideline, so the two never disagree the other way: there
    # is no geometry where the guideline escapes but the full rule set does not.
    for row in [
        {"nt_to_last_junction": float(LAST_JUNCTION_NT)},
        {"in_last_exon": True, "nt_to_last_junction": float("nan")},
        {"nt_from_start": START_PROXIMAL_NT - 3, "nt_to_last_junction": 800.0},
        {"ptc_exon_length": LONG_EXON_NT + 1, "nt_to_last_junction": 800.0},
    ]:
        r = _one(row)
        assert not (r["escape_guideline"] and not r["escape_full_rules"])
