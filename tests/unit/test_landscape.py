"""The landscape joins three conditions without combining them into one number."""

import pandas as pd
import pytest

from riborescue.variants.landscape import LAST_JUNCTION_RULE_NT, Thresholds, landscape, summarise


@pytest.fixture
def contexts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "variant_id": ["v1", "v2", "v3", "v4"],
            "scoreable": [True, True, True, False],
            "gene_symbol": ["AAA", "BBB", "CCC", "DDD"],
            "protein_position": [10, 20, 30, 40],
            "stop_type": ["uga", "uag", "uga", "uaa"],
            "original_aa": ["R", "Q", "W", "E"],
            "review_stars": [2, 1, 3, 0],
            "in_last_exon": [True, False, False, False],
            "nt_to_last_junction": [None, 10, 5000, 900],
        }
    )


@pytest.fixture
def scores() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "variant_id": ["v1", "v1", "v2", "v3"],
            "therapy_id": ["G418", "DAP", "G418", "G418"],
            "readthrough_predicted": [0.01, 0.03, 0.02, 0.001],
            "readthrough_low": [0.004, 0.012, 0.008, 0.0004],
            "readthrough_high": [0.02, 0.05, 0.04, 0.002],
            "status": ["present", "present", "present", "present"],
        }
    )


def test_only_placed_variants_appear(contexts, scores):
    assert list(landscape(contexts, scores)["variant_id"]) == ["v1", "v2", "v3"]


def test_the_best_scoring_therapy_is_the_one_carried(contexts, scores):
    table = landscape(contexts, scores).set_index("variant_id")
    assert table.loc["v1", "best_therapy"] == "DAP"
    assert table.loc["v1", "best_readthrough"] == 0.03


def test_a_stop_in_the_last_exon_escapes_decay(contexts, scores):
    table = landscape(contexts, scores).set_index("variant_id")
    assert bool(table.loc["v1", "escapes_decay_by_rule"])


def test_a_stop_close_to_the_last_junction_escapes_decay(contexts, scores):
    table = landscape(contexts, scores).set_index("variant_id")
    assert contexts.set_index("variant_id").loc["v2", "nt_to_last_junction"] < LAST_JUNCTION_RULE_NT
    assert bool(table.loc["v2", "escapes_decay_by_rule"])


def test_a_stop_far_upstream_does_not(contexts, scores):
    table = landscape(contexts, scores).set_index("variant_id")
    assert not bool(table.loc["v3", "escapes_decay_by_rule"])


def test_the_tolerable_share_is_a_fraction_of_what_the_stop_allows(contexts, scores):
    """Arginine is one of six residues a UGA stop can be read as, and substitutes for itself."""

    table = landscape(contexts, scores).set_index("variant_id")
    share = float(table.loc["v1", "tolerable_insertion_share"])  # type: ignore[arg-type]
    assert 0.0 < share <= 1.0


def test_the_conjunction_never_exceeds_its_weakest_condition(contexts, scores):
    table = landscape(contexts, scores)
    counts = summarise(table, Thresholds(readthrough=(0.005,)))
    row = counts.iloc[0]
    assert row["all_conditions"] <= row["escapes_decay"]
    assert row["all_conditions"] <= row["reaches_threshold"]
    assert row["all_conditions"] <= row["most_insertions_tolerable"]


def test_a_higher_threshold_can_only_narrow_the_set(contexts, scores):
    counts = summarise(landscape(contexts, scores), Thresholds(readthrough=(0.005, 0.01, 0.02)))
    assert list(counts["reaches_threshold"]) == sorted(counts["reaches_threshold"], reverse=True)
    assert list(counts["all_conditions"]) == sorted(counts["all_conditions"], reverse=True)


def test_the_lower_bound_is_never_more_generous_than_the_point_estimate(contexts, scores):
    counts = summarise(landscape(contexts, scores), Thresholds())
    assert (counts["all_conditions_lower_bound"] <= counts["all_conditions"]).all()
