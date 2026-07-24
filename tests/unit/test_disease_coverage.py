import pandas as pd

from riborescue.variants.disease_coverage import disease_coverage, disease_reach_frontier


def _diseases(rows: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """Disease rows from (variant_id, gene, medgen, status), with placeholder xrefs filled in."""

    return pd.DataFrame(
        [
            {
                "variant_id": vid,
                "gene_symbol": gene,
                "condition_name": f"disease-{medgen}",
                "medgen": medgen,
                "omim": "1" if status == "mapped" else "",
                "orphanet": "",
                "mondo": "",
                "mesh": "",
                "mapping_status": status,
                "reason": "",
            }
            for vid, gene, medgen, status in rows
        ]
    )


def _contexts(rows: list[tuple[str, bool, str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"variant_id": vid, "scoreable": ok, "stop_type": stop, "original_aa": aa}
            for vid, ok, stop, aa in rows
        ]
    )


def test_a_disease_fraction_is_covered_over_eligible_not_a_reached_flag():
    diseases = _diseases([("V1", "A", "C1", "mapped"), ("V2", "A", "C1", "mapped")])
    # One of the disease's two variants is scoreable, the other is not.
    contexts = _contexts([("V1", True, "uga", "R"), ("V2", False, "", "")])
    coverage = disease_coverage(diseases, contexts)
    row = coverage.iloc[0]
    assert row["eligible_variants"] == 2
    assert row["model_covered"] == 1
    assert row["covered_fraction"] == 0.5
    assert row["reach"]  # reached, but reach is not the coverage number
    assert not row["complete"]  # and reach is not complete coverage either


def test_placeholder_and_unmapped_conditions_are_not_diseases():
    diseases = _diseases([("V1", "A", "CN169374", "placeholder"), ("V2", "A", "", "unmapped")])
    contexts = _contexts([("V1", True, "uga", "R"), ("V2", True, "uga", "R")])
    assert disease_coverage(diseases, contexts).empty


def test_designs_contributing_are_the_covering_designs():
    diseases = _diseases([("V1", "A", "C1", "mapped"), ("V2", "B", "C1", "mapped")])
    contexts = _contexts([("V1", True, "uga", "R"), ("V2", True, "uaa", "K")])
    row = disease_coverage(diseases, contexts).iloc[0]
    assert row["designs"] == "UAA-K;UGA-R"
    assert row["designs_contributing"] == 2
    assert row["genes"] == 2
    assert row["complete"]  # both eligible variants covered


def test_mapping_completeness_reflects_the_cross_references():
    diseases = _diseases([("V1", "A", "C1", "mapped"), ("V2", "A", "C2", "medgen_only")])
    contexts = _contexts([("V1", True, "uga", "R"), ("V2", True, "uga", "R")])
    coverage = disease_coverage(diseases, contexts).set_index("medgen")
    assert coverage.loc["C1", "mapping_completeness"] == "mapped"
    assert coverage.loc["C2", "mapping_completeness"] == "medgen_only"


def test_a_variant_missing_from_contexts_counts_as_eligible_but_uncovered():
    diseases = _diseases([("V1", "A", "C1", "mapped")])
    contexts = _contexts([("V2", True, "uga", "R")])  # V1 absent
    row = disease_coverage(diseases, contexts).iloc[0]
    assert row["eligible_variants"] == 1
    assert row["model_covered"] == 0
    assert not row["reach"]


def test_the_reach_frontier_reuses_the_greedy_engine_and_reaches_the_most_diseases():
    # One design (UGA-R) reaches two diseases; another (UAA-K) reaches one. Greedy takes the first.
    diseases = _diseases(
        [
            ("V1", "A", "C1", "mapped"),
            ("V2", "B", "C2", "mapped"),
            ("V3", "C", "C3", "mapped"),
        ]
    )
    contexts = _contexts(
        [("V1", True, "uga", "R"), ("V2", True, "uga", "R"), ("V3", True, "uaa", "K")]
    )
    frontier = disease_reach_frontier(diseases, contexts)
    assert list(frontier["design_id"]) == ["UGA-R", "UAA-K"]
    assert list(frontier["marginal"]) == [2, 1]
    assert frontier["cumulative"].iloc[-1] == 3
