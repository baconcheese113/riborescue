import pandas as pd
import pytest

from riborescue.variants.suppressor_panels import (
    brute_force_optimal,
    coverage_frontier,
    covered_sets,
    greedy_panel,
    is_partition,
)


def _contexts(rows: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """Contexts from (variant_id, gene, stop, residue) tuples, plus one unscoreable row."""

    scoreable = [
        {
            "variant_id": vid,
            "gene_symbol": gene,
            "stop_type": stop,
            "original_aa": aa,
            "scoreable": True,
        }
        for vid, gene, stop, aa in rows
    ]
    dropped = {
        "variant_id": "X",
        "gene_symbol": "ZZZ",
        "stop_type": "uga",
        "original_aa": "R",
        "scoreable": False,
    }
    return pd.DataFrame([*scoreable, dropped])


def test_a_design_covers_a_variant_only_when_stop_and_residue_both_match():
    contexts = _contexts(
        [("V1", "A", "uga", "R"), ("V2", "A", "uga", "K"), ("V3", "B", "uaa", "R")]
    )
    sets = covered_sets(contexts)
    assert sets == {
        "UGA-R": frozenset({"V1"}),
        "UGA-K": frozenset({"V2"}),
        "UAA-R": frozenset({"V3"}),
    }


def test_unscoreable_variants_are_never_covered():
    contexts = _contexts([("V1", "A", "uga", "R")])
    covered = frozenset().union(*covered_sets(contexts).values())
    assert "X" not in covered


def test_exact_restoration_partitions_the_variants():
    contexts = _contexts(
        [("V1", "A", "uga", "R"), ("V2", "A", "uga", "K"), ("V3", "B", "uga", "R")]
    )
    # V1 and V3 share a design; no variant belongs to two designs, so the sets are disjoint.
    assert is_partition(covered_sets(contexts))


def test_greedy_takes_the_largest_design_first_and_reports_the_marginal():
    contexts = _contexts(
        [("V1", "A", "uga", "R"), ("V2", "B", "uga", "R"), ("V3", "C", "uaa", "K")]
    )
    frontier = coverage_frontier(contexts, objective="variants")
    assert list(frontier["design_id"]) == ["UGA-R", "UAA-K"]
    assert list(frontier["marginal"]) == [2, 1]
    assert list(frontier["cumulative"]) == [2, 3]
    assert list(frontier["uncovered"]) == [1, 0]
    assert frontier["cumulative_fraction"].iloc[-1] == pytest.approx(1.0)


def test_genes_objective_counts_distinct_genes_not_variants():
    # Two variants in gene A share one design; covering that design covers gene A once.
    contexts = _contexts(
        [("V1", "A", "uga", "R"), ("V2", "A", "uga", "R"), ("V3", "B", "uaa", "K")]
    )
    frontier = coverage_frontier(contexts, objective="genes")
    assert frontier["cumulative"].iloc[0] == 1
    assert frontier["cumulative"].iloc[-1] == 2


def test_ties_break_on_design_id_so_the_panel_is_deterministic():
    # Two designs each cover one variant; the lexicographically smaller id is chosen first.
    contexts = _contexts([("V1", "A", "uga", "R"), ("V2", "B", "uaa", "K")])
    frontier = coverage_frontier(contexts, objective="variants")
    assert list(frontier["design_id"]) == ["UAA-K", "UGA-R"]


def test_greedy_matches_the_brute_force_optimum_on_the_exact_predicate():
    # Disjoint sets make greedy provably optimal; the brute force is the independent check.
    contexts = _contexts(
        [
            ("V1", "A", "uga", "R"),
            ("V2", "B", "uga", "R"),
            ("V3", "C", "uga", "K"),
            ("V4", "D", "uaa", "L"),
            ("V5", "E", "uaa", "L"),
        ]
    )
    sets = covered_sets(contexts)
    frontier = coverage_frontier(contexts, objective="variants")
    for k in range(1, len(sets) + 1):
        assert int(frontier["cumulative"].iloc[k - 1]) == brute_force_optimal(sets, k)


def test_greedy_stays_within_the_submodular_bound_on_overlapping_sets():
    # A hand-built overlapping instance — the case a later conservative-substitution predicate makes
    # real. Greedy takes the big set A first, which blocks the optimal B+C pair at k=2, so it is
    # genuinely suboptimal here (5 of 6) yet never below the (1 - 1/e) guarantee.
    sets = {
        "A": frozenset({"1", "2", "3", "4"}),
        "B": frozenset({"1", "2", "5"}),
        "C": frozenset({"3", "4", "6"}),
    }
    steps = greedy_panel(sets)
    assert steps[1].cumulative == 5 and brute_force_optimal(sets, 2) == 6
    for step in steps:
        assert step.cumulative >= (1 - 1 / 2.718281828) * brute_force_optimal(sets, step.rank)


def test_marginals_sum_to_the_total_covered():
    contexts = _contexts(
        [("V1", "A", "uga", "R"), ("V2", "B", "uga", "K"), ("V3", "C", "uaa", "L")]
    )
    frontier = coverage_frontier(contexts, objective="variants")
    assert frontier["marginal"].sum() == frontier["cumulative"].iloc[-1]


def test_an_unknown_objective_is_refused():
    contexts = _contexts([("V1", "A", "uga", "R")])
    with pytest.raises(ValueError, match="objective must be"):
        coverage_frontier(contexts, objective="patients")
