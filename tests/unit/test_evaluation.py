import pandas as pd
import pytest

from riborescue.core.contracts import EvalConfig
from riborescue.variants.evaluation import (
    ROUNDS,
    BootstrapCI,
    UnsupportedEvalConfigError,
    bootstrap_ci,
    split,
)


@pytest.fixture
def features() -> pd.DataFrame:
    """Sixty variants over twelve genes, five variants each — the real data's shape in miniature."""

    genes = [f"GENE{i:02d}" for i in range(12) for _ in range(5)]
    return pd.DataFrame(
        {"gene": genes, "RT_binomial": [0.01 + i / 1000 for i in range(60)]},
        index=pd.RangeIndex(1, 61, name="row"),
    )


@pytest.mark.parametrize("config", [EvalConfig.published_random_cv, EvalConfig.grouped_by_gene])
def test_a_split_holds_out_a_tenth_of_the_rows_each_round(
    features: pd.DataFrame, config: EvalConfig
):
    rounds = split(features, config)
    assert sorted(rounds["round"].unique()) == list(range(1, ROUNDS + 1))
    for _, group in rounds.groupby("round"):
        assert group["row"].is_unique
        assert set(group["row"]) <= set(features.index)
        assert 0 < len(features) - len(group) < len(features)


def test_a_gene_grouped_split_never_puts_one_gene_on_both_sides(features: pd.DataFrame):
    for _, group in split(features, EvalConfig.grouped_by_gene).groupby("round"):
        held_out = features.index.difference(pd.Index(group["row"]))
        trained_genes = set(features.loc[group["row"], "gene"])
        held_out_genes = set(features.loc[held_out, "gene"])
        assert trained_genes & held_out_genes == set()


def test_a_random_split_does_put_genes_on_both_sides(features: pd.DataFrame):
    """The leak the grouped split exists to remove, shown rather than assumed."""

    shared = 0
    for _, group in split(features, EvalConfig.published_random_cv).groupby("round"):
        held_out = features.index.difference(pd.Index(group["row"]))
        trained_genes = set(features.loc[group["row"], "gene"])
        shared += len(trained_genes & set(features.loc[held_out, "gene"]))
    assert shared > 0


def test_the_same_seed_gives_the_same_split(features: pd.DataFrame):
    first = split(features, EvalConfig.grouped_by_gene, seed=7)
    assert first.equals(split(features, EvalConfig.grouped_by_gene, seed=7))
    assert not first.equals(split(features, EvalConfig.grouped_by_gene, seed=8))


def test_a_config_without_a_splitter_says_which_ones_exist(features: pd.DataFrame):
    with pytest.raises(UnsupportedEvalConfigError, match="grouped_by_gene"):
        split(features, EvalConfig.grouped_by_sequence_cluster)


def test_a_bootstrap_interval_brackets_the_mean():
    interval = bootstrap_ci(pd.Series([0.70, 0.72, 0.74, 0.76, 0.78]))
    assert interval.low < interval.point < interval.high
    assert interval.point == pytest.approx(0.74)
    assert not interval.includes_zero()


def test_an_interval_spanning_zero_is_reported_as_such():
    assert BootstrapCI(low=-0.01, point=0.002, high=0.015).includes_zero()
