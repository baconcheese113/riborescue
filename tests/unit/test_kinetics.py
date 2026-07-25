"""The kinetics comparison, and the structural fact that decides what it can be.

The first test here is the one the whole design rests on: a per-codon score added to the published
model is a column its design matrix already contains, so the fit does not change. If that ever stops
being true the four models in ADR-0020 are the wrong four, and this test is how anyone finds out.
"""

import numpy as np
import pandas as pd
import pytest

from riborescue.core.contracts import EvalConfig
from riborescue.riboseq.codon_occupancy import GENETIC_CODE, SENSE_CODONS, SYNONYMOUS
from riborescue.variants.evaluation import ShuffleKind, split
from riborescue.variants.kinetics import (
    MODELS,
    attach,
    codon_scores,
    head_to_head,
    improvement,
    permute_scores,
    shuffle,
)
from riborescue.variants.readthrough_model import FORMULA, fit_round

_ORACLE = "tests/fixtures/oracle"


@pytest.fixture(scope="module")
def features() -> pd.DataFrame:
    return pd.read_csv(f"{_ORACLE}/features_G418.tsv.gz", sep="\t")


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    """A codon table with 61 distinct scores, so a permutation is always visible."""

    return pd.DataFrame(
        {
            "codon": SENSE_CODONS,
            "amino_acid": [GENETIC_CODE[c] for c in SENSE_CODONS],
            "site": "a",
            "occupancy": np.linspace(0.5, 2.0, len(SENSE_CODONS)),
        }
    )


@pytest.fixture(scope="module")
def scores(table: pd.DataFrame) -> pd.Series:
    return codon_scores(table)


@pytest.fixture(scope="module")
def rounds(features: pd.DataFrame) -> pd.DataFrame:
    return split(features, EvalConfig.published_random_cv, rounds=2)


class TestTheSpanClaim:
    def test_a_codon_score_added_to_the_baseline_changes_nothing_at_all(self, features, scores):
        # The claim ADR-0020 rests on. `up_123nt` is a 64-level factor, so k_up is a linear
        # combination of columns already present and the added column is dropped as rank-deficient.
        attached = attach(features, scores)
        rows = pd.Series(attached.index[: int(0.9 * len(attached))])
        plain = fit_round(attached, rows, 1, FORMULA)
        augmented = fit_round(attached, rows, 1, f"{FORMULA} + k_up + k_down")
        pd.testing.assert_series_equal(plain.predictions, augmented.predictions)
        assert plain.r2 == augmented.r2

    def test_the_declared_kinetic_terms_are_not_spanned_and_do_change_the_fit(
        self, features, scores
    ):
        attached = attach(features, scores)
        rows = pd.Series(attached.index[: int(0.9 * len(attached))])
        plain = fit_round(attached, rows, 1, MODELS["B"])
        structured = fit_round(attached, rows, 1, MODELS["K1"])
        assert not np.allclose(plain.predictions, structured.predictions)
        # Six design columns, five of them identifiable beyond the baseline. The three
        # stop-type-specific k_up slopes sum to the marginal k_up effect, which the test above
        # shows the baseline already contains, so they carry two new dimensions rather than three.
        # Only how decoding demand varies with stop type is new, never its overall level.
        assert len(structured.coefficients) == len(plain.coefficients) + 5


class TestAttach:
    def test_both_triplets_are_scored_from_the_same_table(self, features, scores):
        attached = attach(features, scores)
        assert attached["k_up"].notna().all()
        assert attached["k_down"].notna().all()
        row = attached.iloc[0]
        assert row["k_up"] == scores[row["up_123nt"].upper().replace("U", "T")]

    def test_a_triplet_the_table_cannot_score_is_refused_not_left_empty(self, features, scores):
        thin = scores.drop("TTC")
        with pytest.raises(ValueError, match="scores no TTC"):
            attach(features, thin)

    def test_the_site_is_chosen_rather_than_taken_from_whichever_row_came_first(self, table):
        both = pd.concat([table, table.assign(site="p", occupancy=1.0)], ignore_index=True)
        assert codon_scores(both, site="p").eq(1.0).all()
        assert not codon_scores(both, site="a").eq(1.0).all()
        with pytest.raises(ValueError, match="no rows for the e site"):
            codon_scores(both, site="e")


class TestShuffles:
    def test_the_global_shuffle_keeps_the_scores_and_moves_them(self, scores):
        permuted = permute_scores(scores, ShuffleKind.global_, seed=1)
        assert sorted(permuted) == pytest.approx(sorted(scores))
        assert not permuted.equals(scores)

    def test_the_context_matched_shuffle_never_moves_a_score_out_of_its_family(self, scores):
        permuted = permute_scores(scores, ShuffleKind.context_matched, seed=1)
        for family in SYNONYMOUS.values():
            assert sorted(permuted[list(family)]) == pytest.approx(sorted(scores[list(family)]))
        assert not permuted.equals(scores)

    def test_the_context_matched_shuffle_leaves_every_amino_acid_mean_untouched(self, scores):
        permuted = permute_scores(scores, ShuffleKind.context_matched, seed=1)
        for family in SYNONYMOUS.values():
            assert permuted[list(family)].mean() == pytest.approx(scores[list(family)].mean())

    def test_the_single_codon_families_cannot_move(self, scores):
        permuted = permute_scores(scores, ShuffleKind.context_matched, seed=3)
        for codon in ("ATG", "TGG"):
            assert permuted[codon] == scores[codon]

    def test_within_gene_keeps_each_genes_kinetics_and_reassigns_them(self, features, scores):
        plain = attach(features, scores)
        permuted = shuffle(features, scores, ShuffleKind.within_gene, seed=1)
        for _gene, rows in plain.groupby("gene"):
            moved = permuted.loc[rows.index]
            assert sorted(moved["k_up"]) == pytest.approx(sorted(rows["k_up"]))
        assert not permuted["k_up"].equals(plain["k_up"])

    def test_within_gene_moves_the_two_columns_together(self, features, scores):
        plain = attach(features, scores)
        permuted = shuffle(features, scores, ShuffleKind.within_gene, seed=1)
        # Every (k_up, k_down) pair still exists somewhere; only which variant holds it changed.
        before = sorted(zip(plain["k_up"], plain["k_down"], strict=True))
        after = sorted(zip(permuted["k_up"], permuted["k_down"], strict=True))
        assert before == pytest.approx(after)

    def test_no_shuffle_touches_a_sequence_column(self, features, scores):
        columns = ["gene", "stop_type", "up_123nt", "down_123nt", "RT_binomial"]
        for kind in ShuffleKind:
            permuted = shuffle(features, scores, kind, seed=1)
            pd.testing.assert_frame_equal(permuted[columns], features[columns])

    def test_within_gene_is_not_a_table_permutation(self, scores):
        with pytest.raises(ValueError, match="permutes variants"):
            permute_scores(scores, ShuffleKind.within_gene, seed=1)


class TestClusterSplit:
    def test_a_triplet_pair_never_sits_on_both_sides_of_the_split(self, features):
        rounds = split(features, EvalConfig.grouped_by_sequence_cluster, rounds=3)
        pairs = features["up_123nt"] + "|" + features["down_123nt"]
        for round_ in rounds["round"].unique():
            held_in = features.index.isin(rounds.loc[rounds["round"] == round_, "row"])
            assert not set(pairs[held_in]) & set(pairs[~held_in])

    def test_the_gene_split_still_holds_out_whole_genes(self, features):
        rounds = split(features, EvalConfig.grouped_by_gene, rounds=3)
        for round_ in rounds["round"].unique():
            held_in = features.index.isin(rounds.loc[rounds["round"] == round_, "row"])
            assert not set(features["gene"][held_in]) & set(features["gene"][~held_in])


class TestHeadToHead:
    def test_every_model_reports_every_round(self, features, scores, rounds):
        models = {name: MODELS[name] for name in ("B", "K1", "K2")}
        scored = head_to_head(features, rounds, scores, models=models)
        assert len(scored) == 3 * 2
        assert set(scored["model"]) == {"B", "K1", "K2"}

    def test_the_baseline_is_identical_under_every_shuffle(self, features, scores, rounds):
        models = {"B": MODELS["B"], "K1": MODELS["K1"]}
        plain = head_to_head(features, rounds, scores, models=models)
        for kind in ShuffleKind:
            shuffled = head_to_head(features, rounds, scores, models=models, kind=kind)
            baseline = shuffled.loc[shuffled["model"] == "B", "r2"].to_numpy()
            assert baseline == pytest.approx(plain.loc[plain["model"] == "B", "r2"].to_numpy())

    def test_a_shuffle_moves_the_kinetic_model_and_not_the_baseline(self, features, scores, rounds):
        models = {"B": MODELS["B"], "K1": MODELS["K1"]}
        plain = head_to_head(features, rounds, scores, models=models)
        shuffled = head_to_head(features, rounds, scores, models=models, kind=ShuffleKind.global_)
        assert not np.allclose(
            plain.loc[plain["model"] == "K1", "r2"].to_numpy(),
            shuffled.loc[shuffled["model"] == "K1", "r2"].to_numpy(),
        )

    def test_the_gain_is_paired_within_round(self, features, scores, rounds):
        models = {"B": MODELS["B"], "K1": MODELS["K1"]}
        scored = head_to_head(features, rounds, scores, models=models)
        gains = improvement(scored)
        assert set(gains["model"]) == {"K1"}
        wide = scored.pivot(index="round", columns="model", values="r2")
        expected = (wide["K1"] - wide["B"]).to_numpy()
        assert gains.sort_values("round")["gain"].to_numpy() == pytest.approx(expected)

    def test_a_run_without_the_baseline_is_refused(self, features, scores, rounds):
        scored = head_to_head(features, rounds, scores, models={"K1": MODELS["K1"]})
        with pytest.raises(ValueError, match="no 'B' rows"):
            improvement(scored)
