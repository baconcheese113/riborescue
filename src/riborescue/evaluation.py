"""The evaluation seam: named splits, their bootstrap intervals, and the shuffle controls.

A split is named in configuration and generated from a seed, never improvised at a call site, so
that every reported number says which protocol produced it. The published protocol holds out random
variants; the grouped protocols hold out whole genes, because a gene contributes many variants with
near-identical local context, and a random draw puts those on both sides of the split.

The controls are the most important tests in the project and the ones most likely to be quietly
dropped. Their signatures live here so the protected control suite can assert against them before
the kinetics comparison itself exists.
"""

from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from sklearn.model_selection import GroupShuffleSplit, ShuffleSplit

from riborescue.baseline import cross_validate
from riborescue.contracts import EvalConfig

__all__ = [
    "ROUNDS",
    "SEED",
    "TRAIN_FRACTION",
    "BootstrapCI",
    "ShuffleKind",
    "UnsupportedEvalConfigError",
    "bootstrap_ci",
    "evaluate",
    "grouped_split_leakage",
    "run_shuffle_control",
    "split",
]

SEED = 721
ROUNDS = 10
TRAIN_FRACTION = 0.9


class ShuffleKind(StrEnum):
    global_ = "global"
    within_gene = "within_gene"
    context_matched = "context_matched"


class UnsupportedEvalConfigError(NotImplementedError):
    pass


@dataclass(frozen=True)
class BootstrapCI:
    low: float
    point: float
    high: float

    def includes_zero(self) -> bool:
        return self.low <= 0.0 <= self.high


def _random_rounds(features: pd.DataFrame, rounds: int, train_fraction: float, seed: int):
    splitter = ShuffleSplit(n_splits=rounds, train_size=train_fraction, random_state=seed)
    return splitter.split(features)


def _gene_grouped_rounds(features: pd.DataFrame, rounds: int, train_fraction: float, seed: int):
    splitter = GroupShuffleSplit(n_splits=rounds, train_size=train_fraction, random_state=seed)
    return splitter.split(features, groups=features["gene"])


_SPLITTERS: dict[EvalConfig, Callable] = {
    EvalConfig.published_random_cv: _random_rounds,
    EvalConfig.grouped_by_gene: _gene_grouped_rounds,
}


def split(
    features: pd.DataFrame,
    config: EvalConfig,
    *,
    rounds: int = ROUNDS,
    train_fraction: float = TRAIN_FRACTION,
    seed: int = SEED,
) -> pd.DataFrame:
    """Generate the training rows of each round under a named evaluation config.

    The result carries the same `round` and `row` columns the oracle's own assignments use, so a
    generated split and a published one are consumed identically.
    """

    if (splitter := _SPLITTERS.get(config)) is None:
        supported = ", ".join(sorted(c.value for c in _SPLITTERS))
        raise UnsupportedEvalConfigError(f"{config.value} has no splitter; supported: {supported}")

    rows = features.index.to_numpy()
    return pd.concat(
        pd.DataFrame({"round": round_, "row": rows[train]})
        for round_, (train, _) in enumerate(splitter(features, rounds, train_fraction, seed), 1)
    ).reset_index(drop=True)


def evaluate(features: pd.DataFrame, rounds: pd.DataFrame) -> pd.DataFrame:
    """Fit and score every round of a split, one row per round."""

    return pd.DataFrame(
        {"round": fit.round, "r2": fit.r2, "held_out": len(fit.predictions)}
        for fit in cross_validate(features, rounds)
    )


def bootstrap_ci(values: pd.Series, *, confidence: float = 0.95, seed: int = SEED) -> BootstrapCI:
    """Return a percentile bootstrap interval for the mean of per-round scores."""

    sample = np.asarray(values, dtype=float)
    interval = bootstrap(
        (sample,), np.mean, confidence_level=confidence, method="percentile", rng=seed
    ).confidence_interval
    return BootstrapCI(
        low=float(interval.low), point=float(sample.mean()), high=float(interval.high)
    )


def run_shuffle_control(kind: ShuffleKind) -> BootstrapCI:
    """Return the bootstrap CI of the kinetics improvement after shuffling.

    A genuine kinetic signal collapses under every shuffle, so a passing control has a CI that
    includes zero.
    """

    raise NotImplementedError("shuffle controls are implemented with the head-to-head comparison")


def grouped_split_leakage(config: str) -> BootstrapCI:
    """Return the improvement retained under a gene- or cluster-grouped split.

    Near-identical local contexts leaking across a random split would inflate this; grouping removes
    the leak.
    """

    raise NotImplementedError("grouped-split evaluation lands with the head-to-head comparison")
