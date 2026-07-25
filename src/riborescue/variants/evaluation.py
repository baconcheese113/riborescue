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
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from sklearn.model_selection import GroupShuffleSplit, ShuffleSplit

from riborescue.core.contracts import EvalConfig
from riborescue.variants.readthrough_model import cross_validate

__all__ = [
    "CODON_TABLE",
    "CONTROLLED",
    "ORACLE",
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


def _cluster_grouped_rounds(features: pd.DataFrame, rounds: int, train_fraction: float, seed: int):
    """Hold out whole sequence contexts, so the same triplet pair cannot sit on both sides.

    A scored variant carries nine nucleotides and no more, so its context *is* the pair of triplets
    either side of the stop. Grouping on that exact pair is what the cluster split can mean here:
    58% of variants share their pair with another, and a random draw puts those across the boundary,
    where a model with enough capacity is rewarded for memorising rather than for generalising. Stop
    type is deliberately not part of the group — it has three levels, and grouping on it would hold
    out a third of the library at a time.
    """

    groups = features["up_123nt"].astype(str) + "|" + features["down_123nt"].astype(str)
    splitter = GroupShuffleSplit(n_splits=rounds, train_size=train_fraction, random_state=seed)
    return splitter.split(features, groups=groups)


_SPLITTERS: dict[EvalConfig, Callable] = {
    EvalConfig.published_random_cv: _random_rounds,
    EvalConfig.grouped_by_gene: _gene_grouped_rounds,
    EvalConfig.grouped_by_sequence_cluster: _cluster_grouped_rounds,
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


ORACLE = Path("tests/fixtures/oracle")
"""The reproduction fixtures: per-drug features, the authors' rounds, the reliability ceilings."""

CODON_TABLE = Path("tests/fixtures/kinetics/codon_occupancy.tsv")
"""The committed codon table the controls run against.

Committed for the reason the oracle fixtures are: a control that reads a regenerable results tree
passes or fails according to what was last run there, which is the opposite of what a control is
for. It is 61 rows.
"""

CONTROLLED = "K1"
"""The model the controls stand over — the one whose gain is the claim."""


def _control_inputs():
    """Load what every control needs, deferring the imports that would close a cycle."""

    from riborescue.core.tables import read_table
    from riborescue.variants.kinetics import codon_scores

    if not CODON_TABLE.exists():
        raise FileNotFoundError(
            f"{CODON_TABLE} is absent; the controls run against a committed codon table, "
            "which `pixi run codon-occupancy` produces"
        )
    scores = codon_scores(read_table(CODON_TABLE))
    drugs = sorted(
        path.name.removeprefix("features_").removesuffix(".tsv.gz")
        for path in ORACLE.glob("features_*.tsv.gz")
    )
    return read_table, scores, drugs


def _gains(config: EvalConfig, kind: "ShuffleKind | None") -> dict[str, pd.Series]:
    """Each drug's per-round gain of the controlled model over the baseline, under one split."""

    from riborescue.variants.kinetics import MODELS, head_to_head, improvement

    read_table, scores, drugs = _control_inputs()
    models = {name: MODELS[name] for name in ("B", CONTROLLED)}
    found: dict[str, pd.Series] = {}
    for drug in drugs:
        features = read_table(ORACLE / f"features_{drug}.tsv.gz")
        rounds = (
            read_table(ORACLE / f"rounds_{drug}.tsv.gz")
            if config is EvalConfig.published_random_cv
            else split(features, config)
        )
        scored = head_to_head(features, rounds, scores, models=models, kind=kind)
        found[drug] = improvement(scored).sort_values("round")["gain"].reset_index(drop=True)
    return found


def run_shuffle_control(kind: ShuffleKind) -> BootstrapCI:
    """Return the bootstrap CI of the kinetics improvement after shuffling.

    A genuine kinetic signal collapses under every shuffle, so a passing control has a CI that
    includes zero.

    Run under the grouped-by-gene split, which is the one the claim is made on, and reported for the
    drug whose shuffled gain is largest. Taking the worst case rather than pooling across drugs is
    what stops one leaking drug from being averaged away by five clean ones.
    """

    gains = _gains(EvalConfig.grouped_by_gene, kind)
    worst = max(gains, key=lambda drug: gains[drug].mean())
    return bootstrap_ci(gains[worst])


def grouped_split_leakage(config: str) -> BootstrapCI:
    """Return the improvement the published random split adds over a grouped one.

    Near-identical local contexts sit on both sides of a random split, where a model with enough
    capacity is rewarded for memorising them. The quantity is therefore the *excess* of the random
    split's gain over the grouped split's, round by round — what a gain would lose if the leak were
    closed. A pipeline that is not leaking has an interval including zero, and ADR-0020's rule that
    a gain appearing only under the published split is capacity rather than signal is this
    measurement.

    Reported for the drug whose excess is largest, for the same reason the shuffles are.
    """

    grouped = _gains(EvalConfig(config), None)
    random = _gains(EvalConfig.published_random_cv, None)
    excess = {drug: random[drug] - grouped[drug] for drug in grouped}
    worst = max(excess, key=lambda drug: excess[drug].mean())
    return bootstrap_ci(excess[worst])
