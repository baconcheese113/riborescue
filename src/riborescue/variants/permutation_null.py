"""The kinetics gain against a null built from many shuffles, corrected for testing six drugs.

ADR-0021. ADR-0020 froze the shuffle control as one permutation, with its interval taken over the
ten rounds of that single arrangement. Ten overlapping training folds say how *stably* one shuffled
mapping helps; they say nothing about how often an arbitrary mapping helps at all. A single draw is
one sample from a null whose spread is wide enough to produce either verdict — on this data, seed
721 made SRI look like a leak and SJ6986 look clean, from the draw alone.

The correction is the ordinary permutation test. Many independent shuffles build the null, and the
observed gain is placed in it. Two details make it a *familywise* test rather than six separate
ones:

- **The permutation is synchronised across drugs.** One shuffled mapping is applied to every drug in
  a permutation, so the null preserves the dependence between them — six drugs measured on one
  library over the same genes and the same contexts are not six independent experiments.
- **The statistic is the maximum over drugs.** Each drug's observed gain is compared against the
  distribution of the largest gain any drug achieved under a shuffle. That is a single-step
  Westfall-Young procedure, chosen here rather than a correction picked after the p-values existed.

The baseline reads no kinetic column, so its held-out score is identical under every permutation and
is fitted once per drug and round. Only the kinetic model is refitted, which is what makes a
thousand permutations a few hours rather than a few days.
"""

from __future__ import annotations

from collections.abc import Collection, Iterator
from dataclasses import dataclass

import numpy as np
import pandas as pd

from riborescue.core.contracts import EvalConfig
from riborescue.core.tables import read_table
from riborescue.variants.evaluation import ORACLE, SEED, ShuffleKind, split
from riborescue.variants.kinetics import MODELS, attach, permute_scores, shuffle
from riborescue.variants.readthrough_model import fit_round

__all__ = [
    "DrugFolds",
    "familywise_pvalues",
    "load_folds",
    "null_gains",
    "observed_gains",
]


@dataclass(frozen=True)
class DrugFolds:
    """One drug's features, its rounds, and the baseline scores that never change.

    The baseline is fitted once. It reads no kinetic column, so refitting it inside the permutation
    loop would compute the same numbers a thousand times.
    """

    drug: str
    features: pd.DataFrame
    rounds: pd.DataFrame
    baseline: np.ndarray

    @property
    def numbers(self) -> list[int]:
        return sorted(int(r) for r in self.rounds["round"].unique())


def load_folds(
    drugs: tuple[str, ...],
    scores: pd.Series,
    *,
    config: EvalConfig = EvalConfig.grouped_by_gene,
    oracle=ORACLE,
) -> list[DrugFolds]:
    """Each drug's folds and its baseline held-out scores, computed once.

    The folds are generated from the fixed seed, so the observed fits and every permuted fit see the
    same held-out rows. A null built on different splits would carry the split variation as well.
    """

    loaded: list[DrugFolds] = []
    for drug in drugs:
        features = attach(read_table(oracle / f"features_{drug}.tsv.gz"), scores)
        rounds = (
            read_table(oracle / f"rounds_{drug}.tsv.gz")
            if config is EvalConfig.published_random_cv
            else split(features, config)
        )
        numbers = sorted(int(r) for r in rounds["round"].unique())
        baseline = np.array(
            [
                fit_round(features, rounds.loc[rounds["round"] == n, "row"], n, MODELS["B"]).r2
                for n in numbers
            ]
        )
        loaded.append(DrugFolds(drug=drug, features=features, rounds=rounds, baseline=baseline))
    return loaded


def _gain(folds: DrugFolds, features: pd.DataFrame) -> np.ndarray:
    """The kinetic model's held-out gain over the cached baseline, round by round."""

    kinetic = np.array(
        [
            fit_round(
                features, folds.rounds.loc[folds.rounds["round"] == n, "row"], n, MODELS["K1"]
            ).r2
            for n in folds.numbers
        ]
    )
    return kinetic - folds.baseline


def observed_gains(folds: list[DrugFolds]) -> pd.DataFrame:
    """The measured table's gain, per drug and round. Nothing is shuffled here."""

    return pd.DataFrame(
        {"drug": one.drug, "round": n, "gain": g}
        for one in folds
        for n, g in zip(one.numbers, _gain(one, one.features), strict=True)
    )


def null_gains(
    folds: list[DrugFolds],
    scores: pd.Series,
    kind: ShuffleKind,
    *,
    permutations: int,
    offset: int = 0,
    seed: int = SEED,
    skip: Collection[int] = (),
) -> Iterator[pd.DataFrame]:
    """Yield each shuffle's per-drug mean gain, one permutation at a time.

    One permutation index gives one shuffled mapping, applied to every drug, so the drugs stay as
    dependent under the null as they are in the data. `offset` lets a long run be sharded across
    processes without two shards drawing the same permutation, and `skip` lets a killed shard resume
    without recomputing what it already wrote.

    Yielded rather than returned, because a run of a thousand permutations is hours long and a
    function that returns at the end is one that cannot be watched, cannot be resumed, and loses
    everything if it is interrupted. All three of those happened.
    """

    skipped = set(skip)
    for index in range(offset, offset + permutations):
        if index in skipped:
            continue
        draw = seed + 1 + index
        permuted = (
            None if kind is ShuffleKind.within_gene else permute_scores(scores, kind, seed=draw)
        )
        yield pd.DataFrame(
            {
                "permutation": index,
                "shuffle": kind.value,
                "drug": one.drug,
                "gain": float(
                    np.mean(
                        _gain(
                            one,
                            shuffle(one.features, scores, kind, seed=draw)
                            if permuted is None
                            else attach(one.features, permuted),
                        )
                    )
                ),
            }
            for one in folds
        )


def familywise_pvalues(observed: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    """Each drug's gain against the null of the largest gain any drug reached under a shuffle.

    Single-step Westfall-Young. The maximum is taken within a permutation, across drugs, so a drug
    is asked to beat not its own null but the best any of the six managed on the same shuffled
    mapping. That is what controls the family error rate without a correction chosen afterwards.

    The empirical p-value adds the observed to its own null, which is what keeps it from ever being
    reported as zero: with `n` permutations the smallest attainable value is `1 / (n + 1)`.
    """

    per_drug = observed.groupby("drug")["gain"].mean()
    maxima = null.groupby("permutation")["gain"].max().to_numpy()
    draws = len(maxima)
    rows = []
    for drug, gain in per_drug.items():
        beaten = int(np.sum(maxima >= gain))
        rows.append(
            {
                "drug": drug,
                "gain": float(gain),
                "permutations": draws,
                "null_mean": float(maxima.mean()),
                "null_sd": float(maxima.std(ddof=1)),
                "null_max": float(maxima.max()),
                "p_familywise": (1 + beaten) / (draws + 1),
                "resolution": 1 / (draws + 1),
            }
        )
    return pd.DataFrame(rows).sort_values("p_familywise", ignore_index=True)
