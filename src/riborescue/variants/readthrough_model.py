"""The published drug-specific readthrough model: a binomial GLM over stop type and the triplets
either side of the stop.

Performance is squared Pearson correlation on held-out data, the authors' definition, not the
deviance-based pseudo-r².
"""

from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy.highlevel import build_design_matrices, dmatrices

__all__ = [
    "FORMULA",
    "FittedModel",
    "RoundPrediction",
    "cross_validate",
    "fit",
    "fit_round",
    "r_squared",
    "relative_error_quantile",
]

FORMULA = "RT_binomial ~ 0 + stop_type + down_123nt + up_123nt + stop_type:down_123nt"

# An empty triplet-by-stop cell leaves a dependent column. R aliases the later one to NA and
# predicts as if it were zero; matching its predictions means matching that choice.
_RANK_TOLERANCE = 1e-7


@dataclass(frozen=True)
class RoundPrediction:
    """One cross-validation round: the fitted coefficients and the held-out predictions."""

    round: int
    coefficients: pd.Series
    predictions: pd.Series
    observed: pd.Series

    @property
    def r2(self) -> float:
        return r_squared(self.predictions, self.observed)


def relative_error_quantile(held_out: pd.DataFrame, quantile: float = 0.95) -> float:
    """The relative error within which the given share of held-out predictions fall.

    The model's own standard errors describe a single Bernoulli trial rather than the assay's
    measurement error, so the interval comes from what its held-out predictions actually did. Error
    grows with the prediction, so it is a fraction of the value rather than a fixed band.
    """

    error = (held_out["predicted"] - held_out["observed"]).abs() / held_out["predicted"]
    return float(error.quantile(quantile))


def r_squared(predicted: pd.Series, observed: pd.Series) -> float:
    """Squared Pearson correlation, the authors' measure of held-out performance."""

    return float(np.corrcoef(predicted, observed)[0, 1] ** 2)


def identifiable_columns(design: pd.DataFrame) -> list[str]:
    """The design columns that are not linear combinations of the ones before them.

    An unpivoted QR puts each column's norm, once the earlier ones are projected out, on the R
    diagonal, so the first of an aliased pair wins — the rule R applies.
    """

    matrix = design.to_numpy(dtype=float)
    residual = np.abs(np.diag(np.linalg.qr(matrix, mode="r")))
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1.0)
    identifiable = residual > _RANK_TOLERANCE * scale
    return [str(name) for name, keep in zip(design.columns, identifiable, strict=True) if keep]


@dataclass(frozen=True)
class FittedModel:
    """A model fitted on every measured variant, and the feature levels it was fitted over.

    A context carrying a level the library never held has no coefficient and cannot be scored.
    """

    identifiable: tuple[str, ...]
    encoding: Any
    fitted: Any
    levels: dict[str, frozenset[str]]

    def known_levels(self, features: pd.DataFrame) -> pd.Series:
        """Which rows carry only feature levels the training library contained."""

        known = pd.Series(True, index=features.index)
        for name, levels in self.levels.items():
            known &= features[name].isin(levels)
        return known

    def predict(self, features: pd.DataFrame) -> pd.Series:
        """Predict readthrough for rows whose levels the fit knows."""

        design = cast(
            pd.DataFrame,
            build_design_matrices([self.encoding], features, return_type="dataframe")[0],
        )
        return pd.Series(
            self.fitted.predict(design[list(self.identifiable)]),
            index=features.index,
            name="readthrough_predicted",
        )


def fit(features: pd.DataFrame) -> FittedModel:
    """Fit the published model on every measured variant."""

    response, built = dmatrices(FORMULA, features, return_type="dataframe")
    design = cast(pd.DataFrame, built)
    identifiable = identifiable_columns(design)
    fitted = sm.GLM(response, design[identifiable], family=sm.families.Binomial()).fit()
    return FittedModel(
        identifiable=tuple(identifiable),
        encoding=cast(Any, built).design_info,
        fitted=fitted,
        levels={
            name: frozenset(features[name].unique())
            for name in ("stop_type", "up_123nt", "down_123nt")
        },
    )


def fit_round(
    features: pd.DataFrame, train_rows: pd.Series, round_: int, formula: str = FORMULA
) -> RoundPrediction:
    """Fit on the given training rows and predict the rows left out of them.

    A grouped split can leave a whole triplet level out of training, and patsy cannot code a level
    it never saw. Such rows are dropped from the held-out set and counted, rather than crashing the
    round or being scored against a coefficient that does not exist. Under the published random
    rounds nothing is ever dropped, which is why parity is unaffected.
    """

    held_in = features.index.isin(train_rows)
    train, test = features[held_in], features[~held_in]
    for column in ("stop_type", "up_123nt", "down_123nt"):
        if column in features.columns:
            test = test[test[column].isin(set(train[column]))]

    # patsy ships no type information.
    response, built = dmatrices(formula, train, return_type="dataframe")
    design = cast(pd.DataFrame, built)
    identifiable = identifiable_columns(design)
    fitted = sm.GLM(response, design[identifiable], family=sm.families.Binomial()).fit()

    encoding = cast(Any, built).design_info
    held_out = cast(
        pd.DataFrame, build_design_matrices([encoding], test, return_type="dataframe")[0]
    )
    predicted = fitted.predict(held_out[identifiable])

    return RoundPrediction(
        round=round_,
        coefficients=fitted.params,
        predictions=pd.Series(predicted, index=test.index, name="predicted"),
        observed=test["RT_binomial"].rename("observed"),
    )


def cross_validate(
    features: pd.DataFrame, rounds: pd.DataFrame, formula: str = FORMULA
) -> list[RoundPrediction]:
    """Fit every round the oracle defined, in its order."""

    numbered = sorted(int(r) for r in rounds["round"].unique())
    return [
        fit_round(features, rounds.loc[rounds["round"] == round_, "row"], round_, formula)
        for round_ in numbered
    ]
