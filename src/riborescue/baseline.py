"""The published drug-specific readthrough model, refit in Python.

Toledano's model is a binomial GLM with a logit link over four terms, where the nucleotides on
either side of the premature stop enter as triplets rather than as separate positional terms. It is
reproduced here on the oracle's own fold assignments, so the Python fit is judged by whether it
returns the authors' predictions rather than by whether it lands near their published r².

Performance is squared Pearson correlation on held-out data, matching the authors' definition; it is
not the deviance-based pseudo-r² a binomial GLM would otherwise report.
"""

from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy.highlevel import build_design_matrices, dmatrices

__all__ = ["FORMULA", "FoldPrediction", "cross_validate", "fit_fold", "r_squared"]

FORMULA = "RT_binomial ~ 0 + stop_type + down_123nt + up_123nt + stop_type:down_123nt"

# A triplet-by-stop-codon cell missing from a training fold leaves its column linearly dependent on
# the others. R's pivoted QR keeps the earlier column and aliases the later one to NA, which
# `predict` then treats as zero; least squares by pseudo-inverse would instead spread the fit across
# both. Reproducing the authors' predictions means reproducing their choice.
_RANK_TOLERANCE = 1e-7


@dataclass(frozen=True)
class FoldPrediction:
    """One cross-validation round: the fitted coefficients and the held-out predictions."""

    round: int
    coefficients: pd.Series
    predictions: pd.Series
    observed: pd.Series

    @property
    def r2(self) -> float:
        return r_squared(self.predictions, self.observed)


def r_squared(predicted: pd.Series, observed: pd.Series) -> float:
    """Squared Pearson correlation, the authors' measure of held-out performance."""

    return float(np.corrcoef(predicted, observed)[0, 1] ** 2)


def identifiable_columns(design: pd.DataFrame) -> list[str]:
    """Return the design columns that are not linear combinations of the ones before them.

    The R factor of an unpivoted QR carries, on its diagonal, the norm of each column once the
    columns before it have been projected out. A negligible diagonal entry therefore marks a column
    that adds nothing the earlier ones do not already span, so the first of an aliased pair wins —
    the rule R applies when it reports a coefficient as NA.
    """

    matrix = design.to_numpy(dtype=float)
    residual = np.abs(np.diag(np.linalg.qr(matrix, mode="r")))
    scale = np.maximum(np.linalg.norm(matrix, axis=0), 1.0)
    identifiable = residual > _RANK_TOLERANCE * scale
    return [str(name) for name, keep in zip(design.columns, identifiable, strict=True) if keep]


def fit_fold(features: pd.DataFrame, train_rows: pd.Series, round_: int) -> FoldPrediction:
    """Fit one round on the given training rows and predict the rows left out of them.

    `features` is indexed by the oracle's `row` numbering, so a fold is named by the rows it trains
    on and the complement is the held-out set.
    """

    held_in = features.index.isin(train_rows)
    train, test = features[held_in], features[~held_in]

    # patsy ships no type information, so the design it hands back is cast at this boundary.
    response, built = dmatrices(FORMULA, train, return_type="dataframe")
    design = cast(pd.DataFrame, built)
    identifiable = identifiable_columns(design)
    fitted = sm.GLM(response, design[identifiable], family=sm.families.Binomial()).fit()

    encoding = cast(Any, built).design_info
    held_out = cast(
        pd.DataFrame, build_design_matrices([encoding], test, return_type="dataframe")[0]
    )
    predicted = fitted.predict(held_out[identifiable])

    return FoldPrediction(
        round=round_,
        coefficients=fitted.params,
        predictions=pd.Series(predicted, index=test.index, name="predicted"),
        observed=test["RT_binomial"].rename("observed"),
    )


def cross_validate(features: pd.DataFrame, folds: pd.DataFrame) -> list[FoldPrediction]:
    """Fit every round the oracle defined, in its order."""

    rounds = sorted(int(r) for r in folds["round"].unique())
    return [
        fit_fold(features, folds.loc[folds["round"] == round_, "row"], round_) for round_ in rounds
    ]
