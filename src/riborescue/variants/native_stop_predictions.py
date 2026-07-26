"""Scoring the amenability model at canonical stop codons — an extrapolation, labelled as one.

The model was fit on premature stops in a reporter. A canonical stop is out of that distribution: it
is a real termination signal the cell uses on purpose, in contexts the training never saw as
stops. Scoring it produces a number: a prediction of how readable-through the stop's
sequence context is — not whether the cell tolerates reading through it. It is kept in its
own column, never blended with the measured occupancy it is compared against.

The comparison is rank-only. The model predicts a readthrough efficiency on a reporter scale; the
ribo-seq measures downstream occupancy. Whether they agree is a question about ordering — does the
model rank native stops the way the cell's ribosomes do — so it is answered with a rank correlation,
not by expecting the two numbers to match.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy.stats import bootstrap, spearmanr

from riborescue.riboseq.readthrough_assay import gencode_sequences

__all__ = [
    "concordance",
    "four_quadrants",
    "native_stop_features",
]


def _transcribe(dna: str) -> str:
    """Lower-case RNA, as the model's training levels are spelled."""

    return dna.lower().replace("t", "u")


def native_stop_features(transcripts: Path, annotation: pd.DataFrame) -> pd.DataFrame:
    """The model's features read at each transcript's own stop codon.

    `stop_type` is the native stop, `up_123nt` the last sense codon, `down_123nt` the first codon of
    the 3' untranslated region — the same three features the premature-stop contexts carry, so the
    fitted model can score them unchanged.
    """

    sequences = gencode_sequences(transcripts)
    coding = annotation.loc[annotation["l_cds"] >= 6]
    rows = []
    for transcript, utr5, cds in zip(
        coding["transcript"], coding["l_utr5"], coding["l_cds"], strict=True
    ):
        sequence = sequences.get(str(transcript))
        if sequence is None:
            continue
        stop_start = int(utr5) + int(cds) - 3
        stop = sequence[stop_start : stop_start + 3]
        up = sequence[stop_start - 3 : stop_start]
        down = sequence[stop_start + 3 : stop_start + 6]
        if len(stop) < 3 or len(up) < 3 or len(down) < 3:
            continue
        if _transcribe(stop) not in ("uaa", "uag", "uga"):
            continue
        rows.append(
            {
                "transcript": str(transcript),
                "stop_type": _transcribe(stop),
                "up_123nt": _transcribe(up),
                "down_123nt": _transcribe(down),
            }
        )
    return pd.DataFrame(rows)


def concordance(
    predicted: pd.Series, measured: pd.Series, seed: int = 20240724, draws: int = 2000
) -> dict:
    """Spearman rank correlation between prediction and measurement, with a bootstrap interval.

    Rank-based on purpose: the two are on different scales, so what is testable is whether the model
    orders native stops the way the ribosomes do. The interval is a percentile bootstrap over the
    shared transcripts — the honest width given a few thousand of them, and it is reported whatever
    it shows, including a correlation indistinguishable from zero.
    """

    paired = pd.DataFrame({"p": predicted, "m": measured}).dropna()
    n = len(paired)
    if n < 20:
        return {"n": n, "rho": float("nan"), "low": float("nan"), "high": float("nan")}

    def rho(p: np.ndarray, m: np.ndarray) -> float:
        # scipy's correlation result ships no type information.
        return float(cast(Any, spearmanr(p, m)).statistic)

    # Resampled as pairs: a prediction and the measurement it is being ranked against have to travel
    # together, or the correlation is taken over a pairing the data never had.
    values = (paired["p"].to_numpy(), paired["m"].to_numpy())
    resampled = bootstrap(
        values,
        rho,
        paired=True,
        vectorized=False,
        n_resamples=draws,
        method="percentile",
        rng=seed,
    ).bootstrap_distribution

    # A resample that draws one rank over and over is constant, and a rank correlation over a
    # constant is undefined rather than zero. Those draws are dropped: the interval is taken over
    # the draws that defined one, which is what keeps a rare rank from erasing the interval
    # entirely instead of widening it.
    defined = resampled[np.isfinite(resampled)]
    low, high = np.percentile(defined, [2.5, 97.5]) if defined.size else (np.nan, np.nan)
    return {
        "n": n,
        "rho": round(rho(*values), 4),
        "low": round(float(low), 4),
        "high": round(float(high), 4),
    }


def four_quadrants(
    predicted: pd.Series, measured: pd.Series, predicted_cut: float, measured_cut: float
) -> pd.Series:
    """Label each native stop by whether prediction and measurement each call it susceptible.

    Four groups, not one score: they agree it opens, only the model says so, only the measurement
    says so, or neither. Where they disagree is the point — a blended number would hide it.
    """

    high_p = predicted >= predicted_cut
    high_m = measured >= measured_cut
    labels = pd.Series("neither", index=predicted.index)
    labels[high_p & high_m] = "both"
    labels[high_p & ~high_m] = "predicted only"
    labels[~high_p & high_m] = "measured only"
    return labels
