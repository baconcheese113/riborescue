"""NMDetective-AI's NMD-efficiency score per variant, and how it tracks the rule tier.

NMDetective-AI (Veiner et al. 2026; Vejni/NMDetectiveAI, MIT) is the model tier's deep member — a
fine-tuned Orthrus/Mamba encoder that reads the transcript sequence around a premature stop and
predicts an NMD-*efficiency* score: higher means stronger decay, lower (and negative) means escape.
It is not a rule and gives no yes/no verdict, so nothing here thresholds it into one — the score is
carried continuously, and its agreement with the rules is measured as a separation of group means,
not a shared label. See `scripts/nmdetective_predict.py` for the GPU inference that produces it.
"""

from __future__ import annotations

import pandas as pd

__all__ = [
    "nmdetective_summary",
    "read_nmdetective",
]


def read_nmdetective(path: str) -> pd.DataFrame:
    """The per-variant NMDetective-AI table: availability, reason, and the efficiency score."""

    table = pd.read_csv(path, sep="\t")
    if "nmd_efficiency" not in table:
        table["nmd_efficiency"] = pd.NA
    return table


def nmdetective_summary(rule_predictors: pd.DataFrame, verdicts: pd.DataFrame) -> dict:
    """Where NMDetective-AI's efficiency lands relative to the rule tier's escape calls.

    The check is threshold-free: over the variants the model scored, it compares the mean efficiency
    of the stops the rules call escape against those they call decay. A model concordant with the
    rules gives escaping stops a *lower* efficiency, so a positive decay-minus-escape separation is
    agreement in direction — without turning the continuous score into a verdict it does not make.
    """

    available = verdicts[verdicts["nmd_available"].fillna(False).astype(bool)]
    merged = rule_predictors.merge(available, on="variant_id", how="inner")
    summary: dict[str, object] = {
        "scoreable": len(rule_predictors),
        "available": len(available),
        "both_available": len(merged),
        "available_fraction": round(float(verdicts["nmd_available"].mean()), 4)
        if len(verdicts)
        else 0.0,
    }
    if merged.empty:
        return summary

    efficiency = merged["nmd_efficiency"].astype(float)
    summary["efficiency_median"] = round(float(efficiency.median()), 4)

    def separation(rule_col: str) -> dict:
        escape = efficiency[merged[rule_col].astype(bool)]
        decay = efficiency[~merged[rule_col].astype(bool)]
        return {
            "escape_mean_efficiency": round(float(escape.mean()), 4) if len(escape) else None,
            "decay_mean_efficiency": round(float(decay.mean()), 4) if len(decay) else None,
            "separation": round(float(decay.mean() - escape.mean()), 4)
            if len(escape) and len(decay)
            else None,
        }

    summary["guideline"] = separation("escape_guideline")
    summary["full_rules"] = separation("escape_full_rules")
    return summary
