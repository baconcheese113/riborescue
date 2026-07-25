"""Which pathogenic nonsense variants are plausibly addressable, and on which conditions.

Surviving transcript, readthrough of the stop and a tolerated insertion are reported as separate
conditions and never multiplied: how they combine is not a calibrated equation. Counts are counts of
variants, never of people — ClinVar carries no allele frequencies.

The transcript-survival rule is a screen rather than a verdict. At least thirty percent of the stops
it expects to be degraded escape decay, and readthrough itself antagonises decay.
"""

from dataclasses import dataclass

import pandas as pd

from riborescue.variants.nmd_rules import LAST_JUNCTION_NT
from riborescue.variants.residue import NEAR_COGNATE, conservative

__all__ = ["LAST_JUNCTION_RULE_NT", "TOLERABLE_SHARE", "Thresholds", "landscape", "summarise"]

LAST_JUNCTION_RULE_NT = LAST_JUNCTION_NT
"""The guideline NMD rule's last-junction window, shared with the ensemble (`nmd.py`) so the single
rule the landscape reports and the ensemble's `guideline` predictor can never diverge."""

TOLERABLE_SHARE = 0.5
"""The share of available insertions that must be conservative for the outcome to be robust."""


@dataclass(frozen=True)
class Thresholds:
    """The readthrough levels a landscape is reported under.

    The question inherits the uncertainty of every layer beneath it, so the ranking across
    thresholds carries the meaning rather than the count at any one.
    """

    readthrough: tuple[float, ...] = (0.005, 0.01, 0.02)


def _escapes_decay(contexts: pd.DataFrame) -> pd.Series:
    # The window is strictly 5' of the junction: a last-exon stop has a negative distance and is
    # carried by in_last_exon, not by this term — matching the ensemble's guideline predictor.
    distance = contexts["nt_to_last_junction"]
    within_rule = (distance > 0) & (distance <= LAST_JUNCTION_RULE_NT)
    return contexts["in_last_exon"] | within_rule.fillna(False)


def _near_cognates(stop_type: object) -> tuple[str, ...]:
    return NEAR_COGNATE[str(stop_type).upper().replace("U", "T")]


def _tolerable_share(contexts: pd.DataFrame) -> pd.Series:
    """The share of the insertions available here that evolution accepts for the residue lost.

    A share rather than a yes: the original residue is always among the stop's near-cognates, since
    the stop came from it by one base, so *whether* a tolerable insertion exists is always yes.
    """

    shares = {
        (original, stop): sum(conservative(original, inserted) for inserted in _near_cognates(stop))
        / len(_near_cognates(stop))
        for original, stop in contexts[["original_aa", "stop_type"]]
        .dropna()
        .drop_duplicates()
        .itertuples(index=False)
    }
    keys = list(zip(contexts["original_aa"], contexts["stop_type"], strict=True))
    return pd.Series([shares.get(key) for key in keys], index=contexts.index, dtype="Float64")


def landscape(contexts: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """One row per variant, each condition separate, alongside its best-scoring therapy."""

    present = scores[scores["status"] == "present"]
    best = present.loc[present.groupby("variant_id")["readthrough_predicted"].idxmax()]
    best = best.set_index("variant_id")[
        ["therapy_id", "readthrough_predicted", "readthrough_low", "readthrough_high"]
    ]
    best.columns = ["best_therapy", "best_readthrough", "best_readthrough_low", "best_high"]

    placed = contexts[contexts["scoreable"]].set_index("variant_id")
    joined = placed.join(best, how="left")
    return pd.DataFrame(
        {
            "gene_symbol": joined["gene_symbol"],
            "protein_position": joined["protein_position"],
            "stop_type": joined["stop_type"],
            "original_aa": joined["original_aa"],
            "review_stars": joined["review_stars"],
            "escapes_decay_by_rule": _escapes_decay(joined),
            "nt_to_last_junction": joined["nt_to_last_junction"],
            "best_therapy": joined["best_therapy"],
            "best_readthrough": joined["best_readthrough"],
            "best_readthrough_low": joined["best_readthrough_low"],
            "tolerable_insertion_share": _tolerable_share(joined),
        },
        index=joined.index,
    ).reset_index()


def summarise(landscape_table: pd.DataFrame, thresholds: Thresholds) -> pd.DataFrame:
    """Count variants meeting each condition, and all of them together, at each threshold.

    `all_conditions` is a conjunction, not a combined score.
    """

    escapes = landscape_table["escapes_decay_by_rule"]
    tolerated = landscape_table["tolerable_insertion_share"] >= TOLERABLE_SHARE
    rows = []
    for threshold in thresholds.readthrough:
        reads = landscape_table["best_readthrough"] >= threshold
        confident = landscape_table["best_readthrough_low"] >= threshold
        rows.append(
            {
                "readthrough_threshold": threshold,
                "variants": len(landscape_table),
                "escapes_decay": int(escapes.sum()),
                "reaches_threshold": int(reads.sum()),
                "reaches_threshold_lower_bound": int(confident.sum()),
                "most_insertions_tolerable": int(tolerated.sum()),
                "all_conditions": int((escapes & reads & tolerated).sum()),
                "all_conditions_lower_bound": int((escapes & confident & tolerated).sum()),
            }
        )
    return pd.DataFrame(rows)
