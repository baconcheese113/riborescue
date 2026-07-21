"""The amenability landscape: which pathogenic nonsense variants are plausibly addressable, and how.

Three factors are brought together here and none is multiplied by another. Rescue requires surviving
transcript, readthrough of the stop, and a tolerated inserted residue, but the relationship between
them is not a calibrated equation, and presenting it as one would invent a precision nobody has
measured. Each condition keeps its own value, its own evidence and its own reason for failing,
and a variant meets a threshold or it does not.

Every count here is a count of *variants*. ClinVar carries no allele frequencies, so nothing in this
module says how many people a therapy or design would reach.

The transcript-survival condition is the canonical positional rule and is the weakest of the three.
At least thirty percent of the stops it expects to be degraded escape decay in practice, readthrough
itself antagonises decay, and the rule rests on positional assumptions the three published
predictors share. It is a screen, not a verdict, which is why it stays its own column rather than
folding into a score.
"""

from dataclasses import dataclass

import pandas as pd

from riborescue.residue import NEAR_COGNATE, conservative

__all__ = ["LAST_JUNCTION_RULE_NT", "TOLERABLE_SHARE", "Thresholds", "landscape", "summarise"]

LAST_JUNCTION_RULE_NT = 55
"""Distance from the last exon junction within which a premature stop escapes decay by the rule."""

TOLERABLE_SHARE = 0.5
"""The share of available insertions that must be conservative for the outcome to be robust."""


@dataclass(frozen=True)
class Thresholds:
    """The readthrough levels a landscape is reported under.

    The coverage question inherits the uncertainty of every layer beneath it, so the answer is given
    at several thresholds and the ranking carries the meaning, not the counts at any one of them.
    """

    readthrough: tuple[float, ...] = (0.005, 0.01, 0.02)


def _escapes_decay(contexts: pd.DataFrame) -> pd.Series:
    within_rule = contexts["nt_to_last_junction"] < LAST_JUNCTION_RULE_NT
    return contexts["in_last_exon"] | within_rule.fillna(False)


def _near_cognates(stop_type: object) -> tuple[str, ...]:
    return NEAR_COGNATE[str(stop_type).upper().replace("U", "T")]


def _tolerable_share(contexts: pd.DataFrame) -> pd.Series:
    """The share of the insertions available here that evolution accepts for the residue lost.

    Asking whether *any* insertion is tolerable answers nothing. A nonsense variant is one base from
    the stop it creates, so reverting that base restores the original residue and the original is
    always among the stop's near-cognates — the answer is trivially yes for every variant.

    What varies, and therefore what discriminates, is how much of the available insertion space is
    tolerable. Which residue a compound actually inserts is a property of the compound that nobody
    has measured across this population, so a variant where most of the possibilities are acceptable
    is a safer bet than one where only the exact original will do.
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

    `all_conditions` is a conjunction, not a combined score: the transcript is expected to survive,
    the best therapy's predicted readthrough clears the threshold, and at least one residue the
    compound could insert there, most are conservative substitutions for the one that was lost.
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
