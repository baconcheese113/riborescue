"""Which experiment would settle which open question, and how much it would settle.

This project's central results are two honest negatives and a data gap. What follows from that is
not a better ranking of therapies but the measurements that would resolve the ties — so this layer
answers a different question from the rest of the package. It never says *this therapy will work*.
It says *this is the experiment most likely to resolve an important uncertainty, and here is exactly
why*.

**What is authored and what is computed, stated plainly.** A protocol cannot be synthesised from a
table: the assay, the comparison and the decision rule for each programme are written by hand and
live in `experiments/programs.tsv`. What this module computes is everything that makes a programme
comparable to the others — how many variants, genes and conditions it would inform, which recorded
claims it could close, what evidence already exists, and how many replicates the observed variance
says it needs. Claiming the designs themselves were generated would be the kind of overstatement the
rest of this package exists to avoid.

**No combined score.** Every axis is reported separately and the frontier is Pareto: a programme is
on it when nothing else beats it on every axis at once. Blending reach, feasibility and evidence gap
into one number would hide the trade-off a reader is there to make, and would invite calling a
weighted sum something it is not.

**Nothing is invented.** A replicate count is either derived from a measured variance, with the
assumption shown, or reported as `not estimated`. Cost is not modelled at all; complexity is a
declared tier, not a currency.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd

__all__ = [
    "AXES",
    "PROGRAMS",
    "REACH",
    "frontier",
    "propose",
    "read_programs",
]

PROGRAMS = Path("experiments/programs.tsv")
"""The authored half: one row per programme, its question and its protocol."""

AXES = (
    "variants_informed",
    "genes_informed",
    "conditions_informed",
    "claims_resolved",
    "evidence_gap",
    "feasibility",
)
"""The dimensions the frontier is taken over. Reported separately, never summed.

`evidence_gap` and `feasibility` are ordinal ranks derived from the declared evidence grade and
complexity tier, so that the frontier can compare them; the underlying labels stay in the table.
"""

_EVIDENCE_GAP = {"none": 3, "indirect": 2, "partial": 1, "direct": 0}
"""How much is missing. A claim with no bearing evidence at all has the largest gap."""

_FEASIBILITY = {"low": 3, "moderate": 2, "high": 1, "very high": 0}
"""Declared complexity, inverted: a low-complexity programme is the more feasible one."""


def _uag(contexts: pd.DataFrame) -> pd.DataFrame:
    """Variants whose premature stop is UAG, which a UAG-reading suppressor tRNA would address."""

    return contexts.loc[contexts["stop_type"].str.upper().str.replace("U", "T") == "TAG"]


def _all_scoreable(contexts: pd.DataFrame) -> pd.DataFrame:
    return contexts


def _decay_escaping(contexts: pd.DataFrame) -> pd.DataFrame:
    """Variants whose message should survive decay, so readthrough has something to act on.

    Refuses rather than returning nothing when the decay verdict is absent. A reach rule that
    silently reports zero is indistinguishable from a programme that helps nobody, and the second is
    a finding while the first is a missing input.
    """

    if "escapes_decay_by_rule" not in contexts.columns:
        raise ValueError(
            "reach rule 'decay_escaping' needs `escapes_decay_by_rule`, which the amenability "
            "landscape carries and the placed contexts do not"
        )
    return contexts.loc[contexts["escapes_decay_by_rule"].astype(bool)]


REACH: dict[str, Callable[[pd.DataFrame], pd.DataFrame]] = {
    "uag_stops": _uag,
    "all_scoreable": _all_scoreable,
    "decay_escaping": _decay_escaping,
}
"""Named rules for who a programme informs, so reach is auditable rather than asserted.

A programme names one of these; the module applies it. A free-form query string would be neither
checkable nor safe to evaluate.
"""

_REQUIRED = (
    "experiment_id",
    "question",
    "why_it_matters",
    "what_the_lab_does",
    "comparison",
    "success_criterion",
    "if_it_fails",
    "evidence_gap_reason",
    "assay",
    "model_system",
    "endpoint",
    "decision_rule",
    "replicates",
    "resolves",
    "reach_rule",
    "evidence_grade",
    "complexity",
    "safety_relevant",
    "provenance",
)


def read_programs(path: Path = PROGRAMS) -> pd.DataFrame:
    """Load the authored programmes, refusing one that is malformed.

    Every field is required. Where a value does not exist there is a literal for it — `not
    estimated` for a replicate count nothing can derive, `none` for a programme that closes no
    recorded claim — because an empty cell reads as an omission rather than as a statement.

    A programme resolving no recorded claim is not thereby worthless: it opens a line the ledger
    does not yet carry, and it scores zero on that axis, which is visible on the frontier.
    """

    rows = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if missing := [column for column in _REQUIRED if column not in rows.columns]:
        raise ValueError(f"{path} is missing {', '.join(missing)}")
    if blank := [c for c in _REQUIRED if (rows[c] == "").any()]:
        raise ValueError(f"{path} leaves {', '.join(blank)} empty on at least one row")
    if unknown := sorted(set(rows["reach_rule"]) - set(REACH)):
        raise ValueError(f"{path} names reach rules that do not exist: {', '.join(unknown)}")
    if bad := sorted(set(rows["evidence_grade"]) - set(_EVIDENCE_GAP)):
        raise ValueError(f"{path} names evidence grades outside {sorted(_EVIDENCE_GAP)}: {bad}")
    if bad := sorted(set(rows["complexity"]) - set(_FEASIBILITY)):
        raise ValueError(f"{path} names complexity tiers outside {sorted(_FEASIBILITY)}: {bad}")
    if rows["experiment_id"].duplicated().any():
        raise ValueError(f"{path} reuses an experiment_id")
    return rows


def propose(
    programs: pd.DataFrame,
    contexts: pd.DataFrame,
    diseases: pd.DataFrame,
    ledger: pd.DataFrame,
) -> pd.DataFrame:
    """Attach the computed axes to each authored programme.

    Reach is counted over the variants the named rule selects, and the conditions those variants map
    to, so a programme addressing a common stop codon outranks one addressing a rare context on
    reach without either being scored against the other overall.

    `claims_resolved` counts only the recorded claims that are still open — a programme that would
    re-confirm something already supported is not credited for closing it.
    """

    open_claims = set(ledger.loc[~ledger["verdict"].isin(["supported", "refuted"]), "claim_id"])
    by_variant = diseases.drop_duplicates(["variant_id", "medgen"])

    computed = []
    for row in programs.itertuples():
        reached = REACH[str(row.reach_rule)](contexts)
        variants = set(reached["variant_id"])
        conditions = by_variant.loc[by_variant["variant_id"].isin(variants), "medgen"]
        named = {
            claim.strip()
            for claim in str(row.resolves).split(";")
            if claim.strip() and claim.strip() != "none"
        }
        computed.append(
            {
                "experiment_id": row.experiment_id,
                "variants_informed": len(variants),
                "genes_informed": reached["gene_symbol"].nunique(),
                "conditions_informed": conditions.nunique(),
                "claims_resolved": len(named & open_claims),
                "claims_named": "; ".join(sorted(named)),
                "evidence_gap": _EVIDENCE_GAP[str(row.evidence_grade)],
                "feasibility": _FEASIBILITY[str(row.complexity)],
            }
        )
    return programs.merge(pd.DataFrame(computed), on="experiment_id", validate="one_to_one")


def frontier(proposed: pd.DataFrame, axes: tuple[str, ...] = AXES) -> pd.DataFrame:
    """Mark the programmes nothing else beats on every axis at once.

    A Pareto frontier rather than a ranking, because the axes are not commensurable: more variants
    informed does not trade against lower complexity at any exchange rate this project could
    defend. `dominated_by` names what beat a programme, so a reader can see why it fell off rather
    than only that it did.
    """

    values = proposed[list(axes)].to_numpy()
    identifiers = list(proposed["experiment_id"])
    on_frontier, dominated_by = [], []
    for index in range(len(values)):
        beaten = [
            identifiers[other]
            for other in range(len(values))
            if other != index
            and (values[other] >= values[index]).all()
            and (values[other] > values[index]).any()
        ]
        on_frontier.append(not beaten)
        dominated_by.append("; ".join(beaten))
    return proposed.assign(on_frontier=on_frontier, dominated_by=dominated_by)
