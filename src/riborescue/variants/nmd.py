"""Whether a premature stop escapes nonsense-mediated decay, by more than one rule.

`escapes_decay` used to be a single rule — the 50-nt rule ACMG applies. This computes two named
rule-based predictors from the transcript geometry the context step already records, and exposes
where they disagree, because the disagreement is the point: the 50-nt rule alone mislabels a large
share of escaping stops, and a curator setting PVS1 strength should see that boundary rather than
one verdict. See ADR-0016; the model predictors (NMDetective-AI, predNMD) are specified there.

Two predictors, deterministic and one-directional:

- `guideline` — the PTC is in the last exon, or within 55 nt of the last exon-exon junction.
- `full_rules` — the guideline conditions, plus start-proximal escape (within 150 nt of the start)
  and long-exon escape (the PTC's exon is longer than 407 nt), after Lindeboom, Supek & Lehner
  (Nat. Genet. 2016).

`full_rules` is a superset of `guideline`, so they disagree only where the fuller set escapes and
the guideline does not — exactly the start-proximal and long-exon stops the 50-nt rule misses.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import pandas as pd

__all__ = [
    "LAST_JUNCTION_NT",
    "LONG_EXON_NT",
    "START_PROXIMAL_NT",
    "NmdRules",
    "disagreement_atlas",
    "nmd_predictors",
]

# The published thresholds (Lindeboom, Supek & Lehner, Nat. Genet. 2016; Nagy & Maquat 1998).
LAST_JUNCTION_NT = 55
START_PROXIMAL_NT = 150
LONG_EXON_NT = 407


@dataclass(frozen=True)
class NmdRules:
    """The individual escape rules for one stop, and the two predictors they compose into."""

    last_exon: bool
    within_last_junction: bool
    start_proximal: bool
    long_exon: bool

    @property
    def guideline(self) -> bool:
        """The 50-nt rule as ACMG applies it: last exon, or within 55 nt of the last junction."""

        return self.last_exon or self.within_last_junction

    @property
    def full_rules(self) -> bool:
        """The guideline plus the start-proximal and long-exon escapes of the fuller rule set."""

        return self.guideline or self.start_proximal or self.long_exon

    @property
    def disagree(self) -> bool:
        """Whether the two predictors split — only ever full-escape over guideline-decay."""

        return self.guideline != self.full_rules


def nmd_predictors(contexts: pd.DataFrame) -> pd.DataFrame:
    """Per scoreable variant, the escape rules and the two predictors, with a disagreement flag.

    `contexts` must carry the geometry the context step records — `in_last_exon`,
    `nt_to_last_junction`, `nt_from_start`, `ptc_exon_length`. Unscoreable rows have no geometry and
    are dropped: a stop that was never placed has no NMD verdict to give.
    """

    scoreable = contexts[contexts["scoreable"]]
    rows = []
    for row in scoreable.itertuples():
        # itertuples hands back a scalar of pandas' own union type; the arithmetic below wants plain
        # numbers, so each numeric field is cast before it is compared.
        junction = cast(Any, row.nt_to_last_junction)
        rules = NmdRules(
            last_exon=bool(row.in_last_exon),
            # A NaN junction distance means the stop is in the last exon (no downstream junction),
            # already carried by the last-exon flag — so an absent distance is just not in-window.
            within_last_junction=bool(pd.notna(junction) and junction < LAST_JUNCTION_NT),
            start_proximal=int(cast(Any, row.nt_from_start)) < START_PROXIMAL_NT,
            long_exon=int(cast(Any, row.ptc_exon_length)) > LONG_EXON_NT,
        )
        rows.append(
            {
                "variant_id": row.variant_id,
                "gene_symbol": row.gene_symbol,
                "rule_last_exon": rules.last_exon,
                "rule_within_last_junction": rules.within_last_junction,
                "rule_start_proximal": rules.start_proximal,
                "rule_long_exon": rules.long_exon,
                "escape_guideline": rules.guideline,
                "escape_full_rules": rules.full_rules,
                "predictors_disagree": rules.disagree,
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "variant_id",
            "gene_symbol",
            "rule_last_exon",
            "rule_within_last_junction",
            "rule_start_proximal",
            "rule_long_exon",
            "escape_guideline",
            "escape_full_rules",
            "predictors_disagree",
        ],
    )


def disagreement_atlas(predictors: pd.DataFrame) -> dict:
    """Where the two predictors split, and which added rule drives each split.

    Every disagreement is the fuller rule set escaping a stop the guideline calls decay, so the
    split is attributable: it is driven by the start-proximal rule, the long-exon rule, or both.
    That attribution is the atlas — it says which blind spot of the 50-nt rule is doing the work.
    """

    total = len(predictors)
    disagree = predictors[predictors["predictors_disagree"]]
    start = disagree["rule_start_proximal"]
    long_exon = disagree["rule_long_exon"]

    def fraction(count: int) -> float:
        return round(count / total, 4) if total else 0.0

    guideline = int(predictors["escape_guideline"].sum())
    full_rules = int(predictors["escape_full_rules"].sum())
    return {
        "scoreable": total,
        "escape_guideline": guideline,
        "escape_full_rules": full_rules,
        "guideline_fraction": fraction(guideline),
        "full_rules_fraction": fraction(full_rules),
        "disagree": len(disagree),
        "disagree_fraction": fraction(len(disagree)),
        "driven_by_start_proximal": int((start & ~long_exon).sum()),
        "driven_by_long_exon": int((~start & long_exon).sum()),
        "driven_by_both": int((start & long_exon).sum()),
    }
