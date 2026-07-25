"""The compact aggregate the researcher dashboard reads — coverage, not per-variant rows.

The dashboard never loads the variant table; it loads this. It carries the three coverage frontiers
(variants, genes, condition entities), the per-entity coverage with its denominators, the mapping
counts, and the provenance every number rests on. Everything is an aggregate with the denominator it
was computed over, so no figure on the page can be read without the count behind it.

A MedGen concept is a ClinVar condition — a disease, but also possibly a finding, susceptibility, or
broad label — so these are condition entities, not verified diseases. Model coverage is the exact-
restoration predicate (ADR-0014); the condition frontier is a *reach* frontier (ADR-0015). Nothing
states unmet need, which would need a treatment-status source this does not have.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import pandas as pd

from riborescue.variants.disease_coverage import disease_coverage, disease_reach_frontier
from riborescue.variants.nmd import disagreement_atlas, nmd_predictors
from riborescue.variants.panels import coverage_frontier

__all__ = [
    "ResearchAggregate",
    "build_research_aggregate",
]


def _frontier(frame: pd.DataFrame) -> list[dict]:
    """A coverage frontier trimmed to what a curve needs: size, design, cumulative, marginal."""

    out = []
    for i in range(len(frame)):
        row = frame.iloc[i]
        out.append(
            {
                "rank": int(row["rank"]),
                "design_id": row["design_id"],
                "cumulative": int(row["cumulative"]),
                "cumulative_fraction": round(float(row["cumulative_fraction"]), 4),
                "marginal": int(row["marginal"]),
            }
        )
    return out


@dataclass(frozen=True)
class ResearchAggregate:
    """The whole dashboard payload: provenance, frontiers, and per-condition-entity coverage."""

    provenance: dict
    mapping_completeness: dict
    reach_denominator: dict
    nmd: dict
    frontiers: dict
    condition_coverage_top: list[dict]
    caveats: dict

    def to_json(self) -> str:
        # allow_nan=False so a stray NaN fails here rather than becoming a token no browser parses.
        return json.dumps(
            {
                "provenance": self.provenance,
                "mapping_completeness": self.mapping_completeness,
                "reach_denominator": self.reach_denominator,
                "nmd": self.nmd,
                "frontiers": self.frontiers,
                "condition_coverage_top": self.condition_coverage_top,
                "caveats": self.caveats,
            },
            indent=2,
            allow_nan=False,
        )


def build_research_aggregate(
    diseases: pd.DataFrame,
    contexts: pd.DataFrame,
    *,
    clinvar_release: str,
    commit: str = "",
    top: int = 50,
) -> ResearchAggregate:
    """Assemble the researcher aggregate from the disease and context tables.

    `top` bounds the per-entity list to the largest by eligible denominator, so the payload stays
    small; the frontiers and completeness counts are whole-set. `clinvar_release` and `commit` are
    the provenance the whole page rests on.
    """

    coverage = disease_coverage(diseases, contexts)
    completeness = {
        status: int(count) for status, count in diseases["mapping_status"].value_counts().items()
    }
    scoreable = int(contexts["scoreable"].fillna(False).astype(bool).sum())

    eligible = diseases[diseases["mapping_status"].isin(("mapped", "medgen_only"))]
    diseases_frontier = _frontier(disease_reach_frontier(diseases, contexts))
    reachable = diseases_frontier[-1]["cumulative"] if diseases_frontier else 0

    head = coverage.head(top)
    top_rows = []
    for i in range(len(head)):
        row = head.iloc[i]
        top_rows.append(
            {
                "condition_name": row["condition_name"],
                "medgen": row["medgen"],
                "omim": row["omim"],
                "orphanet": row["orphanet"],
                "eligible_variants": int(row["eligible_variants"]),
                "model_covered": int(row["model_covered"]),
                "covered_fraction": round(float(row["covered_fraction"]), 4),
                "reach": bool(row["reach"]),
                "complete": bool(row["complete"]),
                "genes": int(row["genes"]),
                "designs_contributing": int(row["designs_contributing"]),
                "mapping_completeness": row["mapping_completeness"],
            }
        )

    return ResearchAggregate(
        provenance={
            "clinvar_release": clinvar_release,
            "commit": commit,
            "qualifying_variants": int(contexts["variant_id"].nunique()),
            "scoreable_variants": scoreable,
            "condition_entities": len(coverage),
        },
        mapping_completeness=completeness,
        nmd=disagreement_atlas(nmd_predictors(contexts)),
        reach_denominator={
            "eligible_condition_entities": int(eligible["medgen"].nunique()),
            "unique_variant_condition_pairs": len(eligible),
            "excluded_placeholder_or_unmapped_rows": int(len(diseases) - len(eligible)),
            "reachable_entities": int(reachable),
            "note": (
                "Entities are deduplicated to distinct MedGen concepts; reaching every reachable "
                "entity is partly a closure property of the design universe, not a therapeutic "
                "result."
            ),
        },
        frontiers={
            "variants": _frontier(coverage_frontier(contexts, "variants")),
            "genes": _frontier(coverage_frontier(contexts, "genes")),
            "conditions": diseases_frontier,
        },
        condition_coverage_top=top_rows,
        caveats={
            "entities": (
                "A MedGen concept is a ClinVar condition — a disease, but possibly a finding, "
                "susceptibility, or broad umbrella label. These are condition entities, not "
                "verified diseases."
            ),
            "coverage": (
                "Model coverage is exact restoration — a suppressor design decodes the stop and "
                "reinserts the native residue. Not a therapeutic or clinical claim."
            ),
            "three_metrics": (
                "reach = at least one eligible variant covered; covered_fraction = covered over "
                "all eligible variants for the entity; complete = every eligible variant covered."
            ),
            "cross_references": (
                "OMIM and Orphanet are ClinVar-provided cross-reference identifiers only — they "
                "add no prevalence, treatment availability, or unmet-need evidence."
            ),
            "unmet_need": (
                "No unmet-need claim is made; that needs a treatment-status source not used here."
            ),
            "nmd": (
                "NMD escape is two rule-based predictors — the 50-nt guideline and the fuller "
                "Lindeboom rule set — not yet the ML models (NMDetective-AI, predNMD). A "
                "disagreement is the guideline's blind spot, not model uncertainty."
            ),
        },
    )
