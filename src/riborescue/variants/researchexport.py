"""The compact aggregate the researcher dashboard reads — coverage, not per-variant rows.

The dashboard never loads the variant table; it loads this. It carries the three coverage frontiers
(variants, genes, diseases), the per-disease coverage with its denominators, the mapping counts, and
the provenance every number rests on. Everything is an aggregate with the denominator it was
computed over, so no figure on the page can be read without the count behind it.

Model coverage here is the exact-restoration predicate (ADR-0014); the disease frontier is a *reach*
frontier (ADR-0015 coverage semantics), and both are labelled as such. Nothing states unmet
therapeutic need, which would need a treatment-status source this does not have.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import pandas as pd

from riborescue.variants.disease_coverage import disease_coverage, disease_reach_frontier
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
    """The whole dashboard payload: provenance, frontiers, and per-disease coverage."""

    provenance: dict
    mapping_completeness: dict
    frontiers: dict
    disease_coverage_top: list[dict]
    caveats: dict

    def to_json(self) -> str:
        # allow_nan=False so a stray NaN fails here rather than becoming a token no browser parses.
        return json.dumps(
            {
                "provenance": self.provenance,
                "mapping_completeness": self.mapping_completeness,
                "frontiers": self.frontiers,
                "disease_coverage_top": self.disease_coverage_top,
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

    `top` bounds the per-disease list to the largest by eligible denominator, so the payload stays
    small; the frontiers and completeness counts are whole-set. `clinvar_release` and `commit` are
    the provenance the whole page rests on.
    """

    coverage = disease_coverage(diseases, contexts)
    completeness = {
        status: int(count) for status, count in diseases["mapping_status"].value_counts().items()
    }
    scoreable = int(contexts["scoreable"].fillna(False).astype(bool).sum())

    head = coverage.head(top)
    top_rows = []
    for i in range(len(head)):
        row = head.iloc[i]
        top_rows.append(
            {
                "disease_name": row["disease_name"],
                "medgen": row["medgen"],
                "omim": row["omim"],
                "orphanet": row["orphanet"],
                "eligible_variants": int(row["eligible_variants"]),
                "model_covered": int(row["model_covered"]),
                "covered_fraction": round(float(row["covered_fraction"]), 4),
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
            "diseases": len(coverage),
        },
        mapping_completeness=completeness,
        frontiers={
            "variants": _frontier(coverage_frontier(contexts, "variants")),
            "genes": _frontier(coverage_frontier(contexts, "genes")),
            "diseases": _frontier(disease_reach_frontier(diseases, contexts)),
        },
        disease_coverage_top=top_rows,
        caveats={
            "coverage": (
                "Model coverage is exact restoration — a suppressor design decodes the stop and "
                "reinserts the native residue. Not a therapeutic or clinical claim."
            ),
            "disease_frontier": (
                "A reach frontier: a design reaches a disease when it restores at least one of its "
                "variants. Not a claim the disease is fully covered."
            ),
            "unmet_need": (
                "No unmet-need claim is made; that needs a treatment-status source not used here."
            ),
        },
    )
