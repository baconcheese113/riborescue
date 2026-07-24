"""Coverage counted over diseases, with the denominator it is a fraction of.

A disease is not covered because one of its variants is. For each MedGen disease this reports the
eligible nonsense-variant denominator, how many are model-covered and the fraction that is, the
genes and designs that carry it, and how complete its cross-references are. "At least one variant
covered" is kept as *reach* — a separate, weaker flag — never as the coverage number.

Model coverage is the exact-restoration predicate of ADR-0014: a suppressor design decodes the
variant's stop and reinserts the native residue. Under it a variant is covered exactly when it is
scoreable — a premature stop with a determinable residue — so a disease's covered count is the
number of its eligible variants that reach a scoreable context. It is not a therapeutic or clinical
claim, and makes no statement about unmet need, which would need a treatment-status source it lacks.
"""

from __future__ import annotations

from dataclasses import asdict

import pandas as pd

from riborescue.variants.panels import greedy_panel

__all__ = [
    "disease_coverage",
    "disease_reach_frontier",
]

# Conditions that carry a real MedGen concept and can key a disease; placeholders and the
# no-MedGen rows are excluded here because they are not diseases (ADR-0015).
_DISEASE_STATUS = ("mapped", "medgen_only")

_COVERAGE_COLUMNS = [
    "medgen",
    "disease_name",
    "omim",
    "orphanet",
    "eligible_variants",
    "model_covered",
    "covered_fraction",
    "genes",
    "designs_contributing",
    "designs",
    "reach",
    "mapping_completeness",
]


def _designs(contexts: pd.DataFrame) -> pd.DataFrame:
    """variant_id, whether it is covered, and the design covering it — from the contexts table."""

    covered = contexts["scoreable"].fillna(False).astype(bool)
    design = [
        f"{str(stop).upper()}-{aa}" if is_covered else ""
        for stop, aa, is_covered in zip(
            contexts["stop_type"], contexts["original_aa"], covered, strict=True
        )
    ]
    return pd.DataFrame(
        {"variant_id": contexts["variant_id"], "covered": covered, "design_id": design}
    )


def disease_coverage(diseases: pd.DataFrame, contexts: pd.DataFrame) -> pd.DataFrame:
    """Per-disease coverage with its eligible denominator, one row per MedGen concept.

    `diseases` is the normalized table (ADR-0015); `contexts` says which variants are scoreable and
    at which stop and residue. Only real diseases are kept — placeholder and no-MedGen conditions
    are not diseases. The fraction is model-covered over eligible, never a bare "reached" flag.
    """

    keyed = diseases[diseases["mapping_status"].isin(_DISEASE_STATUS)]
    merged = keyed.merge(_designs(contexts), on="variant_id", how="left")
    merged["covered"] = merged["covered"].fillna(False).astype(bool)

    rows = []
    for medgen, group in merged.groupby("medgen"):
        hit = group[group["covered"]]
        designs = sorted({d for d in hit["design_id"] if d})
        omim = next((x for x in group["omim"] if x), "")
        orphanet = next((x for x in group["orphanet"] if x), "")
        eligible = group["variant_id"].nunique()
        model_covered = hit["variant_id"].nunique()
        names = group["condition_name"].mode()
        rows.append(
            {
                "medgen": medgen,
                "disease_name": names.iloc[0] if not names.empty else "",
                "omim": omim,
                "orphanet": orphanet,
                "eligible_variants": eligible,
                "model_covered": model_covered,
                "covered_fraction": round(model_covered / eligible, 4) if eligible else 0.0,
                "genes": group["gene_symbol"].nunique(),
                "designs_contributing": len(designs),
                "designs": ";".join(designs),
                "reach": bool(model_covered >= 1),
                "mapping_completeness": "mapped" if (omim or orphanet) else "medgen_only",
            }
        )
    frame = pd.DataFrame(rows, columns=_COVERAGE_COLUMNS)
    return frame.sort_values("eligible_variants", ascending=False, ignore_index=True)


def disease_reach_frontier(diseases: pd.DataFrame, contexts: pd.DataFrame) -> pd.DataFrame:
    """The greedy panel that *reaches* the most diseases — a disease with any covered variant.

    This is the reach frontier, and it is labelled that way: a design reaches a disease when it
    restores at least one of the disease's variants. It answers "which designs touch the most
    distinct diseases", not "which diseases are fully covered", and reuses the ADR-0014 engine.
    """

    keyed = diseases[diseases["mapping_status"].isin(_DISEASE_STATUS)]
    scoreable = contexts[contexts["scoreable"].fillna(False).astype(bool)]
    reached = keyed.merge(_designs(scoreable), on="variant_id", how="inner")

    sets: dict[str, set[str]] = {}
    for design, medgen in zip(reached["design_id"], reached["medgen"], strict=True):
        if design:
            sets.setdefault(design, set()).add(medgen)
    frozen = {design: frozenset(members) for design, members in sets.items()}
    return pd.DataFrame([asdict(step) for step in greedy_panel(frozen)])
