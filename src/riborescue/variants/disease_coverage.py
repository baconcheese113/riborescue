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

from riborescue.variants.suppressor_panels import greedy_panel

__all__ = [
    "disease_coverage",
    "disease_reach_frontier",
]

# Rows that carry a real MedGen concept and can key a condition entity; placeholders and the
# no-MedGen rows are excluded because they name no entity (ADR-0015). A MedGen concept is a ClinVar
# *condition* — a disease, but also possibly a finding, a susceptibility, or a broad label — so
# the counts here are condition entities, not verified diseases.
_ENTITY_STATUS = ("mapped", "medgen_only")

_COVERAGE_COLUMNS = [
    "medgen",
    "condition_name",
    "omim",
    "orphanet",
    "eligible_variants",
    "model_covered",
    "covered_fraction",
    "reach",
    "complete",
    "genes",
    "designs_contributing",
    "designs",
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
    """Per-condition-entity coverage with its eligible denominator, one row per MedGen concept.

    `diseases` is the normalized table (ADR-0015); `contexts` says which variants are scoreable and
    at which stop and residue. Only rows naming a real MedGen concept are kept — placeholder and
    no-MedGen rows name no entity. Three distinct metrics ride together and are never conflated:
    `reach` (at least one eligible variant covered), `covered_fraction` (the variant fraction:
    covered over eligible), and `complete` (every eligible variant covered).
    """

    keyed = diseases[diseases["mapping_status"].isin(_ENTITY_STATUS)].copy()
    # A TSV round-trip reads empty cells back as NaN, and NaN is truthy — it would slip through the
    # "first non-empty" pick below and reach the JSON as a token no browser parses. Coerce it to "".
    for column in ("condition_name", "medgen", "omim", "orphanet"):
        keyed[column] = keyed[column].fillna("")
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
                "condition_name": names.iloc[0] if not names.empty else "",
                "omim": omim,
                "orphanet": orphanet,
                "eligible_variants": eligible,
                "model_covered": model_covered,
                "covered_fraction": round(model_covered / eligible, 4) if eligible else 0.0,
                "reach": bool(model_covered >= 1),
                "complete": bool(eligible > 0 and model_covered == eligible),
                "genes": group["gene_symbol"].nunique(),
                "designs_contributing": len(designs),
                "designs": ";".join(designs),
                "mapping_completeness": "mapped" if (omim or orphanet) else "medgen_only",
            }
        )
    frame = pd.DataFrame(rows, columns=_COVERAGE_COLUMNS)
    return frame.sort_values("eligible_variants", ascending=False, ignore_index=True)


def disease_reach_frontier(diseases: pd.DataFrame, contexts: pd.DataFrame) -> pd.DataFrame:
    """The greedy panel that *reaches* the most condition entities — one with any covered variant.

    This is the reach frontier, and it is labelled that way: a design reaches an entity when it
    restores at least one of its variants. It answers "which designs touch the most distinct
    condition entities", not "which are fully covered", and reuses the ADR-0014 engine. The panel
    reaching every entity is partly a closure property of the design universe — every scoreable
    variant defines a restoring design — not a therapeutic result.
    """

    keyed = diseases[diseases["mapping_status"].isin(_ENTITY_STATUS)]
    scoreable = contexts[contexts["scoreable"].fillna(False).astype(bool)]
    reached = keyed.merge(_designs(scoreable), on="variant_id", how="inner")

    sets: dict[str, set[str]] = {}
    for design, medgen in zip(reached["design_id"], reached["medgen"], strict=True):
        if design:
            sets.setdefault(design, set()).add(medgen)
    frozen = {design: frozenset(members) for design, members in sets.items()}
    return pd.DataFrame([asdict(step) for step in greedy_panel(frozen)])
