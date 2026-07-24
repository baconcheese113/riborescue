"""Normalizing each scored variant's ClinVar conditions to disease identifiers.

ClinVar states a variant's conditions as names (`CLNDN`) with a parallel list of cross-references
(`CLNDISDB`), one entry per condition. Each entry is a comma-separated set of `Source:ID` pairs —
MedGen, OMIM, Orphanet, MONDO, MeSH. This turns those two parallel lists into one row per
variant-condition, keyed on the MedGen concept, with the other sources kept as attached references.

Nothing licensed is imported: these are the identifiers ClinVar itself publishes, from the same
pinned release the variant set is read from. A condition that is a placeholder ("not provided",
"not specified") or that resolves to MedGen alone is kept and labelled, not dropped — its
completeness is reported, never silently filtered. See ADR-0015.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

__all__ = [
    "Condition",
    "normalize_conditions",
    "parse_conditions",
]

# MedGen concept ids for ClinVar's non-disease placeholders — a real concept, but not a disease.
# C3661900 is the trap: a plain CUI, not a CN* id, that ClinVar uses for "not provided", so it slips
# past a filter that only knows the CN* placeholders and becomes a 25,000-variant phantom disease.
_PLACEHOLDER_MEDGEN = {
    "CN517202": "not provided",
    "CN169374": "not specified",
    "C3661900": "not provided",
}
_SOURCES = ("MedGen", "OMIM", "Orphanet", "MONDO", "MeSH")


@dataclass(frozen=True)
class Condition:
    """One asserted condition and its cross-references, drawn from a single CLNDN/CLNDISDB pair."""

    name: str
    medgen: str
    omim: str
    orphanet: str
    mondo: str
    mesh: str
    mapping_status: str
    reason: str


def _xrefs(entry: str) -> dict[str, list[str]]:
    """The `Source:ID` pairs in one CLNDISDB entry, grouped by source in first-seen order.

    A source can appear more than once — OMIM lists a phenotype and its phenotypic series — so every
    id is kept. Orphanet ids come as `Orphanet:791`; MONDO doubles its prefix
    (`MONDO:MONDO:0019200`), so the split is on the first colon only.
    """

    found: dict[str, list[str]] = {source: [] for source in _SOURCES}
    for token in entry.split(","):
        source, _, ident = token.partition(":")
        if source in found and ident and ident not in found[source]:
            found[source].append(ident)
    return found


def parse_conditions(clndn: str, clndisdb: str) -> list[Condition]:
    """The conditions for one variant, from its parallel CLNDN and CLNDISDB fields.

    The two are pipe-separated and index-aligned. A condition with no MedGen id, or one whose MedGen
    id is a ClinVar placeholder, is kept with a `mapping_status` that says why it is not a disease,
    so a downstream denominator can exclude it deliberately rather than by accident.
    """

    names = clndn.split("|") if clndn else []
    entries = clndisdb.split("|") if clndisdb else []
    conditions: list[Condition] = []
    for index, name in enumerate(names):
        entry = entries[index] if index < len(entries) else ""
        refs = _xrefs(entry)
        medgen = refs["MedGen"][0] if refs["MedGen"] else ""
        omim, orphanet = ";".join(refs["OMIM"]), ";".join(refs["Orphanet"])
        if medgen in _PLACEHOLDER_MEDGEN:
            status, reason = "placeholder", f"ClinVar placeholder ({_PLACEHOLDER_MEDGEN[medgen]})"
        elif not medgen:
            status, reason = "unmapped", "no MedGen concept on the condition"
        elif not omim and not orphanet:
            status, reason = (
                "medgen_only",
                "MedGen concept with no OMIM or Orphanet cross-reference",
            )
        else:
            status, reason = "mapped", ""
        conditions.append(
            Condition(
                name=name,
                medgen=medgen,
                omim=omim,
                orphanet=orphanet,
                mondo=";".join(refs["MONDO"]),
                mesh=";".join(refs["MeSH"]),
                mapping_status=status,
                reason=reason,
            )
        )
    return conditions


def normalize_conditions(variants: pd.DataFrame) -> pd.DataFrame:
    """One row per variant-condition, keyed on MedGen, from the scored nonsense variants.

    `variants` must carry `variant_id`, `gene_symbol`, `conditions` (CLNDN) and `condition_xrefs`
    (CLNDISDB). A variant asserting several conditions yields several rows; the mapping never picks
    one disease for it. The result is the join key for disease-level coverage.
    """

    rows = []
    for variant in variants.itertuples():
        for condition in parse_conditions(str(variant.conditions), str(variant.condition_xrefs)):
            rows.append(
                {
                    "variant_id": variant.variant_id,
                    "gene_symbol": variant.gene_symbol,
                    "condition_name": condition.name,
                    "medgen": condition.medgen,
                    "omim": condition.omim,
                    "orphanet": condition.orphanet,
                    "mondo": condition.mondo,
                    "mesh": condition.mesh,
                    "mapping_status": condition.mapping_status,
                    "reason": condition.reason,
                }
            )
    return pd.DataFrame(rows, columns=_DISEASE_COLUMNS)


_DISEASE_COLUMNS = [
    "variant_id",
    "gene_symbol",
    "condition_name",
    "medgen",
    "omim",
    "orphanet",
    "mondo",
    "mesh",
    "mapping_status",
    "reason",
]
