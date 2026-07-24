"""ClinVar's pathogenic nonsense variants — the population RiboRescue scores.

A variant qualifies on three counts, each read from the record rather than inferred: it is
classified pathogenic or likely pathogenic, its molecular consequence includes a premature stop,
and it is a single-nucleotide substitution. Everything else — a nonsense change of uncertain
significance, a pathogenic frameshift, a deletion that happens to remove a stop — is excluded here
rather than filtered downstream, so the population being scored is defined in one place.

Review status is carried through as ClinVar's own star rating. A one-star assertion and a
guideline-backed one are both pathogenic, and the difference belongs in the output rather than in a
threshold applied here.
"""

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from cyvcf2 import VCF

__all__ = [
    "NONSENSE_SO",
    "PATHOGENIC_CLASSES",
    "REVIEW_STARS",
    "NonsenseVariants",
    "pathogenic_nonsense",
]

BASES = frozenset("ACGT")

NONSENSE_SO = "SO:0001587"
"""Sequence Ontology term for a stop-gained (nonsense) consequence."""

PATHOGENIC_CLASSES = frozenset({"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"})
"""Classifications RiboRescue treats as disease-causing.

Conflicting and uncertain classifications are excluded: a therapy match for a variant that may not
cause disease is a claim the evidence does not support.
"""

REVIEW_STARS: dict[str, int] = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
}
"""ClinVar's review status as its published star rating; anything unlisted carries no assertion."""


@dataclass(frozen=True)
class NonsenseVariants:
    """The scoreable variants, and a count of those an ambiguous allele put out of reach.

    A few ClinVar records give the alternate allele as an IUPAC ambiguity code, so which base
    substitutes is unknown and the resulting stop codon cannot be determined. They are excluded and
    counted rather than dropped quietly, because a shrinking population that nobody notices is how a
    coverage claim becomes wrong.
    """

    variants: pd.DataFrame
    ambiguous_alleles: int


def _records(path: Path) -> Iterator[dict[str, object]]:
    for record in VCF(str(path)):
        info = record.INFO
        if info.get("CLNVC") != "single_nucleotide_variant":
            continue
        if NONSENSE_SO not in str(info.get("MC", "")):
            continue
        significance = str(info.get("CLNSIG", ""))
        if significance.split("|")[0] not in PATHOGENIC_CLASSES:
            continue
        if not record.ALT or len(record.ALT) != 1:
            continue

        symbol, _, gene_id = str(info.get("GENEINFO", "")).split("|")[0].partition(":")
        if not gene_id.isdigit():
            continue
        review = str(info.get("CLNREVSTAT", ""))
        yield {
            "variant_id": str(info.get("CLNHGVS", "")),
            "allele_id": int(info.get("ALLELEID", 0)),
            "chrom": record.CHROM,
            "pos": record.POS,
            "ref": record.REF,
            "alt": record.ALT[0],
            "gene_symbol": symbol,
            "gene_id": int(gene_id),
            "clinical_significance": significance,
            "review_status": review,
            "review_stars": REVIEW_STARS.get(review, 0),
            "conditions": str(info.get("CLNDN", "")),
        }


def pathogenic_nonsense(path: Path) -> NonsenseVariants:
    """Read the pathogenic nonsense substitutions out of a ClinVar VCF, in file order."""

    found = pd.DataFrame(_records(path))
    unambiguous = found["ref"].isin(BASES) & found["alt"].isin(BASES)
    return NonsenseVariants(
        variants=found[unambiguous].reset_index(drop=True),
        ambiguous_alleles=int((~unambiguous).sum()),
    )
