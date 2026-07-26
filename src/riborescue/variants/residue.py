"""Which amino acid a readthrough event inserts, and how far it is from the one it replaces.

A small molecule promotes an endogenous near-cognate tRNA, so the genetic code fixes which residues
are available; an engineered suppressor carries whichever of the twenty it was charged with.

Compatibility is BLOSUM62, a hypothesis about tolerance that knows nothing of a position's
structural or catalytic role. It ranks candidates; it does not establish that any restores function.
"""

from dataclasses import dataclass
from typing import Any, cast

import pandas as pd
from Bio.Align import substitution_matrices
from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

from riborescue.core.sequences import as_dna

__all__ = [
    "AMINO_ACIDS",
    "NEAR_COGNATE",
    "SuppressorDesign",
    "conservative",
    "coverage_by_design",
    "near_cognate_residues",
    "similarity",
]

_BLOSUM62 = substitution_matrices.load("BLOSUM62")
_BASES = "ACGT"

AMINO_ACIDS = tuple(sorted(set(standard_dna_table.forward_table.values())))
"""The twenty amino acids a suppressor tRNA can be charged with."""


def near_cognate_residues(stop_codon: str) -> tuple[str, ...]:
    """Amino acids encoded by codons one base from the stop, in DNA to match the codon table."""

    codon = as_dna(stop_codon)
    found = set()
    for index, original in enumerate(codon):
        for base in _BASES:
            if base == original:
                continue
            candidate = codon[:index] + base + codon[index + 1 :]
            residue = str(Seq(candidate).translate())
            if residue != "*":
                found.add(residue)
    return tuple(sorted(found))


NEAR_COGNATE: dict[str, tuple[str, ...]] = {
    stop: near_cognate_residues(stop) for stop in ("TAA", "TAG", "TGA")
}
"""What each stop codon's near-cognate tRNAs can insert."""


def similarity(original: str, inserted: str) -> int:
    """BLOSUM62 score for replacing one residue with another; higher is more interchangeable."""

    # Biopython's substitution matrix is indexed by a residue pair; it ships no type information.
    return int(cast(Any, _BLOSUM62)[original, inserted])


def conservative(original: str, inserted: str) -> bool:
    """Whether the substitution is one evolution accepts more often than chance."""

    return similarity(original, inserted) >= 0


@dataclass(frozen=True)
class SuppressorDesign:
    """A suppressor tRNA is defined by the stop it decodes and the residue it inserts."""

    recognized_stop: str
    inserted_aa: str

    @property
    def design_id(self) -> str:
        return f"{self.recognized_stop}-{self.inserted_aa}"


def coverage_by_design(contexts: pd.DataFrame) -> pd.DataFrame:
    """How many variants each suppressor design reaches, one row per stop and inserted residue.

    `restores_exactly` counts variants whose original residue the design puts back; `conservative`
    counts those where the insertion is one evolution accepts more often than chance. Variants, not
    patients. Says nothing about the native stops the design would also read through.
    """

    rows = []
    for stop, at_stop in contexts.groupby("stop_type"):
        originals = at_stop["original_aa"].value_counts()
        for inserted in AMINO_ACIDS:
            rows.append(
                {
                    "design_id": f"{str(stop).upper()}-{inserted}",
                    "recognized_stop": str(stop).upper(),
                    "inserted_aa": inserted,
                    "variants_at_stop": len(at_stop),
                    "restores_exactly": int(originals.get(inserted, 0)),
                    "conservative": int(
                        sum(
                            n
                            for original, n in originals.items()
                            if conservative(str(original), inserted)
                        )
                    ),
                }
            )
    return pd.DataFrame(rows).sort_values("conservative", ascending=False, ignore_index=True)
