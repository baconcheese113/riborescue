"""Tabular contracts — the schemas every table crossing a pipeline boundary is checked against.

Schemas are strict and never coerce: a column of the wrong type is a defect in whatever produced it,
not something to repair on read. Measured quantities that can be absent are nullable, so a missing
replicate standard deviation stays missing rather than arriving as zero.
"""

from pathlib import Path

import pandas as pd
import pandera.pandas as pa
from pandera.typing import Series

from riborescue.contracts import Consequence, StopCodon, TriageClass

__all__ = [
    "PathogenicNonsense",
    "ReadthroughLabels",
    "TriageInput",
    "TriageOutput",
    "read_table",
    "write_table",
]

_STOP_CODONS = [c.value for c in StopCodon]
_CONSEQUENCES = [c.value for c in Consequence]
_TRIAGE_CLASSES = [c.value for c in TriageClass]


class _Strict(pa.DataFrameModel):
    class Config(pa.DataFrameModel.Config):
        strict = True
        coerce = False


class ReadthroughLabels(_Strict):
    """Measured readthrough efficiency per variant × therapy — the training and parity labels.

    `censored` marks measurements at the assay's readthrough ceiling, which are bounds rather than
    values and must be modelled as such.
    """

    variant_id: Series[str]
    gene_id: Series[str]
    therapy_id: Series[str]
    stop_codon: Series[str] = pa.Field(isin=_STOP_CODONS)
    readthrough_efficiency: Series[float] = pa.Field(ge=0.0, le=1.0)
    replicate_sd: Series[float] = pa.Field(ge=0.0, nullable=True)
    censored: Series[bool]


class PathogenicNonsense(_Strict):
    """Pathogenic nonsense substitutions drawn from ClinVar.

    `review_stars` is ClinVar's own confidence rating, carried through so a one-star assertion and a
    guideline-backed one stay distinguishable rather than being flattened into one population.
    """

    variant_id: Series[str]
    allele_id: Series[int]
    chrom: Series[str]
    pos: Series[int] = pa.Field(gt=0)
    ref: Series[str] = pa.Field(isin=list("ACGT"))
    alt: Series[str] = pa.Field(isin=list("ACGT"))
    gene_symbol: Series[str]
    gene_id: Series[int]
    clinical_significance: Series[str]
    review_status: Series[str]
    review_stars: Series[int] = pa.Field(ge=0, le=4)
    conditions: Series[str]


class TriageInput(_Strict):
    """Variants submitted for triage."""

    variant_id: Series[str]
    consequence: Series[str] = pa.Field(isin=_CONSEQUENCES)
    transcript_supported: Series[bool]


class TriageOutput(TriageInput):
    """Triaged variants, each carrying the verdict and the reason behind it."""

    triage_class: Series[str] = pa.Field(isin=_TRIAGE_CLASSES)
    applies: Series[bool]
    reason: Series[str]


def read_table(path: Path) -> pd.DataFrame:
    """Read a tab-separated table. Whatever it holds is what the schema then judges.

    Types are inferred over the whole column rather than chunk by chunk: a chromosome column holding
    both 7 and X would otherwise come back as integers in one chunk and strings in another, and a
    table would fail to validate on the way in having passed on the way out.
    """

    return pd.read_csv(path, sep="\t", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)
