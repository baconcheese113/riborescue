"""Data contracts — the validation boundary for every RiboRescue record.

Contracts refuse invalid state rather than coercing it. Four rules are enforced here rather than
left to discipline: the canonical row key carries variant, transcript, biological context and
therapy together; the four rescue factors never collapse into one number; a missing quantity never
defaults to zero; and cross-modality ranking requires an explicit acknowledgement that the scales
are incomparable.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum

from pydantic import BaseModel, ConfigDict, model_validator

from riborescue._version import CONTRACTS_VERSION

__all__ = [
    "CONTRACTS_VERSION",
    "REPORTER_CONTEXT_NT",
    "AenmdResult",
    "BiologicalContext",
    "ConfidenceTier",
    "Consequence",
    "EvalConfig",
    "FieldStatus",
    "IncomparableModalitiesError",
    "Measurement",
    "MissingReason",
    "Modality",
    "Namespace",
    "NmdResult",
    "RowKey",
    "ScoreRecord",
    "StopCodon",
    "TherapyIdentity",
    "TranscriptIdentity",
    "TriageClass",
    "UnresolvedReporterWindowError",
    "VariantIdentity",
    "WindowSpec",
    "enforce_modality_comparability",
]


class StopCodon(StrEnum):
    uaa = "UAA"
    uag = "UAG"
    uga = "UGA"


class Modality(StrEnum):
    small_molecule = "small_molecule"
    suppressor_trna = "suppressor_tRNA"


class FieldStatus(StrEnum):
    present = "present"
    missing = "missing"


class MissingReason(StrEnum):
    not_applicable = "not_applicable"
    not_available = "not_available"
    tool_failed = "tool_failed"
    transcript_unmapped = "transcript_unmapped"
    insufficient_coverage = "insufficient_coverage"
    unsupported_tissue = "unsupported_tissue"


class EvalConfig(StrEnum):
    published_random_cv = "published_random_cv"
    grouped_by_gene = "grouped_by_gene"
    grouped_by_sequence_cluster = "grouped_by_sequence_cluster"
    external_transfer_test = "external_transfer_test"


class Namespace(StrEnum):
    sequence = "sequence"
    kinetics = "kinetics"
    nmd = "nmd"
    residue = "residue"
    protein = "protein"
    safety = "safety"


class ConfidenceTier(StrEnum):
    mechanistic = "mechanistic"
    experimentally_anchored = "experimentally_anchored"
    uncertain = "uncertain"


class Consequence(StrEnum):
    stop_gained = "stop_gained"
    missense_variant = "missense_variant"
    frameshift_variant = "frameshift_variant"
    splice_acceptor_variant = "splice_acceptor_variant"
    splice_donor_variant = "splice_donor_variant"
    splice_region_variant = "splice_region_variant"
    synonymous_variant = "synonymous_variant"
    stop_lost = "stop_lost"
    inframe_deletion = "inframe_deletion"
    inframe_insertion = "inframe_insertion"


class TriageClass(StrEnum):
    supported_nonsense = "supported_nonsense"
    missense = "missense"
    frameshift = "frameshift"
    splice = "splice"
    synonymous = "synonymous"
    stop_loss = "stop_loss"
    unsupported_transcript = "unsupported_transcript"
    unsupported = "unsupported"


class Measurement(BaseModel):
    """A number that knows whether it is present, and if absent, why.

    A missing measurement carries a reason and never a value, so absence can never be silently read
    as zero.
    """

    model_config = ConfigDict(frozen=True)

    value: float | None = None
    status: FieldStatus
    reason: MissingReason | None = None

    @model_validator(mode="after")
    def _consistent(self) -> "Measurement":
        if self.status is FieldStatus.present:
            if self.value is None:
                raise ValueError("a present measurement requires a value")
            if self.reason is not None:
                raise ValueError("a present measurement cannot carry a missing reason")
        else:
            if self.reason is None:
                raise ValueError("a missing measurement requires a reason")
            if self.value is not None:
                raise ValueError("a missing measurement cannot carry a value")
        return self

    @classmethod
    def present(cls, value: float) -> "Measurement":
        return cls(value=value, status=FieldStatus.present)

    @classmethod
    def absent(cls, reason: MissingReason) -> "Measurement":
        return cls(status=FieldStatus.missing, reason=reason)


class VariantIdentity(BaseModel):
    model_config = ConfigDict(frozen=True)

    genome_build: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene_id: str | None = None
    gene_symbol: str | None = None
    hgvs_g: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    original_aa: str | None = None
    protein_position: int | None = None
    stop_codon: StopCodon | None = None
    exon_number: int | None = None
    mane_select: bool = False


class TranscriptIdentity(BaseModel):
    model_config = ConfigDict(frozen=True)

    transcript_id: str
    version: int
    is_mane_select: bool = False


class BiologicalContext(BaseModel):
    """Cell line or tissue a score was derived in — part of the row key, never metadata.

    A score derived from HEK293 must never surface as a fibroblast or liver score.
    """

    model_config = ConfigDict(frozen=True)

    context_id: str
    description: str
    source: str


class TherapyIdentity(BaseModel):
    model_config = ConfigDict(frozen=True)

    therapy_id: str
    modality: Modality
    name: str
    recognized_stop: StopCodon
    inserted_aa: str | None = None
    insertion_distribution: tuple[tuple[str, float], ...] | None = None
    evidence: str | None = None
    tissue_context: str | None = None

    @model_validator(mode="after")
    def _insertion_defined(self) -> "TherapyIdentity":
        if self.modality is Modality.suppressor_trna and self.inserted_aa is None:
            raise ValueError("a suppressor tRNA must declare its inserted amino acid")
        if self.modality is Modality.small_molecule and self.insertion_distribution is None:
            raise ValueError("a small molecule must declare an insertion distribution")
        return self


class RowKey(BaseModel):
    """The canonical key: variant × transcript × biological context × therapy.

    All four are required; none is optional metadata.
    """

    model_config = ConfigDict(frozen=True)

    variant: VariantIdentity
    transcript: TranscriptIdentity
    context: BiologicalContext
    therapy: TherapyIdentity


class ScoreRecord(BaseModel):
    """One scored row. The four rescue factors are kept separate and never combined here."""

    model_config = ConfigDict(frozen=True)

    key: RowKey
    modality: Modality
    transcript_availability: Measurement
    readthrough_efficiency: Measurement
    residue_compatibility: Measurement
    protein_function: Measurement
    confidence_tier: ConfidenceTier
    model_disagreement: Measurement
    evidence_provenance: str
    limitations: str


class AenmdResult(BaseModel):
    model_config = ConfigDict(frozen=True)

    escape_predicted: bool
    triggered_rules: tuple[str, ...]


class NmdResult(BaseModel):
    """Native outputs of the three NMD predictors — aenmd stays a rule set, never a score."""

    model_config = ConfigDict(frozen=True)

    aenmd: AenmdResult
    prednmd_score: Measurement
    nmdetective_ai_score: Measurement
    concordance: Measurement


REPORTER_CONTEXT_NT: int | None = None
"""Reporter context length in nucleotides: 147 (72 + 3 + 72) versus 150 (75 + 75).

Unresolved until read from the oligo design table. An off-by-three shifts every positional feature,
so `WindowSpec` refuses to construct while this is None.
"""


class UnresolvedReporterWindowError(RuntimeError):
    pass


@dataclass(frozen=True)
class WindowSpec:
    """A feature window around a premature stop codon.

    Cannot be constructed until the reporter context length is resolved, gating all feature work on
    a decision that must be made from the oligo design table rather than assumed.
    """

    upstream_nt: int
    downstream_nt: int

    def __post_init__(self) -> None:
        if REPORTER_CONTEXT_NT is None:
            raise UnresolvedReporterWindowError(
                "reporter context length is unresolved (147 = 72+3+72 vs 150 = 75+75); "
                "resolve it from the oligo design table before constructing feature windows"
            )
        total = self.upstream_nt + 3 + self.downstream_nt
        if total != REPORTER_CONTEXT_NT:
            raise ValueError(
                f"window of {total} nt does not match "
                f"the reporter context of {REPORTER_CONTEXT_NT} nt"
            )


class IncomparableModalitiesError(RuntimeError):
    pass


def enforce_modality_comparability(
    modalities: Iterable[Modality], *, incomparable_scales: bool = False
) -> None:
    """Reject a ranking that spans modalities without acknowledging their scales differ.

    Small-molecule and suppressor-tRNA labels come from different assays on different scales;
    ranking across them implies a calibration that does not exist.
    """

    if len(set(modalities)) > 1 and not incomparable_scales:
        raise IncomparableModalitiesError(
            "ranking spans modalities with no shared functional evidence; "
            "set incomparable_scales to acknowledge that the scores are not comparable"
        )
