import pytest
from pydantic import ValidationError

from riborescue.core.contracts import (
    CONTRACTS_VERSION,
    REPORTER_DOWNSTREAM_NT,
    REPORTER_UPSTREAM_NT,
    AenmdResult,
    BiologicalContext,
    ConfidenceTier,
    FieldStatus,
    IncomparableModalitiesError,
    Measurement,
    Modality,
    NmdResult,
    RowKey,
    ScoreRecord,
    StopCodon,
    TherapyIdentity,
    TranscriptIdentity,
    VariantIdentity,
    WindowSpec,
    enforce_modality_comparability,
)
from riborescue.core.contracts import MissingReason as Reason


def _variant() -> VariantIdentity:
    return VariantIdentity(genome_build="GRCh38", chrom="4", pos=979_000, ref="G", alt="A")


def _therapy() -> TherapyIdentity:
    return TherapyIdentity(
        therapy_id="uga-trp",
        modality=Modality.suppressor_trna,
        name="ACE-tRNA UGA-Trp",
        recognized_stop=StopCodon.uga,
        inserted_aa="W",
    )


def _row_key() -> RowKey:
    return RowKey(
        variant=_variant(),
        transcript=TranscriptIdentity(
            transcript_id="ENST00000000001", version=1, is_mane_select=True
        ),
        context=BiologicalContext(
            context_id="hek293", description="HEK293T", source="reporter assay"
        ),
        therapy=_therapy(),
    )


def test_contracts_version_looks_like_semver():
    assert CONTRACTS_VERSION.count(".") == 2


def test_present_measurement_requires_a_value():
    with pytest.raises(ValidationError):
        Measurement(status=FieldStatus.present)


def test_missing_measurement_requires_a_reason():
    with pytest.raises(ValidationError):
        Measurement(status=FieldStatus.missing)


def test_missing_measurement_cannot_be_zero():
    with pytest.raises(ValidationError):
        Measurement(value=0.0, status=FieldStatus.missing, reason=Reason.not_available)


def test_present_measurement_cannot_carry_a_reason():
    with pytest.raises(ValidationError):
        Measurement(value=0.3, status=FieldStatus.present, reason=Reason.not_available)


def test_measurement_helpers():
    assert Measurement.present(0.0).value == 0.0
    assert Measurement.absent(Reason.tool_failed).reason is Reason.tool_failed


def test_row_key_requires_all_four_components():
    with pytest.raises(ValidationError):
        RowKey(variant=_variant(), transcript=TranscriptIdentity(transcript_id="x", version=1))  # type: ignore[call-arg]


def test_score_record_keeps_the_four_factors_separate():
    fields = ScoreRecord.model_fields
    for factor in (
        "transcript_availability",
        "readthrough_efficiency",
        "residue_compatibility",
        "protein_function",
    ):
        assert factor in fields
    for banned in ("combined_score", "score", "fitness", "combined_efficacy", "overall"):
        assert banned not in fields


def test_score_record_constructs():
    record = ScoreRecord(
        key=_row_key(),
        modality=Modality.suppressor_trna,
        transcript_availability=Measurement.present(0.7),
        readthrough_efficiency=Measurement.present(0.4),
        residue_compatibility=Measurement.present(0.9),
        protein_function=Measurement.absent(Reason.not_available),
        confidence_tier=ConfidenceTier.mechanistic,
        model_disagreement=Measurement.present(0.1),
        evidence_provenance="reporter assay + ACE-tRNA panel",
        limitations="function layer is hypothesis-tier for this gene",
    )
    assert record.modality is Modality.suppressor_trna


def test_small_molecule_requires_an_insertion_distribution():
    with pytest.raises(ValidationError):
        TherapyIdentity(
            therapy_id="g418",
            modality=Modality.small_molecule,
            name="G418",
            recognized_stop=StopCodon.uga,
        )


def test_suppressor_trna_requires_an_inserted_residue():
    with pytest.raises(ValidationError):
        TherapyIdentity(
            therapy_id="uga-x",
            modality=Modality.suppressor_trna,
            name="unnamed",
            recognized_stop=StopCodon.uga,
        )


def test_nmd_keeps_aenmd_as_rules():
    result = NmdResult(
        aenmd=AenmdResult(escape_predicted=True, triggered_rules=("last_exon", "start_proximal")),
        prednmd_score=Measurement.present(0.2),
        nmdetective_ai_score=Measurement.present(0.31),
        concordance=Measurement.present(0.8),
    )
    assert result.aenmd.triggered_rules == ("last_exon", "start_proximal")


def test_a_window_spanning_the_whole_reporter_context_is_accepted():
    window = WindowSpec(upstream_nt=REPORTER_UPSTREAM_NT, downstream_nt=min(REPORTER_DOWNSTREAM_NT))
    assert (window.upstream_nt, window.downstream_nt) == (72, 72)


def test_a_window_reaching_past_the_shorter_oligo_design_is_refused():
    with pytest.raises(ValueError, match="shorter oligo design"):
        WindowSpec(upstream_nt=72, downstream_nt=max(REPORTER_DOWNSTREAM_NT))


def test_a_window_reaching_past_the_upstream_context_is_refused():
    with pytest.raises(ValueError, match="upstream"):
        WindowSpec(upstream_nt=REPORTER_UPSTREAM_NT + 1, downstream_nt=0)


def test_ranking_across_modalities_requires_acknowledgement():
    both = [Modality.small_molecule, Modality.suppressor_trna]
    with pytest.raises(IncomparableModalitiesError):
        enforce_modality_comparability(both)
    enforce_modality_comparability(both, incomparable_scales=True)
    enforce_modality_comparability([Modality.small_molecule, Modality.small_molecule])
