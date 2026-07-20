import pytest

from riborescue.contracts import Consequence, TriageClass
from riborescue.triage import classify


@pytest.mark.parametrize(
    ("consequence", "expected", "applies"),
    [
        (Consequence.stop_gained, TriageClass.supported_nonsense, True),
        (Consequence.missense_variant, TriageClass.missense, False),
        (Consequence.frameshift_variant, TriageClass.frameshift, False),
        (Consequence.splice_donor_variant, TriageClass.splice, False),
        (Consequence.splice_acceptor_variant, TriageClass.splice, False),
        (Consequence.synonymous_variant, TriageClass.synonymous, False),
        (Consequence.stop_lost, TriageClass.stop_loss, False),
    ],
)
def test_classify_maps_consequence_to_verdict(consequence, expected, applies):
    result = classify(consequence)
    assert result.triage_class is expected
    assert result.applies is applies


def test_unsupported_transcript_short_circuits():
    result = classify(Consequence.stop_gained, transcript_supported=False)
    assert result.triage_class is TriageClass.unsupported_transcript
    assert result.applies is False


def test_only_nonsense_variants_are_eligible():
    for consequence in Consequence:
        result = classify(consequence)
        assert result.applies == (result.triage_class is TriageClass.supported_nonsense)
