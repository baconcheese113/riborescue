"""Placing variants on the hand-worked transcripts from tests/fixtures/mane.

Both carry AACCATGAAAGGGCCCTAAT with coding starting at offset 4, so codon 2 is AAA at offsets 7-9.
Substituting its first base gives TAA, a premature stop; the surrounding bases are known, so what a
correct placement must return is known too.
"""

from pathlib import Path

import pytest

from riborescue.context import ContextFailure, context_for
from riborescue.transcripts import load_transcripts

MANE = Path(__file__).parents[1] / "fixtures" / "mane"


@pytest.fixture
def wide_window(monkeypatch: pytest.MonkeyPatch):
    """Narrow the window to what a twenty-base transcript can supply, leaving the logic intact."""

    monkeypatch.setattr("riborescue.context._WINDOW", 3)


@pytest.fixture
def models():
    return load_transcripts(MANE / "sample.gff", MANE / "sample.fna")


@pytest.fixture
def plus(models):
    return models["NM_000001.1"]


@pytest.fixture
def minus(models):
    return models["NM_000002.1"]


def test_a_substitution_making_a_stop_is_read_as_a_premature_stop(plus, wide_window):
    # Offset 7 is genomic 108; A>T turns AAA into TAA.
    context = context_for(plus, 108, "A", "T")
    assert not isinstance(context, ContextFailure)
    assert context.stop_type == "uaa"
    assert context.reference_codon == "AAA"
    assert context.protein_position == 2
    assert context.up_123nt == "aug"
    assert context.down_123nt == "ggg"


def test_the_same_variant_on_the_reverse_strand_reads_the_complementary_bases(minus, wide_window):
    # The transcript runs from 410 down, so offset 7 is genomic 403; forward-strand, that is T>A.
    context = context_for(minus, 403, "T", "A")
    assert not isinstance(context, ContextFailure)
    assert context.stop_type == "uaa"
    assert context.protein_position == 2


def test_a_reference_base_that_does_not_match_the_transcript_is_refused(plus):
    assert context_for(plus, 108, "G", "T") is ContextFailure.reference_mismatch


def test_a_substitution_that_makes_no_stop_is_refused(plus, wide_window):
    assert context_for(plus, 108, "A", "C") is ContextFailure.not_a_premature_stop


def test_a_position_outside_the_transcript_is_refused(plus):
    assert context_for(plus, 150, "A", "T") is ContextFailure.outside_transcript


def test_a_stop_without_the_full_window_either_side_is_refused(plus):
    """The synthetic transcript is twenty bases, far short of the reporter context."""

    assert context_for(plus, 108, "A", "T") is ContextFailure.truncated_context
