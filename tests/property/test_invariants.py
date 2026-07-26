from hypothesis import given
from hypothesis import strategies as st

from riborescue.core.contracts import Consequence
from riborescue.variants.triage import classify


@given(st.sampled_from(Consequence), st.booleans())
def test_classify_is_total(consequence, transcript_supported):
    result = classify(consequence, transcript_supported=transcript_supported)
    assert isinstance(result.reason, str)
    assert result.reason
    assert result.applies in (True, False)
