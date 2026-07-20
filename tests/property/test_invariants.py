import pytest
from hypothesis import given
from hypothesis import strategies as st

from riborescue.contracts import (
    Consequence,
    FieldStatus,
    Measurement,
    MissingReason,
    UnresolvedReporterWindowError,
    WindowSpec,
)
from riborescue.triage import classify


@given(st.sampled_from(Consequence), st.booleans())
def test_classify_is_total(consequence, transcript_supported):
    result = classify(consequence, transcript_supported=transcript_supported)
    assert isinstance(result.reason, str)
    assert result.reason
    assert result.applies in (True, False)


@given(st.floats(allow_nan=False, allow_infinity=False))
def test_present_measurement_round_trips_any_finite_value(value):
    measurement = Measurement.present(value)
    assert measurement.value == value
    assert measurement.status is FieldStatus.present
    assert measurement.reason is None


@given(st.sampled_from(MissingReason))
def test_absent_measurement_never_holds_a_value(reason):
    measurement = Measurement.absent(reason)
    assert measurement.value is None
    assert measurement.status is FieldStatus.missing


@given(
    st.integers(min_value=0, max_value=200),
    st.integers(min_value=0, max_value=200),
)
def test_window_spec_refuses_any_arguments_while_unresolved(upstream, downstream):
    with pytest.raises(UnresolvedReporterWindowError):
        WindowSpec(upstream_nt=upstream, downstream_nt=downstream)
