"""The controls that decide the kinetics claim, and the guard that keeps them honest.

Every negative control is registered in CONTROLS and named after its test. Each is `xfail(strict)`
until the head-to-head comparison implements it: the test encodes the real invariant now, cannot
pass by accident, and cannot silently pass without forcing its marker to be removed. When a control
graduates, its stage joins WAVES_COMPLETE and the guard below fails until the marker is gone.
"""

import sys

import pytest

from riborescue.evaluation import ShuffleKind, grouped_split_leakage, run_shuffle_control

CONTROLS: dict[str, str] = {
    "shuffle_global": "head_to_head",
    "shuffle_within_gene": "head_to_head",
    "shuffle_context_matched": "head_to_head",
    "grouped_split_leakage": "head_to_head",
}
WAVES_COMPLETE: frozenset[str] = frozenset()

_PENDING = "shuffle controls land with the head-to-head comparison"


@pytest.mark.control
@pytest.mark.xfail(strict=True, reason=_PENDING)
def test_shuffle_global():
    assert run_shuffle_control(ShuffleKind.global_).includes_zero()


@pytest.mark.control
@pytest.mark.xfail(strict=True, reason=_PENDING)
def test_shuffle_within_gene():
    assert run_shuffle_control(ShuffleKind.within_gene).includes_zero()


@pytest.mark.control
@pytest.mark.xfail(strict=True, reason=_PENDING)
def test_shuffle_context_matched():
    assert run_shuffle_control(ShuffleKind.context_matched).includes_zero()


@pytest.mark.control
@pytest.mark.xfail(strict=True, reason=_PENDING)
def test_grouped_split_leakage():
    assert grouped_split_leakage("grouped_by_gene").includes_zero()


def test_no_control_is_silently_removed():
    module = sys.modules[__name__]
    for name in CONTROLS:
        assert hasattr(module, f"test_{name}"), f"control '{name}' is registered but has no test"


def test_graduated_controls_drop_their_xfail():
    module = sys.modules[__name__]
    for name, stage in CONTROLS.items():
        if stage in WAVES_COMPLETE:
            marks = {m.name for m in getattr(getattr(module, f"test_{name}"), "pytestmark", [])}
            assert "xfail" not in marks, f"'{name}' has graduated; remove its xfail marker"
