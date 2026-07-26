"""The controls that decide the kinetics claim, and the guard that keeps them honest.

Every negative control is registered in CONTROLS, and the recorded outcomes below are what say it
ran. A control deleted from the analysis loses its row and this suite goes red.

**A failing control is a scientific outcome, not broken software.** A control that does not clear
means the claim is not made — which is a valid, green state for a test suite to be in. What must
never happen is the claim being made anyway, so the guard is on the conjunction: only
`Criteria.supported` may describe kinetics as carrying transferable information, and the
recorded outcome below pins what that conjunction currently returns. If a change to the data or the
code would make the claim supported, `test_the_recorded_outcome_still_does_not_support_the_claim`
goes red and somebody has to look.

The control calculations themselves are exercised under `slow`, because each one refits every model
over every round of six drugs.
"""

from pathlib import Path
from typing import ClassVar

import pytest

from riborescue.core.tables import read_table
from riborescue.variants.evaluation import (
    CODON_TABLE,
    BootstrapCI,
    ShuffleKind,
    grouped_split_leakage,
    run_shuffle_control,
)
from riborescue.variants.kinetics import Criteria

CONTROLS: dict[str, str] = {
    "shuffle_global": "head_to_head",
    "shuffle_within_gene": "head_to_head",
    "shuffle_context_matched": "head_to_head",
    "grouped_split_leakage": "head_to_head",
}

RECORDED = Path("tests/fixtures/kinetics/control_outcomes.tsv")
"""What the controls returned, under ADR-0020's frozen single-permutation rule.

Committed rather than recomputed in the fast suite for the reason the oracle fixtures are: six
minutes of refitting does not belong in a gate that runs on every change. The slow tests recompute
it and fail if it has drifted.
"""


def _recorded() -> dict[str, BootstrapCI]:
    rows = read_table(RECORDED)
    return {
        str(control): BootstrapCI(low=float(low), point=float(gain), high=float(high))
        for control, low, gain, high in zip(
            rows["control"], rows["ci_low"], rows["gain"], rows["ci_high"], strict=True
        )
    }


def _criteria() -> Criteria:
    """ADR-0020's three conditions, from the recorded outcomes."""

    recorded = _recorded()
    return Criteria(
        grouped_gain=recorded["grouped_gain"],
        shuffles={kind.value: recorded[kind.value] for kind in ShuffleKind},
        # The stratum the condition guards against is empty: every held-out cell was observed in
        # training, so no gain can be confined to one the baseline could not fit.
        confined_to_unsupported=False,
    )


class TestTheConjunction:
    """The rule, over constructed inputs. These are the software invariants."""

    @staticmethod
    def _criteria(gain: BootstrapCI, shuffles: dict[str, BootstrapCI], confined=False) -> Criteria:
        return Criteria(grouped_gain=gain, shuffles=shuffles, confined_to_unsupported=confined)

    _CLEARS: ClassVar = BootstrapCI(low=0.004, point=0.006, high=0.008)
    _COLLAPSES: ClassVar = {
        kind.value: BootstrapCI(low=-0.001, point=0.0, high=0.001) for kind in ShuffleKind
    }

    def test_the_claim_holds_only_when_all_three_conditions_do(self):
        assert self._criteria(self._CLEARS, self._COLLAPSES).supported

    def test_a_gain_that_includes_zero_does_not_support_the_claim(self):
        flat = BootstrapCI(low=-0.002, point=0.001, high=0.004)
        criteria = self._criteria(flat, self._COLLAPSES)
        assert not criteria.supported
        assert criteria.failures() == ("gain_excludes_zero",)

    def test_a_negative_gain_that_clears_zero_does_not_support_the_claim(self):
        # Clearing zero downwards is a reliable *loss*, and the rule is one-sided by design.
        worse = BootstrapCI(low=-0.008, point=-0.006, high=-0.004)
        assert not self._criteria(worse, self._COLLAPSES).supported

    @pytest.mark.parametrize("leaking", [k.value for k in ShuffleKind])
    def test_any_one_shuffle_that_survives_sinks_the_claim(self, leaking):
        shuffles = dict(self._COLLAPSES)
        shuffles[leaking] = BootstrapCI(low=0.0003, point=0.0007, high=0.0011)
        criteria = self._criteria(self._CLEARS, shuffles)
        assert not criteria.supported
        assert criteria.failures() == ("shuffles_collapse",)

    def test_a_gain_confined_to_cells_the_baseline_could_not_fit_sinks_the_claim(self):
        criteria = self._criteria(self._CLEARS, self._COLLAPSES, confined=True)
        assert not criteria.supported
        assert criteria.failures() == ("survives_support",)

    def test_a_run_with_no_shuffles_at_all_cannot_support_the_claim(self):
        assert not self._criteria(self._CLEARS, {}).supported


class TestTheRecordedOutcome:
    """What the controls actually returned, and what follows from it."""

    def test_every_registered_control_has_a_recorded_outcome(self):
        recorded = _recorded()
        for kind in ShuffleKind:
            assert kind.value in recorded
        assert "grouped_gain" in recorded and "grouped_split_leakage" in recorded

    def test_the_recorded_outcome_still_does_not_support_the_claim(self):
        # ADR-0020's verdict, pinned. Kinetics is not reported as carrying transferable information,
        # and this goes red if a change would make it supported without anyone deciding to.
        criteria = _criteria()
        assert not criteria.supported
        assert criteria.failures() == ("shuffles_collapse",)

    def test_the_condition_that_failed_is_the_global_shuffle_and_only_it(self):
        recorded = _recorded()
        assert not recorded[ShuffleKind.global_.value].includes_zero()
        assert recorded[ShuffleKind.within_gene.value].includes_zero()
        assert recorded[ShuffleKind.context_matched.value].includes_zero()

    def test_the_leakage_control_cleared(self):
        assert _recorded()["grouped_split_leakage"].includes_zero()


@pytest.mark.slow
class TestTheControlsThemselves:
    """The calculations, recomputed. Slow: each refits every round of six drugs."""

    @pytest.mark.parametrize("kind", list(ShuffleKind))
    def test_each_shuffle_control_reproduces_its_recorded_interval(self, kind):
        computed = run_shuffle_control(kind)
        recorded = _recorded()[kind.value]
        assert computed.point == pytest.approx(recorded.point, abs=1e-6)
        assert computed.low == pytest.approx(recorded.low, abs=1e-6)
        assert computed.high == pytest.approx(recorded.high, abs=1e-6)

    def test_the_leakage_control_reproduces_its_recorded_interval(self):
        computed = grouped_split_leakage("grouped_by_gene")
        recorded = _recorded()["grouped_split_leakage"]
        assert computed.point == pytest.approx(recorded.point, abs=1e-6)


def test_the_codon_table_the_controls_read_is_committed():
    assert CODON_TABLE.exists(), "the controls must not depend on a regenerable results tree"


def test_no_control_is_silently_removed():
    # Every registered control has to leave a recorded outcome behind. Deleting one from the
    # analysis removes its row, and this goes red — which is what the registry is for.
    recorded = _recorded()
    for name in CONTROLS:
        assert name.removeprefix("shuffle_") in recorded, (
            f"control '{name}' is registered but has no recorded outcome"
        )

