"""What a design could have resolved, over counts small enough to work through by hand.

The arithmetic that matters here is the split between the error a deeper run removes and the error
it does not, so the fixtures are built to make one dominate and then the other.
"""

import math

import numpy as np
import pandas as pd
import pytest

from riborescue.riboseq.calibration import LibraryVerdict
from riborescue.riboseq.detectability import (
    Detectability,
    bootstrap_intervals,
    counting_error,
    detectability,
    failure_kind,
)


def _library(**overrides) -> LibraryVerdict:
    fields = {
        "sample": "lib",
        "psites": 2_000_000,
        "frame0_share": 0.65,
        "dominant_length": 30,
        "offset_from_5": 12,
        "failures": (),
    }
    return LibraryVerdict(**(fields | overrides))


def _counts(sample: str, transcripts: int, extension_frame0: int, seed: int) -> pd.DataFrame:
    """One library's qualifying transcripts, with the downstream signal unevenly spread.

    Uneven on purpose: real libraries carry most of their downstream signal on a few transcripts,
    and a fixture that spread it evenly would be exactly Poisson, which is the one case where the
    transcript bootstrap and the counting error agree.
    """

    generator = np.random.default_rng(seed)
    weights = generator.dirichlet(np.full(transcripts, 0.3))
    return pd.DataFrame(
        {
            "transcript": [f"ENST{i:05d}" for i in range(transcripts)],
            "sample": sample,
            "cds_frame0": 40_000,
            "cds_frame1": 15_000,
            "cds_frame2": 15_000,
            "extension_frame0": generator.multinomial(extension_frame0, weights),
            "extension_frame1": 1,
            "extension_frame2": 1,
            "termination": 2_500,
        }
    )


class TestFailureKind:
    def test_a_passing_library_has_no_failure_to_classify(self):
        assert failure_kind(_library()) == "none"

    def test_a_depth_shortfall_alone_is_precision(self):
        thin = _library(psites=659_427, failures=("659,427 P-sites over the selected lengths",))
        assert failure_kind(thin) == "precision"

    @pytest.mark.parametrize("offset", [9, 10, 15, 20])
    def test_an_offset_outside_the_window_is_validity(self, offset):
        displaced = _library(offset_from_5=offset, failures=("5' offset outside 11-14",))
        assert failure_kind(displaced) == "validity"

    def test_a_library_failing_on_both_is_validity(self):
        both = _library(psites=500, offset_from_5=9, failures=("thin", "displaced"))
        assert failure_kind(both) == "validity"

    def test_a_frame0_share_below_the_floor_is_validity(self):
        flat = _library(frame0_share=0.31, failures=("frame-0 share 31.0%, under 40%",))
        assert failure_kind(flat) == "validity"


class TestCountingError:
    def test_ten_times_the_downstream_counts_halves_the_relative_error(self):
        thin = _counts("thin", transcripts=100, extension_frame0=100, seed=1)
        deep = _counts("deep", transcripts=100, extension_frame0=10_000, seed=1)
        errors = counting_error(pd.concat([thin, deep]))
        occupancy = errors["downstream_occupancy"]
        # The Poisson error on a ratio scales as 1/sqrt(count), so a hundredfold count is a tenth
        # the relative error. Both denominators are four million, large enough that the error is
        # the numerator's alone.
        denominator = 100 * 40_000
        thin_relative = occupancy["thin"] / (100 / denominator)
        deep_relative = occupancy["deep"] / (10_000 / denominator)
        assert deep_relative == pytest.approx(thin_relative / 10, rel=0.01)

    def test_a_library_with_no_downstream_signal_still_has_a_finite_error(self):
        empty = _counts("empty", transcripts=50, extension_frame0=0, seed=2)
        errors = counting_error(empty)
        assert math.isfinite(errors["downstream_occupancy"]["empty"])
        assert errors["downstream_occupancy"]["empty"] > 0


class TestBootstrap:
    def test_the_interval_brackets_the_value_the_whole_library_gives(self):
        counts = _counts("lib", transcripts=200, extension_frame0=800, seed=3)
        observed = counts["extension_frame0"].sum() / counts["cds_frame0"].sum()
        interval = bootstrap_intervals(counts, draws=400, seed=0)
        row = interval.set_index("quantity").loc["downstream_occupancy"]
        assert row["bootstrap_low"] < observed < row["bootstrap_high"]

    def test_it_is_wider_than_the_counting_error_because_transcripts_differ(self):
        counts = _counts("lib", transcripts=200, extension_frame0=800, seed=4)
        bootstrap = bootstrap_intervals(counts, draws=800, seed=0)
        row = bootstrap.set_index("quantity").loc["downstream_occupancy"]
        counted = counting_error(counts)["downstream_occupancy"]["lib"]
        assert row["bootstrap_se"] > counted

    def test_the_same_seed_gives_the_same_interval(self):
        counts = _counts("lib", transcripts=80, extension_frame0=300, seed=5)
        first = bootstrap_intervals(counts, draws=200, seed=7)
        second = bootstrap_intervals(counts, draws=200, seed=7)
        pd.testing.assert_frame_equal(first, second)


def _sized(treated, control, treated_counting, control_counting, target) -> Detectability:
    return Detectability(
        quantity="downstream_occupancy",
        treated=treated,
        control=control,
        treated_counting=treated_counting,
        control_counting=control_counting,
        target=target,
    )


class TestDetectability:
    def test_the_reference_effect_sign_does_not_matter_only_its_size(self):
        rising = _sized((0.02, 0.03), (0.01, 0.011), (0.001, 0.001), (0.001, 0.001), 0.01)
        falling = _sized((0.02, 0.03), (0.01, 0.011), (0.001, 0.001), (0.001, 0.001), -0.01)
        assert rising.depth_multiplier == falling.depth_multiplier
        assert rising.replicates_required == falling.replicates_required

    def test_depth_is_unreachable_when_the_libraries_disagree_beyond_their_counting_error(self):
        # The two treated libraries sit a long way apart on an error of 0.0001, so the spread is
        # between libraries rather than in the counts, and sequencing cannot close it.
        spread = _sized((0.05, 0.20), (0.01, 0.02), (1e-4, 1e-4), (1e-4, 1e-4), 0.01)
        assert spread.between_library_share > 0.99
        assert math.isinf(spread.depth_multiplier)
        assert spread.replicates_required is not None

    def test_depth_is_the_lever_when_the_spread_is_all_counting_error(self):
        # Libraries agreeing to within their own counting error leave no residual, so the whole
        # variance falls as depth rises and a deep enough run resolves the target.
        agreeing = _sized((0.0300, 0.0302), (0.0100, 0.0102), (0.02, 0.02), (0.02, 0.02), 0.005)
        assert agreeing.between_library_share == pytest.approx(0.0)
        assert agreeing.irreducible_effect == pytest.approx(0.0)
        assert math.isfinite(agreeing.depth_multiplier)
        assert agreeing.depth_multiplier > 1

    def test_a_design_that_already_resolves_the_target_needs_no_more_of_either(self):
        precise = _sized((0.030, 0.0301), (0.010, 0.0101), (1e-5, 1e-5), (1e-5, 1e-5), 0.01)
        assert precise.minimum_detectable_effect < 0.01
        assert precise.depth_multiplier == 1.0
        assert precise.replicates_required == 2

    def test_the_floor_at_unlimited_depth_is_never_above_what_this_depth_resolves(self):
        effect = _sized((0.05, 0.09), (0.01, 0.02), (0.01, 0.01), (0.01, 0.01), 0.01)
        assert effect.irreducible_effect <= effect.minimum_detectable_effect

    def test_two_libraries_an_arm_gives_the_variance_two_degrees_of_freedom(self):
        effect = _sized((0.05, 0.09), (0.01, 0.02), (0.01, 0.01), (0.01, 0.01), 0.01)
        assert effect.degrees_of_freedom == 2

    def test_more_replicates_resolve_a_smaller_effect_than_the_design_has(self):
        effect = _sized((0.05, 0.09), (0.01, 0.02), (0.01, 0.01), (0.01, 0.01), 0.005)
        assert effect.replicates_required is not None
        assert effect.replicates_required > 2


class TestAssembly:
    def test_every_quantity_needs_a_reference_effect_before_it_is_sized(self):
        ratios = pd.DataFrame(
            {
                "sample": ["a", "b"],
                "downstream_occupancy": [0.02, 0.01],
                "termination_occupancy": [0.05, 0.06],
                "frame_gap": [-0.1, -0.2],
            }
        )
        counting = ratios.set_index("sample")
        arms = pd.DataFrame({"sample": ["a", "b"], "treatment": ["drug", "none"]})
        with pytest.raises(ValueError, match="no reference effect for 'frame_gap'"):
            detectability(
                ratios,
                counting,
                arms,
                "drug",
                "none",
                {"downstream_occupancy": 0.01, "termination_occupancy": -0.01},
            )

    def test_an_arm_with_no_libraries_is_refused_rather_than_sized_on_one_side(self):
        ratios = pd.DataFrame(
            {
                "sample": ["a"],
                "downstream_occupancy": [0.02],
                "termination_occupancy": [0.05],
                "frame_gap": [-0.1],
            }
        )
        arms = pd.DataFrame({"sample": ["a"], "treatment": ["drug"]})
        targets: dict[str, float] = dict.fromkeys(
            ("downstream_occupancy", "termination_occupancy", "frame_gap"), 0.01
        )
        with pytest.raises(ValueError, match="no libraries for 'none'"):
            detectability(ratios, ratios.set_index("sample"), arms, "drug", "none", targets)
