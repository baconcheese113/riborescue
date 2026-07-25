"""How large an effect a design could have resolved, when its libraries are too thin to say more.

ADR-0019. A dataset whose libraries fail calibration on depth alone has not measured the wrong
thing, only imprecisely, and an interval describes that completely. What it cannot do is answer the
question it was asked, so the answer worth having from it is a different one: how large an effect
this design *could* have seen, and what a design that saw a smaller one would need.

That question has two answers, because uncertainty here has two components that buy different
experiments. The **counting** component comes from the finite number of P-sites a ratio is built
from; it falls as 1/sqrt(depth) and it is the only component sequencing removes. The
**between-library** component is whatever variation survives across libraries of the same arm once
counting is accounted for — preparation, biology, batch — and no depth touches it. Only replication
does. Reporting which of the two binds is what turns a refused dataset into a specification for the
experiment that would settle it.

A validity failure is not admitted here. An inferred P-site offset outside its window means reads
are assigned to the wrong codon, which displaces an estimate rather than widening it, and an
interval around a displaced estimate is a precise statement about the wrong quantity.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from riborescue.riboseq.calibration import (
    MINIMUM_FRAME0_SHARE,
    OFFSET_FROM_5,
    LibraryVerdict,
)

__all__ = [
    "CONFIDENCE",
    "POWER",
    "QUANTITIES",
    "Detectability",
    "bootstrap_intervals",
    "counting_error",
    "detectability",
    "failure_kind",
]

CONFIDENCE = 0.95
"""Two-sided significance the minimum detectable effect is stated at."""

POWER = 0.80
"""Probability of detecting the reference effect, at that significance."""

QUANTITIES = ("downstream_occupancy", "termination_occupancy", "frame_gap")
"""The quantities ADR-0010 fixes. Nothing about the assay is redefined here."""

_MAXIMUM_REPLICATES = 1000
"""Where the replicate search gives up. Any answer near it is a statement that this is hopeless."""

_MAXIMUM_DEPTH = 1e9
"""Where the depth search gives up, which is far past any sequencing anyone would do."""


def failure_kind(verdict: LibraryVerdict) -> str:
    """Whether a library's calibration failure is one an interval can represent.

    Derived from the manifest's own fields rather than from its failure text, so a reworded message
    cannot silently reclassify a library. `precision` is a depth shortfall and nothing else;
    everything else — an offset outside the window, a coding frame-0 share below the floor — is
    `validity`, which means the library is not measuring what the quantity names.
    """

    if verdict.passes:
        return "none"
    if not OFFSET_FROM_5[0] <= verdict.offset_from_5 <= OFFSET_FROM_5[1]:
        return "validity"
    if verdict.frame0_share < MINIMUM_FRAME0_SHARE:
        return "validity"
    return "precision"


def _pooled(counts: pd.DataFrame) -> pd.DataFrame:
    """Per library, the sums every quantity is built from."""

    columns = [
        "cds_frame0",
        "cds_frame1",
        "cds_frame2",
        "extension_frame0",
        "extension_frame1",
        "extension_frame2",
        "termination",
    ]
    return counts.groupby("sample")[columns].sum()


def _ratio_error(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """Delta-method error on a ratio of two counts, each taken as Poisson.

    A count of zero has no relative error to speak of, so it is floored at one event: the resulting
    figure is the error the smallest observable count would carry, which is the conservative reading
    and keeps the ratio's error finite.
    """

    lower = numerator.clip(lower=1.0)
    safe = denominator.clip(lower=1.0)
    return (lower / safe) * np.sqrt(1.0 / lower + 1.0 / safe)


def _share_error(part: pd.Series, whole: pd.Series) -> pd.Series:
    """Binomial error on a proportion built from `whole` counts."""

    safe = whole.clip(lower=1.0)
    share = part / safe
    return (share * (1.0 - share) / safe) ** 0.5


def counting_error(counts: pd.DataFrame) -> pd.DataFrame:
    """Per library and quantity, the standard error the finite P-site count alone imposes.

    This is the component depth buys down. It is smaller than the transcript bootstrap below,
    because it holds the transcript universe fixed and asks only how well its counts are resolved.
    """

    pooled = _pooled(counts)
    extension = pooled[["extension_frame0", "extension_frame1", "extension_frame2"]].sum(axis=1)
    coding = pooled[["cds_frame0", "cds_frame1", "cds_frame2"]].sum(axis=1)
    gap = np.sqrt(
        _share_error(pooled["extension_frame0"], extension) ** 2
        + _share_error(pooled["cds_frame0"], coding) ** 2
    )
    return pd.DataFrame(
        {
            "downstream_occupancy": _ratio_error(pooled["extension_frame0"], pooled["cds_frame0"]),
            "termination_occupancy": _ratio_error(pooled["termination"], pooled["cds_frame0"]),
            "frame_gap": gap,
        }
    ).rename_axis("sample")


def bootstrap_intervals(counts: pd.DataFrame, draws: int = 2000, seed: int = 0) -> pd.DataFrame:
    """Per library and quantity, a percentile interval from resampling qualifying transcripts.

    The honest within-library precision, and what a per-library value is quoted with. It is wider
    than the counting error because it also carries the heterogeneity between transcripts, which is
    real uncertainty about the library's value but is not something more reads remove.

    Transcripts are resampled by multinomial weight rather than by index, which is the same
    procedure written as one matrix product per library.
    """

    generator = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for sample, block in counts.groupby("sample", sort=True):
        size = len(block)
        weights = generator.multinomial(size, np.full(size, 1.0 / size), size=draws).astype(float)
        summed = {
            column: weights @ block[column].to_numpy(dtype=float)
            for column in (
                "cds_frame0",
                "cds_frame1",
                "cds_frame2",
                "extension_frame0",
                "extension_frame1",
                "extension_frame2",
                "termination",
            )
        }
        extension = sum(summed[f"extension_frame{frame}"] for frame in (0, 1, 2))
        coding = sum(summed[f"cds_frame{frame}"] for frame in (0, 1, 2))
        with np.errstate(divide="ignore", invalid="ignore"):
            replicated = {
                "downstream_occupancy": summed["extension_frame0"] / summed["cds_frame0"],
                "termination_occupancy": summed["termination"] / summed["cds_frame0"],
                "frame_gap": summed["extension_frame0"] / extension - summed["cds_frame0"] / coding,
            }
        for quantity, values in replicated.items():
            finite = values[np.isfinite(values)]
            low, high = np.percentile(finite, [2.5, 97.5]) if finite.size else (np.nan, np.nan)
            rows.append(
                {
                    "sample": sample,
                    "quantity": quantity,
                    "bootstrap_low": float(low),
                    "bootstrap_high": float(high),
                    "bootstrap_se": float(finite.std(ddof=1)) if finite.size > 1 else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def _detectable(error: float, n_treated: int, n_control: int) -> float:
    """The smallest effect resolvable at CONFIDENCE and POWER, for arms of these sizes.

    The reference distribution takes the design's degrees of freedom, `n_treated + n_control - 2`,
    rather than the Welch approximation ADR-0010 uses to report an observed contrast. The two answer
    different questions and only one of them is stable here: Welch's degrees of freedom fall as the
    two arms' variances diverge, so scaling a hypothetical depth would move the reference
    distribution as well as the error, and a design could come out *less* able to resolve an effect
    at unlimited depth than at the depth it has. How many libraries were sequenced is a property of
    the design and does not change with how deeply they were read.
    """

    if error == 0:
        return 0.0
    df = float(n_treated + n_control - 2)
    significance = float(stats.t.ppf(1.0 - (1.0 - CONFIDENCE) / 2.0, df))
    detectable = float(stats.t.ppf(POWER, df))
    return (significance + detectable) * error


@dataclass(frozen=True)
class Detectability:
    """One quantity, asked how large an effect the design behind it could have resolved.

    `counting` holds each library's counting standard error, in the same order as its arm's values.
    `target` is the reference effect from a dataset that passed and completed — the size the arm is
    asked whether it could have seen — and its sign is irrelevant, only its magnitude.
    """

    quantity: str
    treated: tuple[float, ...]
    control: tuple[float, ...]
    treated_counting: tuple[float, ...]
    control_counting: tuple[float, ...]
    target: float

    @property
    def mean_difference(self) -> float:
        return float(np.mean(self.treated) - np.mean(self.control))

    @property
    def degrees_of_freedom(self) -> int:
        """How many observations the between-library variance rests on, in total across both arms.

        Two libraries per arm gives one degree of freedom each. Every figure below that involves the
        between-library component inherits that, and is a rough figure for exactly this reason.
        """

        return len(self.treated) - 1 + len(self.control) - 1

    def _variances(self, depth: float) -> tuple[float, float]:
        """Per-arm variance of a library's value, with the counting component scaled by `depth`.

        The observed spread across an arm's libraries already contains its counting error, so the
        part that survives at unlimited depth is the residual. A negative residual — ordinary when
        the variance rests on one degree of freedom — is reported as zero, which is the assumption
        that this arm's libraries differ only by counting.

        What is subtracted is the counting error and not the wider transcript bootstrap, and the
        difference matters. Every library of a contrast is measured over the same qualifying
        transcripts, so how those transcripts differ from each other is common to all of them and
        contributes nothing to how the libraries differ. Subtracting the bootstrap would remove a
        term the observed spread never contained, and would credit depth with closing a gap that is
        between libraries.
        """

        variances = []
        for values, counting in (
            (self.treated, self.treated_counting),
            (self.control, self.control_counting),
        ):
            observed = float(np.var(values, ddof=1)) if len(values) > 1 else 0.0
            counted = float(np.mean(np.square(counting)))
            between = max(0.0, observed - counted)
            variances.append(between + counted / depth)
        return variances[0], variances[1]

    def _effect_at(self, depth: float, replicates: int | None = None) -> float:
        n_treated = replicates if replicates is not None else len(self.treated)
        n_control = replicates if replicates is not None else len(self.control)
        variance_treated, variance_control = self._variances(depth)
        error = math.sqrt(variance_treated / n_treated + variance_control / n_control)
        return _detectable(error, n_treated, n_control)

    @property
    def between_library_share(self) -> float:
        """The fraction of the pooled variance that no amount of sequencing removes."""

        total_treated, total_control = self._variances(1.0)
        floor_treated, floor_control = self._variances(math.inf)
        total = total_treated + total_control
        return (floor_treated + floor_control) / total if total else float("nan")

    @property
    def minimum_detectable_effect(self) -> float:
        """The smallest effect this design, at the depth and arm sizes it has, could resolve."""

        return self._effect_at(1.0)

    @property
    def irreducible_effect(self) -> float:
        """What stays detectable at unlimited depth and these arm sizes — the replication floor."""

        return self._effect_at(math.inf)

    @property
    def depth_multiplier(self) -> float:
        """How much deeper these libraries would have to be sequenced to resolve the target.

        A lower bound, and labelled one wherever it is reported: the counting component is scaled as
        1/sqrt(depth), which is right for it, while the transcript heterogeneity inside a library's
        real precision does not shrink with reads at all. Infinite where the between-library
        component alone already exceeds the target, which means depth is not the binding constraint
        and no sequencing run resolves this.
        """

        target = abs(self.target)
        if self.irreducible_effect >= target:
            return math.inf
        if self.minimum_detectable_effect <= target:
            return 1.0
        low, high = 1.0, _MAXIMUM_DEPTH
        for _ in range(200):
            middle = math.sqrt(low * high)
            if self._effect_at(middle) > target:
                low = middle
            else:
                high = middle
        return high

    @property
    def replicates_required(self) -> int | None:
        """Libraries per arm needed to resolve the target at the depth these libraries have.

        None where even `_MAXIMUM_REPLICATES` per arm does not reach it, which at these variances
        would mean the quantity is not measurable by replication either.
        """

        target = abs(self.target)
        for replicates in range(2, _MAXIMUM_REPLICATES + 1):
            if self._effect_at(1.0, replicates) <= target:
                return replicates
        return None


def detectability(
    ratios: pd.DataFrame,
    counting: pd.DataFrame,
    arms: pd.DataFrame,
    treated: str,
    control: str,
    targets: Mapping[str, float],
) -> dict[str, Detectability]:
    """Every quantity's detectability for one contrast, from per-library values and counting errors.

    `ratios` is what `library_ratios` produced, `counting` what `counting_error` produced, and
    `arms` names each library's treatment. `targets` is the reference effect per quantity, measured
    on a dataset that passed its calibration and completed its analysis.
    """

    merged = ratios.merge(arms, on="sample", validate="one_to_one").merge(
        counting.reset_index(), on="sample", suffixes=("", "_counting"), validate="one_to_one"
    )
    found: dict[str, Detectability] = {}
    for quantity in QUANTITIES:
        if quantity not in targets:
            raise ValueError(f"no reference effect for {quantity!r}")
        groups = {}
        for arm in (treated, control):
            block = merged.loc[merged["treatment"] == arm]
            if block.empty:
                raise ValueError(f"no libraries for {arm!r}")
            groups[arm] = (
                tuple(block[quantity].astype(float)),
                tuple(block[f"{quantity}_counting"].astype(float)),
            )
        found[quantity] = Detectability(
            quantity=quantity,
            treated=groups[treated][0],
            control=groups[control][0],
            treated_counting=groups[treated][1],
            control_counting=groups[control][1],
            target=float(targets[quantity]),
        )
    return found
