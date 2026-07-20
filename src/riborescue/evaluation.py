"""The evaluation seam: split definitions and the shuffle controls that decide the kinetics claim.

The controls are the most important tests in the project and the ones most likely to be quietly
dropped. Their signatures live here so the protected control suite can assert against them before
the kinetics comparison itself exists.
"""

from dataclasses import dataclass
from enum import StrEnum

__all__ = ["BootstrapCI", "ShuffleKind", "grouped_split_leakage", "run_shuffle_control"]


class ShuffleKind(StrEnum):
    global_ = "global"
    within_gene = "within_gene"
    context_matched = "context_matched"


@dataclass(frozen=True)
class BootstrapCI:
    low: float
    point: float
    high: float

    def includes_zero(self) -> bool:
        return self.low <= 0.0 <= self.high


def run_shuffle_control(kind: ShuffleKind) -> BootstrapCI:
    """Return the bootstrap CI of the kinetics improvement after shuffling.

    A genuine kinetic signal collapses under every shuffle, so a passing control has a CI that
    includes zero.
    """

    raise NotImplementedError("shuffle controls are implemented with the head-to-head comparison")


def grouped_split_leakage(config: str) -> BootstrapCI:
    """Return the improvement retained under a gene- or cluster-grouped split.

    Near-identical local contexts leaking across a random split would inflate this; grouping removes
    the leak.
    """

    raise NotImplementedError("grouped-split evaluation lands with the head-to-head comparison")
