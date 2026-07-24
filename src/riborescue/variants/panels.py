"""Choosing a panel of suppressor tRNAs that covers the most pathogenic nonsense variants.

A suppressor design is the stop codon it decodes and the residue it inserts. Under the coverage
predicate frozen for this first version — exact stop-codon recognition and restoration of the
residue the protein natively carries — a design covers a variant when both match. This is *model
coverage*: the design would put the correct residue back at that stop. It is not a claim the protein
folds or functions, and not a claim any tRNA is safe or clinical.

The objectives are kept separate and explicit — variants covered, genes covered, and the marginal
coverage each added design buys. There is deliberately no safety axis: the native-stop atlas
measures G418 in one cell line and does not transfer to engineered tRNAs.

Under exact restoration each scoreable variant maps to exactly one design — its own stop and its own
residue — so the covered sets are disjoint: they partition the variants. Greedy maximum coverage is
then provably optimal with no submodular gap, and the frontier is just the designs sorted by size.
The greedy engine is written generically anyway, so a later overlapping predicate — a conservative
substitution covers a variant many ways — drops in without a rewrite, and its parity against the
brute-force optimum is what the tests check.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from itertools import combinations

import pandas as pd

__all__ = [
    "PanelStep",
    "brute_force_optimal",
    "coverage_frontier",
    "covered_sets",
    "greedy_panel",
    "is_partition",
]

_ELEMENT = {"variants": "variant_id", "genes": "gene_symbol"}


@dataclass(frozen=True)
class PanelStep:
    """One design added to the panel, with what it buys and what still is not covered.

    `marginal` is the coverage this design adds over those before it; `cumulative` is the running
    total; `uncovered` is what remains. Across ranks they are the whole auditable frontier — every
    panel of size k is the first k rows.
    """

    rank: int
    design_id: str
    recognized_stop: str
    inserted_aa: str
    marginal: int
    cumulative: int
    cumulative_fraction: float
    uncovered: int


def covered_sets(
    contexts: pd.DataFrame, *, element: str = "variant_id"
) -> dict[str, frozenset[str]]:
    """Map each suppressor design to the set of elements it restores exactly.

    `element` is the column whose distinct values are covered — `variant_id` for variants,
    `gene_symbol` for genes. A design's id is the stop it decodes and the residue it inserts, and a
    scoreable variant is restored by the design built from its own stop and original residue.
    """

    scoreable = contexts[contexts["scoreable"]]
    sets: dict[str, set[str]] = {}
    for row in scoreable.itertuples():
        design = f"{str(row.stop_type).upper()}-{row.original_aa}"
        sets.setdefault(design, set()).add(str(getattr(row, element)))
    return {design: frozenset(members) for design, members in sets.items()}


def is_partition(sets: Mapping[str, frozenset[str]]) -> bool:
    """Whether the covered sets are pairwise disjoint — the property that makes greedy optimal."""

    total = sum(len(members) for members in sets.values())
    union = len(frozenset().union(*sets.values())) if sets else 0
    return total == union


def greedy_panel(sets: Mapping[str, frozenset[str]]) -> list[PanelStep]:
    """The maximum-coverage greedy order: at each step, the design that adds the most new coverage.

    Ties break on design id so the panel is deterministic and reproducible. Selection stops when
    no remaining design adds anything, so the frontier never lists a design that buys zero coverage.
    """

    total = len(frozenset().union(*sets.values())) if sets else 0
    covered: set[str] = set()
    available = dict(sets)
    steps: list[PanelStep] = []
    while available:
        # Largest new coverage first; the lexicographic design id is the deterministic tie-break.
        best = min(available, key=lambda design: (-len(available[design] - covered), design))
        gain = len(available[best] - covered)
        if gain == 0:
            break
        covered |= available.pop(best)
        # partition, not split, so a generic set id with no "-" leaves the residue field empty
        # rather than raising — the engine does not depend on the design-id format.
        stop, _, inserted = best.partition("-")
        steps.append(
            PanelStep(
                rank=len(steps) + 1,
                design_id=best,
                recognized_stop=stop,
                inserted_aa=inserted,
                marginal=gain,
                cumulative=len(covered),
                cumulative_fraction=len(covered) / total if total else 0.0,
                uncovered=total - len(covered),
            )
        )
    return steps


def brute_force_optimal(sets: Mapping[str, frozenset[str]], k: int) -> int:
    """The most elements any k designs can cover — the exact optimum greedy is checked against.

    Enumerates every k-subset, so it is only for small fixtures. It proves the greedy engine right
    where the answer can be computed exactly, not to run on the genome.
    """

    designs = list(sets)
    best = 0
    for combo in combinations(designs, min(k, len(designs))):
        reached = len(frozenset().union(*(sets[design] for design in combo))) if combo else 0
        best = max(best, reached)
    return best


def coverage_frontier(contexts: pd.DataFrame, objective: str = "variants") -> pd.DataFrame:
    """The full auditable frontier for one objective, one row per added design.

    `objective` is `variants` or `genes`. The result holds every panel size at once: the panel of
    size k is the first k rows, with its cumulative coverage, the marginal each design added, and
    the count still uncovered at that budget.
    """

    if objective not in _ELEMENT:
        raise ValueError(f"objective must be one of {sorted(_ELEMENT)}, not {objective!r}")
    sets = covered_sets(contexts, element=_ELEMENT[objective])
    steps = greedy_panel(sets)
    return pd.DataFrame([asdict(step) for step in steps])
