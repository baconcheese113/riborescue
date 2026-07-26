#!/usr/bin/env python
"""Record what the controls returned, under the rule ADR-0020 froze.

Committed as a fixture for the reason the oracle fixtures are: the fast gate cannot spend six
minutes refitting six drugs on every change, and a control that reads a regenerable results tree
passes or fails according to what was last run there. The slow suite recomputes these and fails if
they have drifted, so the fixture is a cache rather than a claim.

The grouped gain is taken on the drug whose gain is largest, matching the worst-case logic the
shuffle controls already use: the claim is asked to survive on its best case, and the controls are
asked to clear on theirs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from riborescue.core.contracts import EvalConfig
from riborescue.core.tables import write_table
from riborescue.variants.evaluation import (
    ShuffleKind,
    _gains,
    bootstrap_ci,
    grouped_split_leakage,
    run_shuffle_control,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out", type=Path, default=Path("tests/fixtures/kinetics/control_outcomes.tsv")
    )
    arguments = parser.parse_args()

    rows: list[dict[str, object]] = []

    unshuffled = _gains(EvalConfig.grouped_by_gene, None)
    best = max(unshuffled, key=lambda drug: unshuffled[drug].mean())
    interval = bootstrap_ci(unshuffled[best])
    rows.append(
        {
            "control": "grouped_gain",
            "drug": best,
            "gain": interval.point,
            "ci_low": interval.low,
            "ci_high": interval.high,
            "includes_zero": interval.includes_zero(),
        }
    )

    for kind in ShuffleKind:
        interval = run_shuffle_control(kind)
        rows.append(
            {
                "control": kind.value,
                "drug": "worst case",
                "gain": interval.point,
                "ci_low": interval.low,
                "ci_high": interval.high,
                "includes_zero": interval.includes_zero(),
            }
        )

    interval = grouped_split_leakage("grouped_by_gene")
    rows.append(
        {
            "control": "grouped_split_leakage",
            "drug": "worst case",
            "gain": interval.point,
            "ci_low": interval.low,
            "ci_high": interval.high,
            "includes_zero": interval.includes_zero(),
        }
    )

    table = pd.DataFrame(rows)
    write_table(table, arguments.out)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
