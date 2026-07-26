"""The ledger of what this project predicted and what happened, misses included.

A model that is only ever compared against the data it was fitted on has not been tested. The ledger
records each claim the project made, the observation that bears on it, and what followed — including
the rows where nothing followed, because a prediction quietly dropped when it stopped looking good
is the single easiest way to make a project look better than it is.

**Four verdicts, and they are not interchangeable.** The distinction that matters most is between
the three ways of not being supported, which a single "failed" would collapse:

- `supported` — the observation met the criteria fixed before it was made.
- `refuted` — the observation is evidence *against* the prediction, not merely absent.
- `inconclusive` — the comparison ran and did not clear its criteria, on data that could in
  principle have cleared them. The claim is unsupported; the world is not thereby settled.
- `untestable` — the comparison could not have resolved the prediction whatever it showed, so its
  result carries no information about the claim. A row here is a **data gap**, and it carries what
  would close it.

`untestable` is a verdict rather than a gap in the table. A row whose observation cannot be compared
to a prediction records why, and stays.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

__all__ = [
    "LEDGER",
    "VERDICTS",
    "read_ledger",
    "summarise",
]

LEDGER = Path("ledger/falsification.tsv")
"""Where the ledger lives. Authored rather than computed, and committed."""

VERDICTS = ("supported", "refuted", "inconclusive", "untestable")
"""Every verdict a row may carry. A row outside this set is refused rather than counted."""

_REQUIRED = (
    "claim_id",
    "claim",
    "assay",
    "endpoint",
    "comparability",
    "verdict",
    "detail",
    "what_would_settle_it",
    "source",
)
"""Every column a row must carry.

`comparability` is the honest part: whether the observation measures the thing the claim is about at
all. `what_would_settle_it` is required even where the verdict is `supported`, because a supported
claim still has a next experiment and writing it down is what stops one dataset from becoming proof.
"""


def read_ledger(path: Path = LEDGER) -> pd.DataFrame:
    """Load the ledger, refusing one that is malformed rather than reporting from it.

    A ledger with a missing column or an invented verdict is worse than none: it reads as a complete
    record while quietly omitting the thing a reader is looking for.
    """

    rows = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if missing := [column for column in _REQUIRED if column not in rows.columns]:
        raise ValueError(f"{path} is missing {', '.join(missing)}")
    if unknown := sorted(set(rows["verdict"]) - set(VERDICTS)):
        raise ValueError(f"{path} carries verdicts outside {VERDICTS}: {', '.join(unknown)}")
    if blank := [c for c in _REQUIRED if (rows[c] == "").any()]:
        raise ValueError(f"{path} leaves {', '.join(blank)} empty on at least one row")
    if rows["claim_id"].duplicated().any():
        raise ValueError(f"{path} reuses a claim_id")
    return rows


def summarise(ledger: pd.DataFrame) -> pd.DataFrame:
    """How many rows sit under each verdict, in the order the verdicts are declared.

    Every verdict appears, including the ones with no rows, so a reader can see that a category is
    empty rather than having to notice that it is absent.
    """

    counts = ledger["verdict"].value_counts()
    return pd.DataFrame(
        {"verdict": VERDICTS, "rows": [int(counts.get(verdict, 0)) for verdict in VERDICTS]}
    )
