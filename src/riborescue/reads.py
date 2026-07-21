"""What trimming and alignment did to each library.

Cutadapt and STAR each report their own counts. Summarising them gives the raw-against-cleaned
comparison and the alignment metrics, and the trimming summary gives the check that matters more:
an adapter that is not in the reads is trimmed from almost none of them, and every read then
carries linker sequence into alignment. That failure is silent — the reads still align, softly
clipped — so the rate at which the adapter was found is asserted rather than reported.
"""

import json
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

__all__ = [
    "ADAPTER_REACHED_BY",
    "MINIMUM_ADAPTER_RATE",
    "AdapterNotFoundError",
    "TrimSummary",
    "summarise_alignment",
    "summarise_trimming",
]

ADAPTER_REACHED_BY = frozenset({"riboseq"})
"""Assays whose insert is short enough that every read runs into the adapter."""

MINIMUM_ADAPTER_RATE = 0.5
"""Below this share of reads carrying the declared adapter, the adapter is presumed wrong.

Only libraries whose insert is shorter than the read can be held to this. A footprint library sits
near 1.0, because every read runs off the end of a 30 nt fragment into the linker. A transcriptome
library sits far lower for a legitimate reason — most fragments are longer than the read, so the
adapter is never reached — and asserting the same floor there would reject sound data.
"""


class AdapterNotFoundError(ValueError):
    pass


@dataclass(frozen=True)
class TrimSummary:
    """One library's read counts either side of trimming."""

    sample: str
    reads_raw: int
    reads_cleaned: int
    reads_with_adapter: int
    basepairs_raw: int
    basepairs_cleaned: int

    @property
    def adapter_rate(self) -> float:
        return self.reads_with_adapter / self.reads_raw

    @property
    def retained(self) -> float:
        return self.reads_cleaned / self.reads_raw


def _read_report(path: Path) -> TrimSummary:
    report = json.loads(path.read_text())
    counts = report["read_counts"]
    basepairs = report["basepair_counts"]
    return TrimSummary(
        sample=path.stem.removesuffix(".cutadapt"),
        reads_raw=counts["input"],
        reads_cleaned=counts["output"],
        reads_with_adapter=counts["read1_with_adapter"],
        basepairs_raw=basepairs["input"],
        basepairs_cleaned=basepairs["output"],
    )


def summarise_trimming(
    reports: list[Path],
    adapter_expected: Collection[str] = (),
    minimum_rate: float = MINIMUM_ADAPTER_RATE,
) -> pd.DataFrame:
    """Read cutadapt's JSON reports into one table.

    A library named in `adapter_expected` is refused if its declared adapter was found in too few
    reads, which is what a wrong adapter looks like.
    """

    summaries = sorted((_read_report(path) for path in reports), key=lambda s: s.sample)
    checked = [s for s in summaries if s.sample in adapter_expected]
    if failed := [s for s in checked if s.adapter_rate < minimum_rate]:
        detail = ", ".join(f"{s.sample} {s.adapter_rate:.1%}" for s in failed)
        raise AdapterNotFoundError(
            f"the declared adapter was found in under {minimum_rate:.0%} of reads: {detail}"
        )
    return pd.DataFrame(
        {
            "sample": [s.sample for s in summaries],
            "reads_raw": [s.reads_raw for s in summaries],
            "reads_cleaned": [s.reads_cleaned for s in summaries],
            "reads_retained": [round(s.retained, 5) for s in summaries],
            "reads_with_adapter": [s.reads_with_adapter for s in summaries],
            "adapter_rate": [round(s.adapter_rate, 5) for s in summaries],
            "basepairs_raw": [s.basepairs_raw for s in summaries],
            "basepairs_cleaned": [s.basepairs_cleaned for s in summaries],
        }
    )


# STAR writes its final log as "label |\tvalue", with percentages carrying a trailing sign.
_STAR_FIELDS = {
    "Number of input reads": "reads_input",
    "Uniquely mapped reads number": "reads_unique",
    "Uniquely mapped reads %": "unique_rate",
    "% of reads mapped to multiple loci": "multimapped_rate",
    "% of reads mapped to too many loci": "over_multimapped_rate",
    "% of reads unmapped: too short": "unmapped_too_short_rate",
    "Average mapped length": "mapped_length_mean",
}


def _read_star_log(path: Path) -> dict[str, float | str]:
    found: dict[str, float | str] = {"sample": path.name.split(".")[0]}
    for line in path.read_text().splitlines():
        label, _, value = line.partition("|")
        if (field := _STAR_FIELDS.get(label.strip())) is not None:
            found[field] = float(value.strip().rstrip("%"))
    return found


def summarise_alignment(logs: list[Path]) -> pd.DataFrame:
    """Read STAR's final logs into the alignment metrics table, one row per library."""

    summary = pd.DataFrame(
        sorted((_read_star_log(path) for path in logs), key=lambda r: r["sample"])
    )
    rates = [c for c in summary.columns if c.endswith("_rate")]
    summary["mapped_rate"] = summary[rates].drop(columns="unmapped_too_short_rate").sum(axis=1)
    return summary
