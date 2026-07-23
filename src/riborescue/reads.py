"""What trimming and alignment did to each library.

Cutadapt and STAR each report their own counts. Summarising them gives the raw-against-cleaned
comparison and the alignment metrics, and the trimming summary gives the check that matters more:
an adapter that is not in the reads is trimmed from almost none of them, and every read then
carries linker sequence into alignment. That failure is silent — the reads still align, softly
clipped — so the rate at which the adapter was found is asserted rather than reported.
"""

import gzip
import json
import statistics
from collections.abc import Collection, Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

__all__ = [
    "ADAPTER_REACHED_BY",
    "MINIMUM_ADAPTER_RATE",
    "SURVEY_READS",
    "AdapterNotFoundError",
    "AdapterSurvey",
    "TrimSummary",
    "summarise_alignment",
    "summarise_trimming",
    "survey_adapter",
    "survey_adapters",
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


SURVEY_READS = 200_000
"""Reads drawn from the head of a library to survey it.

The head of a FASTQ is not a random sample of the flow cell, but an adapter is either in the library
chemistry or it is not, and a fifth of a million reads settles that to well under a percent.
"""


@dataclass(frozen=True)
class AdapterSurvey:
    """What a library's own reads say about the adapter it was declared with."""

    sample: str
    adapter: str
    reads: int
    reads_with_adapter: int
    insert_median: int | None
    insert_p10: int | None
    insert_p90: int | None

    @property
    def adapter_rate(self) -> float:
        return self.reads_with_adapter / self.reads


def _sequences(fastq: Path, limit: int) -> Iterator[str]:
    with gzip.open(fastq, "rt") as handle:
        for index, line in enumerate(handle):
            if index >= limit * 4:
                return
            if index % 4 == 1:
                yield line.rstrip("\n")


def survey_adapter(
    sample: str, fastq: Path, adapter: str, limit: int = SURVEY_READS
) -> AdapterSurvey:
    """Look for the declared adapter in a library's reads, and measure what sits in front of it.

    A degenerate prefix is not searched for — it is random by construction, so only the fixed part
    of the linker identifies it, and the insert length is measured to where that fixed part begins.
    """

    fixed = adapter.lstrip("N")
    reads = inserts = 0
    lengths: list[int] = []
    for sequence in _sequences(fastq, limit):
        reads += 1
        found = sequence.find(fixed)
        if found >= 0:
            inserts += 1
            lengths.append(found)
    if reads == 0:
        raise ValueError(f"{fastq} holds no reads")
    quantiles = statistics.quantiles(lengths, n=10) if len(lengths) > 10 else None
    return AdapterSurvey(
        sample=sample,
        adapter=adapter,
        reads=reads,
        reads_with_adapter=inserts,
        insert_median=round(statistics.median(lengths)) if lengths else None,
        insert_p10=round(quantiles[0]) if quantiles else None,
        insert_p90=round(quantiles[-1]) if quantiles else None,
    )


def survey_adapters(runs: pd.DataFrame, minimum_rate: float = MINIMUM_ADAPTER_RATE) -> pd.DataFrame:
    """Survey every library whose insert is short enough to reach its adapter.

    This runs on the staged reads, before trimming and hours of alignment. Three of the four series
    this project uses describe their adapter wrongly or not at all, and a wrong adapter is silent
    downstream: the reads still align, softly clipped, carrying linker into every P-site.
    """

    footprints = runs.loc[runs["assay"].isin(ADAPTER_REACHED_BY)]
    surveys = [
        survey_adapter(str(row["sample"]), Path(str(row["fastq_1"])), str(row["adapter_3p"]))
        for _, row in footprints.iterrows()
    ]
    if failed := [s for s in surveys if s.adapter_rate < minimum_rate]:
        detail = ", ".join(f"{s.sample} {s.adapter_rate:.1%}" for s in failed)
        raise AdapterNotFoundError(
            f"the declared adapter is in under {minimum_rate:.0%} of the reads: {detail}"
        )
    return pd.DataFrame(
        {
            "sample": [s.sample for s in surveys],
            "adapter_3p": [s.adapter for s in surveys],
            "reads_surveyed": [s.reads for s in surveys],
            "adapter_rate": [round(s.adapter_rate, 5) for s in surveys],
            "insert_p10": [s.insert_p10 for s in surveys],
            "insert_median": [s.insert_median for s in surveys],
            "insert_p90": [s.insert_p90 for s in surveys],
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


def _read_star_log(path: Path) -> dict[str, float | int | str]:
    found: dict[str, float | int | str] = {"sample": path.name.split(".")[0]}
    for line in path.read_text().splitlines():
        label, _, value = line.partition("|")
        if (field := _STAR_FIELDS.get(label.strip())) is not None:
            number = float(value.strip().rstrip("%"))
            found[field] = int(number) if field.startswith("reads_") else number
    return found


def summarise_alignment(logs: list[Path]) -> pd.DataFrame:
    """Read STAR's final logs into the alignment metrics table, one row per library."""

    summary = pd.DataFrame(
        sorted((_read_star_log(path) for path in logs), key=lambda r: r["sample"])
    )
    mapped = ["unique_rate", "multimapped_rate", "over_multimapped_rate"]
    summary["mapped_rate"] = summary[mapped].sum(axis=1).round(2)
    return summary
