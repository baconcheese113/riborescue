"""Whether a compound makes ribosomes carry on past a stop codon.

Readthrough is a redistribution and both halves are required: occupancy at the termination site
falls while in-frame occupancy downstream of it rises. Out-of-frame downstream occupancy is carried
alongside as the control, because readthrough continues the reading frame and must raise the
in-frame share specifically.

Every quantity is a ratio taken within one library against that library's own coding sequence, so
library depth cancels. Transcript abundance does not: counts are pooled before dividing, so a highly
expressed transcript carries more weight than a lightly expressed one. That is deliberate for a
proportion, and it is why the per-transcript median and the share of transcripts contributing
anything are reported beside every pooled figure. These libraries have no matched RNA-seq, so
abundance is unknown and cannot be adjusted for.
"""

import gzip
import math
import re
from collections import defaultdict
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "MINIMUM_CDS_PSITES",
    "MINIMUM_UTR3_LENGTH",
    "PROGRAMMED_READTHROUGH",
    "Effect",
    "PairedEffect",
    "UnpairedEffect",
    "library_ratios",
    "overlapping_downstream_cds",
    "paired_effect",
    "qualifying",
    "signature",
    "stalling",
    "transcript_genes",
    "unpaired_effect",
]

_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
_BIN = 100_000

# 95% two-sided t, by degrees of freedom, rather than pulling in scipy for a handful of points.
_T_CRITICAL = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306}

CDS_FRAMES = ["cds_frame0", "cds_frame1", "cds_frame2"]
EXTENSION_FRAMES = ["extension_frame0", "extension_frame1", "extension_frame2"]

MINIMUM_UTR3_LENGTH = 50
"""Shorter than this, the termination peak spills into the window being measured."""

MINIMUM_CDS_PSITES = 100
"""Below this the per-transcript ratio has no stable denominator."""

PROGRAMMED_READTHROUGH = frozenset(
    {
        "AQP4",
        "MAPK10",
        "OPRK1",
        "OPRL1",
        "VDR",
        "MTCH2",
        "AGO1",
        "LDHB",
        "SACM1L",
        "MDH1",
        "EEF1B2",
        "TMED2",
        "CDKN2A",
        "BRI3",
        "SLC35A4",
    }
)
"""Genes reported to read through their stop codon natively, which do so with or without a drug."""


STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


def _sequences(transcripts: Path) -> dict[str, str]:
    """Transcript sequences keyed by accession, from a GENCODE FASTA."""

    found: dict[str, str] = {}
    name, parts = None, []
    with gzip.open(transcripts, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    found[name] = "".join(parts)
                name, parts = line[1:].split("|", 1)[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        found[name] = "".join(parts)
    return found


def next_in_frame_stop(sequence: str, stop_start: int) -> int | None:
    """How far past the native stop the next stop in the same frame begins.

    `stop_start` indexes the first base of the native stop codon. The answer is a multiple of three
    and is the width of the window a readthrough ribosome can occupy. None where the frame runs off
    the end of the transcript without meeting another stop, and none where the annotated coding
    sequence does not in fact end in a stop codon: an incomplete annotation gives a length that is
    not a multiple of three, and every frame measured from it would be off by one or two bases.
    """

    if sequence[stop_start : stop_start + 3].upper() not in STOP_CODONS:
        return None
    offset = 3
    position = stop_start + 3
    while position + 3 <= len(sequence):
        if sequence[position : position + 3].upper() in STOP_CODONS:
            return offset
        position += 3
        offset += 3
    return None


def extension_windows(transcripts: Path, annotation: pd.DataFrame) -> pd.DataFrame:
    """The extension window of every coding transcript, one row each.

    Looked up by dictionary rather than by scanning the sequence names once per transcript, which
    turns the whole build from quadratic into a single pass.
    """

    sequences = _sequences(transcripts)
    coding = annotation.loc[annotation["l_cds"] > 0]
    rows: list[tuple[str, int | None]] = []
    for transcript, utr5, cds in zip(
        coding["transcript"], coding["l_utr5"], coding["l_cds"], strict=True
    ):
        sequence = sequences.get(transcript)
        if sequence is None:
            continue
        # The coding sequence ends with its stop codon, so the stop begins three bases back.
        rows.append((transcript, next_in_frame_stop(sequence, int(utr5) + int(cds) - 3)))
    return pd.DataFrame(rows, columns=["transcript", "extension"])


def _features(gtf: Path, wanted: frozenset[str]):
    """Yield the GTF rows of the requested feature types, already split."""

    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            field = line.split("\t")
            if field[2] in wanted:
                yield field


def transcript_genes(gtf: Path) -> pd.Series:
    """Map each transcript to the gene symbol it belongs to."""

    found: dict[str, str] = {}
    for field in _features(gtf, frozenset({"transcript"})):
        attributes = field[8]
        transcript = _TRANSCRIPT_ID.search(attributes)
        name = _GENE_NAME.search(attributes)
        if transcript and name:
            found[transcript.group(1)] = name.group(1)
    return pd.Series(found, name="gene_name")


def overlapping_downstream_cds(gtf: Path) -> frozenset[str]:
    """Transcripts whose 3' untranslated region runs into another gene's coding sequence.

    Ribosomes on the neighbouring gene land in this transcript's downstream window and read as
    readthrough that never happened. Only a different gene counts: a transcript's own coding
    sequence overlapping its own untranslated region is ordinary isoform structure.
    """

    coding: dict[tuple[str, str, int], list[tuple[int, int, str]]] = defaultdict(list)
    for field in _features(gtf, frozenset({"CDS"})):
        gene = _GENE_ID.search(field[8])
        if gene is None:
            continue
        start, end = int(field[3]), int(field[4])
        for window in range(start // _BIN, end // _BIN + 1):
            coding[(field[0], field[6], window)].append((start, end, gene.group(1)))

    contaminated: set[str] = set()
    for field in _features(gtf, frozenset({"three_prime_utr"})):
        gene = _GENE_ID.search(field[8])
        transcript = _TRANSCRIPT_ID.search(field[8])
        if gene is None or transcript is None:
            continue
        start, end = int(field[3]), int(field[4])
        for window in range(start // _BIN, end // _BIN + 1):
            for other_start, other_end, other_gene in coding.get((field[0], field[6], window), ()):
                if other_gene != gene.group(1) and other_start <= end and start <= other_end:
                    contaminated.add(transcript.group(1))
                    break
    return frozenset(contaminated)


def qualifying(
    counts: pd.DataFrame,
    genes: pd.Series | None = None,
    excluded_transcripts: frozenset[str] = frozenset(),
    samples: Collection[str] | None = None,
) -> pd.DataFrame:
    """The transcripts the comparison is allowed to use.

    Every exclusion is a property of the transcript rather than of the treatment, so the same rows
    are dropped from both conditions. `samples` narrows the table to the libraries the comparison
    actually uses, before the shared universe is worked out: a library that takes no part must not
    decide which transcripts the ones that do are allowed to keep.
    """

    if samples is not None:
        counts = counts.loc[counts["sample"].isin(samples)]
    keep = (
        (counts["l_utr3"] >= MINIMUM_UTR3_LENGTH)
        & (counts["cds_frame0"] >= MINIMUM_CDS_PSITES)
        & counts["extension"].notna()
    )
    if genes is not None:
        symbol = counts["transcript"].map(genes)
        keep &= ~symbol.isin(PROGRAMMED_READTHROUGH)
    if excluded_transcripts:
        keep &= ~counts["transcript"].isin(excluded_transcripts)
    passing = counts.loc[keep]

    # One universe for every library. Coverage is a property of the library, so a threshold applied
    # library by library would let the treated and untreated medians be taken over different
    # transcripts, and the comparison would be partly about which transcripts cleared the bar.
    libraries = len(set(samples)) if samples is not None else counts["sample"].nunique()
    everywhere = passing.groupby("transcript")["sample"].nunique() == libraries
    shared = everywhere.index[everywhere]
    return passing.loc[passing["transcript"].isin(shared)].copy()


def library_ratios(counts: pd.DataFrame) -> pd.DataFrame:
    """Per library, the quantities the control is judged on.

    Counts are summed over the library's qualifying transcripts before any ratio is taken. A frame
    composition is a proportion, and the per-transcript median of one sits at zero whenever most
    transcripts carry no downstream signal, which is the usual case. The median and the share of
    transcripts carrying anything are reported beside the pooled figures so that a proportion
    resting on few transcripts stays visible.
    """

    per_transcript = counts.assign(
        _rt=counts["extension_frame0"] / counts["cds_frame0"],
        _any=(counts[EXTENSION_FRAMES].sum(axis=1) > 0),
    )
    pooled = counts.groupby("sample")[CDS_FRAMES + EXTENSION_FRAMES + ["termination"]].sum()
    grouped = per_transcript.groupby("sample")

    downstream_total = pooled[EXTENSION_FRAMES].sum(axis=1)
    cds_total = pooled[CDS_FRAMES].sum(axis=1)
    downstream_share0 = pooled["extension_frame0"] / downstream_total
    cds_share0 = pooled["cds_frame0"] / cds_total
    return pd.DataFrame(
        {
            "transcripts": grouped.size(),
            "with_downstream": grouped["_any"].mean(),
            "cds_frame0_total": pooled["cds_frame0"],
            "downstream_total": downstream_total,
            "termination_total": pooled["termination"],
            "downstream_occupancy": pooled["extension_frame0"] / pooled["cds_frame0"],
            "termination_occupancy": pooled["termination"] / pooled["cds_frame0"],
            "downstream_share0": downstream_share0,
            "cds_share0": cds_share0,
            "frame_gap": downstream_share0 - cds_share0,
            "median_readthrough": grouped["_rt"].median(),
        }
    ).reset_index()


@dataclass(frozen=True)
class PairedEffect:
    """One quantity compared treated against untreated, paired within replicate."""

    quantity: str
    differences: tuple[float, ...]
    ratios: tuple[float, ...]

    @property
    def mean_difference(self) -> float:
        return float(np.mean(self.differences))

    @property
    def interval(self) -> tuple[float, float]:
        """A t interval on the paired differences, at n = 3 a statement of direction not size."""

        values = np.asarray(self.differences, dtype=float)
        n = len(values)
        if n < 2:
            return (float("nan"), float("nan"))
        # 95% two-sided for n-1 degrees of freedom, read from the t table rather than pulled in
        # from scipy for three points.
        critical = _T_CRITICAL.get(n - 1, 1.96)
        half = critical * values.std(ddof=1) / np.sqrt(n)
        return (float(values.mean() - half), float(values.mean() + half))

    @property
    def consistent(self) -> bool:
        """Whether every replicate moved the same way."""

        return bool(
            np.all(np.asarray(self.differences) > 0) or np.all(np.asarray(self.differences) < 0)
        )


@dataclass(frozen=True)
class UnpairedEffect:
    """One quantity compared between two groups of libraries that are not matched.

    The primary shape for a dataset whose arms are separate libraries, where a treated library has
    no partner it was prepared beside.
    """

    quantity: str
    treated: tuple[float, ...]
    control: tuple[float, ...]

    @property
    def mean_difference(self) -> float:
        return float(np.mean(self.treated) - np.mean(self.control))

    @property
    def interval(self) -> tuple[float, float]:
        """A Welch interval, which does not assume the two groups vary alike."""

        a, b = np.asarray(self.treated, dtype=float), np.asarray(self.control, dtype=float)
        if len(a) < 2 or len(b) < 2:
            return (float("nan"), float("nan"))
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        error = np.sqrt(va + vb)
        if error == 0:
            return (self.mean_difference, self.mean_difference)
        df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
        # Rounding down is the conservative direction: fewer degrees of freedom, wider interval.
        critical = _T_CRITICAL.get(max(1, math.floor(df)), 1.96)
        half = critical * error
        return (self.mean_difference - half, self.mean_difference + half)

    @property
    def consistent(self) -> bool:
        """Whether the groups separate completely, every treated library past every control one."""

        a, b = np.asarray(self.treated), np.asarray(self.control)
        return bool(a.min() > b.max() or a.max() < b.min())


Effect = PairedEffect | UnpairedEffect


def signature(effects: Mapping[str, Effect]) -> dict[str, bool]:
    """The three conditions required together, each reported separately.

    The frame condition alone must clear an interval that excludes zero: it is the claim the control
    exists to make. The other two are directional, where agreement across replicates is the evidence
    and an interval from three points is too wide to ask more of.
    """

    downstream = effects["downstream_occupancy"]
    termination = effects["termination_occupancy"]
    gap = effects["frame_gap"]
    low, _ = gap.interval
    return {
        "downstream_rose": bool(downstream.mean_difference > 0 and downstream.consistent),
        "termination_fell": bool(termination.mean_difference < 0 and termination.consistent),
        "frame_moved_to_coding": bool(gap.mean_difference > 0 and gap.consistent and low > 0),
    }


def stalling(effects: Mapping[str, Effect]) -> bool:
    """Whether a compound behaves as a stalling agent rather than a readthrough one.

    The negative control has to fail the signature in the direction its mechanism predicts —
    occupancy at the stop does not fall, and occupancy beyond it does not rise — rather than fail
    the frame condition, which is unstable when there is little downstream signal to compose.
    """

    termination = effects["termination_occupancy"]
    downstream = effects["downstream_occupancy"]
    return bool(
        termination.mean_difference >= 0
        and downstream.mean_difference <= 0
        and termination.consistent
    )


def paired_effect(
    ratios: pd.DataFrame,
    quantity: str,
    replicates: pd.DataFrame,
    treated: str,
    control: str,
) -> PairedEffect:
    """Compare treated against control within each replicate, for one quantity.

    Only for a design that documents which libraries were prepared beside each other. Where the
    metadata does not say, matching replicate numbers across arms mean nothing and the unpaired
    comparison is the honest one.
    """

    merged = ratios.merge(replicates, on="sample", validate="one_to_one")
    wide = merged.pivot(index="replicate", columns="treatment", values=quantity)
    for arm in (treated, control):
        if arm not in wide.columns:
            raise ValueError(f"no libraries for {arm!r}")
    missing = wide[[treated, control]].isna().any(axis=1)
    if missing.any():
        raise ValueError(f"replicates without both conditions: {list(wide.index[missing])}")
    differences = (wide[treated] - wide[control]).to_numpy(dtype=float)
    ratios_ = (wide[treated] / wide[control]).to_numpy(dtype=float)
    return PairedEffect(quantity=quantity, differences=tuple(differences), ratios=tuple(ratios_))


def unpaired_effect(
    ratios: pd.DataFrame,
    quantity: str,
    replicates: pd.DataFrame,
    treated: str,
    control: str,
) -> UnpairedEffect:
    """Compare two groups of libraries that were not prepared as matched pairs."""

    merged = ratios.merge(replicates, on="sample", validate="one_to_one")
    groups = {arm: merged.loc[merged["treatment"] == arm, quantity] for arm in (treated, control)}
    for arm, values in groups.items():
        if values.empty:
            raise ValueError(f"no libraries for {arm!r}")
    return UnpairedEffect(
        quantity=quantity,
        treated=tuple(groups[treated].astype(float)),
        control=tuple(groups[control].astype(float)),
    )
