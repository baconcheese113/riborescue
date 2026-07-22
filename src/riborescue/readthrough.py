"""Whether a compound makes ribosomes carry on past a stop codon.

Readthrough is a redistribution and both halves are required: occupancy at the termination site
falls while in-frame occupancy downstream of it rises. Out-of-frame downstream occupancy is carried
alongside as the control, because readthrough continues the reading frame and must raise the
in-frame share specifically.

Ratios are formed within a transcript, against that transcript's own coding sequence, so neither
library depth nor transcript abundance enters, which matters here, because the footprint libraries
have no matched RNA-seq and abundance is therefore unknown.
"""

import gzip
import re
from collections import defaultdict
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "MINIMUM_CDS_PSITES",
    "MINIMUM_UTR3_LENGTH",
    "PROGRAMMED_READTHROUGH",
    "PairedEffect",
    "frame_specific",
    "library_ratios",
    "overlapping_downstream_cds",
    "paired_effect",
    "qualifying",
    "signature",
    "transcript_genes",
]

_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
_BIN = 100_000

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
    the end of the transcript without meeting another stop.
    """

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
            for other_start, other_end, other_gene in coding.get(
                (field[0], field[6], window), ()
            ):
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
        & (counts["cds_inframe"] >= MINIMUM_CDS_PSITES)
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
    libraries = counts["sample"].nunique()
    everywhere = passing.groupby("transcript")["sample"].nunique() == libraries
    shared = everywhere.index[everywhere]
    return passing.loc[passing["transcript"].isin(shared)].copy()


def library_ratios(counts: pd.DataFrame) -> pd.DataFrame:
    """Per library, the median transcript's readthrough, termination and out-of-frame ratios.

    A median over transcripts rather than a pooled sum, so that a handful of highly expressed genes
    do not decide the answer for the library.
    """

    per_transcript = counts.assign(
        readthrough=counts["extension_inframe"] / counts["cds_inframe"],
        out_of_frame=counts["extension_outframe"] / counts["cds_inframe"],
        termination=counts["termination"] / counts["cds_inframe"],
    )
    grouped = per_transcript.groupby("sample")
    return pd.DataFrame(
        {
            "transcripts": grouped.size(),
            "readthrough": grouped["readthrough"].median(),
            "out_of_frame": grouped["out_of_frame"].median(),
            "termination": grouped["termination"].median(),
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
        critical = {2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571}.get(n - 1, 1.96)
        half = critical * values.std(ddof=1) / np.sqrt(n)
        return (float(values.mean() - half), float(values.mean() + half))

    @property
    def consistent(self) -> bool:
        """Whether every replicate moved the same way."""

        return bool(np.all(np.asarray(self.differences) > 0) or
                    np.all(np.asarray(self.differences) < 0))


def frame_specific(in_frame: PairedEffect, out_of_frame: PairedEffect) -> bool:
    """Whether the rise is specific to the reading frame, replicate by replicate.

    Readthrough continues the frame, so it must raise in-frame occupancy by more than out-of-frame
    occupancy in every replicate. A treatment that lifts both equally has changed something that is
    not decoding — degradation, contamination or mapping — and the signature is not readthrough.
    """

    gains = np.asarray(in_frame.differences) - np.asarray(out_of_frame.differences)
    return bool(np.all(gains > 0))


def signature(effects: dict[str, PairedEffect]) -> dict[str, bool]:
    """The three conditions required together, each reported separately."""

    downstream = effects["readthrough"]
    termination = effects["termination"]
    return {
        "downstream_rose": downstream.mean_difference > 0 and downstream.consistent,
        "termination_fell": termination.mean_difference < 0 and termination.consistent,
        "frame_specific": frame_specific(downstream, effects["out_of_frame"]),
    }


def paired_effect(ratios: pd.DataFrame, quantity: str, replicates: pd.DataFrame) -> PairedEffect:
    """Compare treated against untreated within each replicate, for one quantity.

    The replicate is the unit of inference. The three experiments were prepared differently and
    preparation moves these ratios more than treatment plausibly does, so a treated library is only
    ever compared with the untreated library prepared beside it.
    """

    merged = ratios.merge(replicates, on="sample", validate="one_to_one")
    wide = merged.pivot(index="replicate", columns="treatment", values=quantity)
    missing = wide.isna().any(axis=1)
    if missing.any():
        raise ValueError(f"replicates without both conditions: {list(wide.index[missing])}")
    differences = (wide["g418"] - wide["untreated"]).to_numpy(dtype=float)
    ratios_ = (wide["g418"] / wide["untreated"]).to_numpy(dtype=float)
    return PairedEffect(quantity=quantity, differences=tuple(differences), ratios=tuple(ratios_))
