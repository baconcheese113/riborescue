"""Codon occupancy: how long ribosomes sit over each codon, measured transcriptome-wide.

ADR-0020. The quantity is **occupancy**, also called a pause score — the number of P-sites assigned
to a codon position, normalised within its own transcript. It is not a dwell time, a decoding rate
or a velocity: those are inferred from occupancy under assumptions about initiation and flux that
this project neither makes nor checks.

Two conventions decide what the number means, and both are named rather than assumed. The codon
being *decoded* sits in the **A site**, which is one codon 3' of the P site riboWaltz reports, and
the A site is what the kinetic features are built from; the P site is exported beside it. The
internal-CDS window belongs to whichever site is being scored, so a position enters on its own
offsets rather than on the ones the export happened to record.

Normalisation is within transcript before anything is pooled, which is what makes a codon score a
statement about codons rather than about which genes were most expressed. A position with no P-sites
is a measured zero and counts: dropping it would let a codon that ribosomes never pause on look
ordinary.
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

__all__ = [
    "GENETIC_CODE",
    "MINIMUM_CDS_PSITES",
    "MINIMUM_SCORABLE_CODONS",
    "SENSE_CODONS",
    "SITE_SHIFT",
    "STOP_CODONS",
    "SYNONYMOUS",
    "aggregate_libraries",
    "as_dna",
    "library_table",
    "scorable_bounds",
    "transcript_codons",
]

MINIMUM_CDS_PSITES = 100
"""In-frame P-sites a transcript needs over the scored window, matching the readthrough assay.

One convention across the project is worth more here than a separately optimal one.
"""

MINIMUM_SCORABLE_CODONS = 64
"""Scorable positions a transcript needs, so it cannot enter the table on a single window."""

# The internal-CDS window of ADR-0011, in codons. A position's first base must lie at least 18 nt
# past the start codon and at most -16 nt from the last base of the stop, which excludes the
# initiation and termination peaks from this quantity as it excludes them from every other
# denominator in the project.
_WINDOW_FROM_START = 18
_WINDOW_FROM_STOP = -16
_FIRST_SCORABLE = -(-_WINDOW_FROM_START // 3)

SITE_SHIFT = {"a": 1, "p": 0}
"""Codons to add to the P-site index to reach the scored site. The A site is one codon 3'."""

GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}  # fmt: skip
"""The standard genetic code. The synonymous families below are read from it, not listed twice."""

STOP_CODONS = frozenset(codon for codon, residue in GENETIC_CODE.items() if residue == "*")
SENSE_CODONS = tuple(sorted(set(GENETIC_CODE) - STOP_CODONS))

SYNONYMOUS: dict[str, tuple[str, ...]] = {
    residue: tuple(sorted(c for c in SENSE_CODONS if GENETIC_CODE[c] == residue))
    for residue in sorted({GENETIC_CODE[c] for c in SENSE_CODONS})
}
"""Each amino acid's codons. Methionine and tryptophan have one each and are invariant under the
context-matched shuffle, which is reported rather than hidden."""


def as_dna(triplet: str) -> str:
    """`uuc` and `TTC` name the same codon. The table is keyed on the DNA spelling."""

    return triplet.upper().replace("U", "T")


def scorable_bounds(l_cds: int) -> tuple[int, int]:
    """The inclusive codon-index range a transcript's scored site may occupy.

    Indices count codons from the start codon, which is index 0, so codon `i` begins `3i` bases into
    the coding sequence and its distance from the last base of the stop is `3i - l_cds + 1`.
    """

    return _FIRST_SCORABLE, (l_cds - 1 + _WINDOW_FROM_STOP) // 3


def transcript_codons(sequences: dict[str, str], annotation: pd.DataFrame) -> dict[str, str]:
    """Each coding transcript's codon string over its scorable window, one transcript per entry.

    Returned as the spliced substring rather than a list of codons: slicing it by index is how both
    the composition and the per-position lookup below are taken, and it holds a transcriptome of
    codons in one string each.
    """

    found: dict[str, str] = {}
    coding = annotation.loc[annotation["l_cds"] > 0]
    for transcript, utr5, cds in zip(
        coding["transcript"], coding["l_utr5"], coding["l_cds"], strict=True
    ):
        sequence = sequences.get(transcript)
        if sequence is None:
            continue
        first, last = scorable_bounds(int(cds))
        if last < first:
            continue
        start = int(utr5) + 3 * first
        window = sequence[start : int(utr5) + 3 * (last + 1)].upper()
        # A transcript whose recorded coding length runs past the sequence it was annotated on would
        # otherwise contribute a truncated final codon.
        found[transcript] = window[: len(window) - len(window) % 3]
    return found


def _composition(codons: str) -> Counter[str]:
    return Counter(codons[i : i + 3] for i in range(0, len(codons), 3))


def library_table(
    counts: pd.DataFrame,
    annotation: pd.DataFrame,
    codons: dict[str, str],
    *,
    site: str = "a",
    lengths: tuple[int, ...] | None = None,
) -> pd.DataFrame:
    """One library's codon occupancy, 64 rows, from its per-codon-position P-site counts.

    `counts` is what `run_psite.R` exported: in-frame P-sites per transcript, footprint length and
    P-site codon index. `lengths` narrows it to a dataset's calibrated set before anything is
    summed; passing none keeps every length, right only where the window is already fixed.

    A transcript's positions are normalised by that transcript's own mean occupancy, so abundance
    cancels, and a codon's score is the mean of those normalised values over every position of it in
    every qualifying transcript. Positions carrying no P-sites are in that mean.
    """

    if (shift := SITE_SHIFT.get(site)) is None:
        raise ValueError(f"site must be one of {sorted(SITE_SHIFT)}, not {site!r}")
    if lengths is not None:
        counts = counts.loc[counts["length"].isin(list(lengths))]
        if counts.empty:
            raise ValueError(f"no counts at lengths {sorted(lengths)}")

    # The window is taken from the codon strings rather than recomputed from `l_cds`, so a
    # transcript whose annotation outruns its sequence cannot be scored at a position the string
    # does not hold. `annotation` is still required, because a transcript absent from it is not a
    # coding transcript and takes no part.
    positions = pd.Series({t: len(c) // 3 for t, c in codons.items()}, name="positions")
    positions = positions.loc[positions.index.isin(set(annotation["transcript"]))]
    last = positions + (_FIRST_SCORABLE - 1)

    scored = counts.assign(index=counts["codon_index"] + shift)
    scored = scored.groupby(["transcript", "index"], as_index=False, sort=False)["n"].sum()
    scored = scored.join(last.rename("last"), on="transcript", how="inner")
    scored = scored.loc[(scored["index"] >= _FIRST_SCORABLE) & (scored["index"] <= scored["last"])]

    # A transcript enters on the P-sites it carries inside the window and on how much window it has,
    # never on its expression. Both floors are properties of the transcript, so the same transcripts
    # are dropped from every library of a comparison.
    per_transcript = scored.groupby("transcript")["n"].sum().to_frame("psites")
    per_transcript = per_transcript.join(positions, how="inner")
    qualifying = per_transcript.loc[
        (per_transcript["psites"] >= MINIMUM_CDS_PSITES)
        & (per_transcript["positions"] >= MINIMUM_SCORABLE_CODONS)
    ]
    if qualifying.empty:
        raise ValueError("no transcript clears the coverage floors")

    # The mean occupancy of a position in this transcript, which every one of its positions is
    # divided by. Taken over the whole window rather than over the positions that happen to carry a
    # P-site, so an empty position is a zero rather than an absence.
    mean_occupancy = qualifying["psites"] / qualifying["positions"]

    kept = scored.loc[scored["transcript"].isin(qualifying.index)]
    identity = [
        codons[transcript][3 * (index - _FIRST_SCORABLE) : 3 * (index - _FIRST_SCORABLE) + 3]
        for transcript, index in zip(kept["transcript"], kept["index"], strict=True)
    ]
    observed = (
        kept.assign(
            codon=identity,
            normalised=kept["n"].to_numpy() / mean_occupancy.loc[kept["transcript"]].to_numpy(),
        )
        .groupby("codon")["normalised"]
        .sum()
    )

    # Every scorable position of every qualifying transcript, whether it carried a P-site or not.
    # This is the denominator, and it is why a codon ribosomes never pause on scores near zero
    # instead of dropping out of the table.
    available: Counter[str] = Counter()
    for transcript in qualifying.index:
        available.update(_composition(codons[transcript]))

    table = pd.DataFrame(
        {
            "codon": SENSE_CODONS,
            "amino_acid": [GENETIC_CODE[c] for c in SENSE_CODONS],
            "positions": [available.get(c, 0) for c in SENSE_CODONS],
            "psites": [float(observed.get(c, 0.0)) for c in SENSE_CODONS],
        }
    )
    table["occupancy"] = table["psites"] / table["positions"].replace(0, np.nan)
    table["site"] = site
    table["transcripts"] = len(qualifying)
    return table[["codon", "amino_acid", "site", "occupancy", "positions", "transcripts"]]


def aggregate_libraries(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """The codon table across an arm's libraries, with the spread between them beside each score.

    The mean of the per-library scores, not a score over pooled counts: a library sequenced deeper
    would otherwise become the table. The standard deviation across libraries is what says whether a
    codon's score is a property of the cell or of one preparation.
    """

    stacked = pd.concat(
        [table.assign(sample=sample) for sample, table in tables.items()], ignore_index=True
    )
    grouped = stacked.groupby(["codon", "amino_acid", "site"], as_index=False)
    summarised = grouped.agg(
        occupancy=("occupancy", "mean"),
        occupancy_sd=("occupancy", lambda values: values.std(ddof=1)),
        positions=("positions", "min"),
        libraries=("occupancy", "size"),
    )
    return summarised.sort_values(["site", "codon"], ignore_index=True)
