"""Gene expression from the matched RNA-seq arm: counts to TPM, and what is most expressed.

The footprint libraries say where ribosomes sit; they do not say how much message there is. That is
what the matched RNA-seq libraries are for, and the two must not be confused — occupancy is not
abundance, and a gene can carry many footprints because it is heavily translated or merely because
it is heavily transcribed.

Counts come from featureCounts over union exons, so the length used as the denominator is the union
exon length of the gene rather than any one transcript's. **TPM, not FPKM**: dividing by length
before normalising to the library total makes the values sum to a million in every library, which is
what allows one library's number to be compared with another's.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

__all__ = [
    "composition",
    "gene_class",
    "gene_symbols",
    "library_depth",
    "read_counts",
    "top_expressed",
    "tpm",
]

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_NAME = re.compile(r'gene_name "([^"]+)"')

# Depletion removes cytoplasmic ribosomal RNA; it does not touch the mitochondrial rRNAs or the
# small structural RNAs, which therefore survive in whatever quantity the protocol left behind. They
# are named here because TPM normalises to the library total: a library carrying more of them
# deflates every other gene in it, so their share has to be reported before any figure is read.
_MITOCHONDRIAL = re.compile(r"^MT-")
_STRUCTURAL = re.compile(r"^(RN7SL|RN7SK|RNU\d|RNA5S|RNA5-8S|RNA18S|RNA28S|RMRP|RPPH1|RNY\d)")


def gene_class(symbol: str) -> str:
    """Whether a gene is mitochondrial, a small structural RNA, or ordinary nuclear message."""

    if _MITOCHONDRIAL.match(symbol):
        return "mitochondrial"
    if _STRUCTURAL.match(symbol):
        return "structural"
    return "nuclear"


def read_counts(path: Path) -> pd.DataFrame:
    """featureCounts output as gene by library counts, with the union exon length alongside.

    The libraries arrive named by BAM path; they are renamed to the sample the path was built from,
    because a column called `results/reads/alignments/x.sorted.bam` is not a library name.
    """

    counts = pd.read_csv(path, sep="\t", comment="#")
    renamed = {
        column: Path(column).name.replace(".sorted.bam", "")
        for column in counts.columns
        if column.endswith(".bam")
    }
    counts = counts.rename(columns=renamed)
    keep = ["Geneid", "Length", *renamed.values()]
    return counts[keep].rename(columns={"Geneid": "gene_id", "Length": "length"})


def gene_symbols(gtf: Path) -> dict[str, str]:
    """Gene id to gene symbol, read from the annotation's gene records only.

    The counts are keyed by the stable gene id, which is what makes them joinable, but a table of
    the most expressed genes is unreadable without symbols.
    """

    symbols: dict[str, str] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            if len(fields) > 2 and fields[2] == "gene":
                identifier = _GENE_ID.search(line)
                name = _GENE_NAME.search(line)
                if identifier is not None and name is not None:
                    symbols[identifier.group(1)] = name.group(1)
    return symbols


def library_depth(counts: pd.DataFrame) -> pd.Series:
    """Assigned reads per library — the denominator behind every TPM, reported not hidden."""

    libraries = [column for column in counts.columns if column not in ("gene_id", "length")]
    return counts[libraries].sum()


def tpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Transcripts per million per gene per library.

    Reads per kilobase first, then normalised to that library's total — the order that makes the
    columns comparable across libraries of different depth.
    """

    libraries = [column for column in counts.columns if column not in ("gene_id", "length")]
    per_kilobase = counts[libraries].div(counts["length"] / 1000, axis=0)
    scaled = per_kilobase.div(per_kilobase.sum(), axis=1) * 1e6
    return pd.concat([counts[["gene_id"]], scaled], axis=1)


def composition(expression: pd.DataFrame, symbols: dict[str, str]) -> pd.DataFrame:
    """The share of each library that is mitochondrial, structural, or nuclear message.

    This is the figure to read before the expression figure. Every TPM is a fraction of its own
    library, so two libraries whose surviving structural RNA differs are not measuring the same
    denominator, and a gene can appear to move between them without changing at all.
    """

    libraries = [column for column in expression.columns if column != "gene_id"]
    classes = pd.Series(
        [gene_class(symbols.get(str(gene), "")) for gene in expression["gene_id"]],
        index=expression.index,
    )
    shares = expression.groupby(classes)[libraries].sum() / 1e4
    return shares.rename_axis("gene_class").reset_index()


def top_expressed(
    expression: pd.DataFrame,
    symbols: dict[str, str],
    top: int = 20,
    exclude: tuple[str, ...] = (),
) -> pd.DataFrame:
    """The most expressed genes by mean TPM across the libraries, with each library's own value.

    Ranking on the mean rather than on any one library keeps the table from being a description of
    whichever library was sequenced deepest, and every library's value stays in the row so a gene
    carried by a single library is visible as one.

    `exclude` drops whole gene classes — naming `mitochondrial` and `structural` gives the ranking
    of ordinary nuclear message, which is what the cell's expression programme looks like once the
    surviving contaminants of the preparation are set aside. Both rankings are worth having: the
    unfiltered one is what the libraries contain, the filtered one is what the cell was doing.
    """

    libraries = [column for column in expression.columns if column != "gene_id"]
    named = expression.assign(
        gene_symbol=[
            symbols.get(str(gene), str(gene).split(".")[0]) for gene in expression["gene_id"]
        ]
    )
    if exclude:
        classes = named["gene_symbol"].map(gene_class)
        named = named[~classes.isin(exclude)]
    ranked = named.assign(mean_tpm=named[libraries].mean(axis=1)).nlargest(top, "mean_tpm")
    return ranked[["gene_id", "gene_symbol", "mean_tpm", *libraries]].reset_index(drop=True)
