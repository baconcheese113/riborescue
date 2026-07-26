"""MANE Select transcripts: their exon structure, their sequence, and the map between the two.

One transcript per gene, so every layer scores a variant on the same molecule.

Coordinates follow their sources: genomic positions are 1-based inclusive, as GFF and VCF give them,
and transcript offsets are 0-based, as Python indexes a string.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from riborescue.core.sequences import STOP_CODONS, read_fasta

__all__ = [
    "PRIMARY_CHROMOSOMES",
    "STOP_CODONS",
    "TranscriptModel",
    "load_transcripts",
    "read_annotation",
    "read_sequences",
]

PRIMARY_CHROMOSOMES = frozenset(f"chr{name}" for name in [*map(str, range(1, 23)), "X", "Y"])
"""The primary assembly. Alternate and fix patches carry duplicate genes on other coordinates."""

# Exons carry transcript_id, coding features only Parent, so the parent names both.
_TRANSCRIPT_ID = re.compile(r"Parent=rna-([^;]+)")
_GENE_ID = re.compile(r"GeneID:(\d+)")
_PROTEIN_ID = re.compile(r"protein_id=([^;]+)")


@dataclass(frozen=True)
class TranscriptModel:
    """One MANE Select transcript: where its exons sit, where its coding sequence starts, and its
    spliced sequence."""

    transcript_id: str
    gene_id: int
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    cds_start: int
    cds_end: int
    sequence: str
    protein_id: str

    def offset_of(self, position: int) -> int | None:
        """Return the 0-based offset of a genomic position in the spliced transcript.

        A position in an intron, or outside the transcript, has no offset and returns None.
        """

        ordered = self.exons if self.strand == "+" else tuple(reversed(self.exons))
        travelled = 0
        for start, end in ordered:
            if start <= position <= end:
                within = position - start if self.strand == "+" else end - position
                return travelled + within
            travelled += end - start + 1
        return None

    @property
    def exon_lengths(self) -> tuple[int, ...]:
        """Exon lengths in transcript order, reversing the genomic order on the minus strand."""

        ordered = self.exons if self.strand == "+" else tuple(reversed(self.exons))
        return tuple(end - start + 1 for start, end in ordered)

    @property
    def junction_offsets(self) -> tuple[int, ...]:
        """Transcript offsets where one exon ends and the next begins.

        The exon junction complex is deposited near these, and how far a premature stop sits from
        the last of them is what every nonsense-mediated decay rule turns on.
        """

        offsets, running = [], 0
        for length in self.exon_lengths[:-1]:
            running += length
            offsets.append(running)
        return tuple(offsets)

    @property
    def coding_offset(self) -> int | None:
        """The 0-based offset at which the coding sequence starts."""

        return self.offset_of(self.cds_start if self.strand == "+" else self.cds_end)

    @property
    def coding_end_offset(self) -> int | None:
        """The 0-based offset of the last coding base — the 3' base of the natural stop codon."""

        return self.offset_of(self.cds_end if self.strand == "+" else self.cds_start)

    def codon_covering(self, position: int) -> tuple[int, str] | None:
        """The codon containing a genomic position: its offset in the transcript, and its bases.

        A position 3' of the natural stop — in the 3' UTR — has no reading-frame codon: translation
        has already terminated, so a substitution there is not a premature stop however the bases
        fall. Such positions return None rather than a spurious codon past the coding sequence.
        """

        offset, coding = self.offset_of(position), self.coding_offset
        coding_end = self.coding_end_offset
        if offset is None or coding is None or offset < coding:
            return None
        start = coding + ((offset - coding) // 3) * 3
        if coding_end is not None and start + 2 > coding_end:
            return None
        codon = self.sequence[start : start + 3]
        return (start, codon) if len(codon) == 3 else None


def read_annotation(gff: Path) -> pd.DataFrame:
    """Read the exon and coding features of every MANE Select transcript from a RefSeq GFF.

    MANE Plus Clinical transcripts and the alternate assembly are left out: the transcript policy is
    one Select transcript per gene on the primary assembly, so no variant is scored twice on two
    different molecules.
    """

    frame = pd.read_csv(
        gff,
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attrs"],
        dtype={"chrom": str},
    )
    frame = frame[
        frame["feature"].isin({"exon", "CDS"})
        & frame["attrs"].str.contains("tag=MANE Select", regex=False)
        & frame["chrom"].isin(PRIMARY_CHROMOSOMES)
    ].copy()
    frame["transcript_id"] = frame["attrs"].str.extract(_TRANSCRIPT_ID, expand=False)
    frame["gene_id"] = frame["attrs"].str.extract(_GENE_ID, expand=False).astype("Int64")
    frame["protein_id"] = frame["attrs"].str.extract(_PROTEIN_ID, expand=False)
    return frame.dropna(subset=["transcript_id", "gene_id"]).drop(
        columns=["attrs", "source", "score"]
    )


def read_sequences(fasta: Path) -> dict[str, str]:
    """Read transcript sequences, keyed by accession, from a RefSeq FASTA."""

    return read_fasta(fasta)


def load_transcripts(gff: Path, fasta: Path) -> dict[str, TranscriptModel]:
    """Build the transcript models, keyed by accession, from an annotation and its sequences."""

    annotation = read_annotation(gff)
    sequences = read_sequences(fasta)
    models: dict[str, TranscriptModel] = {}

    for transcript_id, features in annotation.groupby("transcript_id", sort=False):
        sequence = sequences.get(str(transcript_id))
        coding = features[features["feature"] == "CDS"]
        exons = features[features["feature"] == "exon"]
        if sequence is None or coding.empty or exons.empty:
            continue
        models[str(transcript_id)] = TranscriptModel(
            transcript_id=str(transcript_id),
            gene_id=int(features["gene_id"].iloc[0]),
            chrom=str(features["chrom"].iloc[0]),
            strand=str(features["strand"].iloc[0]),
            exons=tuple(sorted(zip(exons["start"], exons["end"], strict=True))),
            cds_start=int(coding["start"].min()),
            cds_end=int(coding["end"].max()),
            sequence=sequence,
            protein_id=str(coding["protein_id"].iloc[0]),
        )
    return models
