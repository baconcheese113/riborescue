"""Reading sequence files, and the codon facts every layer reads them against.

The stop codons and the genetic code come from Biopython's standard table rather than being
transcribed, so a codon assignment has one definition and no layer can drift from another. Parsing
is Biopython's too: a FASTA differs between sources only in how its header names the record, which
is what `key` supplies.

The two spellings live here for the same reason. A codon is written `TTC` where it indexes the
genetic code and `uuc` where it names a level of the readthrough model, and every layer that crosses
between them does so through `as_dna` and `as_rna` rather than its own case-and-substitute.
"""

import gzip
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import IO

from Bio import SeqIO
from Bio.Data.CodonTable import standard_dna_table
from Bio.SeqRecord import SeqRecord

__all__ = [
    "GENETIC_CODE",
    "SENSE_CODONS",
    "STOP_CODONS",
    "SYNONYMOUS",
    "accession",
    "as_dna",
    "as_rna",
    "gencode_accession",
    "gencode_sequences",
    "open_text",
    "read_fasta",
    "records",
]

STOP_CODONS = frozenset(standard_dna_table.stop_codons)
SENSE_CODONS = tuple(sorted(standard_dna_table.forward_table))

GENETIC_CODE: dict[str, str] = standard_dna_table.forward_table | dict.fromkeys(STOP_CODONS, "*")
"""The standard genetic code, stops as `*`, keyed on the DNA spelling."""

SYNONYMOUS: dict[str, tuple[str, ...]] = {
    residue: tuple(sorted(c for c in SENSE_CODONS if GENETIC_CODE[c] == residue))
    for residue in sorted({GENETIC_CODE[c] for c in SENSE_CODONS})
}
"""Each amino acid's codons. Methionine and tryptophan have one each and are invariant under the
context-matched shuffle, which is reported rather than hidden."""

_TO_DNA = str.maketrans("acgtuACGTU", "ACGTTACGTT")
_TO_RNA = str.maketrans("acgtuACGTU", "acguuacguu")


def as_dna(sequence: str) -> str:
    """Upper-case DNA. `uuc` and `TTC` name the same codon; the tables are keyed on the latter."""

    return sequence.translate(_TO_DNA)


def as_rna(sequence: str) -> str:
    """Lower-case RNA — how the readthrough assay spelled the features the model was fitted on."""

    return sequence.translate(_TO_RNA)


def open_text(path: Path) -> IO[str]:
    """Open a text file, decompressing it where the name says it is compressed."""

    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def accession(record: SeqRecord) -> str:
    """The first whitespace-delimited token of the header — how RefSeq names a record."""

    return str(record.id)


def gencode_accession(record: SeqRecord) -> str:
    """The first pipe-delimited field — how GENCODE names a record."""

    return str(record.id).split("|", 1)[0]


def records(path: Path) -> Iterator[SeqRecord]:
    """Every record of a FASTA, compressed or not."""

    with open_text(path) as handle:
        yield from SeqIO.parse(handle, "fasta")


def read_fasta(path: Path, key: Callable[[SeqRecord], str] = accession) -> dict[str, str]:
    """A FASTA as upper-case sequences, keyed by whatever `key` reads out of each header.

    Upper case because the callers compare against codons: a soft-masked reference would otherwise
    match no stop codon and read as a transcript with no termination signal at all.
    """

    return {key(record): str(record.seq).upper() for record in records(path)}


def gencode_sequences(transcripts: Path) -> dict[str, str]:
    """Transcript sequences keyed by accession, from a GENCODE FASTA."""

    return read_fasta(transcripts, key=gencode_accession)
