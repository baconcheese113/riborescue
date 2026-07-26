"""Reading sequence files, and the codon facts every layer reads them against.

The stop codons and the genetic code come from Biopython's standard table rather than being
transcribed, so a codon assignment has one definition and no layer can drift from another. Parsing
is Biopython's too: a FASTA differs between sources only in how its header names the record, which
is what `key` supplies.
"""

import gzip
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import IO

from Bio import SeqIO
from Bio.Data.CodonTable import standard_dna_table
from Bio.SeqRecord import SeqRecord

__all__ = [
    "STOP_CODONS",
    "accession",
    "gencode_accession",
    "open_text",
    "read_fasta",
    "records",
]

STOP_CODONS = frozenset(standard_dna_table.stop_codons)


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
