"""The non-coding sequence a ribosome profiling library is mostly made of.

Two thirds of the footprints in an undepleted library are ribosomal RNA, and rDNA repeats and rRNA
pseudogenes sit in many places in the genome, so those reads align to many loci rather than to none.
Left in, they crowd out the coding signal and drag the unique-mapping rate to a fifth of the reads.

Depletion is an alignment against these sequences first, keeping only what fails to match. GENCODE
writes the biotype into each FASTA header, so most of the contaminant set is a selection from a file
already pinned rather than a separately curated reference.

The exception is ribosomal RNA itself. GENCODE annotates almost none, because the rDNA arrays are
not assembled into the primary reference, so the repeating unit is carried separately. Selecting on
biotype alone removes a tenth of a footprint library; with the repeat unit it removes three fifths.
"""

import gzip
from collections.abc import Iterator, Sequence
from pathlib import Path

__all__ = ["CONTAMINANT_BIOTYPES", "write_contaminants"]

CONTAMINANT_BIOTYPES = frozenset(
    {
        "rRNA",
        "rRNA_pseudogene",
        "Mt_rRNA",
        "Mt_tRNA",
        "snRNA",
        "snoRNA",
        "scaRNA",
        "misc_RNA",
        "vault_RNA",
        "sRNA",
        "scRNA",
    }
)
"""What a footprint library contains besides messenger RNA.

Structural and small non-coding RNA, not the long non-coding transcripts: those are translated
often enough that discarding them would remove signal rather than contamination.
"""

# GENCODE headers are pipe-separated, with the transcript biotype in the eighth field.
_BIOTYPE_FIELD = 7


def _records(handle: Iterator[str]) -> Iterator[tuple[str, list[str]]]:
    header, sequence = None, []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, sequence
            header, sequence = line, []
        elif header is not None:
            sequence.append(line)
    if header is not None:
        yield header, sequence


def write_contaminants(transcripts: Path, out: Path, include: Sequence[Path] = ()) -> int:
    """Write the structural RNA from a GENCODE transcript FASTA, plus any FASTA named in `include`.

    Returns how many sequences were written.
    """

    selected = 0
    with gzip.open(transcripts, "rt") as source, out.open("w") as sink:
        for header, sequence in _records(source):
            fields = header.rstrip("\n").split("|")
            if len(fields) > _BIOTYPE_FIELD and fields[_BIOTYPE_FIELD] in CONTAMINANT_BIOTYPES:
                sink.write(header)
                sink.writelines(sequence)
                selected += 1
        if selected == 0:
            raise ValueError(f"{transcripts} yielded no structural RNA; is it GENCODE FASTA?")

        added = 0
        for extra in include:
            text = extra.read_text()
            sink.write(text if text.endswith("\n") else text + "\n")
            added += text.count(">")
    return selected + added
