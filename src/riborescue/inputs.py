"""Public inputs — declared once, fetched from their source, verified by checksum.

Nothing large is committed. Reproducibility rests on regeneration, so every input names where it
comes from, what it should hash to, and what licence it carries. Fetching is `pooch`'s job: it
downloads to a temporary file and only moves it into place once the digest matches, so a corrupt
or truncated file is never mistaken for the real input.
"""

import os
from dataclasses import dataclass
from pathlib import Path

import pooch

__all__ = ["INPUTS", "Input", "UnknownInputError", "data_root", "fetch"]


class UnknownInputError(KeyError):
    pass


@dataclass(frozen=True)
class Input:
    """One fetchable public file: where it lives, what it must hash to, and who published it."""

    name: str
    url: str
    md5: str
    path: Path
    source: str
    licence: str

    def resolve(self, root: Path) -> Path:
        return root / self.path


INPUTS: dict[str, Input] = {
    "toledano_treated_samples": Input(
        name="toledano_treated_samples",
        url="https://ndownloader.figshare.com/files/44388752",
        md5="bb3474efa01f074466912a36964827c2",
        path=Path("toledano/treated_samples.rds"),
        source="doi:10.6084/m9.figshare.25138712.v6",
        licence="CC BY 4.0",
    ),
    "clinvar_grch38": Input(
        name="clinvar_grch38",
        url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260715.vcf.gz",
        md5="5151e8d3ca7ffb26ecfc296040bc4f48",
        path=Path("clinvar/clinvar_20260715.vcf.gz"),
        source="NCBI ClinVar weekly release 2026-07-15, GRCh38",
        licence="public domain (NCBI)",
    ),
    "mane_annotation": Input(
        name="mane_annotation",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
            "MANE.GRCh38.v1.5.refseq_genomic.gff.gz"
        ),
        md5="011190055c489f332fb8995fe00822c7",
        path=Path("mane/MANE.GRCh38.v1.5.refseq_genomic.gff.gz"),
        source="NCBI MANE release 1.5, GRCh38",
        licence="public domain (NCBI)",
    ),
    "mane_transcripts": Input(
        name="mane_transcripts",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
            "MANE.GRCh38.v1.5.refseq_rna.fna.gz"
        ),
        md5="f4e69d96b9c9d2698541c41452bc16b9",
        path=Path("mane/MANE.GRCh38.v1.5.refseq_rna.fna.gz"),
        source="NCBI MANE release 1.5, GRCh38",
        licence="public domain (NCBI)",
    ),
    "mane_proteins": Input(
        name="mane_proteins",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
            "MANE.GRCh38.v1.5.refseq_protein.faa.gz"
        ),
        md5="efadce6e19c3c7acfd356cf7915c31bd",
        path=Path("mane/MANE.GRCh38.v1.5.refseq_protein.faa.gz"),
        source="NCBI MANE release 1.5, GRCh38",
        licence="public domain (NCBI)",
    ),
}


def data_root() -> Path:
    """The working directory for fetched inputs, overridable for a machine with data elsewhere."""

    return Path(os.environ.get("RIBORESCUE_DATA", "data"))


def fetch(name: str, root: Path | None = None, *, force: bool = False) -> Path:
    """Fetch a declared input into the data root and return its path.

    A file already present with the right digest is left alone; one with the wrong digest is
    fetched again.
    """

    if (declared := INPUTS.get(name)) is None:
        raise UnknownInputError(f"{name!r} is not a declared input; known: {', '.join(INPUTS)}")

    destination = declared.resolve(root if root is not None else data_root())
    if force:
        destination.unlink(missing_ok=True)
    pooch.retrieve(
        declared.url,
        known_hash=f"md5:{declared.md5}",
        fname=destination.name,
        path=destination.parent,
    )
    return destination
