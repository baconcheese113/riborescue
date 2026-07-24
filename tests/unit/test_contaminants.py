import gzip
from pathlib import Path

import pytest

from riborescue.riboseq.contaminants import CONTAMINANT_BIOTYPES, write_contaminants

# GENCODE headers are pipe-separated with the transcript biotype eighth.
GENCODE = """\
>ENST01.1|ENSG01.1|OTTHUMG01|OTTHUMT01|RNA5-8S1-201|RNA5-8S1|157|rRNA|
GGGTCGATGAAGAACGCAGC
>ENST02.1|ENSG02.1|OTTHUMG02|OTTHUMT02|ACTB-201|ACTB|1812|protein_coding|
ATGGATGATGATATCGCCGC
>ENST03.1|ENSG03.1|OTTHUMG03|OTTHUMT03|SNORD3A-201|SNORD3A|217|snoRNA|
AAGACTATACTTTCAGGGAT
>ENST04.1|ENSG04.1|OTTHUMG04|OTTHUMT04|MALAT1-201|MALAT1|8779|lncRNA|
AGGCGGAGTTTGAGACCAGC
"""


@pytest.fixture
def transcripts(tmp_path: Path) -> Path:
    path = tmp_path / "gencode.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GENCODE)
    return path


def test_only_structural_rna_is_selected(transcripts: Path, tmp_path: Path):
    out = tmp_path / "contaminants.fa"
    assert write_contaminants(transcripts, out) == 2
    written = out.read_text()
    assert "RNA5-8S1" in written and "SNORD3A" in written
    assert "ACTB" not in written


def test_long_non_coding_rna_is_kept_out_of_the_contaminant_set(transcripts: Path, tmp_path: Path):
    """It is translated often enough that depleting it would remove signal."""

    assert "lncRNA" not in CONTAMINANT_BIOTYPES
    out = tmp_path / "contaminants.fa"
    write_contaminants(transcripts, out)
    assert "MALAT1" not in out.read_text()


def test_the_rdna_unit_is_appended(transcripts: Path, tmp_path: Path):
    """GENCODE annotates almost no rRNA, so the repeating unit comes from its own accession."""

    rdna = tmp_path / "U13369.1.fasta"
    rdna.write_text(">U13369.1 Human ribosomal DNA complete repeating unit\nGGCCGGCGGGT\n")
    out = tmp_path / "contaminants.fa"
    assert write_contaminants(transcripts, out, [rdna]) == 3
    assert "U13369.1" in out.read_text()


def test_a_fasta_without_biotypes_is_refused(tmp_path: Path):
    plain = tmp_path / "plain.fa.gz"
    with gzip.open(plain, "wt") as handle:
        handle.write(">chr1\nACGT\n")
    with pytest.raises(ValueError, match="no structural RNA"):
        write_contaminants(plain, tmp_path / "out.fa")


def test_sequences_survive_the_round_trip_unwrapped(transcripts: Path, tmp_path: Path):
    out = tmp_path / "contaminants.fa"
    write_contaminants(transcripts, out)
    lines = out.read_text().splitlines()
    assert lines[1] == "GGGTCGATGAAGAACGCAGC"
    assert sum(1 for line in lines if line.startswith(">")) == 2
