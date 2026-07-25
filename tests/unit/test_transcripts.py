"""The genomic-to-transcript map, checked against transcripts whose answers are worked by hand.

Both fixtures carry the same sequence, AACCATGAAAGGGCCCTAAT, on opposite strands: two exons of ten
bases, a coding sequence starting at offset 4 and ending in TAA. A strand or off-by-one error moves
an answer that would otherwise still look like a sequence.
"""

from pathlib import Path

import pytest

from riborescue.variants.transcripts import load_transcripts, read_annotation, reverse_complement

MANE = Path(__file__).parents[1] / "fixtures" / "mane"


@pytest.fixture
def models():
    return load_transcripts(MANE / "sample.gff", MANE / "sample.fna")


@pytest.fixture
def plus(models):
    return models["NM_000001.1"]


@pytest.fixture
def minus(models):
    return models["NM_000002.1"]


def test_a_plus_strand_transcript_runs_from_its_lowest_coordinate(plus):
    assert [plus.offset_of(p) for p in (101, 110, 201, 210)] == [0, 9, 10, 19]


def test_a_minus_strand_transcript_runs_from_its_highest_coordinate(minus):
    assert [minus.offset_of(p) for p in (410, 401, 310, 301)] == [0, 9, 10, 19]


@pytest.mark.parametrize("position", [150, 199])
def test_an_intronic_position_has_no_transcript_offset(plus, position):
    assert plus.offset_of(position) is None


@pytest.mark.parametrize("position", [100, 211])
def test_a_position_outside_the_transcript_has_no_offset(plus, position):
    assert plus.offset_of(position) is None


def test_the_spliced_exons_account_for_the_whole_sequence(plus, minus):
    assert plus.spliced_length == len(plus.sequence) == 20
    assert minus.spliced_length == len(minus.sequence) == 20


def test_coding_starts_at_the_first_coding_base_on_either_strand(plus, minus):
    assert plus.coding_offset == minus.coding_offset == 4
    assert plus.sequence[plus.coding_offset : plus.coding_offset + 3] == "ATG"
    assert minus.sequence[minus.coding_offset : minus.coding_offset + 3] == "ATG"


def test_a_codon_is_read_in_frame_from_the_coding_start(plus):
    assert plus.codon_covering(105) == (4, "ATG")
    assert plus.codon_covering(108) == (7, "AAA")
    assert plus.codon_covering(203) == (10, "GGG")


def test_a_position_before_the_coding_start_has_no_codon(plus):
    assert plus.codon_covering(102) is None


def test_the_transcript_base_is_read_on_the_transcripts_own_strand(plus, minus):
    assert plus.base_at(105) == "A"
    assert minus.base_at(406) == "A"


def test_each_transcript_carries_its_gene_and_placement(plus, minus):
    assert (plus.gene_id, plus.chrom, plus.strand) == (9001, "chr9", "+")
    assert (minus.gene_id, minus.chrom, minus.strand) == (9002, "chr9", "-")


def test_only_exon_and_coding_features_are_read():
    annotation = read_annotation(MANE / "sample.gff")
    assert set(annotation["feature"]) == {"exon", "CDS"}
    assert set(annotation["transcript_id"]) == {"NM_000001.1", "NM_000002.1"}


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [("ATGC", "GCAT"), ("AAAA", "TTTT"), ("", ""), ("ACGTN", "NACGT")],
)
def test_reverse_complement_reverses_and_complements(sequence, expected):
    assert reverse_complement(sequence) == expected


def test_a_mane_plus_clinical_transcript_is_left_out(models):
    """The policy is one Select transcript per gene, so a second transcript never competes."""

    assert "NM_000003.1" not in models
    assert [model.gene_id for model in models.values()].count(9001) == 1


def test_a_transcript_on_an_alternate_contig_is_left_out(models):
    assert "NM_000004.1" not in models


def test_each_transcript_names_the_protein_it_encodes(plus, minus):
    assert (plus.protein_id, minus.protein_id) == ("NP_000001.1", "NP_000002.1")


def test_exon_lengths_follow_transcript_order_not_genomic_order(plus, minus):
    assert plus.exon_lengths == (10, 10)
    assert minus.exon_lengths == (10, 10)


def test_a_junction_sits_where_one_exon_ends_and_the_next_begins(plus, minus):
    assert plus.junction_offsets == (10,)
    assert minus.junction_offsets == (10,)


def test_a_single_exon_transcript_has_no_junctions():
    from riborescue.variants.transcripts import TranscriptModel

    single = TranscriptModel(
        transcript_id="NM_single.1",
        gene_id=1,
        chrom="chr1",
        strand="+",
        exons=((100, 199),),
        cds_start=100,
        cds_end=199,
        sequence="A" * 100,
        protein_id="NP_single.1",
    )
    assert single.junction_offsets == ()


def _utr_model():
    """ATG · AAA · TAA(stop) · AAAAAA(3' UTR): coding ends at offset 8, six UTR bases follow."""

    from riborescue.variants.transcripts import TranscriptModel

    return TranscriptModel(
        transcript_id="NM_utr.1",
        gene_id=2,
        chrom="chr1",
        strand="+",
        exons=((100, 114),),
        cds_start=100,
        cds_end=108,
        sequence="ATGAAATAAAAAAAA",
        protein_id="NP_utr.1",
    )


def test_coding_end_offset_is_the_last_base_of_the_natural_stop():
    model = _utr_model()
    assert model.coding_end_offset == 8
    assert model.sequence[6:9] == "TAA"  # the natural stop, ending at the coding-end offset


def test_a_codon_in_the_three_prime_utr_past_the_natural_stop_has_no_codon():
    # Genomic 110 is offset 10, inside the 3' UTR: translation has already terminated, so there is
    # no reading-frame codon there however the bases fall — a substitution is not a premature stop.
    model = _utr_model()
    assert model.codon_covering(106) == (6, "TAA")  # the natural stop is still covered
    assert model.codon_covering(110) is None  # one codon into the 3' UTR
    assert model.codon_covering(113) is None
