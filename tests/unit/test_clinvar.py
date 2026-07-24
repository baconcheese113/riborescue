from pathlib import Path

import pytest

from riborescue.core.tables import PathogenicNonsense
from riborescue.variants.clinvar import REVIEW_STARS, pathogenic_nonsense

SAMPLE = Path(__file__).parents[1] / "fixtures" / "clinvar" / "sample.vcf"


@pytest.fixture
def found():
    return pathogenic_nonsense(SAMPLE)


@pytest.fixture
def variants(found):
    return PathogenicNonsense.validate(found.variants)


def test_only_pathogenic_nonsense_substitutions_survive(variants):
    assert list(variants["allele_id"]) == [100001, 100002, 100007]


def test_a_benign_nonsense_variant_is_excluded(variants):
    """Nonsense alone is not enough: the variant must also be disease-causing."""

    assert 100003 not in set(variants["allele_id"])


def test_a_pathogenic_missense_variant_is_excluded(variants):
    assert 100004 not in set(variants["allele_id"])


def test_a_pathogenic_deletion_is_excluded(variants):
    """Readthrough acts on a substituted stop codon, not on a deletion that removes sequence."""

    assert 100005 not in set(variants["allele_id"])


def test_a_conflicting_classification_is_excluded(variants):
    assert 100006 not in set(variants["allele_id"])


def test_review_status_is_carried_through_as_stars(variants):
    stars = dict(zip(variants["allele_id"], variants["review_stars"], strict=True))
    assert stars == {100001: 3, 100002: 1, 100007: 2}


def test_an_unrecognised_review_status_carries_no_stars():
    assert REVIEW_STARS.get("no_assertion_criteria_provided", 0) == 0


def test_the_first_gene_of_an_overlapping_pair_names_the_variant(variants):
    """A variant listed against a gene and its antisense partner belongs to the coding gene."""

    overlapping = variants[variants["allele_id"] == 100007].iloc[0]
    assert overlapping["gene_symbol"] == "CFTR"
    assert overlapping["gene_id"] == 1080


def test_an_ambiguous_alternate_allele_is_excluded_and_counted(found):
    """An IUPAC ambiguity code leaves the stop codon unknown, so the variant is dropped."""

    assert found.ambiguous_alleles == 1
    assert 100008 not in set(found.variants["allele_id"])


def test_the_genomic_coordinate_is_preserved_verbatim(variants):
    cftr = variants[variants["allele_id"] == 100001].iloc[0]
    assert (cftr["chrom"], cftr["pos"], cftr["ref"], cftr["alt"]) == ("7", 117530975, "C", "T")
    assert cftr["variant_id"] == "NC_000007.14:g.117530975C>T"
