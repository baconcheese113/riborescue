import pandas as pd
import pytest

from riborescue.residue import (
    AMINO_ACIDS,
    NEAR_COGNATE,
    SuppressorDesign,
    conservative,
    coverage_by_design,
    near_cognate_residues,
    similarity,
)


def test_a_suppressor_trna_can_carry_any_of_the_twenty_amino_acids():
    assert len(AMINO_ACIDS) == 20
    assert "*" not in AMINO_ACIDS


@pytest.mark.parametrize("stop", ["TAA", "TAG", "TGA"])
def test_a_near_cognate_set_never_includes_a_stop(stop: str):
    assert all(residue != "*" for residue in NEAR_COGNATE[stop])


def test_the_near_cognate_residues_are_those_one_base_away():
    """TGG encodes tryptophan and differs from TGA at the last base, so UGA can read to W."""

    assert "W" in near_cognate_residues("TGA")
    assert "R" in near_cognate_residues("TGA")  # CGA and AGA
    assert "Q" in near_cognate_residues("TAG")  # CAG
    assert "W" not in near_cognate_residues("TAA")  # TGG differs at two bases


def test_a_stop_codon_written_as_rna_reads_the_same():
    assert near_cognate_residues("uga") == near_cognate_residues("TGA")


def test_a_residue_substitutes_for_itself_most_readily():
    for residue in AMINO_ACIDS:
        assert similarity(residue, residue) >= 4
        assert conservative(residue, residue)


def test_similarity_does_not_depend_on_direction():
    assert similarity("W", "Y") == similarity("Y", "W")


def test_an_aromatic_swap_is_conservative_and_a_charged_one_is_not():
    """The IDUA case: a tyrosine suppressor outperformed others at a tryptophan stop."""

    assert conservative("W", "Y")
    assert not conservative("W", "R")


@pytest.fixture
def contexts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "stop_type": ["uag", "uag", "uag", "uga"],
            "original_aa": ["Q", "Q", "W", "R"],
        }
    )


def test_coverage_counts_every_design_against_every_stop(contexts: pd.DataFrame):
    coverage = coverage_by_design(contexts)
    assert len(coverage) == 2 * len(AMINO_ACIDS)
    assert set(coverage["recognized_stop"]) == {"UAG", "UGA"}


def test_a_design_restoring_the_original_residue_counts_it_exactly(contexts: pd.DataFrame):
    coverage = coverage_by_design(contexts).set_index("design_id")
    assert coverage.loc["UAG-Q", "restores_exactly"] == 2
    assert coverage.loc["UAG-W", "restores_exactly"] == 1
    assert coverage.loc["UGA-R", "restores_exactly"] == 1


def test_conservative_coverage_is_never_smaller_than_exact_restoration(contexts: pd.DataFrame):
    coverage = coverage_by_design(contexts)
    assert (coverage["conservative"] >= coverage["restores_exactly"]).all()


def test_a_design_is_named_by_the_stop_it_reads_and_the_residue_it_inserts():
    assert SuppressorDesign(recognized_stop="UGA", inserted_aa="W").design_id == "UGA-W"
