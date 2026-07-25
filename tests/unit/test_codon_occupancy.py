"""The codon table over a transcript short enough to count by hand.

Every number below is arithmetic a reader can redo. The transcript is built so that one codon is
hit and another is not, because the position that carries no P-site is the one an implementation
most easily drops, and dropping it makes a codon ribosomes never pause on look ordinary.
"""

import pandas as pd
import pytest

from riborescue.riboseq.codon_occupancy import (
    GENETIC_CODE,
    SENSE_CODONS,
    STOP_CODONS,
    SYNONYMOUS,
    aggregate_libraries,
    as_dna,
    library_table,
    scorable_bounds,
    transcript_codons,
)

# Six codons of padding, then exactly the 64 scorable ones the floor demands, then the tail the
# window excludes. The window runs from codon 6 to the last codon whose first base is at least 16
# nucleotides clear of the stop, which for this coding length is codon 69.
_UTR5 = "A" * 9
_HEAD = "ATG" + "GGG" * 5  # codons 0-5, before the window
_BODY = "AAA" * 32 + "CCC" * 32  # codons 6-37 then 38-69
_TAIL = "GGG" * 4 + "TAA"  # codons 70-74, after the window
_CDS = _HEAD + _BODY + _TAIL
_FIRST_CCC = 38
_LAST = 70  # one past the last scorable codon

_SEQUENCES = {"ENST1": _UTR5 + _CDS + "TTTT"}
_ANNOTATION = pd.DataFrame(
    {"transcript": ["ENST1"], "l_utr5": [len(_UTR5)], "l_cds": [len(_CDS)], "l_utr3": [4]}
)


@pytest.fixture
def codons() -> dict[str, str]:
    return transcript_codons(_SEQUENCES, _ANNOTATION)


def _counts(rows: list[tuple[int, int]], sample: str = "lib") -> pd.DataFrame:
    """P-site codon index and count, at one footprint length unless a row says otherwise."""

    return pd.DataFrame(
        {
            "transcript": "ENST1",
            "length": 30,
            "codon_index": [index for index, _ in rows],
            "sample": sample,
            "n": [n for _, n in rows],
        }
    )


class TestGeneticCode:
    def test_the_code_is_complete_and_has_three_stops(self):
        assert len(GENETIC_CODE) == 64
        assert set(STOP_CODONS) == {"TAA", "TAG", "TGA"}
        assert len(SENSE_CODONS) == 61

    def test_synonymous_families_partition_the_sense_codons(self):
        assigned = [codon for family in SYNONYMOUS.values() for codon in family]
        assert sorted(assigned) == sorted(SENSE_CODONS)
        assert len(assigned) == len(set(assigned))

    def test_the_two_single_codon_families_are_the_ones_the_shuffle_cannot_move(self):
        single = {residue for residue, family in SYNONYMOUS.items() if len(family) == 1}
        assert single == {"M", "W"}

    def test_a_triplet_reads_the_same_spelled_as_rna(self):
        assert as_dna("uuc") == "TTC"
        assert as_dna("TTC") == "TTC"


class TestWindow:
    def test_the_window_starts_at_codon_six_whatever_the_transcript(self):
        assert scorable_bounds(300)[0] == 6
        assert scorable_bounds(9_000)[0] == 6

    def test_the_last_scorable_codon_keeps_sixteen_bases_clear_of_the_stop(self):
        _, last = scorable_bounds(len(_CDS))
        # Codon `last` begins 3*last into the CDS; its distance from the stop's last base is
        # 3*last - l_cds + 1, which must be at most -16.
        assert 3 * last - len(_CDS) + 1 <= -16
        assert 3 * (last + 1) - len(_CDS) + 1 > -16

    def test_the_codon_string_holds_exactly_the_scorable_window(self, codons):
        first, last = scorable_bounds(len(_CDS))
        assert len(codons["ENST1"]) == 3 * (last - first + 1)
        assert codons["ENST1"].startswith("AAA")

    def test_a_transcript_with_no_room_for_a_window_is_absent(self):
        annotation = pd.DataFrame(
            {"transcript": ["SHORT"], "l_utr5": [0], "l_cds": [30], "l_utr3": [0]}
        )
        assert transcript_codons({"SHORT": "A" * 30}, annotation) == {}


class TestLibraryTable:
    def test_a_position_with_no_psites_is_a_zero_not_an_absence(self, codons):
        # Every P-site lands on the AAA codons; the CCC codons are scorable and empty.
        counts = _counts([(index, 25) for index in range(6, _FIRST_CCC)])
        table = library_table(counts, _ANNOTATION, codons, site="p").set_index("codon")
        assert table.loc["CCC", "positions"] == 32
        assert table.loc["CCC", "occupancy"] == 0.0

    def test_occupancy_is_relative_to_the_transcripts_own_mean(self, codons):
        # 800 P-sites over 64 scorable positions is a mean of 12.5 per position. The 32 AAA
        # positions carry 25 each, so each normalises to 25 / 12.5 = 2.
        counts = _counts([(index, 25) for index in range(6, _FIRST_CCC)])
        table = library_table(counts, _ANNOTATION, codons, site="p").set_index("codon")
        assert table.loc["AAA", "occupancy"] == pytest.approx(2.0)

    def test_a_uniform_transcript_gives_every_occupied_codon_one(self, codons):
        counts = _counts([(index, 20) for index in range(6, _LAST)])
        table = library_table(counts, _ANNOTATION, codons, site="p").set_index("codon")
        assert table.loc["AAA", "occupancy"] == pytest.approx(1.0)
        assert table.loc["CCC", "occupancy"] == pytest.approx(1.0)

    def test_the_a_site_is_one_codon_past_the_p_site(self, codons):
        # P-sites recorded at codons 5-36 score the A site at codons 6-37, which are the AAA run.
        counts = _counts([(index, 25) for index in range(5, _FIRST_CCC - 1)])
        table = library_table(counts, _ANNOTATION, codons, site="a").set_index("codon")
        assert table.loc["AAA", "occupancy"] == pytest.approx(2.0)
        assert table.loc["CCC", "occupancy"] == 0.0

    def test_psites_outside_the_window_take_no_part(self, codons):
        inside = _counts([(index, 25) for index in range(6, _FIRST_CCC)])
        outside = pd.concat([inside, _counts([(2, 400), (_LAST + 4, 400)])], ignore_index=True)
        assert library_table(inside, _ANNOTATION, codons, site="p").equals(
            library_table(outside, _ANNOTATION, codons, site="p")
        )

    def test_only_the_named_lengths_are_summed(self, codons):
        counts = pd.concat(
            [
                _counts([(index, 25) for index in range(6, _FIRST_CCC)]),
                _counts([(index, 25) for index in range(_FIRST_CCC, _LAST)]).assign(length=21),
            ],
            ignore_index=True,
        )
        kept = library_table(counts, _ANNOTATION, codons, site="p", lengths=(30,))
        assert kept.set_index("codon").loc["CCC", "occupancy"] == 0.0
        both = library_table(counts, _ANNOTATION, codons, site="p", lengths=(21, 30))
        assert both.set_index("codon").loc["CCC", "occupancy"] == pytest.approx(1.0)

    def test_a_transcript_under_the_psite_floor_carries_nothing(self, codons):
        with pytest.raises(ValueError, match="no transcript clears the coverage floors"):
            library_table(_counts([(6, 99)]), _ANNOTATION, codons, site="p")

    def test_every_sense_codon_has_a_row_even_where_the_transcriptome_lacks_it(self, codons):
        counts = _counts([(index, 25) for index in range(6, _FIRST_CCC)])
        table = library_table(counts, _ANNOTATION, codons)
        assert len(table) == 61
        assert table.loc[table["codon"] == "TGG", "positions"].item() == 0

    def test_an_unknown_site_is_refused_rather_than_guessed(self, codons):
        with pytest.raises(ValueError, match="site must be one of"):
            library_table(_counts([(6, 200)]), _ANNOTATION, codons, site="e")


class TestAggregate:
    def test_libraries_are_averaged_not_pooled_so_depth_does_not_weight_the_table(self, codons):
        shallow = library_table(
            _counts([(index, 25) for index in range(6, _FIRST_CCC)]), _ANNOTATION, codons, site="p"
        )
        deep = library_table(
            _counts([(index, 20) for index in range(6, _LAST)]), _ANNOTATION, codons, site="p"
        )
        table = aggregate_libraries({"shallow": shallow, "deep": deep}).set_index("codon")
        # 2.0 and 1.0 averaged, not 800 P-sites against 1,280 pooled.
        assert table.loc["AAA", "occupancy"] == pytest.approx(1.5)
        assert table.loc["AAA", "libraries"] == 2

    def test_the_spread_between_libraries_is_reported_beside_the_score(self, codons):
        shallow = library_table(
            _counts([(index, 25) for index in range(6, _FIRST_CCC)]), _ANNOTATION, codons, site="p"
        )
        deep = library_table(
            _counts([(index, 20) for index in range(6, _LAST)]), _ANNOTATION, codons, site="p"
        )
        table = aggregate_libraries({"shallow": shallow, "deep": deep}).set_index("codon")
        assert table.loc["AAA", "occupancy_sd"] == pytest.approx(2.0**0.5 / 2, rel=1e-6)
