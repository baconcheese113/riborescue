"""Support counted per term, over a library small enough to count by hand.

The point of the module is that support is a property of terms rather than of contexts, so the
fixtures are built so that a variant is well supported on every term while its complete context was
never measured — the case a single "was this observed?" flag gets wrong.
"""

import pandas as pd
import pytest

from riborescue.variants.support_atlas import context_analogues, support_atlas, term_support

_LIBRARY = pd.DataFrame(
    {
        "gene": ["AAA1", "AAA1", "BBB2", "CCC3", "DDD4", "EEE5"],
        "stop_type": ["uga", "uga", "uga", "uaa", "uaa", "uga"],
        "up_123nt": ["uuc", "gaa", "uuc", "gaa", "uuc", "gaa"],
        "down_123nt": ["gaa", "gaa", "ccc", "gaa", "ccc", "ggg"],
        "RT_binomial": [0.05, 0.04, 0.03, 0.02, 0.01, 0.06],
    }
)
"""Six measured variants, arranged so that every kind of support gap has an example.

`uga`+`uuc`+`gaa` was measured exactly. `uga`+`gaa`+`ccc` never was, though each of its three levels
is observed on its own. And `ggg` downstream occurs only with `uga`, so the `uaa`+`ggg` cell is
empty while both of its levels are supported."""


def _scored(**overrides) -> pd.DataFrame:
    row = {
        "variant_id": "V1",
        "gene_symbol": "AAA1",
        "stop_type": "TGA",
        "up_123nt": "GAA",
        "down_123nt": "CCC",
    }
    return pd.DataFrame([row | overrides])


class TestTermSupport:
    def test_every_term_is_counted_separately(self):
        counts = term_support(_LIBRARY)
        assert set(counts) == {
            "stop",
            "upstream",
            "downstream",
            "interaction_cell",
            "complete_context",
        }
        assert counts["stop"]["TGA"] == 4
        assert counts["upstream"]["TTC"] == 3
        assert counts["interaction_cell"][("TGA", "GAA")] == 2

    def test_rna_and_dna_spellings_are_the_same_level(self):
        counts = term_support(_LIBRARY)
        # The library spells triplets as RNA; the scored table spells them as DNA.
        assert "TTC" in counts["upstream"].index
        assert "UUC" not in counts["upstream"].index


class TestSupportAtlas:
    def test_a_variant_can_be_supported_on_every_term_and_never_measured(self):
        atlas = support_atlas(_scored(), _LIBRARY)
        row = atlas.iloc[0]
        # TGA, GAA upstream, CCC downstream: each level is observed, the combination is not.
        assert row["support_stop"] == 4
        assert row["support_upstream"] == 3
        assert row["support_downstream"] == 2
        assert row["support_complete_context"] == 0
        assert not row["measured_exactly"]
        assert not row["aliased_interaction"]

    def test_a_cell_the_library_never_observed_is_flagged_as_aliased(self):
        atlas = support_atlas(_scored(stop_type="TAA", down_123nt="GGG"), _LIBRARY)
        # TAA is observed twice and GGG once, but never together, so the interacted cell is empty
        # and the prediction falls back to the marginal.
        assert atlas.iloc[0]["support_interaction_cell"] == 0
        assert atlas.iloc[0]["aliased_interaction"]

    def test_a_level_the_library_never_held_leaves_the_term_at_zero(self):
        # A stop codon immediately downstream: the reporter library contains none by design,
        # because such a variant terminates again the moment it reads through.
        atlas = support_atlas(_scored(down_123nt="TAA"), _LIBRARY)
        assert atlas.iloc[0]["support_downstream"] == 0
        assert atlas.iloc[0]["weakest_term"] == 0

    def test_the_weakest_term_is_what_the_prediction_rests_on(self):
        atlas = support_atlas(_scored(), _LIBRARY)
        row = atlas.iloc[0]
        assert row["weakest_term"] == min(
            row["support_stop"],
            row["support_upstream"],
            row["support_downstream"],
            row["support_interaction_cell"],
        )

    def test_the_complete_context_is_counted_when_it_was_measured(self):
        atlas = support_atlas(_scored(up_123nt="TTC", down_123nt="GAA"), _LIBRARY)
        assert atlas.iloc[0]["support_complete_context"] == 1
        assert atlas.iloc[0]["measured_exactly"]


class TestContextAnalogues:
    def test_only_variants_sharing_the_complete_context_are_returned(self):
        paired = context_analogues(
            _scored(up_123nt="TTC", down_123nt="GAA"), _LIBRARY, label="Toledano G418"
        )
        assert len(paired) == 1
        assert paired.iloc[0]["analogue_gene"] == "AAA1"
        assert paired.iloc[0]["analogue_readthrough"] == pytest.approx(0.05)

    def test_the_source_of_every_measurement_is_recorded(self):
        paired = context_analogues(
            _scored(up_123nt="TTC", down_123nt="GAA"), _LIBRARY, label="Toledano G418"
        )
        # The label is what keeps a reporter reading from being mistaken for clinical evidence.
        assert set(paired["source"]) == {"Toledano G418"}

    def test_a_variant_with_no_measured_context_returns_nothing_rather_than_a_near_miss(self):
        assert context_analogues(_scored(), _LIBRARY, label="Toledano G418").empty

    def test_a_measurement_in_the_same_gene_is_listed_first(self):
        library = pd.concat(
            [_LIBRARY, _LIBRARY.iloc[[0]].assign(gene="ZZZ9", RT_binomial=0.09)], ignore_index=True
        )
        paired = context_analogues(
            _scored(up_123nt="TTC", down_123nt="GAA"), library, label="Toledano G418"
        )
        # ZZZ9 has the larger measurement, but AAA1 is the scored variant's own gene.
        assert list(paired["analogue_gene"]) == ["AAA1", "ZZZ9"]
        assert list(paired["same_gene"]) == [True, False]

    def test_the_list_is_capped_and_the_uncapped_count_lives_in_the_atlas(self):
        many = pd.concat([_LIBRARY.iloc[[0]]] * 25, ignore_index=True)
        scored = _scored(up_123nt="TTC", down_123nt="GAA")
        paired = context_analogues(scored, many, label="Toledano G418", limit=10)
        assert len(paired) == 10
        assert support_atlas(scored, many).iloc[0]["support_complete_context"] == 25
