"""Base-editing reachability on hand-built transcripts, where every PAM is placed on purpose.

The panel is BE4max (C→T, window 4-8) and ABE7.10 (A→G, window 4-7), SpCas9 NGG, both strands.
Sequences are single-exon unless a test needs a junction, coding from offset 0, so a genomic
position is its transcript offset plus one.
"""

from pathlib import Path

import pandas as pd
import pytest

from riborescue.variants.base_editing import (
    PRIMARY_PANEL,
    SENSITIVITY_PANEL,
    ReachClass,
    reachability_for,
    reachability_table,
)
from riborescue.variants.transcripts import TranscriptModel


def _seq(overrides: dict[int, str], length: int = 45, fill: str = "C") -> str:
    bases = [fill] * length
    for index, base in overrides.items():
        bases[index] = base
    return "".join(bases)


def _model(sequence: str, *, exons=None, strand="+") -> TranscriptModel:
    exons = exons or ((1, len(sequence)),)
    # On the minus strand coding runs from the high genomic coordinate down, so cds_start/end swap.
    cds_start, cds_end = (exons[0][0], exons[-1][1])
    return TranscriptModel(
        transcript_id="NM_TEST.1",
        gene_id=42,
        chrom="chr1",
        strand=strand,
        exons=exons,
        cds_start=cds_start,
        cds_end=cds_end,
        sequence=sequence,
        protein_id="NP_TEST.1",
    )


# A Trp codon TGG at offsets 15-17, with GG at 34-35 so an NGG PAM sits where a protospacer starting
# at 13 places offset 17 at window position 5. Coding runs to a natural TAA at 42-44.
_TGG = {
    0: "A",
    1: "T",
    2: "G",
    15: "T",
    16: "G",
    17: "G",
    34: "G",
    35: "G",
    42: "T",
    43: "A",
    44: "A",
}


def test_an_adenine_editor_restores_a_tryptophan_stop_exactly():
    model = _model(_seq(_TGG))
    # Genomic 18 is offset 17; G>A turns the wildtype TGG into TGA.
    reach = reachability_for(model, 18, "G", "A")
    assert reach is not None
    assert reach.stop_type == "TGA"
    assert reach.reference_codon == "TGG"
    assert reach.reach_class is ReachClass.exact
    exact = [
        g
        for g in reach.guides
        if g.editor == "ABE7.10" and g.strand == "sense" and g.restores == "exact_wildtype"
    ]
    assert exact, "the sense ABE7.10 placement that reverts TGA→TGG should be found"
    assert exact[0].window_position == 5
    assert exact[0].pam.endswith("GG")
    assert exact[0].bystander_free


def test_a_cytosine_editor_cannot_correct_a_nonsense_stop():
    """No stop codon contains a C, and editing either G to A leaves another stop — a fact of the
    code, not of this sequence. BE4max should never produce a stop-removing guide."""

    model = _model(_seq(_TGG))
    reach = reachability_for(model, 18, "G", "A")
    assert reach is not None
    assert reach.guides
    assert not any(g.editor == "BE4max" for g in reach.guides)


def test_a_placement_whose_guide_crosses_an_exon_junction_is_refused():
    """Same sequence, split so a junction at offset 25 falls inside the protospacer 13-35. A guide
    is genomic and is not designed across it, so the sense ABE placement disappears."""

    spliced = _model(_seq(_TGG), exons=((1, 25), (100, 119)))
    assert spliced.junction_offsets == (25,)
    reach = reachability_for(spliced, 18, "G", "A")
    assert reach is not None
    assert not any(
        g.editor == "ABE7.10" and g.strand == "sense" and g.window_position == 5
        for g in reach.guides
    )


def test_a_coding_bystander_is_reported_with_its_amino_acid_change():
    # Put an A at offset 19, inside the editing window 16-19, in codon CAC (His) at 18-20.
    overrides = _TGG | {18: "C", 19: "A", 20: "C"}
    reach = reachability_for(_model(_seq(overrides)), 18, "G", "A")
    assert reach is not None
    guide = next(
        g
        for g in reach.guides
        if g.editor == "ABE7.10" and g.strand == "sense" and g.window_position == 5
    )
    bystander = next(b for b in guide.bystanders if b.offset == 19)
    assert bystander.coding
    assert (bystander.aa_from, bystander.aa_to) == ("H", "R")
    assert not bystander.silent
    assert not guide.bystander_free


def test_a_substitution_that_makes_no_stop_has_no_reachability():
    assert reachability_for(_model(_seq(_TGG)), 18, "G", "C") is None


def test_the_summary_table_carries_one_row_per_variant():
    variants = pd.DataFrame(
        [
            {
                "variant_id": "v1",
                "gene_id": 42,
                "gene_symbol": "TEST",
                "pos": 18,
                "ref": "G",
                "alt": "A",
            },
            {
                "variant_id": "v2",
                "gene_id": 99,
                "gene_symbol": "OTHER",
                "pos": 18,
                "ref": "G",
                "alt": "A",
            },
        ]
    )
    table = reachability_table(variants, {42: _model(_seq(_TGG))})
    assert list(table["variant_id"]) == ["v1", "v2"]
    reached = table.set_index("variant_id")
    assert reached.loc["v1", "reachable"]
    assert reached.loc["v1", "reach_class"] == ReachClass.exact.value
    assert reached.loc["v1", "editor"] == "ABE7.10"
    # A variant whose gene has no model is not scoreable, not a crash.
    assert not reached.loc["v2", "scoreable"]


def test_the_primary_panel_is_the_two_named_editors():
    assert {e.name for e in PRIMARY_PANEL} == {"BE4max", "ABE7.10"}
    assert all(e.pam == "NGG" and e.arm == "primary" for e in PRIMARY_PANEL)


def test_a_minus_strand_transcript_is_scored_on_its_own_reading_frame():
    """The same reachability through the minus-strand path: offsets count down from the high
    coordinate and the forward-strand alleles are complemented before the search."""

    minus = _model(_seq(_TGG), strand="-")
    # Offset 17 is genomic 45-17=28 on the minus strand; forward-strand C>T there is transcript G>A.
    reach = reachability_for(minus, 28, "C", "T")
    assert reach is not None
    assert reach.stop_type == "TGA"
    assert reach.reach_class is ReachClass.exact
    assert any(g.editor == "ABE7.10" and g.restores == "exact_wildtype" for g in reach.guides)


def test_a_relaxed_pam_reaches_a_stop_the_ngg_panel_cannot():
    """With a GT where the NGG panel needs GG, no NGG guide places the target in the window, but the
    SpCas9-NG sensitivity arm — which reads only the G — does."""

    model = _model(_seq(_TGG | {35: "T"}))  # PAM at 33-35 becomes C,G,T: NG but not NGG

    def sense_window5(reach):
        return [g for g in reach.guides if g.strand == "sense" and g.window_position == 5]

    primary = reachability_for(model, 18, "G", "A", panel=PRIMARY_PANEL)
    assert primary is not None
    assert not sense_window5(primary), "NGG should not place the target in the window here"

    sensitive = reachability_for(model, 18, "G", "A", panel=SENSITIVITY_PANEL)
    assert sensitive is not None
    reached = sense_window5(sensitive)
    assert reached and all(g.arm == "sensitivity" for g in reached)


def test_a_failure_names_no_pam_when_no_motif_is_in_reach():
    # An all-A context: no GG on either strand near the stop, so no PAM reaches any target.
    model = _model(
        _seq(
            {0: "A", 1: "T", 2: "G", 15: "T", 16: "G", 17: "G", 42: "T", 43: "A", 44: "A"}, fill="A"
        )
    )
    reach = reachability_for(model, 18, "G", "A")
    assert reach is not None
    assert not reach.guides
    assert reach.reach_class is ReachClass.none
    assert reach.reason == "no_pam"


def test_a_failure_names_the_window_when_a_pam_sits_out_of_reach():
    # A PAM at 24-25 places the target at protospacer position 15 — inside the guide, outside the
    # editing window — so the failure is the window, not the PAM.
    overrides = {
        0: "A",
        1: "T",
        2: "G",
        15: "T",
        16: "G",
        17: "G",
        24: "G",
        25: "G",
        42: "T",
        43: "A",
        44: "A",
    }
    reach = reachability_for(_model(_seq(overrides, fill="A")), 18, "G", "A")
    assert reach is not None
    assert not reach.guides
    assert reach.reason == "target_outside_window"


def test_an_unscoreable_variant_carries_that_reason_in_the_table():
    variants = pd.DataFrame(
        [
            {
                "variant_id": "v2",
                "gene_id": 99,
                "gene_symbol": "OTHER",
                "pos": 18,
                "ref": "G",
                "alt": "A",
            }
        ]
    )
    table = reachability_table(variants, {42: _model(_seq(_TGG))}).set_index("variant_id")
    assert not table.loc["v2", "scoreable"]
    assert table.loc["v2", "reason"] == "unscoreable_context"


def test_the_sensitivity_arm_is_named_enzymes_with_declared_pam_rules():
    names = {e.name for e in SENSITIVITY_PANEL}
    assert {"BE4max+SpCas9-NG", "ABE7.10+SpRY"} <= names
    assert all(e.arm == "sensitivity" for e in SENSITIVITY_PANEL)


# Real GRCh38 / MANE Select transcripts, three variants read out of the ClinVar run and confirmed by
# hand against the published reachability call. Skipped where the fetched annotation is absent (CI);
# runs wherever the pipeline's own inputs are present.
_MANE_GFF = Path("data/mane/MANE.GRCh38.v1.5.refseq_genomic.gff.gz")
_MANE_FNA = Path("data/mane/MANE.GRCh38.v1.5.refseq_rna.fna.gz")


@pytest.mark.skipif(
    not (_MANE_GFF.exists() and _MANE_FNA.exists()),
    reason="MANE annotation is fetched, not committed; the real-variant panel runs where it exists",
)
def test_hand_checked_real_variants_reproduce_their_reachability_call():
    from riborescue.variants.transcripts import load_transcripts

    by_gene = {m.gene_id: m for m in load_transcripts(_MANE_GFF, _MANE_FNA).values()}
    # gene_id, pos, ref, alt, expected class, editor, strand
    panel = [
        (148398, 943995, "C", "T", ReachClass.exact, "ABE7.10", "antisense"),  # SAMD11
        (9636, 1014359, "G", "T", ReachClass.alternative, "ABE7.10", "sense"),  # ISG15
        (375790, 1022368, "C", "A", ReachClass.none, None, None),  # AGRN
    ]
    for gene_id, pos, ref, alt, cls, editor, strand in panel:
        reach = reachability_for(by_gene[gene_id], pos, ref, alt)
        assert reach is not None, f"{gene_id}:{pos} should place on a stop"
        assert reach.reach_class is cls
        if editor is not None:
            assert any(g.editor == editor and g.strand == strand for g in reach.guides)
