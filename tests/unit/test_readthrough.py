import gzip
from pathlib import Path

import pandas as pd
import pytest

from riborescue.readthrough import (
    MINIMUM_CDS_PSITES,
    MINIMUM_UTR3_LENGTH,
    PROGRAMMED_READTHROUGH,
    frame_specific,
    library_ratios,
    overlapping_downstream_cds,
    paired_effect,
    qualifying,
    signature,
    transcript_genes,
)

REPLICATES = pd.DataFrame(
    {
        "sample": [
            "hek293t_untreated_rep1_riboseq",
            "hek293t_untreated_rep2_riboseq",
            "hek293t_untreated_rep3_riboseq",
            "hek293t_g418_rep1_riboseq",
            "hek293t_g418_rep2_riboseq",
            "hek293t_g418_rep3_riboseq",
        ],
        "replicate": [1, 2, 3, 1, 2, 3],
        "treatment": ["untreated"] * 3 + ["g418"] * 3,
    }
)


def counts(**overrides) -> pd.DataFrame:
    base = {
        "transcript": "ENST1",
        "sample": "s1",
        "cds_inframe": 1000,
        "cds_total": 1500,
        "termination": 50,
        "extension": 300,
        "extension_inframe": 10,
        "extension_outframe": 5,
        "l_cds": 900,
        "l_utr3": 500,
    }
    return pd.DataFrame([base | overrides])


def test_a_short_untranslated_region_is_excluded():
    """The termination peak spills into a window this short."""

    assert qualifying(counts(l_utr3=MINIMUM_UTR3_LENGTH - 1)).empty
    assert len(qualifying(counts(l_utr3=MINIMUM_UTR3_LENGTH))) == 1


def test_a_transcript_without_enough_coding_coverage_is_excluded():
    assert qualifying(counts(cds_inframe=MINIMUM_CDS_PSITES - 1)).empty
    assert len(qualifying(counts(cds_inframe=MINIMUM_CDS_PSITES))) == 1


def test_a_programmed_readthrough_gene_is_excluded():
    """It reads through with or without a drug, so it cannot evidence one."""

    genes = pd.Series({"ENST1": "AQP4"})
    assert "AQP4" in PROGRAMMED_READTHROUGH
    assert qualifying(counts(), genes=genes).empty
    assert len(qualifying(counts(), genes=pd.Series({"ENST1": "ACTB"}))) == 1


def test_a_transcript_running_into_another_gene_is_excluded():
    assert qualifying(counts(), excluded_transcripts=frozenset({"ENST1"})).empty


def test_exclusions_do_not_depend_on_treatment():
    """The same transcripts must drop from both conditions, or the comparison is not paired."""

    both = pd.concat([counts(sample="untreated"), counts(sample="g418", extension_inframe=99)])
    kept = qualifying(both)
    assert set(kept["sample"]) == {"untreated", "g418"}


def test_library_ratios_take_the_median_transcript():
    """One enormous transcript must not decide the library."""

    frame = pd.concat(
        [
            counts(transcript="a", extension_inframe=10, cds_inframe=1000),
            counts(transcript="b", extension_inframe=10, cds_inframe=1000),
            counts(transcript="c", extension_inframe=9000, cds_inframe=1000),
        ]
    )
    ratios = library_ratios(frame)
    assert ratios.loc[0, "readthrough"] == pytest.approx(0.01)
    assert ratios.loc[0, "transcripts"] == 3


def test_paired_effect_compares_within_replicate():
    ratios = pd.DataFrame(
        {
            "sample": REPLICATES["sample"],
            "readthrough": [0.010, 0.020, 0.030, 0.015, 0.026, 0.037],
        }
    )
    effect = paired_effect(ratios, "readthrough", REPLICATES)
    assert effect.differences == pytest.approx((0.005, 0.006, 0.007))
    assert effect.consistent is True
    assert effect.mean_difference == pytest.approx(0.006)


def test_a_replicate_missing_a_condition_is_refused():
    ratios = pd.DataFrame({"sample": REPLICATES["sample"][:5], "readthrough": [0.1] * 5})
    with pytest.raises(ValueError, match="without both conditions"):
        paired_effect(ratios, "readthrough", REPLICATES.iloc[:5])


def test_an_inconsistent_effect_is_reported_as_such():
    ratios = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.02, 0.03, 0.02, 0.01, 0.04]}
    )
    effect = paired_effect(ratios, "readthrough", REPLICATES)
    assert effect.consistent is False


def test_the_interval_widens_with_disagreement():
    tight = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.01, 0.01, 0.02, 0.02, 0.02]}
    )
    loose = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.01, 0.01, 0.05, 0.02, 0.001]}
    )
    narrow = paired_effect(tight, "readthrough", REPLICATES).interval
    wide = paired_effect(loose, "readthrough", REPLICATES).interval
    assert (wide[1] - wide[0]) > (narrow[1] - narrow[0])


GTF = """\
chr1\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\tCDS\t100\t500\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\tthree_prime_utr\t501\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\ttranscript\t600\t1200\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "NEIGH";
chr1\tX\tCDS\t600\t1200\t.\t+\t0\tgene_id "G2"; transcript_id "T2"; gene_name "NEIGH";
chr2\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_name "CLEAN";
chr2\tX\tCDS\t100\t500\t.\t+\t0\tgene_id "G3"; transcript_id "T3"; gene_name "CLEAN";
chr2\tX\tthree_prime_utr\t501\t900\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_name "CLEAN";
"""


@pytest.fixture
def gtf(tmp_path: Path) -> Path:
    path = tmp_path / "test.gtf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GTF)
    return path


def test_transcript_genes_maps_every_transcript(gtf: Path):
    genes = transcript_genes(gtf)
    assert genes["T1"] == "ACTB"
    assert genes["T3"] == "CLEAN"


def test_a_neighbouring_gene_in_the_downstream_window_is_found(gtf: Path):
    """T1's 3' UTR runs into G2's coding sequence; T3 sits alone on another chromosome."""

    overlapping = overlapping_downstream_cds(gtf)
    assert "T1" in overlapping
    assert "T3" not in overlapping


def test_a_transcript_overlapping_only_itself_is_kept(gtf: Path):
    """Its own coding sequence crossing its own untranslated region is isoform structure."""

    assert "T3" not in overlapping_downstream_cds(gtf)


def test_a_transcript_qualifying_in_only_one_condition_is_dropped_from_both():
    """Otherwise the treated and untreated medians are taken over different transcripts, and the
    comparison is partly about which transcripts cleared the coverage bar."""

    both = pd.concat(
        [
            counts(transcript="shared", sample="untreated"),
            counts(transcript="shared", sample="g418"),
            counts(transcript="thin", sample="untreated", cds_inframe=MINIMUM_CDS_PSITES - 1),
            counts(transcript="thin", sample="g418"),
        ]
    )
    kept = qualifying(both)
    assert set(kept["transcript"]) == {"shared"}


def test_a_transcript_without_a_usable_extension_window_is_excluded():
    """No next in-frame stop means no window in which readthrough could be measured."""

    assert qualifying(counts(extension=float("nan"))).empty


def test_frame_specificity_requires_in_frame_to_outpace_out_of_frame():
    ratios = pd.DataFrame(
        {
            "sample": REPLICATES["sample"],
            "in_frame": [0.01, 0.01, 0.01, 0.03, 0.03, 0.03],
            "flat": [0.01, 0.01, 0.01, 0.011, 0.011, 0.011],
            "equal": [0.01, 0.01, 0.01, 0.03, 0.03, 0.03],
        }
    )
    rise = paired_effect(ratios, "in_frame", REPLICATES)
    assert frame_specific(rise, paired_effect(ratios, "flat", REPLICATES)) is True
    # Both frames lifted by the same amount is not decoding, it is something else.
    assert frame_specific(rise, paired_effect(ratios, "equal", REPLICATES)) is False


def test_the_signature_requires_all_three_conditions():
    ratios = pd.DataFrame(
        {
            "sample": REPLICATES["sample"],
            "readthrough": [0.01, 0.01, 0.01, 0.03, 0.03, 0.03],
            "termination": [0.50, 0.50, 0.50, 0.30, 0.30, 0.30],
            "out_of_frame": [0.01, 0.01, 0.01, 0.011, 0.011, 0.011],
        }
    )
    effects = {q: paired_effect(ratios, q, REPLICATES) for q in ratios.columns[1:]}
    assert signature(effects) == {
        "downstream_rose": True,
        "termination_fell": True,
        "frame_specific": True,
    }
    # Termination rising instead of falling is stalling, not readthrough.
    stalling = ratios.assign(termination=[0.30, 0.30, 0.30, 0.50, 0.50, 0.50])
    effects = {q: paired_effect(stalling, q, REPLICATES) for q in stalling.columns[1:]}
    assert signature(effects)["termination_fell"] is False


def test_a_library_outside_the_comparison_cannot_narrow_the_universe():
    """Calu-6 sits in the same counts table but takes no part in the HEK293T comparison, so its
    coverage must not decide which transcripts HEK293T is allowed to keep."""

    hek = ["hek293t_untreated_rep1_riboseq", "hek293t_g418_rep1_riboseq"]
    table = pd.concat(
        [counts(transcript="t", sample=name) for name in hek]
        + [counts(transcript="t", sample="calu6_untreated_rep1_riboseq", cds_inframe=1)]
    )
    assert qualifying(table).empty
    assert set(qualifying(table, samples=hek)["sample"]) == set(hek)
