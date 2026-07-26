import gzip
from pathlib import Path

import pandas as pd
import pytest

from riborescue.riboseq.readthrough_assay import (
    MINIMUM_CDS_PSITES,
    MINIMUM_UTR3_LENGTH,
    PROGRAMMED_READTHROUGH,
    extension_windows,
    library_ratios,
    next_in_frame_stop,
    overlapping_downstream_cds,
    paired_effect,
    qualifying,
    signature,
    stalling,
    termination_arms_separate,
    transcript_genes,
    unpaired_effect,
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
        "cds_frame0": 1000,
        "cds_frame1": 250,
        "cds_frame2": 250,
        "termination": 50,
        "extension": 300,
        "extension_frame0": 10,
        "extension_frame1": 3,
        "extension_frame2": 2,
        "l_cds": 900,
        "l_utr3": 500,
    }
    return pd.DataFrame([base | overrides])


def test_a_short_untranslated_region_is_excluded():
    """The termination peak spills into a window this short."""

    assert qualifying(counts(l_utr3=MINIMUM_UTR3_LENGTH - 1)).empty
    assert len(qualifying(counts(l_utr3=MINIMUM_UTR3_LENGTH))) == 1


def test_a_transcript_without_enough_coding_coverage_is_excluded():
    assert qualifying(counts(cds_frame0=MINIMUM_CDS_PSITES - 1)).empty
    assert len(qualifying(counts(cds_frame0=MINIMUM_CDS_PSITES))) == 1


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

    both = pd.concat([counts(sample="untreated"), counts(sample="g418", extension_frame0=99)])
    kept = qualifying(both)
    assert set(kept["sample"]) == {"untreated", "g418"}


def test_pooled_quantities_sum_counts_before_dividing():
    """A frame composition is a proportion, so it pools; the median is reported beside it."""

    frame = pd.concat(
        [
            counts(transcript="a", extension_frame0=10, cds_frame0=1000),
            counts(transcript="b", extension_frame0=10, cds_frame0=1000),
            counts(transcript="c", extension_frame0=9000, cds_frame0=1000),
        ]
    )
    ratios = library_ratios(frame)
    # pooled: (10+10+9000) / (1000+1000+1000)
    assert ratios.loc[0, "downstream_occupancy"] == pytest.approx(9020 / 3000)
    # the median transcript is untouched by the outlier
    assert ratios.loc[0, "median_readthrough"] == pytest.approx(0.01)
    assert ratios.loc[0, "transcripts"] == 3


def test_the_frame_gap_is_measured_against_the_library_own_coding_frame():
    """A fixed third would be the wrong null: coding-frame occupancy varies library to library."""

    frame = counts(
        extension_frame0=50,
        extension_frame1=25,
        extension_frame2=25,
        cds_frame0=600,
        cds_frame1=200,
        cds_frame2=200,
    )
    ratios = library_ratios(frame)
    assert ratios.loc[0, "downstream_share0"] == pytest.approx(0.5)
    assert ratios.loc[0, "cds_share0"] == pytest.approx(0.6)
    assert ratios.loc[0, "frame_gap"] == pytest.approx(-0.1)


def test_the_share_with_any_downstream_signal_is_reported():
    """A pooled proportion resting on a few transcripts must not look like a broad one."""

    frame = pd.concat(
        [
            counts(transcript="a", extension_frame0=10, extension_frame1=0, extension_frame2=0),
            counts(transcript="b", extension_frame0=0, extension_frame1=0, extension_frame2=0),
        ]
    )
    assert library_ratios(frame).loc[0, "with_downstream"] == pytest.approx(0.5)


def test_paired_effect_compares_within_replicate():
    ratios = pd.DataFrame(
        {
            "sample": REPLICATES["sample"],
            "readthrough": [0.010, 0.020, 0.030, 0.015, 0.026, 0.037],
        }
    )
    effect = paired_effect(ratios, "readthrough", REPLICATES, "g418", "untreated")
    assert effect.differences == pytest.approx((0.005, 0.006, 0.007))
    assert effect.consistent is True
    assert effect.mean_difference == pytest.approx(0.006)


def test_a_replicate_missing_a_condition_is_refused():
    ratios = pd.DataFrame({"sample": REPLICATES["sample"][:5], "readthrough": [0.1] * 5})
    with pytest.raises(ValueError, match="without both conditions"):
        paired_effect(ratios, "readthrough", REPLICATES.iloc[:5], "g418", "untreated")


def test_an_inconsistent_effect_is_reported_as_such():
    ratios = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.02, 0.03, 0.02, 0.01, 0.04]}
    )
    effect = paired_effect(ratios, "readthrough", REPLICATES, "g418", "untreated")
    assert effect.consistent is False


def test_the_interval_widens_with_disagreement():
    tight = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.01, 0.01, 0.02, 0.02, 0.02]}
    )
    loose = pd.DataFrame(
        {"sample": REPLICATES["sample"], "readthrough": [0.01, 0.01, 0.01, 0.05, 0.02, 0.001]}
    )
    narrow = paired_effect(tight, "readthrough", REPLICATES, "g418", "untreated").interval
    wide = paired_effect(loose, "readthrough", REPLICATES, "g418", "untreated").interval
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


# GENCODE's GTF writes every untranslated region as a bare `UTR` on both sides of the coding
# sequence and never emits `three_prime_utr`, so a fixture spelling it the Ensembl way exercises a
# feature the real annotation does not contain.
GENCODE_GTF = """\
chr1\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\tCDS\t100\t500\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\tUTR\t1\t99\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\tUTR\t501\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ACTB";
chr1\tX\ttranscript\t600\t1200\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; gene_name "NEIGH";
chr1\tX\tCDS\t600\t1200\t.\t+\t0\tgene_id "G2"; transcript_id "T2"; gene_name "NEIGH";
chr2\tX\ttranscript\t100\t900\t.\t-\t.\tgene_id "G4"; transcript_id "T4"; gene_name "MINUS";
chr2\tX\tCDS\t500\t900\t.\t-\t0\tgene_id "G4"; transcript_id "T4"; gene_name "MINUS";
chr2\tX\tUTR\t100\t499\t.\t-\t.\tgene_id "G4"; transcript_id "T4"; gene_name "MINUS";
chr2\tX\ttranscript\t100\t300\t.\t-\t.\tgene_id "G5"; transcript_id "T5"; gene_name "NEIGH2";
chr2\tX\tCDS\t100\t300\t.\t-\t0\tgene_id "G5"; transcript_id "T5"; gene_name "NEIGH2";
"""


@pytest.fixture
def gencode_gtf(tmp_path: Path) -> Path:
    path = tmp_path / "gencode.gtf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GENCODE_GTF)
    return path


def test_a_bare_utr_is_placed_against_the_coding_sequence(gencode_gtf: Path):
    """The annotation the project actually reads names no 3' UTR, only `UTR` on both sides.

    T1's downstream region runs into G2's coding sequence and its upstream one does not, so the
    exclusion has to tell the two apart by where they sit rather than by what they are called.
    """

    assert "T1" in overlapping_downstream_cds(gencode_gtf)


def test_the_downstream_side_follows_the_strand(gencode_gtf: Path):
    """On the minus strand the 3' region is the one *before* the coding sequence, not after it."""

    assert "T4" in overlapping_downstream_cds(gencode_gtf)


def test_the_upstream_region_is_not_the_downstream_window(gencode_gtf: Path):
    """T1 also carries a 5' UTR; a neighbour there says nothing about readthrough past its stop."""

    upstream_only = """\
chr3\tX\ttranscript\t100\t900\t.\t+\t.\tgene_id "G6"; transcript_id "T6"; gene_name "A";
chr3\tX\tCDS\t500\t900\t.\t+\t0\tgene_id "G6"; transcript_id "T6"; gene_name "A";
chr3\tX\tUTR\t100\t499\t.\t+\t.\tgene_id "G6"; transcript_id "T6"; gene_name "A";
chr3\tX\ttranscript\t100\t300\t.\t+\t.\tgene_id "G7"; transcript_id "T7"; gene_name "B";
chr3\tX\tCDS\t100\t300\t.\t+\t0\tgene_id "G7"; transcript_id "T7"; gene_name "B";
"""
    path = gencode_gtf.with_name("upstream.gtf.gz")
    with gzip.open(path, "wt") as handle:
        handle.write(upstream_only)
    assert "T6" not in overlapping_downstream_cds(path)


def test_a_transcript_qualifying_in_only_one_condition_is_dropped_from_both():
    """Otherwise the treated and untreated medians are taken over different transcripts, and the
    comparison is partly about which transcripts cleared the coverage bar."""

    both = pd.concat(
        [
            counts(transcript="shared", sample="untreated"),
            counts(transcript="shared", sample="g418"),
            counts(transcript="thin", sample="untreated", cds_frame0=MINIMUM_CDS_PSITES - 1),
            counts(transcript="thin", sample="g418"),
        ]
    )
    kept = qualifying(both)
    assert set(kept["transcript"]) == {"shared"}


def test_a_transcript_without_a_usable_extension_window_is_excluded():
    """No next in-frame stop means no window in which readthrough could be measured."""

    assert qualifying(counts(extension=float("nan"))).empty


def test_a_library_outside_the_comparison_cannot_narrow_the_universe():
    """Calu-6 sits in the same counts table but takes no part in the HEK293T comparison, so its
    coverage must not decide which transcripts HEK293T is allowed to keep."""

    hek = ["hek293t_untreated_rep1_riboseq", "hek293t_g418_rep1_riboseq"]
    table = pd.concat(
        [counts(transcript="t", sample=name) for name in hek]
        + [counts(transcript="t", sample="calu6_untreated_rep1_riboseq", cds_frame0=1)]
    )
    assert qualifying(table).empty
    assert set(qualifying(table, samples=hek)["sample"]) == set(hek)


# The stop codon starts at index 6 in each of these; what follows differs.
def test_the_next_stop_immediately_after_is_three_bases_on():
    assert next_in_frame_stop("AAACCCTAATGAAAA", 6) == 3


def test_one_sense_codon_before_the_next_stop_is_six():
    assert next_in_frame_stop("AAACCCTAAGGGTAGAAA", 6) == 6


def test_a_stop_out_of_frame_does_not_count():
    """TGA sits one base off the frame, so the in-frame TAA two codons on is the answer."""

    assert next_in_frame_stop("AAACCCTAAGTGAGGTAA", 6) == 9


def test_no_further_stop_in_frame_gives_nothing():
    assert next_in_frame_stop("AAACCCTAAGGGCCCAAA", 6) is None


def test_a_transcript_ending_at_its_stop_gives_nothing():
    assert next_in_frame_stop("AAACCCTAA", 6) is None


def test_a_trailing_partial_codon_is_not_read():
    assert next_in_frame_stop("AAACCCTAAGGGTA", 6) is None


def test_extension_windows_skips_transcripts_without_a_sequence(tmp_path: Path):
    fasta = tmp_path / "t.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">T1|ENSG1|x|x|N|N|21|protein_coding|\nAAACCCTAAGGGTAGAAA\n")
    annotation = pd.DataFrame({"transcript": ["T1", "MISSING"], "l_utr5": [3, 3], "l_cds": [6, 6]})
    windows = extension_windows(fasta, annotation)
    assert list(windows["transcript"]) == ["T1"]
    assert windows.loc[0, "extension"] == 6


CONFIRM = pd.DataFrame(
    {
        "sample": [f"{arm}_{r}" for arm in ("dmso", "g418", "sri") for r in "ABC"],
        "replicate": list("ABC") * 3,
        "treatment": ["dmso"] * 3 + ["g418"] * 3 + ["sri"] * 3,
    }
)


def _ratios(**arms) -> pd.DataFrame:
    frame = pd.DataFrame({"sample": CONFIRM["sample"]})
    for quantity, values in arms.items():
        frame[quantity] = [v for arm in ("dmso", "g418", "sri") for v in values[arm]]
    return frame


def test_unpaired_effect_separates_groups_without_pairing_them():
    ratios = _ratios(q={"dmso": [0.10, 0.11, 0.12], "g418": [0.20, 0.21, 0.22], "sri": [0.1] * 3})
    effect = unpaired_effect(ratios, "q", CONFIRM, "g418", "dmso")
    assert effect.mean_difference == pytest.approx(0.10)
    assert effect.consistent is True


def test_overlapping_groups_are_not_consistent():
    """Consistency means complete separation, not merely a difference in means."""

    ratios = _ratios(q={"dmso": [0.10, 0.30, 0.12], "g418": [0.20, 0.21, 0.22], "sri": [0.1] * 3})
    assert unpaired_effect(ratios, "q", CONFIRM, "g418", "dmso").consistent is False


def test_an_arm_that_does_not_exist_is_refused():
    ratios = _ratios(q={"dmso": [0.1] * 3, "g418": [0.2] * 3, "sri": [0.1] * 3})
    with pytest.raises(ValueError, match="no libraries for"):
        unpaired_effect(ratios, "q", CONFIRM, "missing", "dmso")


def test_the_signature_needs_all_three_conditions():
    ratios = _ratios(
        downstream_occupancy={"dmso": [0.01] * 3, "g418": [0.03] * 3, "sri": [0.01] * 3},
        termination_occupancy={"dmso": [0.50] * 3, "g418": [0.30] * 3, "sri": [0.70] * 3},
        frame_gap={"dmso": [-0.10, -0.11, -0.12], "g418": [-0.01, -0.02, -0.03], "sri": [-0.1] * 3},
    )
    effects = {
        q: unpaired_effect(ratios, q, CONFIRM, "g418", "dmso")
        for q in ("downstream_occupancy", "termination_occupancy", "frame_gap")
    }
    assert signature(effects) == {
        "downstream_rose": True,
        "termination_fell": True,
        "frame_moved_to_coding": True,
    }


def test_a_frame_gap_interval_spanning_zero_does_not_pass():
    """The frame condition is the claim the control makes, so it alone carries the interval.

    The groups here separate completely and the mean moves the right way, so direction and
    consistency both hold. Only the width of the interval refuses it — remove that clause from
    `signature` and this test fails.
    """

    ratios = _ratios(
        downstream_occupancy={"dmso": [0.01] * 3, "g418": [0.03] * 3, "sri": [0.01] * 3},
        termination_occupancy={"dmso": [0.50] * 3, "g418": [0.30] * 3, "sri": [0.70] * 3},
        # every treated value exceeds every control value, but both groups are so scattered that
        # the difference is not resolved
        frame_gap={"dmso": [-0.50, -0.30, -0.02], "g418": [0.30, 0.02, 0.60], "sri": [-0.1] * 3},
    )
    effects = {
        q: unpaired_effect(ratios, q, CONFIRM, "g418", "dmso")
        for q in ("downstream_occupancy", "termination_occupancy", "frame_gap")
    }
    gap = effects["frame_gap"]
    assert gap.mean_difference > 0
    assert gap.consistent is True
    low, high = gap.interval
    assert low < 0 < high
    assert signature(effects)["frame_moved_to_coding"] is False


def test_a_stalling_compound_is_recognised_by_its_own_direction():
    """Raised occupancy at the stop with none beyond it is stalling, not readthrough."""

    ratios = _ratios(
        downstream_occupancy={"dmso": [0.01] * 3, "g418": [0.03] * 3, "sri": [0.005] * 3},
        termination_occupancy={"dmso": [0.50] * 3, "g418": [0.30] * 3, "sri": [0.70] * 3},
        frame_gap={"dmso": [-0.1] * 3, "g418": [-0.01] * 3, "sri": [-0.1] * 3},
    )
    quantities = ("downstream_occupancy", "termination_occupancy", "frame_gap")
    sri = {q: unpaired_effect(ratios, q, CONFIRM, "sri", "dmso") for q in quantities}
    g418 = {q: unpaired_effect(ratios, q, CONFIRM, "g418", "dmso") for q in quantities}
    assert stalling(sri) is True
    assert stalling(g418) is False


def test_overlapping_arms_qualify_the_stalling_verdict_without_deciding_it():
    """The verdict is on the means, so separation reports beside it rather than inside it.

    The arms here overlap, so the diagnostic fails while the verdict holds. Folding the diagnostic
    into `stalling` would let a stricter condition decide a verdict that is not defined in terms
    of it.
    """

    ratios = _ratios(
        downstream_occupancy={"dmso": [0.01] * 3, "g418": [0.03] * 3, "sri": [0.005] * 3},
        # termination rises on average, but the groups overlap rather than separating
        termination_occupancy={"dmso": [0.50, 0.90, 0.50], "g418": [0.30] * 3, "sri": [0.70] * 3},
        frame_gap={"dmso": [-0.1] * 3, "g418": [-0.01] * 3, "sri": [-0.1] * 3},
    )
    quantities = ("downstream_occupancy", "termination_occupancy", "frame_gap")
    sri = {q: unpaired_effect(ratios, q, CONFIRM, "sri", "dmso") for q in quantities}
    assert sri["termination_occupancy"].mean_difference > 0
    assert sri["termination_occupancy"].consistent is False
    assert stalling(sri) is True
    assert termination_arms_separate(sri) is False


def test_the_welch_interval_rounds_degrees_of_freedom_down():
    """Rounding up would narrow the interval; the conservative direction is fewer degrees."""

    ratios = _ratios(q={"dmso": [0.10, 0.12, 0.20], "g418": [0.30, 0.31, 0.32], "sri": [0.1] * 3})
    low, high = unpaired_effect(ratios, "q", CONFIRM, "g418", "dmso").interval
    # a Welch interval on these groups has df near 2.3; flooring to 2 gives t = 4.303
    assert low < high
    assert (high - low) / 2 > 4.0 * float(
        (
            pd.Series([0.10, 0.12, 0.20]).var(ddof=1) / 3
            + pd.Series([0.30, 0.31, 0.32]).var(ddof=1) / 3
        )
        ** 0.5
    )


def test_an_incomplete_coding_sequence_yields_no_window():
    """A CDS whose length is not a multiple of three ends on an arbitrary triplet, so every frame
    measured from it would be off-register."""

    # position 6 holds GGG, not a stop, so the annotation does not end where it claims
    assert next_in_frame_stop("AAACCCGGGTAAAAA", 6) is None


def test_the_universe_is_counted_over_the_libraries_asked_for():
    """A library named for the contrast but missing from the counts must not relax the bar."""

    present = pd.concat([counts(transcript="t", sample="a"), counts(transcript="t", sample="b")])
    assert len(qualifying(present, samples=["a", "b"])) == 2
    # 'c' was asked for and is not there, so 't' no longer qualifies everywhere
    assert qualifying(present, samples=["a", "b", "c"]).empty


def test_pooled_counts_are_reported_beside_the_proportions():
    """A share resting on a dozen reads must be distinguishable from one resting on thousands."""

    frame = counts(extension_frame0=4, extension_frame1=3, extension_frame2=5, cds_frame0=1000)
    row = library_ratios(frame)
    assert row.loc[0, "downstream_total"] == 12
    assert row.loc[0, "cds_frame0_total"] == 1000
    assert row.loc[0, "termination_total"] == 50


def test_the_termination_ratio_divides_by_in_frame_coding_occupancy():
    frame = counts(termination=90, cds_frame0=900, cds_frame1=100, cds_frame2=100)
    assert library_ratios(frame).loc[0, "termination_occupancy"] == pytest.approx(0.1)


def test_lengths_are_collapsed_by_summing_the_ones_a_dataset_keeps():
    """One pass over the alignments serves the selected set and the published window alike."""

    import pandas as pd

    from riborescue.riboseq.readthrough_assay import collapse_lengths

    counts = pd.DataFrame(
        {
            "transcript": ["T1", "T1", "T1"],
            "sample": ["a", "a", "a"],
            "length": [28, 30, 21],
            "cds_frame0": [100, 200, 7],
            "cds_frame1": [10, 20, 1],
            "cds_frame2": [10, 20, 1],
            "extension_frame0": [5, 6, 0],
            "extension_frame1": [1, 1, 0],
            "extension_frame2": [1, 1, 0],
            "termination": [3, 4, 0],
            "cds_total": [120, 240, 9],
            "extension": [300, 300, 300],
            "l_cds": [900, 900, 900],
            "l_utr3": [400, 400, 400],
        }
    )
    kept = collapse_lengths(counts, [28, 30])
    assert len(kept) == 1
    assert kept.loc[0, "cds_frame0"] == 300
    assert kept.loc[0, "extension_frame0"] == 11
    assert kept.loc[0, "termination"] == 7
    # Carried, not summed: they describe the transcript rather than the reads.
    assert kept.loc[0, "extension"] == 300
    assert kept.loc[0, "l_cds"] == 900


def test_collapsing_refuses_counts_that_are_not_stratified():
    import pandas as pd
    import pytest

    from riborescue.riboseq.readthrough_assay import collapse_lengths

    with pytest.raises(ValueError, match="not stratified"):
        collapse_lengths(pd.DataFrame({"transcript": ["T1"], "sample": ["a"]}), [30])


def test_collapsing_then_pooling_equals_a_direct_pooled_calculation():
    """The point of stratifying: a subset sum must match computing over those reads directly."""

    import pandas as pd

    from riborescue.riboseq.readthrough_assay import collapse_lengths, library_ratios

    # Two transcripts, one library, counts split across three lengths; only two are kept.
    stratified = pd.DataFrame(
        [
            {
                "transcript": t,
                "sample": "a",
                "length": length,
                "cds_frame0": c0,
                "cds_frame1": c1,
                "cds_frame2": c2,
                "extension_frame0": e0,
                "extension_frame1": e1,
                "extension_frame2": e2,
                "termination": term,
                "cds_total": c0 + c1 + c2,
                "extension": 300,
                "l_cds": 900,
                "l_utr3": 500,
            }
            for t, length, c0, c1, c2, e0, e1, e2, term in [
                ("T1", 30, 400, 90, 90, 20, 4, 4, 12),
                ("T1", 31, 300, 70, 70, 10, 2, 2, 8),
                ("T1", 21, 5, 1, 1, 0, 0, 0, 0),  # dropped
                ("T2", 30, 250, 60, 60, 8, 1, 1, 6),
                ("T2", 31, 200, 50, 50, 6, 1, 1, 5),
                ("T2", 21, 9, 2, 2, 0, 0, 0, 0),  # dropped
            ]
        ]
    )

    collapsed = collapse_lengths(stratified, [30, 31])
    from_subset = library_ratios(collapsed)

    # The same reads, summed by hand over the kept lengths, one row per transcript.
    direct = pd.DataFrame(
        [
            {
                "transcript": "T1",
                "sample": "a",
                "cds_frame0": 700,
                "cds_frame1": 160,
                "cds_frame2": 160,
                "extension_frame0": 30,
                "extension_frame1": 6,
                "extension_frame2": 6,
                "termination": 20,
                "extension": 300,
                "l_cds": 900,
                "l_utr3": 500,
            },
            {
                "transcript": "T2",
                "sample": "a",
                "cds_frame0": 450,
                "cds_frame1": 110,
                "cds_frame2": 110,
                "extension_frame0": 14,
                "extension_frame1": 2,
                "extension_frame2": 2,
                "termination": 11,
                "extension": 300,
                "l_cds": 900,
                "l_utr3": 500,
            },
        ]
    )
    from_direct = library_ratios(direct)

    for column in ("downstream_occupancy", "termination_occupancy", "frame_gap"):
        assert from_subset.loc[0, column] == pytest.approx(from_direct.loc[0, column])
