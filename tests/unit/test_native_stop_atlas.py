import gzip
from pathlib import Path

import pandas as pd

from riborescue.riboseq.native_stop_atlas import native_stop_occupancy, translate_extension


def _counts(rows: list[dict]) -> pd.DataFrame:
    base = {"length": 30, "l_utr3": 400, "extension": 300}
    return pd.DataFrame([{**base, **r} for r in rows])


def test_occupancy_pools_replicates_and_reports_depth():
    counts = _counts(
        [
            {"transcript": "T1", "sample": "dmso_a", "extension_frame0": 2, "cds_frame0": 100},
            {"transcript": "T1", "sample": "dmso_b", "extension_frame0": 4, "cds_frame0": 100},
            {"transcript": "T1", "sample": "g418_a", "extension_frame0": 30, "cds_frame0": 100},
            {"transcript": "T1", "sample": "g418_b", "extension_frame0": 30, "cds_frame0": 100},
        ]
    )
    arms = {"dmso": ["dmso_a", "dmso_b"], "g418": ["g418_a", "g418_b"]}
    atlas = native_stop_occupancy(counts, arms, [30], min_depth=100)
    row = atlas.set_index("transcript").loc["T1"].to_dict()
    # pooled: dmso 6/200 = 0.03, g418 60/200 = 0.30
    assert row["dmso_occupancy"] == 0.03
    assert row["g418_occupancy"] == 0.30
    assert row["dmso_depth"] == 200 and row["g418_depth"] == 200
    assert row["included"]


def test_a_transcript_below_depth_in_one_arm_is_not_included():
    counts = _counts(
        [
            {"transcript": "T1", "sample": "dmso_a", "extension_frame0": 1, "cds_frame0": 40},
            {"transcript": "T1", "sample": "g418_a", "extension_frame0": 10, "cds_frame0": 200},
        ]
    )
    atlas = native_stop_occupancy(
        counts, {"dmso": ["dmso_a"], "g418": ["g418_a"]}, [30], min_depth=100
    )
    assert not bool(atlas.set_index("transcript").loc["T1", "included"])


def test_an_overlapping_transcript_is_excluded_like_the_assay_does():
    """A 3'UTR running into another gene fills the window with that gene's ribosomes."""

    counts = _counts(
        [
            {"transcript": "clean", "sample": "g418_a", "extension_frame0": 5, "cds_frame0": 200},
            {
                "transcript": "overlap",
                "sample": "g418_a",
                "extension_frame0": 400,
                "cds_frame0": 200,
            },
        ]
    )
    atlas = native_stop_occupancy(
        counts, {"g418": ["g418_a"]}, [30], excluded_transcripts=frozenset({"overlap"})
    )
    assert set(atlas["transcript"]) == {"clean"}


def test_a_short_three_prime_utr_is_dropped():
    counts = _counts(
        [
            {
                "transcript": "T1",
                "sample": "g418_a",
                "extension_frame0": 5,
                "cds_frame0": 200,
                "l_utr3": 20,
            }
        ]
    )
    atlas = native_stop_occupancy(counts, {"g418": ["g418_a"]}, [30], min_utr3=50)
    assert atlas.empty


def test_the_extension_peptide_translates_from_past_the_native_stop(tmp_path: Path):
    # UTR5 (3) + CDS incl. stop (ATG AAA TAA = 9) + a downstream ORF that ends in TGA.
    fasta = tmp_path / "t.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">ENST1|gene\nCCCATGAAATAAGGGCCCTGATTT\n")
    annotation = pd.DataFrame([{"transcript": "ENST1", "l_utr5": 3, "l_cds": 9}])
    # next in-frame stop is TGA at offset 9 past the native stop's first base.
    extensions = pd.DataFrame([{"transcript": "ENST1", "extension": 9}])
    peptides = translate_extension(fasta, annotation, extensions).set_index("transcript")
    # window is extension-3 = 6 nt after the stop codon: GGG CCC -> G P
    assert peptides.loc["ENST1", "extension_peptide"] == "GP"
    assert peptides.loc["ENST1", "extension_aa"] == 2


def test_a_window_running_past_the_sequence_translates_only_what_is_there(tmp_path: Path):
    """A window wider than the sequence yields the residues that exist, not placeholders for the
    rest.

    `extension_windows` derives the window by finding the next in-frame stop within the sequence, so
    it can never outrun it and this cannot arise from the pipeline. It arises from an extension
    table read off disk that was built against a different reference, and the honest answer there is
    the peptide the sequence supports rather than one padded out to the width that was asked for.
    """

    fasta = tmp_path / "t.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">ENST1|gene\nCCCATGAAATAAGGGCCC\n")
    annotation = pd.DataFrame({"transcript": ["ENST1"], "l_utr5": [3], "l_cds": [9]})
    # 300 nt of window against 6 nt of remaining sequence.
    extensions = pd.DataFrame({"transcript": ["ENST1"], "extension": [303]})

    peptides = translate_extension(fasta, annotation, extensions).set_index("transcript")
    assert peptides.loc["ENST1", "extension_peptide"] == "GP"
    assert peptides.loc["ENST1", "extension_aa"] == 2
