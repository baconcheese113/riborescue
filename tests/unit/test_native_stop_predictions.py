import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from riborescue.variants.native_stop_predictions import (
    concordance,
    four_quadrants,
    native_stop_features,
)


def test_features_are_read_at_the_native_stop(tmp_path: Path):
    # UTR5 (CCC) + CDS: ATG AAA GGG TGA (last sense codon GGG, stop TGA) + 3'UTR CAT...
    fasta = tmp_path / "t.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">ENST1|gene\nCCCATGAAAGGGTGACATTTT\n")
    annotation = pd.DataFrame([{"transcript": "ENST1", "l_utr5": 3, "l_cds": 12}])
    row = native_stop_features(fasta, annotation).set_index("transcript").loc["ENST1"].to_dict()
    assert row["stop_type"] == "uga"
    assert row["up_123nt"] == "ggg"  # last sense codon
    assert row["down_123nt"] == "cau"  # first 3'UTR codon, transcribed


def test_a_transcript_not_ending_in_a_stop_is_skipped(tmp_path: Path):
    fasta = tmp_path / "t.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">ENST1|gene\nCCCATGAAAGGGAAACATTTT\n")  # ends in AAA, not a stop
    annotation = pd.DataFrame([{"transcript": "ENST1", "l_utr5": 3, "l_cds": 12}])
    assert native_stop_features(fasta, annotation).empty


def test_concordance_returns_a_rho_and_a_bootstrap_interval():
    rng = np.random.default_rng(0)
    predicted = pd.Series(rng.random(300))
    measured = predicted * 0.5 + rng.random(300) * 0.5  # correlated
    result = concordance(predicted, measured, draws=200)
    assert result["n"] == 300
    assert result["rho"] > 0.3
    assert result["low"] <= result["rho"] <= result["high"]


def test_concordance_is_undefined_for_too_few_points():
    result = concordance(pd.Series([1.0, 2.0, 3.0]), pd.Series([1.0, 2.0, 3.0]))
    assert result["n"] == 3
    assert np.isnan(result["rho"])


def test_the_four_groups_split_on_both_cuts():
    predicted = pd.Series([0.9, 0.9, 0.1, 0.1], index=["a", "b", "c", "d"])
    measured = pd.Series([0.9, 0.1, 0.9, 0.1], index=["a", "b", "c", "d"])
    groups = four_quadrants(predicted, measured, predicted_cut=0.5, measured_cut=0.5)
    assert groups["a"] == "both"
    assert groups["b"] == "predicted only"
    assert groups["c"] == "measured only"
    assert groups["d"] == "neither"
