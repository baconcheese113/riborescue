"""The numbers the positive control rests on, pinned so a refactor cannot move them quietly.

Every other gate in this project asks whether the code runs and whether the contracts hold. None of
them asks whether a number is the number it was yesterday, and on 2026-07-26 that gap let a refactor
change the G418 result by 1.3% with the whole suite green. It was noticed by accident.

So the six per-library ratios and the three aggregate effects are committed, and this compares
what the pipeline produces now against them — the same idea as the oracle parity fixtures, applied
to the half of the project those fixtures do not cover.

**Nothing here regenerates the expected values.** A scientific correction should be possible and
must never be silent, so changing one means editing the fixture by hand and saying why in the
commit — which is what happened when the neighbouring-gene filter was fixed, and is the only reason
anyone can now see that it moved.

The comparison needs the pipeline to have run, because its input is a 1.5 GB annotation and a
per-transcript count table, neither of which belongs in the repository. It skips with the command
that produces them rather than passing vacuously.
"""

from pathlib import Path

import pandas as pd
import pytest

_PRODUCED = Path("results/readthrough/gse144140")
_EXPECTED = Path("tests/fixtures/readthrough")

TOLERANCE = 1e-6
"""How far a quantity may move before this fails.

The computation is deterministic, so any real difference is a change in behaviour rather than in
arithmetic. This is far tighter than the 0.0024 shift that motivated the fixture and far looser than
float noise.
"""

_QUANTITIES = ["downstream_occupancy", "termination_occupancy", "frame_gap"]


def _pair(produced: str, expected: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    made = _PRODUCED / produced
    if not made.exists():
        pytest.skip(f"{made} is absent; produce it with `pixi run readthrough-gse144140`")
    return pd.read_csv(made, sep="\t"), pd.read_csv(_EXPECTED / expected, sep="\t")


@pytest.mark.control
@pytest.mark.integration
def test_the_positive_control_still_reports_the_effects_it_reported():
    made, expected = _pair("g418_vs_dmso.unpaired.tsv", "gse144140_g418_vs_dmso.tsv")
    merged = expected.merge(made, on="quantity", suffixes=("_expected", "_now"), validate="1:1")
    assert len(merged) == len(expected), "a quantity the fixture pins is no longer reported"
    for row in merged.itertuples():
        assert row.mean_difference_now == pytest.approx(
            row.mean_difference_expected, abs=TOLERANCE
        ), f"{row.quantity} moved; edit the fixture and say why, or find what changed"


@pytest.mark.control
@pytest.mark.integration
def test_every_library_still_reports_the_ratios_it_reported():
    made, expected = _pair(
        "g418_vs_dmso.unpaired_by_library.tsv", "gse144140_g418_vs_dmso_by_library.tsv"
    )
    merged = expected.merge(made, on="sample", suffixes=("_expected", "_now"), validate="1:1")
    assert len(merged) == len(expected), "a library the fixture pins is no longer in the contrast"
    for quantity in _QUANTITIES:
        for row in merged.itertuples():
            assert getattr(row, f"{quantity}_now") == pytest.approx(
                getattr(row, f"{quantity}_expected"), abs=TOLERANCE
            ), f"{row.sample} {quantity} moved"


@pytest.mark.control
@pytest.mark.integration
def test_the_qualifying_universe_is_the_same_size():
    # The count is what changed when the neighbouring-gene filter started working: 495 to 486. It is
    # pinned separately because a shift here explains a shift in every quantity above at once.
    made, expected = _pair(
        "g418_vs_dmso.unpaired_by_library.tsv", "gse144140_g418_vs_dmso_by_library.tsv"
    )
    assert set(made["transcripts"]) == set(expected["transcripts"])
