import json
from pathlib import Path

import pandas as pd
import pytest

from riborescue.riboseq.calibration import (
    MINIMUM_PSITES,
    SURVEY_LENGTHS,
    read_manifest,
    select_lengths,
)


def _frames(rows: list[tuple[str, int, int, int]]) -> pd.DataFrame:
    """One row per sample and length, giving the three frame counts."""

    records = []
    for sample, length, frame0, off in rows:
        records += [
            {"sample": sample, "length": length, "frame": 0, "n": frame0},
            {"sample": sample, "length": length, "frame": 1, "n": off},
            {"sample": sample, "length": length, "frame": 2, "n": off},
        ]
    return pd.DataFrame(records)


def _offsets(samples: list[str], length: int = 30, offset: int = 12) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample": samples,
            "length": [length] * len(samples),
            "corrected_offset_from_5": [offset] * len(samples),
            "total_percentage": [90.0] * len(samples),
        }
    )


def test_a_length_is_kept_only_when_it_is_periodic_in_every_library():
    """One shared set is what makes the libraries of a contrast comparable."""

    frames = _frames(
        [
            ("a", 29, 900_000, 300_000),
            ("a", 30, 900_000, 300_000),
            ("b", 29, 900_000, 300_000),
            ("b", 30, 200_000, 900_000),  # frame 0 loses here
        ]
    )
    manifest = select_lengths("d", frames, _offsets(["a", "b"], length=29))
    assert manifest.lengths == (29,)


def test_a_length_carrying_almost_nothing_is_not_kept():
    """Frame 0 can win on a handful of reads; that is noise, not periodicity."""

    frames = _frames([("a", 30, 900_000, 300_000), ("a", 21, 40, 10)])
    manifest = select_lengths("d", frames, _offsets(["a"]))
    assert manifest.lengths == (30,)


def test_a_shallow_library_fails_and_is_not_dropped():
    frames = _frames([("a", 30, 900_000, 300_000), ("b", 30, 9_000, 3_000)])
    manifest = select_lengths("d", frames, _offsets(["a", "b"]))
    assert not manifest.passes
    assert [lib.sample for lib in manifest.libraries] == ["a", "b"]
    assert "P-sites" in manifest.failures["b"][0]


def test_an_offset_off_the_canonical_value_fails():
    frames = _frames([("a", 30, 900_000, 300_000)])
    manifest = select_lengths("d", frames, _offsets(["a"], offset=4))
    assert not manifest.passes
    assert "offset" in manifest.failures["a"][0]


def test_a_library_without_a_dominant_frame_fails_on_share():
    frames = _frames([("a", 30, 500_000, 480_000)])
    manifest = select_lengths("d", frames, _offsets(["a"]))
    assert not manifest.passes
    assert any("frame-0 share" in failure for failure in manifest.failures["a"])


def test_the_manifest_records_what_was_surveyed_and_what_was_kept(tmp_path: Path):
    frames = _frames([("a", 30, 900_000, 300_000)])
    manifest = select_lengths("d", frames, _offsets(["a"]), script_md5="abc123")
    record = json.loads(manifest.to_json())
    assert record["surveyed"] == list(SURVEY_LENGTHS)
    assert record["lengths"] == [30]
    assert record["script_md5"] == "abc123"
    assert record["passes"] is True


def test_a_failed_manifest_is_refused_when_read(tmp_path: Path):
    """The assay calls this rather than checking the file exists."""

    frames = _frames([("a", 30, 9_000, 3_000)])
    path = tmp_path / "calibration.json"
    path.write_text(select_lengths("d", frames, _offsets(["a"])).to_json())
    with pytest.raises(ValueError, match="did not pass its predeclared calibration"):
        read_manifest(path)


def test_a_passing_manifest_round_trips(tmp_path: Path):
    frames = _frames([("a", 30, 900_000, 300_000), ("a", 31, 900_000, 300_000)])
    path = tmp_path / "calibration.json"
    path.write_text(select_lengths("d", frames, _offsets(["a"])).to_json())
    loaded = read_manifest(path)
    assert loaded.dataset == "d"
    assert loaded.lengths == (30, 31)
    assert loaded.libraries[0].psites >= MINIMUM_PSITES
