import re
from pathlib import Path

import pytest

from riborescue.inputs import INPUTS, Input, UnknownInputError, data_root, digest, fetch

EMPTY_MD5 = "d41d8cd98f00b204e9800998ecf8427e"


def test_every_declared_input_names_its_source_licence_and_checksum():
    for name, declared in INPUTS.items():
        assert declared.name == name
        assert re.fullmatch(r"[0-9a-f]{32}", declared.md5)
        assert declared.url.startswith("https://")
        assert declared.source and declared.licence
        assert not declared.path.is_absolute()


def test_the_data_root_follows_the_environment(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("RIBORESCUE_DATA", "/elsewhere/inputs")
    assert data_root() == Path("/elsewhere/inputs")
    monkeypatch.delenv("RIBORESCUE_DATA")
    assert data_root() == Path("data")


def test_digest_matches_the_content(tmp_path: Path):
    empty = tmp_path / "empty"
    empty.write_bytes(b"")
    assert digest(empty) == EMPTY_MD5


def test_an_undeclared_input_is_refused(tmp_path: Path):
    with pytest.raises(UnknownInputError):
        fetch("nothing_of_the_sort", tmp_path)


def _probe(monkeypatch: pytest.MonkeyPatch) -> Input:
    """Declare an input that hashes to the digest of an empty file, so tests need no network."""

    probe = Input(
        name="probe",
        url="https://example.invalid/probe",
        md5=EMPTY_MD5,
        size=0,
        path=Path("probe/file.bin"),
        source="the test suite",
        licence="none",
    )
    monkeypatch.setitem(INPUTS, "probe", probe)
    return probe


def test_a_download_whose_digest_is_wrong_is_deleted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    probe = _probe(monkeypatch)

    def wrong(_: str, destination: Path) -> None:
        destination.write_bytes(b"not the published file")

    with pytest.raises(OSError, match="expected"):
        fetch(probe.name, tmp_path, retrieve=wrong)
    assert not probe.resolve(tmp_path).exists()


def test_a_verified_file_is_not_fetched_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    probe = _probe(monkeypatch)
    present = probe.resolve(tmp_path)
    present.parent.mkdir(parents=True)
    present.write_bytes(b"")
    calls = 0

    def counted(_: str, destination: Path) -> None:
        nonlocal calls
        calls += 1
        destination.write_bytes(b"")

    assert fetch(probe.name, tmp_path, retrieve=counted) == present
    assert calls == 0
    fetch(probe.name, tmp_path, retrieve=counted, force=True)
    assert calls == 1
