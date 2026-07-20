import re
from pathlib import Path

import pytest

from riborescue import inputs
from riborescue.inputs import INPUTS, Input, UnknownInputError, data_root, fetch


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


def test_an_undeclared_input_is_refused(tmp_path: Path):
    with pytest.raises(UnknownInputError):
        fetch("nothing_of_the_sort", tmp_path)


def _record(monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    """Capture what fetch asks the downloader for, without touching the network."""

    seen: dict[str, object] = {}

    def retrieve(url: str, known_hash: str, fname: str, path: Path) -> str:
        seen.update(url=url, known_hash=known_hash, fname=fname, path=path)
        return str(Path(path) / fname)

    monkeypatch.setattr(inputs.pooch, "retrieve", retrieve)
    return seen


def test_an_input_is_fetched_to_its_declared_path_under_its_published_digest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    seen = _record(monkeypatch)
    declared = INPUTS["toledano_treated_samples"]

    assert fetch(declared.name, tmp_path) == declared.resolve(tmp_path)
    assert seen["url"] == declared.url
    assert seen["known_hash"] == f"md5:{declared.md5}"
    assert seen["fname"] == "treated_samples.rds"
    assert seen["path"] == tmp_path / "toledano"


def test_forcing_a_fetch_discards_the_file_already_there(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    _record(monkeypatch)
    probe = Input(
        name="probe",
        url="https://example.invalid/probe",
        md5="d41d8cd98f00b204e9800998ecf8427e",
        path=Path("probe/file.bin"),
        source="the test suite",
        licence="none",
    )
    monkeypatch.setitem(INPUTS, probe.name, probe)
    present = probe.resolve(tmp_path)
    present.parent.mkdir(parents=True)
    present.write_bytes(b"stale")

    fetch(probe.name, tmp_path)
    assert present.read_bytes() == b"stale"
    fetch(probe.name, tmp_path, force=True)
    assert not present.exists()
