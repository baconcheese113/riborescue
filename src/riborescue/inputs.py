"""Public inputs — declared once, fetched from their source, verified by checksum.

Nothing large is committed. Reproducibility rests on regeneration, so every input names where it
comes from, what it should hash to, and what licence it carries. A download whose digest does not
match is deleted rather than used.
"""

import hashlib
import os
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

__all__ = ["INPUTS", "Input", "UnknownInputError", "data_root", "digest", "fetch"]

_CHUNK = 1 << 20


class UnknownInputError(KeyError):
    pass


@dataclass(frozen=True)
class Input:
    """One fetchable public file: where it lives, what it must hash to, and who published it."""

    name: str
    url: str
    md5: str
    size: int
    path: Path
    source: str
    licence: str

    def resolve(self, root: Path) -> Path:
        return root / self.path


INPUTS: dict[str, Input] = {
    "toledano_treated_samples": Input(
        name="toledano_treated_samples",
        url="https://ndownloader.figshare.com/files/44388752",
        md5="bb3474efa01f074466912a36964827c2",
        size=51_111_391,
        path=Path("toledano/treated_samples.rds"),
        source="doi:10.6084/m9.figshare.25138712.v6",
        licence="CC BY 4.0",
    ),
}


def data_root() -> Path:
    """The working directory for fetched inputs, overridable for a machine with data elsewhere."""

    return Path(os.environ.get("RIBORESCUE_DATA", "data"))


def digest(path: Path) -> str:
    # md5 because that is what figshare publishes; this checks transfer integrity, not trust
    md5 = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            md5.update(chunk)
    return md5.hexdigest()


def fetch(
    name: str,
    root: Path | None = None,
    *,
    force: bool = False,
    retrieve: Callable[[str, Path], None] | None = None,
) -> Path:
    """Fetch an input into the data root and return its path, verifying the published checksum.

    An already-present file with the right digest is left alone; one with the wrong digest is
    re-fetched. A fresh download that fails verification is removed, so a corrupt file can never be
    mistaken for the real input.
    """

    if (declared := INPUTS.get(name)) is None:
        raise UnknownInputError(f"{name!r} is not a declared input; known: {', '.join(INPUTS)}")

    destination = declared.resolve(root if root is not None else data_root())
    if not force and destination.exists() and digest(destination) == declared.md5:
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    (retrieve or _download)(declared.url, destination)

    if (found := digest(destination)) != declared.md5:
        destination.unlink()
        raise OSError(f"{declared.name} downloaded with digest {found}, expected {declared.md5}")
    return destination


def _download(url: str, destination: Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".partial")
    with urllib.request.urlopen(url) as response, partial.open("wb") as handle:
        while chunk := response.read(_CHUNK):
            handle.write(chunk)
    partial.replace(destination)
