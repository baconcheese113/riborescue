"""Public inputs — declared once, fetched from their source, verified by checksum.

Nothing large is committed. Reproducibility rests on regeneration, so every input names where it
comes from, what it should hash to, and what licence it carries. Fetching is `pooch`'s job: it
downloads to a temporary file and only moves it into place once the digest matches, so a corrupt
or truncated file is never mistaken for the real input.
"""

import os
from dataclasses import dataclass
from pathlib import Path

import pooch

__all__ = ["INPUTS", "Input", "UnknownInputError", "data_root", "fetch"]


class UnknownInputError(KeyError):
    pass


@dataclass(frozen=True)
class Input:
    """One fetchable public file: where it lives, what it must hash to, and who published it."""

    name: str
    url: str
    md5: str
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
        path=Path("toledano/treated_samples.rds"),
        source="doi:10.6084/m9.figshare.25138712.v6",
        licence="CC BY 4.0",
    ),
}


def data_root() -> Path:
    """The working directory for fetched inputs, overridable for a machine with data elsewhere."""

    return Path(os.environ.get("RIBORESCUE_DATA", "data"))


def fetch(name: str, root: Path | None = None, *, force: bool = False) -> Path:
    """Fetch a declared input into the data root and return its path.

    A file already present with the right digest is left alone; one with the wrong digest is
    fetched again.
    """

    if (declared := INPUTS.get(name)) is None:
        raise UnknownInputError(f"{name!r} is not a declared input; known: {', '.join(INPUTS)}")

    destination = declared.resolve(root if root is not None else data_root())
    if force:
        destination.unlink(missing_ok=True)
    pooch.retrieve(
        declared.url,
        known_hash=f"md5:{declared.md5}",
        fname=destination.name,
        path=destination.parent,
    )
    return destination
