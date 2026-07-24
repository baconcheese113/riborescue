"""The sequencing runs the Ribo-seq arm consumes, fetched from the archive that published them.

A samplesheet names one run per row with the checksum ENA publishes for each of its FASTQ files, so
a run arrives by the same verified path as every other input: nothing is committed, and a truncated
transfer cannot be mistaken for the data.

Staging turns the declared sheet into one naming local files, which is what the pipeline reads.
"""

import time
from collections.abc import Iterator
from pathlib import Path

import pandas as pd

from riborescue.core.inputs import Input, retrieve

__all__ = ["ATTEMPTS", "FASTQ_SUBDIR", "fastq_inputs", "stage"]

FASTQ_SUBDIR = Path("fastq")

ATTEMPTS = 5
"""How often one file is tried before the staging run gives up.

A public archive drops long transfers, and a run of this size makes that likely rather than
possible. A retry costs nothing: a file already on disk with the right digest is not fetched again.
"""

_SOURCE = "ENA, BioProject PRJNA576648"
_LICENCE = "public domain (INSDC)"


def fastq_inputs(run: pd.Series) -> Iterator[tuple[str, Input]]:
    """The FASTQ files one run publishes, paired with the column each belongs to."""

    for read in ("fastq_1", "fastq_2"):
        url, md5 = run[f"{read}_url"], run[f"{read}_md5"]
        if pd.isna(url):
            continue
        yield (
            read,
            Input(
                name=f"{run['sample']}_{read}",
                url=str(url),
                md5=str(md5),
                path=FASTQ_SUBDIR / str(url).rsplit("/", 1)[-1],
                source=_SOURCE,
                licence=_LICENCE,
            ),
        )


def _retrieve(declared: Input, root: Path | None, *, force: bool) -> Path:
    """Fetch one file, trying again on a dropped or timed-out transfer."""

    for attempt in range(1, ATTEMPTS + 1):
        try:
            return retrieve(declared, root, force=force and attempt == 1)
        except OSError:
            if attempt == ATTEMPTS:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def stage(runs: pd.DataFrame, root: Path | None = None, *, force: bool = False) -> pd.DataFrame:
    """Fetch every run's FASTQ and return the sheet with local paths in place of the URLs."""

    staged = runs.drop(columns=[c for c in runs.columns if c.startswith(("fastq_1", "fastq_2"))])
    located: dict[str, list[str | None]] = {"fastq_1": [], "fastq_2": []}
    for _, run in runs.iterrows():
        paths = {
            read: _retrieve(declared, root, force=force) for read, declared in fastq_inputs(run)
        }
        for read in located:
            path = paths.get(read)
            located[read].append(str(path.resolve()) if path is not None else None)
    return staged.assign(
        **{read: pd.Series(found, index=runs.index) for read, found in located.items()}
    )
