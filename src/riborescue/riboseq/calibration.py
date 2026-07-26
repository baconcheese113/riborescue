"""Which footprint lengths a dataset keeps, and whether its libraries are fit to be compared.

ADR-0011. The lengths are chosen once per dataset from periodicity over internal coding sequence,
before any contrast is run, and the choice is recorded with the evidence behind it. The readthrough
assay reads that record and refuses to run without it, because nothing else in the pipeline asks
whether a library is a ribosome profile at all.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pydantic import BaseModel, ConfigDict, computed_field, field_serializer

__all__ = [
    "MINIMUM_FRAME0_SHARE",
    "MINIMUM_LENGTH_SHARE",
    "MINIMUM_PSITES",
    "OFFSET_FROM_5",
    "SURVEY_LENGTHS",
    "CalibrationManifest",
    "LibraryVerdict",
    "load_manifest",
    "read_manifest",
    "select_lengths",
]

SURVEY_LENGTHS = (18, 40)
"""What the first pass keeps, wide enough to hold every population worth seeing.

The authors of GSE144140 resolve empty-E-site ribosomes at 18-19 nucleotides, empty-A-site at 21-24
and full-length at 28-35. A survey narrower than their union would decide the question it is meant
to answer.
"""

MINIMUM_LENGTH_SHARE = 0.01
"""A length must carry this share of a library's assigned P-sites to be kept."""

MINIMUM_PSITES = 1_000_000
"""Assigned P-sites a library must have over the selected lengths."""

MINIMUM_FRAME0_SHARE = 0.40
"""Coding-frame-0 share a library must reach, pooled over the selected lengths."""

OFFSET_FROM_5 = (11, 14)
"""Inferred 5' offset the dominant length must fall within, inclusive."""


class LibraryVerdict(BaseModel):
    """One library measured against the thresholds fixed before it was looked at."""

    model_config = ConfigDict(frozen=True)

    sample: str
    psites: int
    frame0_share: float
    dominant_length: int
    offset_from_5: int
    failures: tuple[str, ...] = ()

    @field_serializer("frame0_share")
    def _rounded(self, share: float) -> float:
        return round(share, 5)

    @property
    def passes(self) -> bool:
        return not self.failures


class CalibrationManifest(BaseModel):
    """The lengths a dataset keeps, and each library's verdict on them."""

    model_config = ConfigDict(frozen=True)

    dataset: str
    lengths: tuple[int, ...]
    surveyed: tuple[int, int]
    script_md5: str | None = None
    libraries: tuple[LibraryVerdict, ...]

    @computed_field
    @property
    def passes(self) -> bool:
        return bool(self.lengths) and all(library.passes for library in self.libraries)

    @property
    def failures(self) -> dict[str, tuple[str, ...]]:
        return {lib.sample: lib.failures for lib in self.libraries if lib.failures}

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)


def _periodic_lengths(frames: pd.DataFrame) -> set[int]:
    """Lengths where frame 0 beats each off-frame and the length is not a rounding error.

    Judged per library and intersected by the caller: a length kept for the dataset has to hold in
    every library of it, because one shared set is what makes the libraries comparable.
    """

    wide = frames.pivot_table(
        index="length", columns="frame", values="n", aggfunc="sum", fill_value=0
    )
    for frame in (0, 1, 2):
        if frame not in wide.columns:
            wide[frame] = 0
    total = wide[[0, 1, 2]].sum(axis=1)
    share = total / total.sum() if total.sum() else total
    periodic = (wide[0] > wide[1]) & (wide[0] > wide[2]) & (share >= MINIMUM_LENGTH_SHARE)
    return {int(length) for length in wide.index[periodic]}


def select_lengths(
    dataset: str,
    frames: pd.DataFrame,
    offsets: pd.DataFrame,
    script_md5: str | None = None,
) -> CalibrationManifest:
    """Choose one length set for a dataset and judge every library on it.

    `frames` is the per-length three-frame table over internal coding sequence, `offsets` the
    inferred P-site offsets, both written by the first calibration pass.
    """

    samples = sorted(set(frames["sample"]))
    if not samples:
        raise ValueError("the frame table names no libraries")

    per_library = [_periodic_lengths(frames[frames["sample"] == s]) for s in samples]
    lengths = tuple(sorted(set.intersection(*per_library))) if per_library else ()

    verdicts = []
    for sample in samples:
        kept = frames[(frames["sample"] == sample) & (frames["length"].isin(lengths))]
        psites = int(kept["n"].sum())
        frame0 = int(kept.loc[kept["frame"] == 0, "n"].sum())
        share = frame0 / psites if psites else 0.0

        rows = offsets[offsets["sample"] == sample]
        dominant = rows.loc[rows["total_percentage"].idxmax()] if len(rows) else None
        length = int(dominant["length"]) if dominant is not None else 0
        offset = int(dominant["corrected_offset_from_5"]) if dominant is not None else 0

        failures: list[str] = []
        if not lengths:
            failures.append("no length is periodic in every library")
        if psites < MINIMUM_PSITES:
            failures.append(
                f"{psites:,} P-sites over the selected lengths, under {MINIMUM_PSITES:,}"
            )
        if share < MINIMUM_FRAME0_SHARE:
            failures.append(f"frame-0 share {share:.1%}, under {MINIMUM_FRAME0_SHARE:.0%}")
        if not OFFSET_FROM_5[0] <= offset <= OFFSET_FROM_5[1]:
            failures.append(f"5' offset {offset} nt outside {OFFSET_FROM_5[0]}-{OFFSET_FROM_5[1]}")
        verdicts.append(
            LibraryVerdict(
                sample=sample,
                psites=psites,
                frame0_share=share,
                dominant_length=length,
                offset_from_5=offset,
                failures=tuple(failures),
            )
        )

    return CalibrationManifest(
        dataset=dataset,
        lengths=lengths,
        surveyed=SURVEY_LENGTHS,
        libraries=tuple(verdicts),
        script_md5=script_md5,
    )


def load_manifest(path: Path) -> CalibrationManifest:
    """Parse a manifest and say nothing about whether it passed.

    Only the detectability arm of ADR-0019 reads a manifest this way, because its subject is a
    dataset that failed. Every analysis that reports a contrast calls `read_manifest` instead.

    `passes` is recomputed from the verdicts rather than read, so a manifest whose recorded verdict
    was edited to disagree with the libraries it lists cannot be believed on the edited value.
    """

    return CalibrationManifest.model_validate_json(path.read_text())


def read_manifest(path: Path) -> CalibrationManifest:
    """Load a manifest, refusing one that records a failure.

    The assay calls this rather than checking the file exists. A dataset whose libraries did not
    pass their predeclared calibration has no result to report, and the way that goes wrong is a
    contrast quietly computed on libraries nobody looked at.
    """

    manifest = load_manifest(path)
    if not manifest.passes:
        detail = "; ".join(f"{s}: {', '.join(f)}" for s, f in manifest.failures.items())
        raise ValueError(
            f"{manifest.dataset} did not pass its predeclared calibration — "
            f"{detail or 'no length is periodic in every library'}"
        )
    return manifest
