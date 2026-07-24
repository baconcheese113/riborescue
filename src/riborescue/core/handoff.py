"""The upstream handoff manifest — the only door between nf-core/riboseq and RiboRescue.

Downstream code reads the outputs this manifest names and nothing else. Every path is declared
relative to the results root, so no analysis step can reach into an arbitrary location in someone
else's results tree, and the upstream revision is a pin rather than a moving branch.
"""

from collections.abc import Iterator
from pathlib import Path

from pydantic import BaseModel, ConfigDict, field_validator

__all__ = [
    "UPSTREAM_PIPELINE",
    "HandoffOutput",
    "UpstreamHandoff",
]

UPSTREAM_PIPELINE = "nf-core/riboseq"

_MOVING_REVISIONS = frozenset({"dev", "main", "master", "latest"})

type HandoffOutput = tuple[str, Path]


class UpstreamHandoff(BaseModel):
    """The nf-core/riboseq outputs RiboRescue consumes, named once and validated before use."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    pipeline: str
    revision: str
    results_root: Path
    psite_offsets: Path
    codon_coverage: Path
    cds_coverage: Path
    rnaseq_counts: Path
    rnaseq_tpm: Path
    alignments: tuple[Path, ...]
    multiqc: Path

    @field_validator("pipeline")
    @classmethod
    def _is_the_upstream_pipeline(cls, value: str) -> str:
        if value != UPSTREAM_PIPELINE:
            raise ValueError(f"handoff accepts {UPSTREAM_PIPELINE} only, not {value!r}")
        return value

    @field_validator("revision")
    @classmethod
    def _is_pinned(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("revision must name a pinned upstream release")
        if value.lower() in _MOVING_REVISIONS:
            raise ValueError(f"revision {value!r} moves; pin a release tag or commit instead")
        return value

    @field_validator(
        "psite_offsets",
        "codon_coverage",
        "cds_coverage",
        "rnaseq_counts",
        "rnaseq_tpm",
        "alignments",
        "multiqc",
    )
    @classmethod
    def _stays_inside_the_results_tree(
        cls, value: Path | tuple[Path, ...]
    ) -> Path | tuple[Path, ...]:
        for path in value if isinstance(value, tuple) else (value,):
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(
                    f"declared output {str(path)!r} must be relative to the results root"
                )
        return value

    @classmethod
    def from_json(cls, path: Path) -> "UpstreamHandoff":
        return cls.model_validate_json(path.read_text())

    def outputs(self) -> Iterator[HandoffOutput]:
        """Yield each declared output as a field name and a path under the results root."""

        for name in self.__class__.model_fields:
            value = getattr(self, name)
            if name in {"pipeline", "revision", "results_root"}:
                continue
            for path in value if isinstance(value, tuple) else (value,):
                yield name, self.results_root / path

    def missing(self) -> tuple[HandoffOutput, ...]:
        """Return the declared outputs that are not on disk."""

        return tuple((name, path) for name, path in self.outputs() if not path.exists())
