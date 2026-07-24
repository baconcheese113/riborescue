import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from riborescue.core.handoff import UpstreamHandoff

FIXTURE = Path(__file__).parents[2] / "pipeline/tests/data/handoff.json"


def _manifest(**overrides: object) -> dict[str, object]:
    return json.loads(FIXTURE.read_text()) | overrides


def test_the_committed_fixture_manifest_is_valid():
    handoff = UpstreamHandoff.from_json(FIXTURE)
    assert handoff.revision == "1.2.0"
    assert len(tuple(handoff.outputs())) == 7


def test_every_declared_output_is_resolved_against_the_results_root():
    handoff = UpstreamHandoff.model_validate(_manifest(results_root="/somewhere/results"))
    assert all(str(path).startswith("/somewhere/results/") for _, path in handoff.outputs())


def test_the_fixture_tree_is_complete():
    root = FIXTURE.parent / "upstream"
    handoff = UpstreamHandoff.from_json(FIXTURE).model_copy(update={"results_root": root})
    assert handoff.missing() == ()


def test_an_absent_output_is_reported_by_name():
    handoff = UpstreamHandoff.from_json(FIXTURE).model_copy(
        update={"results_root": Path("nowhere")}
    )
    assert {name for name, _ in handoff.missing()} == {
        "psite_offsets",
        "codon_coverage",
        "cds_coverage",
        "rnaseq_counts",
        "rnaseq_tpm",
        "alignments",
        "multiqc",
    }


@pytest.mark.parametrize("revision", ["dev", "main", "master", "latest", "  "])
def test_a_moving_upstream_revision_is_refused(revision: str):
    with pytest.raises(ValidationError):
        UpstreamHandoff.model_validate(_manifest(revision=revision))


def test_another_pipeline_is_refused():
    with pytest.raises(ValidationError):
        UpstreamHandoff.model_validate(_manifest(pipeline="nf-core/rnaseq"))


@pytest.mark.parametrize("escape", ["/etc/passwd", "../../elsewhere/counts.tsv"])
def test_an_output_outside_the_results_tree_is_refused(escape: str):
    with pytest.raises(ValidationError):
        UpstreamHandoff.model_validate(_manifest(rnaseq_counts=escape))


def test_an_unknown_field_is_refused():
    with pytest.raises(ValidationError):
        UpstreamHandoff.model_validate(_manifest(bigwigs="coverage.bw"))
