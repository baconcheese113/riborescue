"""Where a step writes, and where the next one reads.

Every check here exists because a run reported complete success while its outputs went somewhere
nothing downstream looks. The pipeline's own output directory is relative to wherever Nextflow was
launched, and the reads task launches it from `pipeline/`, so a run without an explicit directory
publishes under `pipeline/results` while calibration reads `results/reads`.
"""

import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]
TASKS = tomllib.loads((ROOT / "pyproject.toml").read_text())["tool"]["pixi"]["tasks"]
SAMPLESHEET = ROOT / "pipeline/assets/riboseq_samples.tsv"


def _command(task: str) -> str:
    definition = TASKS[task]
    return definition["cmd"] if isinstance(definition, dict) else definition


def test_the_reads_workflow_names_the_directory_it_publishes_to():
    """It runs from `pipeline/`, where the pipeline's default lands outside the results tree."""

    reads = TASKS["reads"]
    assert reads["cwd"] == "pipeline"
    assert "--outdir" in reads["cmd"]


def test_the_reads_workflow_publishes_where_the_calibration_reads():
    """The one that was wrong: published under `pipeline/results`, read from `results/reads`."""

    command = _command("reads")
    assert "--outdir" in command, "the reads task publishes wherever it was launched from"
    published = command.split("--outdir")[1].split()[0]
    # The task runs one level down, so its directory resolves from there rather than from the root.
    resolved = (ROOT / "pipeline" / published).resolve()

    driver = (ROOT / "scripts/run_psite_all.sh").read_text()
    default = driver.split("alignments=${3:-")[1].split("}")[0]
    assert (resolved / "alignments").resolve() == (ROOT / default).resolve()


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_calibrating_a_dataset_refuses_when_any_of_its_libraries_is_missing(tmp_path: Path):
    """The alignments of another experiment are not a substitute for this one's.

    A directory holding only the previous run's libraries must stop the dataset, by name, before
    anything is calculated — not calibrate whatever happens to be there.
    """

    alignments = tmp_path / "alignments"
    alignments.mkdir()
    (alignments / "hek293t_untreated_rep1_riboseq.Aligned.toTranscriptome.out.bam").touch()

    staged = ROOT / "results/staged_runs.tsv"
    if not staged.exists():
        pytest.skip("no staged samplesheet in this checkout")

    result = subprocess.run(
        ["bash", "scripts/run_psite_all.sh", "gse144140", str(staged), str(alignments)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "no transcriptome alignment for gse144140_" in result.stderr


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_calibrating_a_dataset_names_only_that_dataset(tmp_path: Path):
    """Two studies profile HEK293T, so the cell line cannot separate one experiment from another."""

    alignments = tmp_path / "alignments"
    alignments.mkdir()
    staged = ROOT / "results/staged_runs.tsv"
    if not staged.exists():
        pytest.skip("no staged samplesheet in this checkout")

    result = subprocess.run(
        ["bash", "scripts/run_psite_all.sh", "gse144140", str(staged), str(alignments)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    named = result.stderr.split(" at ")[0].removeprefix("no transcriptome alignment for ").strip()
    assert named.startswith("gse144140_")
