from pathlib import Path

from click.testing import CliRunner

from riborescue.cli import main

DATA = Path(__file__).parents[2] / "pipeline/tests/data"


def test_triage_reports_a_verdict_and_its_reason():
    result = CliRunner().invoke(main, ["triage", "stop_gained"])
    assert result.exit_code == 0
    assert "supported_nonsense: readthrough applies" in result.output


def test_triage_table_writes_a_verdict_per_variant(tmp_path: Path):
    out = tmp_path / "triaged.tsv"
    result = CliRunner().invoke(
        main, ["triage-table", str(DATA / "variants.tsv"), "--out", str(out)]
    )
    assert result.exit_code == 0
    assert "readthrough applies to 2" in result.output
    assert len(out.read_text().splitlines()) == 6


def test_triage_table_refuses_a_table_that_breaks_its_schema(tmp_path: Path):
    bad = tmp_path / "variants.tsv"
    bad.write_text("variant_id\tconsequence\nX:c.1A>T\tstop_gained\n")
    result = CliRunner().invoke(
        main, ["triage-table", str(bad), "--out", str(tmp_path / "out.tsv")]
    )
    assert result.exit_code != 0
    assert "does not match TriageInput" in result.output


def test_validate_labels_counts_the_censored_measurements():
    result = CliRunner().invoke(main, ["validate-labels", str(DATA / "labels.tsv")])
    assert result.exit_code == 0
    assert "3 labels, 1 at the assay ceiling" in result.output


def test_validate_handoff_accepts_the_fixture_manifest():
    result = CliRunner().invoke(main, ["validate-handoff", str(DATA / "handoff.json")])
    assert result.exit_code == 0
    assert "nf-core/riboseq 1.2.0, 7 declared outputs" in result.output


def test_validate_handoff_checks_the_files_it_names():
    args = ["validate-handoff", str(DATA / "handoff.json"), "--check-files"]
    present = CliRunner().invoke(main, [*args, "--results-root", str(DATA / "upstream")])
    absent = CliRunner().invoke(main, [*args, "--results-root", str(DATA)])
    assert present.exit_code == 0
    assert absent.exit_code != 0
    assert "declares outputs that are not present" in absent.output


def test_validate_handoff_refuses_a_malformed_manifest(tmp_path: Path):
    manifest = tmp_path / "handoff.json"
    manifest.write_text('{"pipeline": "nf-core/rnaseq"}')
    result = CliRunner().invoke(main, ["validate-handoff", str(manifest)])
    assert result.exit_code != 0
    assert "not a valid handoff manifest" in result.output
