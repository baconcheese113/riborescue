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


def _counts_row(transcript: str, sample: str) -> dict:
    return {
        "transcript": transcript,
        "sample": sample,
        "cds_frame0": 1000,
        "cds_frame1": 250,
        "cds_frame2": 250,
        "termination": 50,
        "extension": 300,
        "extension_frame0": 10,
        "extension_frame1": 3,
        "extension_frame2": 2,
        "l_cds": 900,
        "l_utr3": 500,
    }


def test_readthrough_refuses_a_contrast_missing_a_library(tmp_path: Path):
    """A named library absent from the counts would shrink the comparison without a word."""

    import gzip

    import pandas as pd

    sheet = tmp_path / "runs.tsv"
    rows = []
    for treatment in ("untreated", "g418"):
        for replicate in (1, 2):
            rows.append(
                {
                    "dataset": "hek293t",
                    "sample": f"{treatment}_rep{replicate}",
                    "run_accession": f"SRR{treatment[:2]}{replicate}",
                    "cell_line": "HEK293T",
                    "variant": "TP53:p.Arg196Ter",
                    "treatment": treatment,
                    "replicate": replicate,
                    "assay": "riboseq",
                    "layout": "single",
                    "read_count": 1000,
                    "adapter_3p": "ACGT",
                    "adapter_3p_2": None,
                    "adapter_overlap": 5,
                    "cut_5p": 0,
                    "fastq_1": "a.fastq.gz",
                    "fastq_2": None,
                }
            )
    pd.DataFrame(rows).to_csv(sheet, sep="\t", index=False)

    # one library is simply not in the counts
    counts = tmp_path / "counts.tsv"
    present = [r["sample"] for r in rows if r["sample"] != "g418_rep2"]
    pd.DataFrame([_counts_row("T1", s) for s in present]).to_csv(counts, sep="\t", index=False)

    gtf = tmp_path / "a.gtf.gz"
    with gzip.open(gtf, "wt") as handle:
        handle.write('chr1\tX\ttranscript\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "T1";\n')

    result = CliRunner().invoke(
        main,
        [
            "readthrough", str(counts), "--gtf", str(gtf), "--samplesheet", str(sheet),
            "--treated", "g418", "--control", "untreated", "--dataset", "hek293t",
            "--out", str(tmp_path / "o.tsv"),
        ],
    )
    assert result.exit_code != 0
    assert "g418_rep2" in result.output
