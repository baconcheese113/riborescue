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
        "cds_total": 1500,
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

    manifest = _passing_manifest(tmp_path, "hek293t", [r["sample"] for r in rows])
    result = CliRunner().invoke(
        main,
        [
            "readthrough", str(counts), "--gtf", str(gtf), "--samplesheet", str(sheet),
            "--treated", "g418", "--control", "untreated", "--dataset", "hek293t",
            "--manifest", str(manifest), "--out", str(tmp_path / "o.tsv"),
        ],
    )
    assert result.exit_code != 0
    assert "g418_rep2" in result.output


def _passing_manifest(tmp_path: Path, dataset: str, samples: list[str]) -> Path:
    import json

    path = tmp_path / f"{dataset}.calibration.json"
    path.write_text(
        json.dumps(
            {
                "dataset": dataset,
                "lengths": [30, 31],
                "surveyed": [18, 40],
                "script_md5": "abc",
                "passes": True,
                "libraries": [
                    {
                        "sample": sample,
                        "psites": 2_000_000,
                        "frame0_share": 0.55,
                        "dominant_length": 30,
                        "offset_from_5": 12,
                        "failures": [],
                    }
                    for sample in samples
                ],
            }
        )
    )
    return path


def test_the_assay_refuses_a_dataset_that_failed_its_calibration(tmp_path: Path):
    """A dataset whose libraries did not pass has no result, and the assay must not produce one."""

    import json

    manifest = tmp_path / "c.json"
    manifest.write_text(
        json.dumps(
            {
                "dataset": "gse144140",
                "lengths": [30],
                "surveyed": [18, 40],
                "script_md5": None,
                "passes": False,
                "libraries": [
                    {
                        "sample": "gse144140_dmso_rep1_riboseq",
                        "psites": 9_000,
                        "frame0_share": 0.51,
                        "dominant_length": 30,
                        "offset_from_5": 12,
                        "failures": ["9,000 P-sites over the selected lengths, under 1,000,000"],
                    }
                ],
            }
        )
    )
    counts = tmp_path / "counts.tsv"
    counts.write_text("sample\n")
    gtf = tmp_path / "a.gtf.gz"
    gtf.write_bytes(b"")
    sheet = tmp_path / "runs.tsv"
    sheet.write_text("dataset\n")

    result = CliRunner().invoke(
        main,
        [
            "readthrough", str(counts), "--gtf", str(gtf), "--samplesheet", str(sheet),
            "--treated", "g418", "--control", "dmso", "--dataset", "gse144140",
            "--manifest", str(manifest), "--out", str(tmp_path / "o.tsv"),
        ],
    )
    assert result.exit_code != 0
    assert "did not pass its predeclared calibration" in result.output


def test_each_arm_reports_the_window_it_actually_used(tmp_path: Path):
    """The sensitivity arm printed the manifest's lengths while summing a different set."""

    import gzip

    import pandas as pd

    sheet = tmp_path / "runs.tsv"
    rows = [
        {
            "dataset": "d", "sample": f"{arm}_rep{n}", "run_accession": f"SRR{arm[:2]}{n}",
            "cell_line": "HEK293T", "variant": "TP53:p.Arg196Ter", "treatment": arm, "replicate": n,
            "assay": "riboseq", "layout": "single", "read_count": 1000, "adapter_3p": "ACGT",
            "adapter_3p_2": None, "adapter_overlap": 5, "cut_5p": 0,
            "fastq_1": "a.fastq.gz", "fastq_2": None,
        }
        for arm in ("dmso", "g418")
        for n in (1, 2)
    ]
    pd.DataFrame(rows).to_csv(sheet, sep="\t", index=False)

    counts = tmp_path / "counts.tsv"
    pd.DataFrame(
        [
            {**_counts_row("T1", r["sample"]), "length": length}
            for r in rows
            for length in (28, 30, 31, 34)
        ]
    ).to_csv(counts, sep="\t", index=False)

    gtf = tmp_path / "a.gtf.gz"
    with gzip.open(gtf, "wt") as handle:
        handle.write('chr1\tX\ttranscript\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "T1";\n')

    manifest = _passing_manifest(tmp_path, "d", [r["sample"] for r in rows])
    common = [
        "readthrough", str(counts), "--gtf", str(gtf), "--samplesheet", str(sheet),
        "--treated", "g418", "--control", "dmso", "--dataset", "d",
        "--manifest", str(manifest), "--out", str(tmp_path / "o.tsv"),
    ]
    selected = CliRunner().invoke(main, common)
    published = CliRunner().invoke(main, [*common, "--published-lengths", "28", "35"])

    assert "selected set: 30, 31 nt" in selected.output
    assert "published window: 28, 29, 30, 31, 32, 33, 34, 35 nt" in published.output
