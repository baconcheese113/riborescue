import json
from pathlib import Path

import pandas as pd
import pytest
from pandera.errors import SchemaErrors

from riborescue.reads import ADAPTER_REACHED_BY, AdapterNotFoundError, summarise_trimming
from riborescue.sequencing import fastq_inputs
from riborescue.tables import SequencingRuns, read_table

SAMPLESHEET = Path(__file__).parents[2] / "pipeline/assets/riboseq_samples.tsv"


def _report(path: Path, *, reads: int, with_adapter: int, output: int | None = None) -> Path:
    path.write_text(
        json.dumps(
            {
                "read_counts": {
                    "input": reads,
                    "output": reads if output is None else output,
                    "read1_with_adapter": with_adapter,
                },
                "basepair_counts": {"input": reads * 50, "output": reads * 30},
            }
        )
    )
    return path


def test_the_declared_runs_validate():
    runs = SequencingRuns.validate(read_table(SAMPLESHEET), lazy=True)
    assert len(runs) == 12
    assert set(runs["assay"]) == {"riboseq", "rnaseq"}
    assert (runs.groupby(["cell_line", "treatment"]).size() > 0).all()


def test_every_paired_run_names_a_second_read_and_every_single_run_does_not():
    runs = read_table(SAMPLESHEET)
    paired = runs["layout"] == "paired"
    assert runs.loc[paired, "fastq_2_url"].notna().all()
    assert runs.loc[~paired, "fastq_2_url"].isna().all()


def test_a_paired_run_missing_its_second_read_is_refused():
    runs = read_table(SAMPLESHEET)
    runs.loc[runs["layout"] == "paired", ["fastq_2_url", "fastq_2_md5"]] = None
    with pytest.raises(SchemaErrors):
        SequencingRuns.validate(runs, lazy=True)


def test_an_overlap_shorter_than_the_degenerate_prefix_is_refused():
    """Six leading Ns match any six bases, so such an adapter is found at every read's end."""

    runs = read_table(SAMPLESHEET)
    runs.loc[runs["adapter_3p"].str.startswith("N"), "adapter_overlap"] = 6
    with pytest.raises(SchemaErrors):
        SequencingRuns.validate(runs, lazy=True)


def test_each_run_yields_one_fastq_input_per_read():
    runs = read_table(SAMPLESHEET)
    for _, run in runs.iterrows():
        expected = 2 if run["layout"] == "paired" else 1
        declared = dict(fastq_inputs(run))
        assert len(declared) == expected
        assert all(i.url.endswith(".fastq.gz") and len(i.md5) == 32 for i in declared.values())


def test_trimming_summary_reports_raw_against_cleaned(tmp_path: Path):
    _report(tmp_path / "a.cutadapt.json", reads=1000, with_adapter=990, output=900)
    summary = summarise_trimming([tmp_path / "a.cutadapt.json"])
    assert summary.loc[0, "reads_raw"] == 1000
    assert summary.loc[0, "reads_cleaned"] == 900
    assert summary.loc[0, "reads_retained"] == pytest.approx(0.9)
    assert summary.loc[0, "adapter_rate"] == pytest.approx(0.99)


def test_an_adapter_absent_from_a_footprint_library_is_refused(tmp_path: Path):
    """The wrong adapter trims nothing and leaves every read carrying linker into alignment."""

    _report(tmp_path / "b.cutadapt.json", reads=1000, with_adapter=3)
    with pytest.raises(AdapterNotFoundError, match="under 50%"):
        summarise_trimming([tmp_path / "b.cutadapt.json"], {"b"})


def test_a_transcriptome_library_may_rarely_reach_its_adapter(tmp_path: Path):
    """Most fragments are longer than the read, so a low rate there is not evidence of a defect."""

    _report(tmp_path / "b.cutadapt.json", reads=1000, with_adapter=32)
    summary = summarise_trimming([tmp_path / "b.cutadapt.json"])
    assert summary.loc[0, "adapter_rate"] == pytest.approx(0.032)


def test_the_assays_whose_reads_must_reach_the_adapter_are_the_footprint_ones():
    runs = read_table(SAMPLESHEET)
    expected = set(runs.loc[runs["assay"].isin(ADAPTER_REACHED_BY), "sample"])
    assert expected == set(runs.loc[runs["assay"] == "riboseq", "sample"])
    assert len(expected) == 8


def test_summaries_are_ordered_by_sample_whatever_order_the_reports_arrive_in(tmp_path: Path):
    reports = [
        _report(tmp_path / f"{name}.cutadapt.json", reads=100, with_adapter=100)
        for name in ("zebra", "alpha", "middle")
    ]
    summary = summarise_trimming(list(reversed(reports)))
    assert list(summary["sample"]) == ["alpha", "middle", "zebra"]


def test_the_declared_read_counts_match_the_archive():
    """ENA's own per-run counts, carried so a truncated transfer is visible before alignment."""

    runs = read_table(SAMPLESHEET)
    assert runs["read_count"].min() > 10_000_000
    assert pd.api.types.is_integer_dtype(runs["read_count"])
