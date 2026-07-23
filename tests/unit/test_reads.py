import json
from pathlib import Path

import pandas as pd
import pytest
from pandera.errors import SchemaErrors

from riborescue.reads import (
    ADAPTER_REACHED_BY,
    AdapterNotFoundError,
    summarise_alignment,
    summarise_trimming,
)
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
    assert len(runs) == 21
    assert set(runs["assay"]) == {"riboseq", "rnaseq"}
    assert (runs.groupby(["dataset", "treatment"]).size() > 0).all()


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
    assert len(expected) == 17


def test_summaries_are_ordered_by_sample_whatever_order_the_reports_arrive_in(tmp_path: Path):
    reports = [
        _report(tmp_path / f"{name}.cutadapt.json", reads=100, with_adapter=100)
        for name in ("zebra", "alpha", "middle")
    ]
    summary = summarise_trimming(list(reversed(reports)))
    assert list(summary["sample"]) == ["alpha", "middle", "zebra"]


STAR_LOG = """\
                                 Started job on |\tJul 21 01:02:03
                             Number of input reads |\t25344634
                         Average input read length |\t31
                                      UNIQUE READS:
                      Uniquely mapped reads number |\t18500000
                           Uniquely mapped reads % |\t73.00%
                             Average mapped length |\t30.55
                                MULTI-MAPPING READS:
        % of reads mapped to multiple loci |\t14.00%
        % of reads mapped to too many loci |\t3.00%
                                     UNMAPPED READS:
                  % of reads unmapped: too short |\t9.50%
"""


def test_alignment_metrics_come_off_the_star_log(tmp_path: Path):
    log = tmp_path / "calu6_untreated_rep1_riboseq.Log.final.out"
    log.write_text(STAR_LOG)
    summary = summarise_alignment([log])
    assert summary.loc[0, "sample"] == "calu6_untreated_rep1_riboseq"
    assert summary.loc[0, "reads_input"] == 25344634
    assert summary.loc[0, "unique_rate"] == pytest.approx(73.0)
    assert summary.loc[0, "multimapped_rate"] == pytest.approx(14.0)
    assert summary.loc[0, "mapped_length_mean"] == pytest.approx(30.55)


def test_the_mapped_share_counts_every_way_a_read_can_map(tmp_path: Path):
    """Uniquely, to several loci, and to too many — but never the reads that did not map."""

    log = tmp_path / "s.Log.final.out"
    log.write_text(STAR_LOG)
    summary = summarise_alignment([log])
    assert summary.loc[0, "mapped_rate"] == pytest.approx(90.0)


def test_the_declared_read_counts_match_the_archive():
    """ENA's own per-run counts, carried so a truncated transfer is visible before alignment."""

    runs = read_table(SAMPLESHEET)
    assert runs["read_count"].min() > 10_000_000
    assert pd.api.types.is_integer_dtype(runs["read_count"])
