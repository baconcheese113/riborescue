import gzip
import json
from pathlib import Path

import pandas as pd
import pytest
from pandera.errors import SchemaErrors

from riborescue.core.tables import SequencingRuns, read_table
from riborescue.riboseq.reads import (
    ADAPTER_REACHED_BY,
    AdapterNotFoundError,
    summarise_alignment,
    summarise_trimming,
    survey_adapter,
    survey_adapters,
)
from riborescue.riboseq.sequencing import fastq_inputs

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
    assert len(runs) == 31
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
    assert len(expected) == 27


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


TRUSEQ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
LINKER = "CTGTAGGCACCATCAAT"


def _library(path: Path, inserts: list[str], trailing: str = LINKER + TRUSEQ) -> Path:
    """A footprint library as it arrives: insert, linker, sequencing adapter, padded to length."""

    with gzip.open(path, "wt") as handle:
        for index, insert in enumerate(inserts):
            read = (insert + trailing)[:80].ljust(80, "A")
            handle.write(f"@r{index}\n{read}\n+\n{'I' * 80}\n")
    return path


def test_the_survey_finds_the_declared_linker_and_measures_the_footprint(tmp_path: Path):
    inserts = ["A" * length for length in (26, 28, 30, 32)] * 25
    survey = survey_adapter("a", _library(tmp_path / "a.fastq.gz", inserts), LINKER)
    assert survey.adapter_rate == 1.0
    assert survey.footprint_median == 29
    assert (survey.footprint_p10, survey.footprint_p90) == (26, 32)


def test_a_degenerate_prefix_is_not_searched_for(tmp_path: Path):
    """Its bases are random by construction, so only the fixed part of the linker identifies it."""

    declared = "NNNNNNCACTCGGGCACCAAGGAC"
    inserts = ["C" * 30] * 100
    library = _library(tmp_path / "a.fastq.gz", inserts, trailing="TTAGCA" + declared.lstrip("N"))
    assert survey_adapter("a", library, declared).adapter_rate == 1.0


def test_two_chemistries_around_the_same_footprint_measure_the_same(tmp_path: Path):
    """One puts ten bases either side of the footprint and one puts none. That is not biology."""

    degenerate = "NNNNNNCACTCGGGCACCAAGGAC"
    bare = _library(tmp_path / "bare.fastq.gz", ["A" * 30] * 100)
    padded = _library(
        tmp_path / "padded.fastq.gz",
        ["TTTT" + "A" * 30] * 100,
        trailing="GACTGA" + degenerate.lstrip("N"),
    )
    assert survey_adapter("bare", bare, LINKER).footprint_median == 30
    assert survey_adapter("padded", padded, degenerate, cut_5p=4).footprint_median == 30


def test_a_linker_the_series_does_not_name_is_caught_before_alignment(tmp_path: Path):
    """GSE179274 records only that Trimmomatic ran, and its reads carry a linker in front of that.

    The sequencing adapter really is there, so its presence proves nothing. Trimming to it would
    leave the linker on every read, and that shows up as a footprint far too long to be a ribosome.
    """

    runs = pd.DataFrame(
        {
            "sample": ["a"],
            "assay": ["riboseq"],
            "adapter_3p": [TRUSEQ],
            "cut_5p": [0],
            "adapter_overlap": [5],
            "fastq_1": [str(_library(tmp_path / "a.fastq.gz", ["G" * 30] * 100))],
        }
    )
    with pytest.raises(AdapterNotFoundError, match="not 20-40 nt"):
        survey_adapters(runs)


def test_an_adapter_absent_from_the_reads_entirely_is_refused(tmp_path: Path):
    runs = pd.DataFrame(
        {
            "sample": ["a"],
            "assay": ["riboseq"],
            "adapter_3p": ["GGGGGGGGGGGGGGGGGG"],
            "cut_5p": [0],
            "adapter_overlap": [8],
            "fastq_1": [str(_library(tmp_path / "a.fastq.gz", ["ACAC" * 8] * 100))],
        }
    )
    with pytest.raises(AdapterNotFoundError, match="under 50%"):
        survey_adapters(runs)


def test_a_transcriptome_library_is_not_surveyed(tmp_path: Path):
    """Most of its fragments are longer than the read, so it never reaches the adapter."""

    runs = pd.DataFrame(
        {
            "sample": ["a"],
            "assay": ["rnaseq"],
            "adapter_3p": [TRUSEQ],
            "cut_5p": [0],
            "adapter_overlap": [5],
            "fastq_1": ["absent.fastq.gz"],
        }
    )
    assert survey_adapters(runs).empty


def test_no_two_declared_runs_are_the_same_bytes():
    """A dataset that confirms another must not be partly the same libraries under new names."""

    runs = read_table(SAMPLESHEET)
    assert runs["fastq_1_md5"].is_unique
    assert runs["run_accession"].is_unique
    second = runs["fastq_2_md5"].dropna()
    assert second.is_unique
    assert not set(second) & set(runs["fastq_1_md5"])


def test_a_read_too_short_to_hold_the_whole_adapter_still_carries_it(tmp_path: Path):
    """A 50 nt read runs out mid-linker, so only its first bases are there — at the read's end."""

    truncated = tmp_path / "short.fastq.gz"
    with gzip.open(truncated, "wt") as handle:
        for index in range(100):
            read = ("A" * 30 + LINKER)[:38]
            handle.write(f"@r{index}\n{read}\n+\n{'I' * len(read)}\n")
    survey = survey_adapter("short", truncated, LINKER, adapter_overlap=5)
    assert survey.adapter_rate == 1.0
    assert survey.footprint_median == 30


def test_too_little_of_the_adapter_to_be_sure_is_not_counted(tmp_path: Path):
    """Four bases of linker at a read's end match by chance often enough to mean nothing."""

    barely = tmp_path / "barely.fastq.gz"
    with gzip.open(barely, "wt") as handle:
        for index in range(100):
            read = "A" * 30 + LINKER[:4]
            handle.write(f"@r{index}\n{read}\n+\n{'I' * len(read)}\n")
    assert survey_adapter("barely", barely, LINKER, adapter_overlap=8).adapter_rate == 0.0
