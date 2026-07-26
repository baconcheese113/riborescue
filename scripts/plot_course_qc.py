#!/usr/bin/env python
"""Course QC, cleaning, and alignment figures for the matched RNA-seq arm."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_INK = "#1c1c1c"
_MUTED = "#5a5a5a"
_BLUE = "#2166ac"
_LIGHT_BLUE = "#92c5de"
_GOLD = "#d8a031"
_GREY = "#d9d9d9"


def _label(sample: str) -> str:
    parts = sample.removeprefix("calu6_").removesuffix("_rnaseq").split("_")
    return f"{parts[0].upper()}\n{parts[1]}" if len(parts) == 2 else sample


def _mean_stat(stats: pd.DataFrame, prefix: str, column: str) -> float:
    rows = stats.loc[stats["Sample"].str.startswith(f"{prefix}_"), column].dropna()
    if rows.empty:
        raise ValueError(f"MultiQC names no {column} rows for {prefix}")
    return float(rows.mean())


def course_summary(
    samples: pd.DataFrame,
    stats: pd.DataFrame,
    trimming: pd.DataFrame,
) -> pd.DataFrame:
    """One row per Calu-6 RNA-seq library, joining raw, cleaned, and FastQC summaries."""

    runs = samples.loc[(samples["dataset"] == "calu6") & (samples["assay"] == "rnaseq")]
    trimmed = trimming.set_index("sample")
    rows: list[dict[str, str | int | float]] = []
    for run in runs.itertuples(index=False):
        sample = str(run.sample)
        accession = str(run.run_accession)
        clean = trimmed.loc[sample]
        rows.append(
            {
                "sample": sample,
                "run_accession": accession,
                "reads_raw": int(clean["reads_raw"]),
                "reads_cleaned": int(clean["reads_cleaned"]),
                "reads_retained_percent": 100 * float(clean["reads_retained"]),
                "adapter_percent": 100 * float(clean["adapter_rate"]),
                "gc_raw_percent": _mean_stat(stats, accession, "fastqc-percent_gc"),
                "gc_cleaned_percent": _mean_stat(stats, sample, "fastqc-percent_gc"),
                "duplicates_raw_percent": _mean_stat(stats, accession, "fastqc-percent_duplicates"),
                "duplicates_cleaned_percent": _mean_stat(
                    stats, sample, "fastqc-percent_duplicates"
                ),
            }
        )
    return pd.DataFrame(rows).round(3)


def _style(axes: plt.Axes, *, percent: bool = False) -> None:
    axes.spines[["top", "right"]].set_visible(False)
    axes.tick_params(length=0, colors=_MUTED)
    axes.grid(axis="y", color="#ececec", linewidth=0.8)
    axes.set_axisbelow(True)
    if percent:
        axes.set_ylim(0, 100)


def plot_qc(summary: pd.DataFrame, out: Path) -> None:
    """Raw GC content and duplication, the two preparation-level QC summaries."""

    labels = [_label(sample) for sample in summary["sample"]]
    x = np.arange(len(labels))
    figure, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharex=True)
    axes[0].bar(x, summary["gc_raw_percent"], color=_BLUE, width=0.7)
    axes[0].set_title("GC content", loc="left", color=_INK)
    axes[0].set_ylabel("Raw reads (%)", color=_MUTED)
    axes[1].bar(x, summary["duplicates_raw_percent"], color=_GOLD, width=0.7)
    axes[1].set_title("Sequence duplication", loc="left", color=_INK)
    for axis in axes:
        axis.set_xticks(x, labels)
        _style(axis, percent=True)
    figure.suptitle("Raw RNA-seq quality control · Calu-6", x=0.06, ha="left", color=_INK)
    figure.tight_layout()
    figure.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_cleaning(summary: pd.DataFrame, out: Path) -> None:
    """Raw against cleaned read counts and duplication."""

    labels = [_label(sample) for sample in summary["sample"]]
    x = np.arange(len(labels))
    width = 0.36
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharex=True)
    axes[0].bar(x - width / 2, summary["reads_raw"] / 1e6, width, label="Raw", color=_GREY)
    axes[0].bar(
        x + width / 2,
        summary["reads_cleaned"] / 1e6,
        width,
        label="Cleaned",
        color=_BLUE,
    )
    axes[0].set_title("Reads retained", loc="left", color=_INK)
    axes[0].set_ylabel("Read pairs (millions)", color=_MUTED)
    axes[1].bar(
        x - width / 2,
        summary["duplicates_raw_percent"],
        width,
        label="Raw",
        color=_GREY,
    )
    axes[1].bar(
        x + width / 2,
        summary["duplicates_cleaned_percent"],
        width,
        label="Cleaned",
        color=_GOLD,
    )
    axes[1].set_title("Sequence duplication", loc="left", color=_INK)
    axes[1].set_ylabel("Reads (%)", color=_MUTED)
    axes[1].set_ylim(0, 100)
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.legend(frameon=False, ncol=2, loc="upper left")
        _style(axis)
    figure.suptitle("RNA-seq cleaning · Calu-6", x=0.06, ha="left", color=_INK)
    figure.tight_layout()
    figure.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_alignment(alignment: pd.DataFrame, out: Path) -> None:
    """Unique, multimapped, and unmapped shares of each RNA-seq library."""

    table = alignment.loc[alignment["sample"].str.match(r"calu6_.*_rnaseq")].copy()
    labels = [_label(sample) for sample in table["sample"]]
    unique = table["unique_rate"].to_numpy()
    multiple = table["multimapped_rate"].to_numpy()
    unmapped = 100 - unique - multiple
    x = np.arange(len(labels))
    figure, axes = plt.subplots(figsize=(7.5, 4.4))
    axes.bar(x, unique, color=_BLUE, label="Unique")
    axes.bar(x, multiple, bottom=unique, color=_LIGHT_BLUE, label="Multimapped")
    axes.bar(x, unmapped, bottom=unique + multiple, color=_GREY, label="Unmapped")
    axes.set_xticks(x, labels)
    axes.set_ylabel("Input reads (%)", color=_MUTED)
    axes.set_title("STAR alignment · Calu-6 RNA-seq", loc="left", color=_INK)
    axes.legend(frameon=False, ncol=3, loc="upper left")
    _style(axes, percent=True)
    figure.tight_layout()
    figure.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--multiqc", required=True, type=Path)
    parser.add_argument("--trimming", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = course_summary(
        pd.read_csv(args.samples, sep="\t"),
        pd.read_csv(args.multiqc, sep="\t"),
        pd.read_csv(args.trimming, sep="\t"),
    )
    summary.to_csv(args.outdir / "course_qc_summary.tsv", sep="\t", index=False)
    plot_qc(summary, args.outdir / "qc_overview.png")
    plot_cleaning(summary, args.outdir / "pre_post_cleaning.png")
    plot_alignment(
        pd.read_csv(args.alignment, sep="\t"),
        args.outdir / "alignment_metrics.png",
    )


if __name__ == "__main__":
    main()
