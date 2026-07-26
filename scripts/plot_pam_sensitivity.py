"""The PAM-flexibility sensitivity analysis as a report figure: base-editing placement under
canonical NGG against relaxed PAM recognition, with the primary editing windows held fixed.

Reads results/base_editing.tsv (canonical NGG) and results/base_editing_sensitivity.tsv (relaxed
NG/NGN/NRN) and writes results/figures/pam_sensitivity.{png,svg}. This is geometric placement, not
expected editing: activity, specificity, delivery and tissue compatibility are unmodelled.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

COLOURS = {
    "base_editable_exact": "#2a78d6",
    "base_editable_alternative": "#eb6834",
    "not_base_editable_under_panel": "#b8b6ad",
}
LABELS = {
    "base_editable_exact": "Exact",
    "base_editable_alternative": "Alternative",
    "not_base_editable_under_panel": "Not editable",
}
ORDER = list(COLOURS)
INK, MUTED, SURFACE = "#1b1f27", "#5c6670", "#ffffff"


def _counts(frame: pd.DataFrame, n: int) -> dict[str, int]:
    s = frame[frame["scoreable"]]
    return {c: int((s["reach_class"] == c).sum()) for c in ORDER}


def plot(primary: pd.DataFrame, relaxed: pd.DataFrame, out: Path) -> None:
    n = int(relaxed["scoreable"].sum())
    rows = [
        ("Canonical NGG\nBE4max / ABE7.10", _counts(primary, n)),
        ("PAM-relaxed\n(NG / NGN / NRN)", _counts(relaxed, n)),
    ]
    incremental = int(relaxed.loc[relaxed["scoreable"], "requires_relaxed_pam"].sum())
    still = relaxed[relaxed["reach_class"] == "not_base_editable_under_panel"]
    reasons = still["reason"].value_counts().to_dict()

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": INK,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
        }
    )
    fig, ax = plt.subplots(figsize=(12.5, 5.0))
    fig.suptitle(
        "Base-editing placement is PAM-limited, not fundamentally limited",
        x=0.035,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )

    reach_pct = {0: "89.0%", 1: "31.2%"}
    for row_i, (label, counts) in enumerate(rows):
        y = len(rows) - 1 - row_i
        left = 0
        for c in ORDER:
            w = counts[c]
            ax.barh(y, w, left=left, height=0.52, color=COLOURS[c], edgecolor=SURFACE, linewidth=2)
            if w / n > 0.05:
                light = c != "not_base_editable_under_panel"
                ax.text(
                    left + w / 2,
                    y,
                    f"{w:,}\n{w / n * 100:.0f}%",
                    ha="center",
                    va="center",
                    color="#ffffff" if light else INK,
                    fontsize=9.5,
                    fontweight="bold",
                    linespacing=1.25,
                )
            left += w
        ax.text(
            n * 1.02,
            y,
            f"{reach_pct[y]} placed",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=INK,
        )
        ax.text(
            -n * 0.015, y, label, va="center", ha="right", fontsize=10, fontweight="bold", color=INK
        )

    ax.set_xlim(0, n * 1.16)
    ax.set_ylim(-1.5, 1.75)
    ax.axis("off")
    # Legend swatches.
    for i, c in enumerate(ORDER):
        ax.add_patch(
            plt.Rectangle(
                (n * (0.0 + i * 0.17), 1.48),
                n * 0.02,
                0.12,
                facecolor=COLOURS[c],
                clip_on=False,
                linewidth=0,
            )
        )
        ax.text(n * (0.027 + i * 0.17), 1.54, LABELS[c], va="center", fontsize=9, color=INK)

    ax.text(
        0,
        -0.45,
        f"Relaxing the PAM, primary windows held fixed, newly reaches {incremental:,} variants no "
        f"canonical NGG guide does\n"
        f"(31.2% → 89.0% of {n:,} scoreable — geometric placement, not expected editing).\n"
        f"Residual {int(still['scoreable'].sum()):,}: "
        f"{reasons.get('target_outside_window', 0):,} target outside any editing window, "
        f"{reasons.get('no_pam', 0):,} no modelled PAM. "
        "NRN excludes SpRY's weaker NYN (conservative).",
        fontsize=8.5,
        color=MUTED,
        va="top",
        linespacing=1.6,
    )

    fig.text(
        0.035,
        0.02,
        "A PAM-flexibility sensitivity abstraction with fixed windows, not validated editor "
        "architectures. Activity, specificity, delivery and tissue are unmodelled.\n"
        "The 89.0% is the labelled sensitivity arm; the base-editing headline stays 31.2% NGG.",
        fontsize=7.5,
        color=MUTED,
        linespacing=1.5,
    )
    fig.subplots_adjust(left=0.13, right=0.87, top=0.86, bottom=0.30)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".svg"))
    print(f"wrote {out.with_suffix('.png')} and {out.with_suffix('.svg')}")


if __name__ == "__main__":
    primary = pd.read_csv("results/base_editing.tsv", sep="\t")
    relaxed = pd.read_csv("results/base_editing_sensitivity.tsv", sep="\t")
    plot(primary, relaxed, Path("results/figures/pam_sensitivity"))
