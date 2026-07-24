"""The compact table the web app reads, built from the results the pipeline already writes.

One row per variant, carrying the best therapy and every therapy's interval, whether the transcript
is expected to escape decay, whether most insertions at the site are tolerated, and which suppressor
tRNA would restore the original residue exactly. The gate statuses ride along at the top, because a
reader has to know what the readthrough control and the safety atlas do and do not cover before
trusting a number on the page.

This is the boundary the frontend depends on. It is small and JSON so the app ships as static files
and never touches a BAM, and it is produced here rather than in the app so the schema has one owner.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

__all__ = [
    "GateStatus",
    "WebTable",
    "build_web_table",
    "safety_summary",
]

# The three-letter stop codons riboWaltz and the triage use, spelled as the suppressor table does.
_STOP_RNA = {"uaa": "UAA", "uag": "UAG", "uga": "UGA", "taa": "UAA", "tag": "UAG", "tga": "UGA"}

# The interval each available therapy carries, keyed by the name it takes in the payload.
_ARM_COLUMNS = {
    "readthrough": "readthrough_predicted",
    "low": "readthrough_low",
    "high": "readthrough_high",
}


@dataclass(frozen=True)
class GateStatus:
    """What the reader must know before trusting a number: which claims are still unbacked."""

    readthrough_control: str = "passed for G418, within one laboratory"
    readthrough_detail: str = (
        "The G418 positive control completes its signature on GSE144140, and SRI-37240 fails it "
        "in the stalling direction — so the assay detects canonical G418 readthrough. This is "
        "confirmation within one laboratory and protocol family, not independent replication, and "
        "it does not license a null about a therapy of different data quality, modality or context."
    )
    safety_atlas: str = "G418 only, in HEK293T"
    safety_detail: str = (
        "Native-stop occupancy is measured for G418 in HEK293T cells only; the other therapies "
        "have no matched empirical safety atlas. It is downstream occupancy past normal stops, not "
        "protein production, toxicity, or a result in any other tissue."
    )


def safety_summary(predicted: pd.DataFrame, therapies: list[str], measured: str = "G418") -> dict:
    """The atlas summary the viewer shows: what is measured, for whom, and how far it reaches.

    Every number here is the measured layer or a rank comparison to it; none is a claim about
    protein, toxicity, or another tissue. The prediction reaches every canonical stop, the
    measurement only the genes the cell line expressed at depth, and the gap between the two counts
    is stated so the small analysed denominator is not mistaken for the whole genome.
    """

    from riborescue.variants.native_stops import concordance

    mane = predicted[predicted["mane_select"]]
    matched = mane.dropna(subset=["predicted_g418", "measured_lift"])
    stats = concordance(matched["predicted_g418"], matched["measured_lift"], draws=1000)
    groups = matched["group"].value_counts()
    return {
        "measured_therapy": measured,
        "cell_line": "HEK293T",
        "canonical_stops_scored": len(mane),
        "analysed": len(matched),
        "unavailable_reason": (
            "the rest were not expressed at footprint depth in this cell line, so they carry a "
            "prediction but no measurement to compare it against"
        ),
        "concordance": {"rho": stats["rho"], "low": stats["low"], "high": stats["high"]},
        "concordance_label": "modest rank concordance, with substantial disagreement",
        "quadrants": {g: int(groups.get(g, 0)) for g in
                      ("both", "predicted only", "measured only", "neither")},
        "per_therapy": {
            t: "native-stop occupancy measured in HEK293T" if t == measured
            else "no matched empirical safety atlas"
            for t in therapies
        },
        "caveat": (
            "Downstream occupancy past normal stops, not protein production, toxicity, or a result "
            "in any other tissue or therapy."
        ),
    }


@dataclass(frozen=True)
class WebTable:
    """The whole payload: the gate statuses, the therapies covered, and the variant rows."""

    status: dict
    therapies: list[str]
    variants: list[dict]
    safety: dict | None = None

    def to_json(self) -> str:
        # allow_nan=False turns a stray NaN or Infinity into an error here rather than the tokens
        # NaN/Infinity, which are not JSON and make the browser's parser reject the whole file.
        return json.dumps(
            {
                "status": self.status,
                "therapies": self.therapies,
                "safety": self.safety,
                "variants": self.variants,
            },
            indent=2,
            allow_nan=False,
        )


def _suppressor(stop_type: str, residue: str) -> dict | None:
    """The suppressor tRNA that reinserts the original residue at this stop, if the stop is known.

    A suppressor charged with the variant's own amino acid restores the protein exactly, which is
    the modality this project ultimately cares about and the one no small molecule can target.
    """

    rna = _STOP_RNA.get(str(stop_type).lower())
    if rna is None or not isinstance(residue, str) or len(residue) != 1:
        return None
    return {"design": f"{rna}-{residue}", "restores_exactly": True}


def build_web_table(
    landscape: pd.DataFrame,
    amenability: pd.DataFrame,
    status: GateStatus | None = None,
    variant_ids: list[str] | None = None,
    safety: dict | None = None,
) -> WebTable:
    """Join the per-variant landscape to the per-therapy intervals into the app's table.

    `variant_ids` narrows the output to a chosen set — the example artifact is a diverse handful,
    the full export is every scoreable variant.
    """

    if variant_ids is not None:
        landscape = landscape[landscape["variant_id"].isin(variant_ids)]
    kept = set(landscape["variant_id"])
    therapy = amenability[amenability["variant_id"].isin(kept)]
    by_variant: dict[object, list[dict]] = {
        vid: rows.to_dict("records") for vid, rows in therapy.groupby("variant_id")
    }
    therapies = sorted(amenability["therapy_id"].unique())

    variants = []
    for row in landscape.to_dict("records"):
        offered = []
        for arm in by_variant.get(row["variant_id"], []):
            present = arm["status"] == "present"
            rounded = (
                {key: round(float(arm[column]), 5) for key, column in _ARM_COLUMNS.items()}
                if present
                else {key: None for key in _ARM_COLUMNS}
            )
            offered.append(
                {
                    "id": arm["therapy_id"],
                    "available": present,
                    **rounded,
                    "reason": None if present else str(arm["reason"]),
                }
            )
        position = row["protein_position"]
        # A variant whose every therapy is unavailable has no best — the landscape leaves it NaN.
        # That must reach the page as null, not as the token NaN, which is not JSON and stops the
        # browser's parser dead on the whole file.
        best = (
            None
            if pd.isna(row["best_readthrough"])
            else {
                "therapy": row["best_therapy"],
                "readthrough": round(float(row["best_readthrough"]), 5),
                "low": round(float(row["best_readthrough_low"]), 5),
            }
        )
        variants.append(
            {
                "id": row["variant_id"],
                "gene": row["gene_symbol"],
                "protein_position": None if pd.isna(position) else int(position),
                "stop": str(row["stop_type"]).upper().replace("T", "U"),
                "residue": row["original_aa"],
                "review_stars": None if pd.isna(row["review_stars"]) else int(row["review_stars"]),
                "escapes_decay": bool(row["escapes_decay_by_rule"]),
                "tolerable_insertion_share": round(float(row["tolerable_insertion_share"]), 3),
                "best": best,
                "therapies": offered,
                "suppressor": _suppressor(str(row["stop_type"]), str(row["original_aa"])),
            }
        )

    return WebTable(
        status=asdict(status or GateStatus()),
        therapies=therapies,
        variants=variants,
        safety=safety,
    )


def diverse_sample(landscape: pd.DataFrame, amenability: pd.DataFrame, size: int = 48) -> list[str]:
    """Variant ids chosen so the example artifact shows every state the app must render.

    A hand-written example drifts from the schema the moment the schema changes. This picks real
    rows instead, and deliberately includes the states a uniform sample would usually miss: a
    therapy that is unavailable for a variant, a site where few insertions are tolerated, one that
    is not expected to escape decay, and a confident readthrough beside an uncertain one.
    """

    wanted: list[str] = []
    missing = amenability.loc[amenability["status"] != "present", "variant_id"].unique()
    wanted += list(missing[:6])

    for column, ascending in [
        ("tolerable_insertion_share", True),
        ("tolerable_insertion_share", False),
        ("best_readthrough", False),
        ("best_readthrough_low", True),
    ]:
        ordered = landscape.sort_values(column, ascending=ascending)["variant_id"]
        wanted += list(ordered.head(6))

    blocked = landscape.loc[~landscape["escapes_decay_by_rule"], "variant_id"]
    wanted += list(blocked.head(6))

    seen: dict[str, None] = {}
    for vid in wanted:
        seen.setdefault(vid, None)
    remaining = [v for v in landscape["variant_id"] if v not in seen]
    for vid in remaining[: max(0, size - len(seen))]:
        seen.setdefault(vid, None)
    return list(seen)


def read_web_inputs(landscape: Path, amenability: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The two result tables the export joins, read with the dtypes the export expects."""

    return (
        pd.read_csv(landscape, sep="\t"),
        pd.read_csv(amenability, sep="\t", dtype={"reason": "string"}),
    )
