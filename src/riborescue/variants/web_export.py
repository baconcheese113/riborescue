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
from typing import Any, cast

import pandas as pd

from riborescue.core.contracts import CONTRACTS_VERSION
from riborescue.core.sequences import STOP_CODONS, as_dna, as_rna

__all__ = [
    "THERAPY_NAMES",
    "GateStatus",
    "WebTable",
    "build_web_table",
    "safety_summary",
    "therapy_name",
]

# The compound each therapy id names. The ids are the readthrough dataset's own labels, and one of
# them is dangerously short: "SRI" there is **SRI-41315**, while the compound in the safety atlas's
# Ribo-seq dataset is **SRI-37240** — two different molecules from the same series, which a reader
# seeing "SRI" beside a sentence about SRI-37240 would otherwise take for one.
THERAPY_NAMES = {
    "CC90009": "CC-90009",
    "Clitocine": "Clitocine",
    "DAP": "2,6-diaminopurine",
    "G418": "G418 (geneticin)",
    "SJ6986": "SJ6986",
    "SRI": "SRI-41315",
}


def therapy_name(therapy_id: str) -> str:
    """The compound a therapy id names, or the id itself where no fuller name is recorded."""

    return THERAPY_NAMES.get(therapy_id, therapy_id)


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
        "in the stalling direction — so the assay detects canonical G418 readthrough. SRI-37240 is "
        "the control compound in that Ribo-seq dataset and is not the SRI-41315 scored here; they "
        "are different molecules. This is confirmation within one laboratory and protocol family, "
        "not independent replication, and it does not license a null about a therapy of different "
        "data quality, modality or context."
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

    from riborescue.variants.native_stop_predictions import concordance

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
        "quadrants": {
            g: int(groups.get(g, 0)) for g in ("both", "predicted only", "measured only", "neither")
        },
        "per_therapy": {
            t: "native-stop occupancy measured in HEK293T"
            if t == measured
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
                "contracts_version": CONTRACTS_VERSION,
                "status": self.status,
                "therapies": self.therapies,
                "therapy_names": {t: therapy_name(t) for t in self.therapies},
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

    if as_dna(str(stop_type)) not in STOP_CODONS:
        return None
    if not isinstance(residue, str) or len(residue) != 1:
        return None
    return {"design": f"{as_rna(str(stop_type)).upper()}-{residue}", "restores_exactly": True}


def _nmd_by_variant(nmd: pd.DataFrame | None) -> dict[object, dict]:
    """Each variant's two NMD verdicts and the rules behind them, keyed for the per-variant join."""

    if nmd is None:
        return {}
    return {
        row.variant_id: {
            "escape_guideline": bool(row.escape_guideline),
            "escape_full_rules": bool(row.escape_full_rules),
            "disagree": bool(row.predictors_disagree),
            "rules": {
                "last_exon": bool(row.rule_last_exon),
                "within_last_junction": bool(row.rule_within_last_junction),
                "start_proximal": bool(row.rule_start_proximal),
                "long_exon": bool(row.rule_long_exon),
            },
        }
        for row in nmd.itertuples()
    }


def _aenmd_by_variant(aenmd: pd.DataFrame | None) -> dict[object, dict]:
    """Each variant's aenmd verdict — available and escape, or the reason it was not scored."""

    if aenmd is None:
        return {}
    return {
        row.variant_id: {
            "available": bool(row.aenmd_available),
            "escape": bool(row.aenmd_escape),
            "reason": str(row.reason) if not bool(row.aenmd_available) else "",
        }
        for row in aenmd.itertuples()
    }


def _nmdetective_by_variant(nmdetective: pd.DataFrame | None) -> dict[object, dict]:
    """Each variant's NMDetective-AI efficiency score, or the reason the model did not place it."""

    if nmdetective is None:
        return {}
    verdicts = {}
    for row in nmdetective.itertuples():
        available = bool(row.nmd_available)
        # itertuples types a scalar as pandas' own union (which includes complex); the score is a
        # plain float, so cast before rounding.
        efficiency = cast(Any, row.nmd_efficiency)
        verdicts[row.variant_id] = {
            "available": available,
            "efficiency": round(float(efficiency), 4)
            if available and pd.notna(efficiency)
            else None,
            "reason": str(row.reason) if not available else "",
        }
    return verdicts


def _with_models(nmd: dict | None, aenmd: dict | None, nmdetective: dict | None) -> dict | None:
    """Attach the model-tier verdicts to a variant's rule-tier NMD object, when it has one."""

    return None if nmd is None else {**nmd, "aenmd": aenmd, "nmdetective": nmdetective}


def build_web_table(
    landscape: pd.DataFrame,
    amenability: pd.DataFrame,
    status: GateStatus | None = None,
    variant_ids: list[str] | None = None,
    safety: dict | None = None,
    nmd: pd.DataFrame | None = None,
    aenmd: pd.DataFrame | None = None,
    nmdetective: pd.DataFrame | None = None,
) -> WebTable:
    """Join the per-variant landscape to the per-therapy intervals into the app's table.

    `variant_ids` narrows the output to a chosen set — the example artifact is a diverse handful,
    the full export is every scoreable variant. `nmd`, when given, adds each variant's two NMD
    predictors and their rules; `aenmd` adds the published tool's verdict beside them, so the
    patient view can show a real second implementation agree or, where unavailable, why it did not.
    """

    if variant_ids is not None:
        landscape = landscape[landscape["variant_id"].isin(variant_ids)]
    kept = set(landscape["variant_id"])
    therapy = amenability[amenability["variant_id"].isin(kept)]
    by_variant: dict[object, list[dict]] = {
        vid: rows.to_dict("records") for vid, rows in therapy.groupby("variant_id")
    }
    nmd_by_variant = _nmd_by_variant(nmd)
    aenmd_by_variant = _aenmd_by_variant(aenmd)
    nmdetective_by_variant = _nmdetective_by_variant(nmdetective)
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
                    "name": therapy_name(str(arm["therapy_id"])),
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
                "stop": as_rna(str(row["stop_type"])).upper(),
                "residue": row["original_aa"],
                "review_stars": None if pd.isna(row["review_stars"]) else int(row["review_stars"]),
                "escapes_decay": bool(row["escapes_decay_by_rule"]),
                "tolerable_insertion_share": round(float(row["tolerable_insertion_share"]), 3),
                "best": best,
                "therapies": offered,
                "suppressor": _suppressor(str(row["stop_type"]), str(row["original_aa"])),
                "nmd": _with_models(
                    nmd_by_variant.get(row["variant_id"]),
                    aenmd_by_variant.get(row["variant_id"]),
                    nmdetective_by_variant.get(row["variant_id"]),
                ),
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
