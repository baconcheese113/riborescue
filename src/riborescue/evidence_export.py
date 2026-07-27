"""The evidence payload: what the project measured, and how it knows the measurement is one.

The patient and researcher payloads answer *what is amenable*. This one answers *why any of it
should be believed* — the positive control and the negative one beside it, the calibration the
contrast rests on, the codon signature the pipeline recovers without being tuned for it, the null
the kinetic claim failed against, and the disagreement between the safety model and the safety
measurement.

Every section is optional and absent rather than empty when its input has not been produced, so a
laptop that ran the dry-lab chain and none of the sequencing still builds a page. Every number is
carried with the spread or interval it was measured with; nothing here is a point estimate standing
on its own.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import pandas as pd

from riborescue.core.contracts import CONTRACTS_VERSION

__all__ = [
    "QUANTITIES",
    "EvidencePayload",
    "build_evidence",
    "codon_signature",
    "kinetics_null",
    "model_parity",
    "periodicity",
    "readthrough_contrast",
    "safety_concordance",
]

QUANTITIES = ("downstream_occupancy", "termination_occupancy", "frame_gap")
"""The three the frame control turns on, in the order they are read."""

PERMUTATIONS_REQUIRED = 999
"""Permutations the shuffle null is declared over. Fewer is a partial run, not a smaller test."""

_METAPROFILE_REGIONS = {
    "Distance from start (nt)": "start",
    "Distance from stop (nt)": "stop",
}


def _int(value: object) -> int:
    """Pandas hands back a scalar of its own union type; every count here is a plain integer."""

    return int(cast(int, value))


def _round(value: object, places: int = 6) -> float | None:
    """A JSON number, or None where the quantity is genuinely absent."""

    if value is None or pd.isna(value):  # type: ignore[arg-type]
        return None
    return round(float(value), places)  # type: ignore[arg-type]


def readthrough_contrast(
    effects: pd.DataFrame, libraries: pd.DataFrame, arms: pd.DataFrame
) -> dict:
    """One contrast: the three quantities with their intervals, and the libraries behind them.

    The per-library values travel with the summary because a mean difference of three against three
    says nothing about whether the arms separated, and the reader can see that only from the points.

    The assay already writes each library's arm beside its ratios; `arms` supplies it only for a
    table that predates that, and joining it over the top would leave two columns of the same name
    and neither of them read.
    """

    labelled = (
        libraries
        if "treatment" in libraries.columns
        else libraries.merge(arms, on="sample", how="left")
    )
    return {
        "quantities": [
            {
                "quantity": row.quantity,
                "mean_difference": _round(row.mean_difference),
                "ci_low": _round(row.ci_low),
                "ci_high": _round(row.ci_high),
                "consistent": bool(row.consistent),
            }
            for row in effects.itertuples()
        ],
        "libraries": [
            {
                "sample": row["sample"],
                "treatment": row.get("treatment"),
                "transcripts": int(row["transcripts"]),
                **{q: _round(row[q]) for q in QUANTITIES if q in row},
            }
            for _, row in labelled.iterrows()
        ],
    }


def periodicity(metaprofile: pd.DataFrame, arms: pd.DataFrame) -> list[dict]:
    """P-site counts either side of the start and stop codons, averaged within a treatment arm.

    The three-nucleotide sawtooth is what says these libraries are ribosome profiles rather than
    fragmented RNA, and it is the one claim no downstream number can substitute for.
    """

    labelled = metaprofile.merge(arms, on="sample", how="left")
    labelled = labelled.assign(region=labelled["region"].map(_METAPROFILE_REGIONS))
    grouped = labelled.dropna(subset=["region"]).groupby(
        ["region", "treatment", "distance"], as_index=False
    )
    averaged = grouped.agg(scaled=("scaled_count", "mean"), libraries=("sample", "nunique"))
    return [
        {
            "region": row.region,
            "treatment": row.treatment,
            "distance": _int(row.distance),
            "scaled": _round(row.scaled, 8),
            "libraries": _int(row.libraries),
        }
        for row in averaged.itertuples()
    ]


def frame_by_length(frames: pd.DataFrame, kept: tuple[int, ...]) -> list[dict]:
    """Per read length, how its P-sites distribute over the three coding frames.

    This is the evidence the length selection was read from: a length is kept when frame 0 leads in
    every library, and which lengths those are is a property of the dataset rather than a default.
    """

    wide = frames.pivot_table(
        index="length", columns="frame", values="n", aggfunc="sum", fill_value=0
    )
    for frame in (0, 1, 2):
        if frame not in wide.columns:
            wide[frame] = 0
    total = wide[[0, 1, 2]].sum(axis=1)
    share = total / total.sum() if total.sum() else total
    return [
        {
            "length": _int(length),
            "frame0": _int(wide.loc[length, 0]),
            "frame1": _int(wide.loc[length, 1]),
            "frame2": _int(wide.loc[length, 2]),
            "frame0_share": _round(wide.loc[length, 0] / total[length] if total[length] else None),
            "library_share": _round(share[length]),
            "kept": _int(length) in kept,
        }
        for length in wide.index
    ]


def codon_signature(tables: list[pd.DataFrame]) -> list[dict]:
    """Occupancy per codon at each scored site, with the spread across libraries beside it.

    Elevated A-site occupancy on glutamate and aspartate, and P-site occupancy on proline, are
    documented human signatures. Nothing here was tuned to reproduce them, so recovering them is a
    statement about the pipeline rather than about the compounds.
    """

    stacked = pd.concat(tables, ignore_index=True)
    return [
        {
            "codon": row.codon,
            "amino_acid": row.amino_acid,
            "site": row.site,
            "occupancy": _round(row.occupancy),
            "occupancy_sd": _round(getattr(row, "occupancy_sd", None)),
            "libraries": int(getattr(row, "libraries", 0) or 0),
        }
        for row in stacked.itertuples()
    ]


def kinetics_null(familywise: pd.DataFrame) -> dict:
    """Each drug's observed gain against the null of the best gain any drug reached under a shuffle.

    Three shuffles break three different things, and only the synonymous one leaves the amino acid
    intact while destroying decoding demand. A gain that clears the first two and not the third is
    attributable to residue identity, which is what these rows are for.

    How many permutations the null was built from travels with them, and a run short of the declared
    count is labelled incomplete. A p-value from a partial run is an estimate on a coarser grid, not
    a bound in either direction, and a page that showed it unlabelled would be reporting a number
    the record does not authorise.
    """

    completed = _int(familywise["permutations"].max()) if len(familywise) else 0
    rows = [
        {
            "drug": row.drug,
            "shuffle": row.shuffle,
            "gain": _round(row.gain),
            "null_mean": _round(row.null_mean),
            "null_sd": _round(row.null_sd),
            "null_max": _round(row.null_max),
            "p_familywise": _round(row.p_familywise),
            "permutations": _int(row.permutations),
        }
        for row in familywise.itertuples()
    ]
    return {
        "permutations_completed": completed,
        "permutations_required": PERMUTATIONS_REQUIRED,
        "analysis_status": "complete" if completed >= PERMUTATIONS_REQUIRED else "incomplete",
        "resolution": _round(1 / (completed + 1)) if completed else None,
        "rows": rows,
    }


def safety_concordance(predicted: pd.DataFrame, points: int = 2000) -> dict:
    """Where the amenability model and the measured native-stop occupancy agree, and where they do
    not.

    The model was fitted around premature stops and a canonical stop is out of that distribution, so
    the comparison is rank-only and the disagreement is the finding. Every stop carrying both a
    prediction and a measurement is returned, capped so the payload stays a page rather than a file.
    """

    from riborescue.variants.native_stop_predictions import concordance

    mane = predicted[predicted["mane_select"].astype(bool)]
    matched = mane.dropna(subset=["predicted_g418", "measured_lift"])
    stats = concordance(matched["predicted_g418"], matched["measured_lift"], draws=1000)
    counts = matched["group"].value_counts()
    shown = matched if len(matched) <= points else matched.sample(points, random_state=0)
    return {
        "rho": stats["rho"],
        "low": stats["low"],
        "high": stats["high"],
        "analysed": len(matched),
        "canonical_stops_scored": len(mane),
        "quadrants": {
            group: int(counts.get(group, 0))
            for group in ("both", "predicted only", "measured only", "neither")
        },
        "points": [
            {
                "gene": row.gene,
                "predicted": _round(row.predicted_g418),
                "measured": _round(row.measured_lift),
                "group": row.group,
            }
            for row in shown.itertuples()
        ],
    }


def model_parity(metrics: pd.DataFrame, reliability: pd.DataFrame) -> list[dict]:
    """Each drug's held-out performance against the ceiling its own replicates set.

    A model at 0.80 against a ceiling of 0.88 is a different claim from the same number against 1.0,
    and the ceiling is per drug rather than one figure for all six.
    """

    rounds = metrics.groupby("drug")["r2"].agg(["mean", "std", "size"])
    ceilings = reliability.set_index("treatment")["ceiling"]
    return [
        {
            "drug": drug,
            "r2_mean": _round(row["mean"]),
            "r2_sd": _round(row["std"]),
            "rounds": int(row["size"]),
            "ceiling": _round(ceilings.get(drug)),
        }
        for drug, row in rounds.iterrows()
    ]


@dataclass(frozen=True)
class EvidencePayload:
    """Every section the page can draw, each absent rather than empty when unmeasured."""

    provenance: dict
    readthrough: dict | None = None
    calibration: dict | None = None
    periodicity: list[dict] | None = None
    codon_occupancy: list[dict] | None = None
    kinetics_null: dict | None = None
    safety: dict | None = None
    model_parity: list[dict] | None = None

    def to_json(self) -> str:
        # allow_nan=False so a stray NaN fails here rather than reaching a browser as a token no
        # JSON parser accepts.
        return json.dumps(
            {
                "contracts_version": CONTRACTS_VERSION,
                "provenance": self.provenance,
                "readthrough": self.readthrough,
                "calibration": self.calibration,
                "periodicity": self.periodicity,
                "codon_occupancy": self.codon_occupancy,
                "kinetics_null": self.kinetics_null,
                "safety": self.safety,
                "model_parity": self.model_parity,
            },
            indent=2,
            allow_nan=False,
        )


def build_evidence(provenance: dict, **sections: object) -> EvidencePayload:
    """Assemble whatever sections were produced, leaving the rest absent."""

    return EvidencePayload(provenance=provenance, **sections)  # type: ignore[arg-type]


def read_arms(samplesheet: Path, dataset: str) -> pd.DataFrame:
    """Each footprint library of a dataset and the treatment arm it belongs to."""

    runs = pd.read_csv(samplesheet, sep="\t")
    footprints = runs.loc[(runs["assay"] == "riboseq") & (runs["dataset"] == dataset)]
    return footprints[["sample", "treatment", "replicate"]].reset_index(drop=True)
