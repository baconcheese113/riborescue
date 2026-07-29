"""The whole scored set as one compact payload the browser can hold.

Columns rather than records, dictionary-encoded strings, and fixed-point integers: the same 70,376
variants that would be 165 MB as objects fit in a few megabytes, which is what lets the lookup view
cover every variant instead of a demonstration slice.

Readthrough values are stored as parts per hundred thousand and residue tolerance per thousand, so
every number is an integer and JSON does not spend bytes on float formatting. A therapy arm that was
not scored is -1 rather than absent, so the six arms of a variant are always at a known offset.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

RESULTS = Path("results")
OUT = Path("frontend/public/riborescue_index.json")
STOP_CODE = {"UAG": 0, "UGA": 1, "UAA": 2, "TAG": 0, "TGA": 1, "TAA": 2}
SCALE_READTHROUGH = 100_000
SCALE_TOLERANCE = 1_000
ABSENT = -1

csv.field_size_limit(sys.maxsize)


def rows(name: str):
    with (RESULTS / name).open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def scaled(value: str | None, scale: int) -> int:
    return round(float(value) * scale) if value else 0


def stop_code(value: str) -> int:
    """The table writes stop codons in lower case. Refusing an unknown one keeps a silent default
    from labelling every variant with the same codon."""
    code = STOP_CODE.get(value.strip().upper())
    if code is None:
        raise ValueError(f"unrecognised stop codon {value!r}")
    return code


def main() -> None:
    landscape = list(rows("amenability_landscape.tsv"))
    order = {row["variant_id"]: i for i, row in enumerate(landscape)}

    # One condition per variant. A placeholder is recorded as absent rather than as a disease name,
    # because "not provided" is what the archive says when nobody entered one.
    condition: dict[str, str] = {}
    extra: dict[str, int] = defaultdict(int)
    for row in rows("diseases.tsv"):
        variant = row["variant_id"]
        if variant not in order:
            continue
        if row["mapping_status"] != "mapped":
            continue
        if variant in condition:
            extra[variant] += 1
            continue
        condition[variant] = row["condition_name"].replace("_", " ")

    # The two NMD rule sets, so the panel can still show where they disagree.
    full_rules: dict[str, int] = {}
    for row in rows("nmd.tsv"):
        if row["variant_id"] in order:
            full_rules[row["variant_id"]] = 1 if row["escape_full_rules"] in {"True", "true", "1"} else 0

    editing: dict[str, str] = {}
    for row in rows("base_editing.tsv"):
        if row["variant_id"] in order and row.get("reach_class"):
            editing[row["variant_id"]] = row["reach_class"]

    genes: dict[str, int] = {}
    conditions: dict[str, int] = {}
    residues: dict[str, int] = {}
    reaches: dict[str, int] = {}

    def code(table: dict[str, int], key: str) -> int:
        if key not in table:
            table[key] = len(table)
        return table[key]

    column: dict[str, list[int]] = defaultdict(list)
    for row in landscape:
        variant = row["variant_id"]
        column["gene"].append(code(genes, row["gene_symbol"]))
        column["position"].append(int(float(row["protein_position"] or 0)))
        column["stop"].append(stop_code(row["stop_type"]))
        column["residue"].append(code(residues, row["original_aa"] or "?"))
        column["escapes"].append(1 if row["escapes_decay_by_rule"] in {"True", "true", "1"} else 0)
        column["escapes_full"].append(full_rules.get(variant, ABSENT))
        column["stars"].append(int(float(row["review_stars"] or 0)))
        column["best"].append(scaled(row["best_readthrough"], SCALE_READTHROUGH))
        column["best_low"].append(scaled(row["best_readthrough_low"], SCALE_READTHROUGH))
        column["tolerance"].append(scaled(row["tolerable_insertion_share"], SCALE_TOLERANCE))
        column["extra_conditions"].append(extra.get(variant, 0))
        name = condition.get(variant)
        column["condition"].append(code(conditions, name) if name else ABSENT)
        reach = editing.get(variant)
        column["reach"].append(code(reaches, reach) if reach else ABSENT)

    therapies: dict[str, int] = {}
    arms: dict[str, dict[int, tuple[int, int, int]]] = defaultdict(dict)
    for row in rows("variant_therapy_scores.tsv"):
        slot = code(therapies, row["therapy_id"])
        if row["variant_id"] not in order or not row["readthrough_predicted"]:
            continue
        arms[row["variant_id"]][slot] = (
            scaled(row["readthrough_predicted"], SCALE_READTHROUGH),
            scaled(row["readthrough_low"], SCALE_READTHROUGH),
            scaled(row["readthrough_high"], SCALE_READTHROUGH),
        )

    width = len(therapies)
    flat: dict[str, list[int]] = {"arm_mid": [], "arm_low": [], "arm_high": []}
    for row in landscape:
        scored = arms.get(row["variant_id"], {})
        for slot in range(width):
            value = scored.get(slot)
            flat["arm_mid"].append(value[0] if value else ABSENT)
            flat["arm_low"].append(value[1] if value else ABSENT)
            flat["arm_high"].append(value[2] if value else ABSENT)

    payload = {
        "scales": {"readthrough": SCALE_READTHROUGH, "tolerance": SCALE_TOLERANCE},
        "therapies": list(therapies),
        "genes": list(genes),
        "conditions": list(conditions),
        "residues": list(residues),
        "reach_classes": list(reaches),
        "ids": [row["variant_id"] for row in landscape],
        **column,
        **flat,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    mapped = sum(1 for value in column["condition"] if value != ABSENT)
    print(
        f"{len(landscape):,} variants · {len(genes):,} genes · {len(conditions):,} conditions "
        f"({mapped / len(landscape):.1%} mapped) · {width} therapy arms · "
        f"{OUT.stat().st_size / 1e6:.1f} MB"
    )


if __name__ == "__main__":
    main()
