import json

import pandas as pd

from riborescue.core.contracts import CONTRACTS_VERSION
from riborescue.variants.research_export import build_research_aggregate


def _diseases() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_id": vid,
                "gene_symbol": gene,
                "condition_name": name,
                "medgen": medgen,
                "omim": "1",
                "orphanet": "",
                "mondo": "",
                "mesh": "",
                "mapping_status": "mapped",
                "reason": "",
            }
            for vid, gene, name, medgen in [
                ("V1", "A", "alpha", "C1"),
                ("V2", "A", "alpha", "C1"),
                ("V3", "B", "beta", "C2"),
            ]
        ]
    )


def _contexts() -> pd.DataFrame:
    geometry = {
        "in_last_exon": False,
        "nt_to_last_junction": 800.0,
        "nt_from_start": 500,
        "ptc_exon_length": 100,
    }
    return pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "gene_symbol": "A",
                "scoreable": True,
                "stop_type": "uga",
                "original_aa": "R",
                **geometry,
            },
            {
                "variant_id": "V2",
                "gene_symbol": "A",
                "scoreable": True,
                "stop_type": "uga",
                "original_aa": "R",
                **geometry,
            },
            {
                "variant_id": "V3",
                "gene_symbol": "B",
                "scoreable": False,
                "stop_type": "",
                "original_aa": "",
                **geometry,
            },
        ]
    )


def test_the_aggregate_carries_provenance_and_denominators():
    aggregate = build_research_aggregate(
        _diseases(), _contexts(), clinvar_release="20260715", commit="abc1234"
    )
    prov = aggregate.provenance
    assert prov["clinvar_release"] == "20260715"
    assert prov["commit"] == "abc1234"
    assert prov["qualifying_variants"] == 3
    assert prov["scoreable_variants"] == 2
    assert prov["condition_entities"] == 2


def test_all_three_frontiers_are_present():
    aggregate = build_research_aggregate(_diseases(), _contexts(), clinvar_release="r")
    assert set(aggregate.frontiers) == {"variants", "genes", "conditions"}
    assert aggregate.frontiers["variants"][0]["design_id"] == "UGA-R"


def test_the_nmd_disagreement_atlas_rides_along():
    aggregate = build_research_aggregate(_diseases(), _contexts(), clinvar_release="r")
    assert aggregate.nmd["scoreable"] == 2
    assert {"escape_guideline", "escape_full_rules", "disagree"} <= set(aggregate.nmd)


def test_the_reach_denominator_is_spelled_out():
    aggregate = build_research_aggregate(_diseases(), _contexts(), clinvar_release="r")
    d = aggregate.reach_denominator
    assert d["eligible_condition_entities"] == 2
    assert d["unique_variant_condition_pairs"] == 3
    assert d["reachable_entities"] == 1  # only C1 has a scoreable variant


def test_the_top_list_is_bounded_and_ordered_by_denominator():
    aggregate = build_research_aggregate(_diseases(), _contexts(), clinvar_release="r", top=1)
    assert len(aggregate.condition_coverage_top) == 1
    assert aggregate.condition_coverage_top[0]["eligible_variants"] == 2
    assert "reach" in aggregate.condition_coverage_top[0]
    assert "complete" in aggregate.condition_coverage_top[0]


def test_the_json_is_strictly_valid_and_round_trips():
    aggregate = build_research_aggregate(_diseases(), _contexts(), clinvar_release="r")
    text = aggregate.to_json()
    assert "NaN" not in text
    parsed = json.loads(text)
    assert parsed["provenance"]["condition_entities"] == 2
    assert parsed["provenance"]["contracts_version"] == CONTRACTS_VERSION
    assert "unmet_need" in parsed["caveats"]
