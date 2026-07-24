import json

import pandas as pd

from riborescue.variants.webexport import GateStatus, build_web_table, diverse_sample


def _landscape() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "gene_symbol": "TP53",
                "protein_position": 196,
                "stop_type": "uga",
                "original_aa": "R",
                "review_stars": 2,
                "escapes_decay_by_rule": True,
                "nt_to_last_junction": -80,
                "best_therapy": "DAP",
                "best_readthrough": 0.02,
                "best_readthrough_low": 0.01,
                "tolerable_insertion_share": 0.6,
            },
            {
                "variant_id": "V2",
                "gene_symbol": "CFTR",
                "protein_position": 1282,
                "stop_type": "uag",
                "original_aa": "W",
                "review_stars": 3,
                "escapes_decay_by_rule": False,
                "nt_to_last_junction": 40,
                "best_therapy": "G418",
                "best_readthrough": 0.005,
                "best_readthrough_low": 0.001,
                "tolerable_insertion_share": 0.1,
            },
        ]
    )


def _amenability() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "therapy_id": "DAP",
                "readthrough_predicted": 0.02,
                "readthrough_low": 0.01,
                "readthrough_high": 0.03,
                "status": "present",
                "reason": None,
            },
            {
                "variant_id": "V1",
                "therapy_id": "G418",
                "readthrough_predicted": 0.015,
                "readthrough_low": 0.008,
                "readthrough_high": 0.022,
                "status": "present",
                "reason": None,
            },
            {
                "variant_id": "V2",
                "therapy_id": "G418",
                "readthrough_predicted": 0.005,
                "readthrough_low": 0.001,
                "readthrough_high": 0.009,
                "status": "present",
                "reason": None,
            },
            {
                "variant_id": "V2",
                "therapy_id": "CC90009",
                "readthrough_predicted": None,
                "readthrough_low": None,
                "readthrough_high": None,
                "status": "missing",
                "reason": "not_available",
            },
        ]
    )


def test_the_payload_carries_the_gate_statuses():
    table = build_web_table(_landscape(), _amenability())
    assert "G418" in table.status["readthrough_control"]
    assert "HEK293T" in table.status["safety_atlas"]


def test_a_variant_carries_every_therapy_and_its_interval():
    table = build_web_table(_landscape(), _amenability())
    v1 = next(v for v in table.variants if v["id"] == "V1")
    assert {t["id"] for t in v1["therapies"]} == {"DAP", "G418"}
    dap = next(t for t in v1["therapies"] if t["id"] == "DAP")
    assert dap["available"] and dap["low"] == 0.01 and dap["high"] == 0.03


def test_an_unavailable_therapy_keeps_its_reason_and_no_number():
    table = build_web_table(_landscape(), _amenability())
    v2 = next(v for v in table.variants if v["id"] == "V2")
    cc = next(t for t in v2["therapies"] if t["id"] == "CC90009")
    assert cc["available"] is False
    assert cc["readthrough"] is None
    assert cc["reason"] == "not_available"


def test_the_suppressor_reinserts_the_original_residue():
    table = build_web_table(_landscape(), _amenability())
    v1 = next(v for v in table.variants if v["id"] == "V1")
    assert v1["suppressor"] == {"design": "UGA-R", "restores_exactly": True}
    assert v1["stop"] == "UGA"


def test_the_status_can_be_overridden_for_a_dataset_that_has_been_confirmed():
    confirmed = GateStatus(readthrough_control="confirmed", readthrough_detail="…")
    table = build_web_table(_landscape(), _amenability(), status=confirmed)
    assert table.status["readthrough_control"] == "confirmed"


def test_the_payload_is_json_and_round_trips():
    table = build_web_table(_landscape(), _amenability())
    parsed = json.loads(table.to_json())
    assert parsed["therapies"] == ["CC90009", "DAP", "G418"]
    assert len(parsed["variants"]) == 2


def test_the_sample_prefers_the_states_a_uniform_draw_would_miss():
    landscape = _landscape()
    amenability = _amenability()
    chosen = diverse_sample(landscape, amenability, size=2)
    # V2 carries the missing therapy and does not escape decay; it must not be dropped.
    assert "V2" in chosen
