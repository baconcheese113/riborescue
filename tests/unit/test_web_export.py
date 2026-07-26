import json

import pandas as pd

from riborescue.core.contracts import CONTRACTS_VERSION
from riborescue.variants.web_export import (
    GateStatus,
    build_web_table,
    diverse_sample,
    escape_summary,
)


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


def test_the_two_sri_compounds_are_never_left_to_be_confused():
    """`SRI` in the readthrough data is SRI-41315; the safety control's compound is SRI-37240.

    They are different molecules from one series, and the page shows both, so the short id alone
    would read as the control compound. The payload therefore carries the full name, and the gate
    text that mentions SRI-37240 says outright that it is not the therapy being scored.
    """

    from riborescue.variants.web_export import therapy_name

    assert therapy_name("SRI") == "SRI-41315"
    detail = build_web_table(_landscape(), _amenability()).status["readthrough_detail"]
    assert "SRI-37240" in detail and "SRI-41315" in detail
    assert "different molecules" in detail


def test_the_payload_stamps_the_contract_version_it_was_written_under():
    """A consumer reading a stored payload has to know which schema produced it."""

    parsed = json.loads(build_web_table(_landscape(), _amenability()).to_json())
    assert parsed["contracts_version"] == CONTRACTS_VERSION


def test_every_therapy_carries_the_compound_it_names():
    parsed = json.loads(build_web_table(_landscape(), _amenability()).to_json())
    assert parsed["therapy_names"]["G418"] == "G418 (geneticin)"
    v1 = next(v for v in parsed["variants"] if v["id"] == "V1")
    assert next(t for t in v1["therapies"] if t["id"] == "DAP")["name"] == "2,6-diaminopurine"


def test_nmd_verdicts_ride_along_when_given_and_are_null_otherwise():
    nmd = pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "rule_last_exon": False,
                "rule_within_last_junction": False,
                "rule_start_proximal": True,
                "rule_long_exon": False,
                "escape_guideline": False,
                "escape_full_rules": True,
                "predictors_disagree": True,
            }
        ]
    )
    with_nmd = build_web_table(_landscape(), _amenability(), nmd=nmd)
    by_id = {v["id"]: v for v in with_nmd.variants}
    assert by_id["V1"]["nmd"]["disagree"] is True
    assert by_id["V1"]["nmd"]["rules"]["start_proximal"] is True
    assert by_id["V2"]["nmd"] is None  # not in the nmd table

    without = build_web_table(_landscape(), _amenability())
    assert all(v["nmd"] is None for v in without.variants)


def test_the_aenmd_verdict_rides_along_inside_the_nmd_object():
    nmd = pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "rule_last_exon": False,
                "rule_within_last_junction": False,
                "rule_start_proximal": True,
                "rule_long_exon": False,
                "escape_guideline": False,
                "escape_full_rules": True,
                "predictors_disagree": True,
            }
        ]
    )
    aenmd = pd.DataFrame(
        [
            {"variant_id": "V1", "aenmd_available": True, "aenmd_escape": True, "reason": ""},
            {
                "variant_id": "V2",
                "aenmd_available": False,
                "aenmd_escape": False,
                "reason": "transcript_absent",
            },
        ]
    )
    table = build_web_table(_landscape(), _amenability(), nmd=nmd, aenmd=aenmd)
    by_id = {v["id"]: v for v in table.variants}
    assert by_id["V1"]["nmd"]["aenmd"] == {"available": True, "escape": True, "reason": ""}
    # V2 has no rule-tier nmd here, so it has no nmd object to hang the aenmd verdict on.
    assert by_id["V2"]["nmd"] is None


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


def test_safety_summary_reports_the_measured_layer_and_its_reach():
    from riborescue.variants.web_export import safety_summary

    predicted = pd.DataFrame(
        {
            "transcript": [f"T{i}" for i in range(60)],
            "mane_select": [True] * 50 + [False] * 10,
            "predicted_g418": [0.01 * (i % 10) for i in range(60)],
            "measured_lift": [0.01 * (i % 7) for i in range(40)] + [None] * 20,
            "group": (["both"] * 10 + ["predicted only"] * 15 + ["neither"] * 15) + [""] * 20,
        }
    )
    summary = safety_summary(predicted, ["G418", "Clitocine"], measured="G418")
    assert summary["measured_therapy"] == "G418"
    assert summary["canonical_stops_scored"] == 50
    assert summary["analysed"] == 40  # MANE with both predicted and measured
    assert summary["per_therapy"]["G418"].startswith("native-stop occupancy measured")
    assert summary["per_therapy"]["Clitocine"] == "no matched empirical safety atlas"
    assert "toxicity" in summary["caveat"]


def _base_editing() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_id": "V1",
                "scoreable": True,
                "reachable": True,
                "reach_class": "base_editable_exact",
                "n_guides": 2,
                "bystander_free": True,
                "editor": "ABE7.10",
                "strand": "antisense",
                "restores": "exact_wildtype",
                "window_position": 5,
            },
            {
                "variant_id": "V2",
                "scoreable": True,
                "reachable": False,
                "reach_class": "not_base_editable_under_panel",
                "n_guides": 0,
                "bystander_free": False,
                "editor": "",
                "strand": "",
                "restores": "",
                "window_position": None,
            },
            {
                "variant_id": "V3",
                "scoreable": False,
                "reachable": False,
                "reach_class": "",
            },
        ]
    )


def test_the_escape_summary_accounts_for_every_row():
    summary = escape_summary(_base_editing())
    assert summary["total"] == 3
    assert summary["scoreable"] == 2
    assert summary["unscoreable"] == 1
    # exact + alternative + not_editable equals scoreable — nothing dropped silently.
    assert (
        summary["exact"] + summary["alternative"] + summary["not_editable"] == summary["scoreable"]
    )
    assert summary["exact"] == 1 and summary["reachable"] == 1
    assert summary["exact_bystander_free"] == 1
    assert "off-target" in summary["caveat"]


def test_base_editing_rides_along_per_variant_and_drives_the_escape_aggregate():
    table = build_web_table(_landscape(), _amenability(), base_editing=_base_editing())
    by_id = {v["id"]: v for v in table.variants}
    assert by_id["V1"]["editing"]["reach_class"] == "base_editable_exact"
    assert by_id["V1"]["editing"]["editor"] == "ABE7.10"
    assert by_id["V1"]["editing"]["window_position"] == 5
    assert by_id["V2"]["editing"]["reach_class"] == "not_base_editable_under_panel"
    parsed = json.loads(table.to_json())
    assert parsed["escape"]["scoreable"] == 2
    # Without the table, the field is absent and every variant's editing is null.
    plain = build_web_table(_landscape(), _amenability())
    assert plain.escape is None
    assert all(v["editing"] is None for v in plain.variants)


def test_a_variant_with_no_available_therapy_has_null_best_not_nan():
    """NaN is not JSON; a browser's parser rejects the whole file on it."""

    landscape = pd.DataFrame(
        [
            {
                "variant_id": "V0",
                "gene_symbol": "STAMBP",
                "protein_position": 424,
                "stop_type": "uga",
                "original_aa": "R",
                "review_stars": 2,
                "escapes_decay_by_rule": True,
                "nt_to_last_junction": -10,
                "best_therapy": float("nan"),
                "best_readthrough": float("nan"),
                "best_readthrough_low": float("nan"),
                "tolerable_insertion_share": 0.167,
            }
        ]
    )
    amenability = pd.DataFrame(
        [
            {
                "variant_id": "V0",
                "therapy_id": "G418",
                "readthrough_predicted": None,
                "readthrough_low": None,
                "readthrough_high": None,
                "status": "missing",
                "reason": "not_available",
            }
        ]
    )
    table = build_web_table(landscape, amenability)
    assert table.variants[0]["best"] is None
    # to_json must not raise and must not contain the NaN token
    text = table.to_json()
    assert "NaN" not in text
    assert json.loads(text)["variants"][0]["best"] is None
