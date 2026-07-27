import json

import pandas as pd
import pytest

from riborescue.evidence_export import (
    build_evidence,
    codon_signature,
    frame_by_length,
    kinetics_null,
    model_parity,
    periodicity,
    readthrough_contrast,
)

ARMS = pd.DataFrame(
    {
        "sample": ["lib_a", "lib_b"],
        "treatment": ["dmso", "g418"],
        "replicate": [1, 1],
    }
)


def _effects() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "quantity": ["downstream_occupancy", "termination_occupancy", "frame_gap"],
            "mean_difference": [0.0096, -0.0116, 0.1806],
            "ci_low": [0.0050, -0.0262, 0.1196],
            "ci_high": [0.0142, 0.0030, 0.2415],
            "consistent": [True, True, True],
        }
    )


def _libraries(**overrides) -> pd.DataFrame:
    base = pd.DataFrame(
        {
            "sample": ["lib_a", "lib_b"],
            "transcripts": [486, 486],
            "downstream_occupancy": [0.0028, 0.0121],
            "termination_occupancy": [0.0707, 0.0591],
            "frame_gap": [-0.2557, -0.0751],
        }
    )
    return base.assign(**overrides)


def test_the_libraries_behind_a_mean_difference_travel_with_it():
    """A three-against-three difference says nothing about separation; the points show it."""

    contrast = readthrough_contrast(_effects(), _libraries(), ARMS)
    assert [q["quantity"] for q in contrast["quantities"]] == [
        "downstream_occupancy",
        "termination_occupancy",
        "frame_gap",
    ]
    assert len(contrast["libraries"]) == 2
    assert contrast["libraries"][0]["treatment"] == "dmso"
    assert contrast["libraries"][1]["downstream_occupancy"] == pytest.approx(0.0121)


def test_an_arm_already_recorded_beside_the_ratios_is_not_joined_over():
    """The assay writes the arm itself; joining again would leave two columns and read neither."""

    already = _libraries(treatment=["untreated", "treated"])
    contrast = readthrough_contrast(_effects(), already, ARMS)
    assert [row["treatment"] for row in contrast["libraries"]] == ["untreated", "treated"]


def test_the_kept_lengths_are_marked_against_the_ones_surveyed():
    """The plot has to show which lengths were chosen, not only how each one behaved."""

    frames = pd.DataFrame(
        {
            "length": [21, 21, 21, 28, 28, 28],
            "frame": [0, 1, 2, 0, 1, 2],
            "sample": ["lib_a"] * 6,
            "n": [800, 100, 100, 100, 400, 500],
        }
    )
    rows = {row["length"]: row for row in frame_by_length(frames, kept=(21,))}
    assert rows[21]["kept"] is True
    assert rows[28]["kept"] is False
    assert rows[21]["frame0_share"] == pytest.approx(0.8)


def test_periodicity_is_averaged_within_an_arm_not_across_the_dataset():
    """Two arms that differ at the stop codon must not be averaged into one curve."""

    meta = pd.DataFrame(
        {
            "sample": ["lib_a", "lib_a", "lib_b", "lib_b"],
            "region": ["Distance from start (nt)"] * 2 + ["Distance from stop (nt)"] * 2,
            "distance": [0, 3, 0, 3],
            "count": [10, 20, 30, 40],
            "scaled_count": [0.1, 0.2, 0.3, 0.4],
        }
    )
    rows = periodicity(meta, ARMS)
    assert {row["region"] for row in rows} == {"start", "stop"}
    assert {row["treatment"] for row in rows} == {"dmso", "g418"}


def test_every_sense_codon_is_carried_for_each_site():
    """The A and P conventions answer different questions and are never merged."""

    def table(site: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "codon": ["GAA", "CCA"],
                "amino_acid": ["E", "P"],
                "site": [site, site],
                "occupancy": [1.4, 0.9],
                "occupancy_sd": [0.01, 0.02],
                "libraries": [3, 3],
            }
        )

    rows = codon_signature([table("a"), table("p")])
    assert len(rows) == 4
    assert {row["site"] for row in rows} == {"a", "p"}


def test_a_null_row_carries_which_shuffle_produced_it():
    """A p-value without its shuffle is unreadable: the three break different things."""

    familywise = pd.DataFrame(
        {
            "drug": ["SJ6986", "SJ6986"],
            "gain": [0.00356, 0.00356],
            "permutations": [199, 199],
            "null_mean": [0.00299, 0.00073],
            "null_sd": [0.00067, 0.00086],
            "null_max": [0.00515, 0.00623],
            "p_familywise": [0.21, 0.02],
            "resolution": [0.005, 0.005],
            "shuffle": ["context_matched", "global"],
        }
    )
    section = kinetics_null(familywise)
    rows = section["rows"]
    assert {row["shuffle"] for row in rows} == {"context_matched", "global"}
    assert rows[0]["p_familywise"] == pytest.approx(0.21)


def test_a_null_short_of_its_declared_count_is_labelled_incomplete():
    """199 permutations is a partial run, and a page must not show it as the declared result."""

    familywise = pd.DataFrame(
        {
            "drug": ["SJ6986"],
            "gain": [0.00356],
            "permutations": [199],
            "null_mean": [0.00299],
            "null_sd": [0.00067],
            "null_max": [0.00515],
            "p_familywise": [0.21],
            "resolution": [0.005],
            "shuffle": ["context_matched"],
        }
    )
    section = kinetics_null(familywise)
    assert section["permutations_completed"] == 199
    assert section["permutations_required"] == 999
    assert section["analysis_status"] == "incomplete"
    assert section["resolution"] == pytest.approx(0.005)


def test_an_analysis_that_ran_and_found_nothing_is_empty_not_absent():
    """A page must tell 'not run' (null) from 'ran, no rows' (empty), so the two differ here."""

    ran = build_evidence({"dataset": "d"}, codon_occupancy=[])
    payload = json.loads(ran.to_json())
    assert payload["codon_occupancy"] == []
    assert payload["readthrough"] is None


def test_a_score_is_reported_against_its_own_drug_ceiling():
    """0.80 against 0.88 is a different claim from 0.80 against 1.0, and the ceiling is per drug."""

    metrics = pd.DataFrame({"drug": ["G418"] * 2, "round": [1, 2], "r2": [0.80, 0.82]})
    reliability = pd.DataFrame({"treatment": ["G418"], "ceiling": [0.88]})
    rows = model_parity(metrics, reliability)
    assert rows[0]["r2_mean"] == pytest.approx(0.81)
    assert rows[0]["ceiling"] == pytest.approx(0.88)
    assert rows[0]["rounds"] == 2


def test_an_unmeasured_section_is_absent_rather_than_empty():
    """A page must be able to tell 'not run' from 'run and found nothing'."""

    payload = json.loads(build_evidence({"dataset": "d", "commit": "abc"}).to_json())
    assert payload["readthrough"] is None
    assert payload["codon_occupancy"] is None
    assert payload["provenance"]["dataset"] == "d"


def test_the_payload_refuses_a_value_no_browser_could_parse():
    """allow_nan=False, for the same reason the other two exports set it."""

    parity = [{"drug": "X", "ceiling": float("nan")}]
    payload = build_evidence({"dataset": "d"}, model_parity=parity)
    with pytest.raises(ValueError):
        payload.to_json()
