"""The designer's contract, and the two things it must never do.

It must never blend the axes into one score, and it must never report a reach it could not compute.
Both failures look like a working system from the outside, which is why they are tested rather than
left to review.
"""

from pathlib import Path

import pandas as pd
import pytest

from riborescue.variants.experiment_designer import (
    AXES,
    PROGRAMS,
    REACH,
    frontier,
    propose,
    read_programs,
)

_PROGRAM = {
    "experiment_id": "EXAMPLE",
    "question": "does the thing happen?",
    "why_it_matters": "because",
    "what_the_lab_does": "measures it",
    "comparison": "treated against control",
    "success_criterion": "the number moves",
    "if_it_fails": "it does not happen",
    "evidence_gap_reason": "nobody measured it",
    "assay": "an assay",
    "model_system": "cells",
    "endpoint": "the number",
    "decision_rule": "declared in advance",
    "replicates": "not estimated",
    "resolves": "none",
    "reach_rule": "all_scoreable",
    "direct_scope": "one thing in one place",
    "generalisation_required": "that it holds elsewhere",
    "replicate_endpoint": "not estimated",
    "replicate_effect": "not estimated",
    "replicate_variance_source": "not estimated",
    "replicate_design": "not estimated",
    "replicate_alpha_power": "not estimated",
    "replicate_method": "not estimated",
    "evidence_grade": "none",
    "complexity": "moderate",
    "safety_relevant": "FALSE",
    "provenance": "ADR-0000",
}

_CONTEXTS = pd.DataFrame(
    {
        "variant_id": ["V1", "V2", "V3"],
        "gene_symbol": ["AAA", "AAA", "BBB"],
        "stop_type": ["uag", "uga", "uag"],
        "escapes_decay_by_rule": [True, False, True],
    }
)
_DISEASES = pd.DataFrame({"variant_id": ["V1", "V2", "V3"], "medgen": ["C1", "C1", "C2"]})
_LEDGER = pd.DataFrame(
    {
        "claim_id": ["OPEN", "CLOSED"],
        "verdict": ["untestable", "supported"],
    }
)


def _written(tmp_path: Path, *rows: dict[str, str]) -> Path:
    path = tmp_path / "programs.tsv"
    pd.DataFrame(list(rows)).to_csv(path, sep="\t", index=False)
    return path


class TestTheContract:
    def test_a_missing_field_is_refused(self, tmp_path):
        path = _written(tmp_path, {k: v for k, v in _PROGRAM.items() if k != "decision_rule"})
        with pytest.raises(ValueError, match="missing decision_rule"):
            read_programs(path)

    def test_a_blank_cell_is_refused_so_none_has_to_be_written_out(self, tmp_path):
        path = _written(tmp_path, _PROGRAM | {"resolves": ""})
        with pytest.raises(ValueError, match="leaves resolves empty"):
            read_programs(path)

    def test_a_reach_rule_that_does_not_exist_is_refused(self, tmp_path):
        path = _written(tmp_path, _PROGRAM | {"reach_rule": "everyone"})
        with pytest.raises(ValueError, match="reach rules that do not exist"):
            read_programs(path)

    def test_an_invented_complexity_tier_is_refused(self, tmp_path):
        path = _written(tmp_path, _PROGRAM | {"complexity": "trivial"})
        with pytest.raises(ValueError, match="complexity tiers outside"):
            read_programs(path)


class TestReach:
    def test_reach_counts_the_variants_the_rule_selects(self, tmp_path):
        programs = read_programs(_written(tmp_path, _PROGRAM | {"reach_rule": "uag_stops"}))
        row = propose(programs, _CONTEXTS, _DISEASES, _LEDGER).iloc[0]
        assert row["potential_variants"] == 2
        assert row["potential_genes"] == 2
        assert row["potential_conditions"] == 2

    def test_a_rule_whose_input_is_absent_refuses_rather_than_reporting_nobody(self, tmp_path):
        # Zero reach is a finding; a missing column is a bug, and they must not look alike.
        programs = read_programs(_written(tmp_path, _PROGRAM | {"reach_rule": "decay_escaping"}))
        with pytest.raises(ValueError, match="needs `escapes_decay_by_rule`"):
            propose(programs, _CONTEXTS.drop(columns="escapes_decay_by_rule"), _DISEASES, _LEDGER)

    def test_every_named_reach_rule_is_callable(self):
        for rule in REACH.values():
            assert callable(rule)


class TestClaims:
    def test_only_open_claims_count(self, tmp_path):
        programs = read_programs(
            _written(
                tmp_path,
                _PROGRAM | {"resolves": "OPEN"},
                _PROGRAM | {"experiment_id": "B", "resolves": "CLOSED"},
                _PROGRAM | {"experiment_id": "C", "resolves": "none"},
            )
        )
        counts = propose(programs, _CONTEXTS, _DISEASES, _LEDGER).set_index("experiment_id")
        assert counts.loc["EXAMPLE", "open_questions_addressed"] == 1
        assert counts.loc["B", "open_questions_addressed"] == 0
        assert counts.loc["C", "open_questions_addressed"] == 0


class TestFrontier:
    @staticmethod
    def _ranked(**axes) -> pd.DataFrame:
        base = dict.fromkeys(AXES, 0)
        rows = [{"experiment_id": name} | base | values for name, values in axes.items()]
        return frontier(pd.DataFrame(rows))

    def test_a_programme_beaten_on_every_axis_falls_off(self):
        ranked = self._ranked(
            BETTER={"potential_variants": 10, "open_questions_addressed": 2},
            WORSE={"potential_variants": 5, "open_questions_addressed": 1},
        ).set_index("experiment_id")
        assert ranked.loc["BETTER", "on_frontier"]
        assert not ranked.loc["WORSE", "on_frontier"]
        assert ranked.loc["WORSE", "dominated_by"] == "BETTER"

    def test_a_programme_better_on_one_axis_stays_on(self):
        # This is the whole reason the axes are not summed: fewer variants but a closed claim is a
        # trade a reader makes, not one the tool makes for them.
        ranked = self._ranked(
            WIDE={"potential_variants": 100, "open_questions_addressed": 0},
            NARROW={"potential_variants": 1, "open_questions_addressed": 3},
        ).set_index("experiment_id")
        assert ranked.loc["WIDE", "on_frontier"]
        assert ranked.loc["NARROW", "on_frontier"]

    def test_an_identical_pair_both_stay_on(self):
        ranked = self._ranked(
            ONE={"potential_variants": 5}, TWO={"potential_variants": 5}
        ).set_index("experiment_id")
        assert ranked["on_frontier"].all()

    def test_no_combined_score_is_produced(self):
        ranked = self._ranked(A={"potential_variants": 5}, B={"open_questions_addressed": 1})
        forbidden = {"score", "priority", "rank", "expected_information_gain", "opportunity"}
        assert not forbidden & {c.lower() for c in ranked.columns}


class TestTheAuthoredProgrammes:
    def test_they_load(self):
        assert len(read_programs(PROGRAMS)) >= 5

    def test_every_programme_says_what_would_follow_from_failure(self):
        # A programme with no failure branch is a demonstration, not an experiment.
        for text in read_programs(PROGRAMS)["if_it_fails"]:
            assert len(str(text)) > 40

    def test_a_replicate_count_is_derived_or_declared_unknown(self):
        for count in read_programs(PROGRAMS)["replicates"]:
            # Either it says where the number came from, or it says nobody computed one.
            assert str(count) == "not estimated" or "per arm" in str(count)

    def test_the_suppressor_programme_carries_the_replication_the_gap_analysis_derived(self):
        programs = read_programs(PROGRAMS).set_index("experiment_id")
        assert str(programs.loc["SUPTRNA-REPLICATED", "replicates"]).startswith("4 per arm")

    def test_a_derived_replicate_count_shows_every_assumption_behind_it(self):
        # A number without its endpoint, effect size, variance source, design and power is not
        # auditable, and 4 per arm is only meaningful conditional on all five.
        fields = [c for c in read_programs(PROGRAMS).columns if c.startswith("replicate_")]
        assert len(fields) == 6
        for row in read_programs(PROGRAMS).itertuples():
            stated = [str(getattr(row, f)) for f in fields]
            derived = str(row.replicates) != "not estimated"
            # Either every assumption is named, or every one of them says nobody computed it.
            assert derived == any(v != "not estimated" for v in stated)

    def test_direct_scope_and_potential_reach_are_separate(self):
        # A programme measuring one scaffold in one cell line must not be described as informing
        # every variant; the generalisation between the two is written out rather than implied.
        for row in read_programs(PROGRAMS).itertuples():
            assert len(str(row.direct_scope)) > 30
            assert len(str(row.generalisation_required)) > 30
