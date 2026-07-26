"""The ledger's shape, and the distinctions it exists to keep.

Most of these guard against the ledger degrading into a list of successes: a row that loses its
verdict, a verdict invented to avoid an awkward one, a claim quietly recorded twice under two names.
"""

from pathlib import Path

import pandas as pd
import pytest

from riborescue.core.falsification import LEDGER, VERDICTS, read_ledger, summarise

_ROW = {
    "claim_id": "EXAMPLE",
    "claim": "a thing was predicted",
    "assay": "some measurement",
    "endpoint": "what was compared",
    "comparability": "direct",
    "verdict": "supported",
    "detail": "what happened",
    "what_would_settle_it": "the next experiment",
    "source": "ADR-0000",
}


def _written(tmp_path: Path, *rows: dict[str, str]) -> Path:
    path = tmp_path / "ledger.tsv"
    pd.DataFrame(list(rows)).to_csv(path, sep="\t", index=False)
    return path


class TestTheContract:
    def test_a_row_missing_a_column_is_refused(self, tmp_path):
        path = _written(tmp_path, {k: v for k, v in _ROW.items() if k != "what_would_settle_it"})
        with pytest.raises(ValueError, match="missing what_would_settle_it"):
            read_ledger(path)

    def test_a_verdict_outside_the_four_is_refused(self, tmp_path):
        path = _written(tmp_path, _ROW | {"verdict": "promising"})
        with pytest.raises(ValueError, match="verdicts outside"):
            read_ledger(path)

    def test_a_blank_cell_is_refused_rather_than_read_as_absent(self, tmp_path):
        # A ledger with holes reads as complete while omitting the thing a reader came for.
        path = _written(tmp_path, _ROW | {"what_would_settle_it": ""})
        with pytest.raises(ValueError, match="leaves what_would_settle_it empty"):
            read_ledger(path)

    def test_the_same_claim_cannot_be_recorded_twice(self, tmp_path):
        path = _written(tmp_path, _ROW, _ROW | {"verdict": "refuted"})
        with pytest.raises(ValueError, match="reuses a claim_id"):
            read_ledger(path)

    def test_a_supported_claim_still_has_to_name_its_next_experiment(self, tmp_path):
        # The column is required on every verdict, so one dataset cannot become proof by omission.
        assert read_ledger(_written(tmp_path, _ROW))["what_would_settle_it"].iloc[0] != ""


class TestTheSummary:
    def test_every_verdict_appears_even_with_no_rows_under_it(self, tmp_path):
        counts = summarise(read_ledger(_written(tmp_path, _ROW)))
        assert list(counts["verdict"]) == list(VERDICTS)
        assert counts.set_index("verdict").loc["refuted", "rows"] == 0

    def test_the_counts_are_the_rows(self, tmp_path):
        rows = [_ROW, _ROW | {"claim_id": "B", "verdict": "untestable"}]
        counts = summarise(read_ledger(_written(tmp_path, *rows))).set_index("verdict")["rows"]
        assert counts["supported"] == 1
        assert counts["untestable"] == 1


class TestTheCommittedLedger:
    def test_it_loads_and_is_well_formed(self):
        assert len(read_ledger(LEDGER)) > 0

    def test_the_untestable_row_carries_what_would_make_it_testable(self):
        ledger = read_ledger(LEDGER)
        for settle in ledger.loc[ledger["verdict"] == "untestable", "what_would_settle_it"]:
            # `untestable` is a data gap, so the row is only useful if it says what closes it.
            assert len(str(settle)) > 40

    def test_the_suppressor_trna_claim_is_untestable_and_not_refuted(self):
        # The design could not have resolved it either way, which is a different statement from
        # evidence that the suppressor tRNA does not work.
        ledger = read_ledger(LEDGER).set_index("claim_id")
        assert ledger.loc["IDUA-W402X-SUPTRNA", "verdict"] == "untestable"

    def test_the_negative_control_is_recorded_beside_the_positive_one(self):
        # A positive control alone shows a response; the pair shows a discrimination.
        ledger = read_ledger(LEDGER).set_index("claim_id")
        assert ledger.loc["G418-READTHROUGH", "verdict"] == "supported"
        assert ledger.loc["SRI37240-STALLING", "verdict"] == "supported"
