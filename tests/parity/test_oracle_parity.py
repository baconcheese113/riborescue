"""Reproduction parity against the R oracle.

The Python baseline is judged against the authors' own fold assignments and predictions, exported as
committed golden fixtures, not against numbers copied from the paper. Matching a summary statistic
is weak evidence; producing identical predictions on identical folds is strong evidence.

The check runs unconditionally and fails — never skips — if a fixture is absent, because a skipped
check is a check that has quietly vanished. Regenerate the fixtures with `pixi run oracle`.
"""

from pathlib import Path

import pandas as pd
import pytest

from riborescue.baseline import cross_validate, fit_fold
from riborescue.tables import read_table

ORACLE = Path(__file__).resolve().parents[1] / "fixtures" / "oracle"
DRUGS = ("CC90009", "Clitocine", "DAP", "G418", "SJ6986", "SRI")

# The fits agree to floating-point noise, so the tolerance is tight enough that a real difference in
# encoding, fold membership or link function cannot hide beneath it.
TOLERANCE = 1e-8


def _fixture(kind: str, drug: str) -> pd.DataFrame:
    path = ORACLE / f"{kind}_{drug}.tsv.gz"
    assert path.exists(), f"{path.name} is absent — regenerate the fixtures with `pixi run oracle`"
    return read_table(path)


def _features(drug: str) -> pd.DataFrame:
    return _fixture("features", drug).set_index("row")


def _oracle_round(drug: str, round_: int) -> pd.DataFrame:
    predictions = _fixture("predictions", drug)
    return predictions[predictions["round"] == round_].set_index("row")


def _metric(drug: str, round_: int) -> float:
    metrics = read_table(ORACLE / "metrics.tsv")
    row = metrics[(metrics["drug"] == drug) & (metrics["round"] == round_)]
    return float(row["r2"].iloc[0])


@pytest.mark.parity
@pytest.mark.parametrize("drug", DRUGS)
def test_python_baseline_reproduces_the_oracle_predictions(drug: str):
    folds = _fixture("folds", drug)
    first = folds[folds["round"] == 1]["row"]
    fit = fit_fold(_features(drug), first, 1)
    oracle = _oracle_round(drug, 1)

    assert list(fit.predictions.index) == list(oracle.index)
    assert (fit.predictions - oracle["predicted"]).abs().max() < TOLERANCE
    assert (fit.observed - oracle["observed"]).abs().max() < TOLERANCE
    assert abs(fit.r2 - _metric(drug, 1)) < TOLERANCE


@pytest.mark.parity
def test_the_oracle_records_the_provenance_of_its_fixtures():
    provenance = ORACLE / "provenance.json"
    assert provenance.exists(), "provenance.json is absent — regenerate with `pixi run oracle`"
    recorded = provenance.read_text()
    assert "lehner-lab/Stop_codon_readthrough" in recorded
    assert "bb3474efa01f074466912a36964827c2" in recorded


@pytest.mark.parity
@pytest.mark.slow
@pytest.mark.parametrize("drug", DRUGS)
def test_python_baseline_reproduces_every_cross_validation_round(drug: str):
    features = _features(drug)
    for fit in cross_validate(features, _fixture("folds", drug)):
        oracle = _oracle_round(drug, fit.round)
        assert (fit.predictions - oracle["predicted"]).abs().max() < TOLERANCE
        assert abs(fit.r2 - _metric(drug, fit.round)) < TOLERANCE
