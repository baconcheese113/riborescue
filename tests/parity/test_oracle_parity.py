"""Reproduction parity against the R oracle.

The Python baseline is judged against the authors' own fold assignments and predictions, exported as
committed golden fixtures, not against numbers copied from the paper. The check runs unconditionally
and stays `xfail(strict)` until the baseline reproduction exists; it never skips, because a skipped
control is a control that has quietly vanished.
"""

from pathlib import Path

import pytest

ORACLE = Path(__file__).resolve().parents[1] / "fixtures" / "oracle"


@pytest.mark.parity
@pytest.mark.xfail(strict=True, reason="baseline reproduction lands with the R oracle parity check")
def test_python_baseline_matches_oracle_folds():
    predictions = ORACLE / "fold_predictions.parquet"
    assert predictions.exists(), "oracle fixtures missing — run `pixi run oracle`"
    raise NotImplementedError("parity comparison lands with the baseline reproduction")
