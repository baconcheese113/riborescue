import numpy as np
import pandas as pd
import pytest

from riborescue.baseline import identifiable_columns, r_squared


def test_independent_columns_are_all_kept():
    design = pd.DataFrame({"a": [1.0, 0.0, 0.0], "b": [0.0, 1.0, 0.0], "c": [0.0, 0.0, 1.0]})
    assert identifiable_columns(design) == ["a", "b", "c"]


def test_an_aliased_column_is_dropped_and_the_earlier_one_wins():
    design = pd.DataFrame({"a": [1.0, 0.0, 1.0], "b": [0.0, 1.0, 0.0], "sum": [1.0, 1.0, 1.0]})
    assert identifiable_columns(design) == ["a", "b"]


def test_a_duplicated_column_is_dropped():
    design = pd.DataFrame({"a": [1.0, 2.0, 3.0], "copy": [1.0, 2.0, 3.0]})
    assert identifiable_columns(design) == ["a"]


def test_an_empty_cell_leaves_no_identifiable_column():
    design = pd.DataFrame({"a": [1.0, 1.0], "absent": [0.0, 0.0]})
    assert identifiable_columns(design) == ["a"]


def test_r_squared_is_the_squared_pearson_correlation():
    predicted = pd.Series([0.1, 0.2, 0.3, 0.45])
    observed = pd.Series([0.12, 0.18, 0.35, 0.4])
    assert r_squared(predicted, observed) == float(np.corrcoef(predicted, observed)[0, 1] ** 2)


def test_a_perfect_prediction_scores_one():
    values = pd.Series([0.01, 0.02, 0.03])
    assert r_squared(values, values) == pytest.approx(1.0)
