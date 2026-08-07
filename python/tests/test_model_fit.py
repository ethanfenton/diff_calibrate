import warnings

import numpy as np
import pytest

from diff_calibrate.model_fit import fit_deg_model, predict_degs, quasi_f_test
from diff_calibrate.model_fit import test_nonlinearity as run_nonlinearity_test
from diff_calibrate.model_fit import test_pooling as run_pooling_test

pytestmark = pytest.mark.filterwarnings("ignore")


def test_nonlinearity_picks_linear_for_linear_data():
    rng = np.random.default_rng(1)
    n = np.tile(np.arange(20, 220, 20), 8)
    degs = 0.6 * n + rng.normal(0, 3, size=len(n))
    out = run_nonlinearity_test(n, degs, sig_level=0.05)
    assert out["decision"] == "linear"
    assert out["p_value"] > 0.05


def test_nonlinearity_picks_spline_for_saturating_data():
    rng = np.random.default_rng(2)
    n = np.tile(np.arange(20, 220, 20), 8)
    degs = 80 * np.log1p(n) + rng.normal(0, 2, size=len(n))
    out = run_nonlinearity_test(n, degs, sig_level=0.05)
    assert out["decision"] == "spline"
    assert out["p_value"] < 0.05


def test_fit_deg_model_forced_model_type():
    rng = np.random.default_rng(3)
    n = np.tile(np.arange(20, 220, 20), 8)
    degs = 80 * np.log1p(n) + rng.normal(0, 2, size=len(n))

    forced_linear = fit_deg_model(n, degs, model_type="linear")
    assert forced_linear.decision == "linear"
    # the nonlinearity test is always computed even when forced
    assert forced_linear.nonlinearity_test["decision"] == "spline"

    forced_spline = fit_deg_model(n, degs, model_type="spline")
    assert forced_spline.decision == "spline"


def test_predict_degs_returns_fit_and_se():
    rng = np.random.default_rng(4)
    n = np.tile(np.arange(20, 220, 20), 8)
    degs = 0.6 * n + rng.normal(0, 3, size=len(n))
    model = fit_deg_model(n, degs, model_type="linear")
    out = predict_degs(model, [100], se=True)
    assert out["fit"].iloc[0] == pytest.approx(60, abs=15)
    assert out["se"].iloc[0] > 0

    point = predict_degs(model, [50, 150])
    assert len(point) == 2


def test_pooling_picks_pooled_for_shared_relationship():
    rng = np.random.default_rng(7)
    n = np.tile(np.arange(20, 220, 20), 10)
    n_all = np.concatenate([n, n])
    degs_all = 0.6 * n_all + rng.normal(0, 3, size=len(n_all))
    cell_type = np.array(["a"] * len(n) + ["b"] * len(n))

    out = run_pooling_test(n_all, degs_all, cell_type, sig_level=0.05)
    assert out["recommendation"] == "pooled"
    assert out["offset_test"]["p_value"] > 0.05


def test_pooling_picks_stratified_for_different_intercepts_and_shapes():
    rng = np.random.default_rng(8)
    n = np.tile(np.arange(20, 220, 20), 10)
    degs_a = 0.3 * n + rng.normal(0, 2, size=len(n))
    degs_b = 40 * np.log1p(n) + 30 + rng.normal(0, 2, size=len(n))
    n_all = np.concatenate([n, n])
    degs_all = np.concatenate([degs_a, degs_b])
    cell_type = np.array(["a"] * len(n) + ["b"] * len(n))

    out = run_pooling_test(n_all, degs_all, cell_type, sig_level=0.05)
    assert out["recommendation"] == "stratified"


def test_quasi_f_test_degenerate_df_short_circuits():
    rng = np.random.default_rng(9)
    n = np.tile(np.arange(20, 220, 20), 5)
    degs = 0.5 * n + rng.normal(0, 2, size=len(n))
    fit = fit_deg_model(n, degs, model_type="linear")
    # comparing a model against itself: df diff is 0, below min_df
    out = quasi_f_test(fit.model, fit.model)
    assert out["p_value"] == 1.0
