import numpy as np
import pandas as pd
import pytest

from burden.adaptive_alpha import apply_adaptive_alpha, calibrate_alpha, predict_alpha
from burden.downsample import run_downsampling


def test_calibrate_alpha_requires_permuted_arm(sim):
    d = sim(n_genes=40, n_per_group=30, seed=1)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=0, seed=1
    )
    with pytest.raises(ValueError):
        calibrate_alpha(res)


def test_calibrate_alpha_shape_and_bounds(sim):
    d = sim(n_genes=60, n_per_group=40, seed=2)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=3, seed=1
    )
    curve = calibrate_alpha(res, target_fpr=0.05)
    assert set(curve.columns) == {"n", "alpha", "n_null_tests"}
    assert (curve["alpha"] >= 0).all()
    assert list(curve["n"]) == sorted(curve["n"])


def test_calibrate_alpha_cap_at_target(sim):
    d = sim(n_genes=60, n_per_group=40, seed=3)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=3, seed=1
    )
    curve = calibrate_alpha(res, target_fpr=0.05, cap_at_target=True)
    assert (curve["alpha"] <= 0.05).all()


def test_predict_alpha_interpolates_and_clamps():
    curve = pd.DataFrame({"n": [20, 50, 100], "alpha": [0.01, 0.03, 0.05]})
    out = predict_alpha(curve, [20, 35, 100, 500, 1])
    assert out[0] == pytest.approx(0.01)
    assert out[1] == pytest.approx(0.02)
    assert out[3] == pytest.approx(0.05)  # clamped above range
    assert out[4] == pytest.approx(0.01)  # clamped below range


def test_apply_adaptive_alpha_adds_column(sim):
    d = sim(n_genes=60, n_per_group=40, seed=4)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=3, seed=1
    )
    curve = calibrate_alpha(res)
    out = apply_adaptive_alpha(res, curve)
    assert "alpha_used" in out.columns
    assert len(out) == len(res["downsampled"])
