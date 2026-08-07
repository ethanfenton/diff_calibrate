import numpy as np
import pandas as pd
import pytest

from diff_calibrate.downsample import (
    downsample_counts,
    downsample_grid,
    n_degs,
    run_downsampling,
    run_downsampling_dataset,
)


def test_downsample_grid_includes_max_and_is_sorted():
    grid = downsample_grid(n_max=200, n_min=20, n_points=6)
    assert grid[-1] == 200
    assert list(grid) == sorted(grid)
    assert len(np.unique(grid)) == len(grid)


def test_downsample_grid_below_min_returns_max_only():
    assert list(downsample_grid(n_max=10, n_min=20)) == [10]


def test_downsample_counts_respects_group_sizes():
    counts = np.arange(40).reshape(4, 10).astype(float)
    group = np.array(["A"] * 5 + ["B"] * 5)
    out = downsample_counts(counts, group, n_per_group=3, rng=np.random.default_rng(1))
    assert out["counts"].shape[1] == 6
    vc = pd.Series(np.asarray(out["group"])).value_counts()
    assert vc["A"] == 3 and vc["B"] == 3


def test_downsample_counts_warns_and_keeps_all_when_too_few():
    counts = np.arange(20).reshape(2, 10).astype(float)
    group = np.array(["A"] * 5 + ["B"] * 5)
    with pytest.warns(UserWarning):
        out = downsample_counts(counts, group, n_per_group=8, rng=np.random.default_rng(1))
    vc = pd.Series(np.asarray(out["group"])).value_counts()
    assert vc["A"] == 5 and vc["B"] == 5


def test_n_degs_counts_below_alpha():
    df = pd.DataFrame({"padj": [0.001, 0.2, np.nan, 0.04]})
    assert n_degs(df, alpha=0.05) == 2


def test_run_downsampling_shapes(sim):
    d = sim(n_genes=80, n_per_group=40, seed=2)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 4}, k=2, k_permute=2, seed=1
    )
    assert set(res.keys()) == {"full", "downsampled", "permuted"}
    assert len(res["full"]) == 80
    assert {"n", "rep", "n_degs", "de_res"} <= set(res["downsampled"].columns)
    assert res["downsampled"]["n"].max() == 40
    assert res["permuted"] is not None
    assert len(res["permuted"]) == len(res["downsampled"])


def test_run_downsampling_no_permute(sim):
    d = sim(n_genes=50, n_per_group=30, seed=3)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=0, seed=1
    )
    assert res["permuted"] is None


def test_run_downsampling_dataset(sim):
    data = {
        "ct1": sim(n_genes=50, n_per_group=30, seed=1),
        "ct2": sim(n_genes=50, n_per_group=50, seed=2),
    }
    res = run_downsampling_dataset(
        data, "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=0, seed=1
    )
    assert set(res.keys()) == {"ct1", "ct2"}


def test_run_downsampling_dataset_requires_data():
    with pytest.raises(ValueError):
        run_downsampling_dataset({}, "wilcoxon")
