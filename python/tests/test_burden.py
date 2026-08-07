import numpy as np
import pytest

from burden import calculate_burden

pytestmark = pytest.mark.filterwarnings("ignore")


@pytest.fixture
def data(sim):
    return {
        "small": sim(n_genes=200, n_per_group=25, n_de=15, effect=2.5, seed=1),
        "medium": sim(n_genes=200, n_per_group=100, n_de=15, effect=2.5, seed=2),
        "large": sim(n_genes=200, n_per_group=300, n_de=15, effect=2.5, seed=3),
    }


def _run(data, **kwargs):
    defaults = dict(de_fun="wilcoxon", k=3, k_permute=3, grid_args={"n_points": 5}, seed=11)
    defaults.update(kwargs)
    return calculate_burden(data, **defaults)


def test_default_pooling_is_pooled(data):
    res = _run(data)
    assert res.pooling_used == "pooled"
    # pooled model -> identical burden across cell types at the same reference_n
    assert res.burden["burden"].nunique() == 1


def test_burden_output_columns(data):
    res = _run(data)
    expected = {
        "cell_type",
        "n_full",
        "reference_n",
        "n_degs_observed",
        "burden",
        "burden_se",
        "extrapolated",
    }
    assert expected <= set(res.burden.columns)
    assert len(res.burden) == 3


def test_reference_n_defaults_to_median_full_n(data):
    res = _run(data)
    assert res.burden["reference_n"].iloc[0] == pytest.approx(100)


def test_reference_n_override(data):
    res = _run(data, reference_n=50)
    assert (res.burden["reference_n"] == 50).all()


def test_extrapolation_flag_and_warning(data):
    with pytest.warns(UserWarning, match="extrapolated"):
        res = _run(data, reference_n=1000)
    assert res.burden["extrapolated"].all()


def test_pooling_test_always_computed_even_when_pooled_forced(data):
    res = _run(data, pooling="pooled")
    assert res.pooling_used == "pooled"
    assert res.pooling_test["recommendation"] in ("pooled", "shared_shape", "stratified")


def test_stratified_pooling_gives_per_cell_type_models(data):
    res = _run(data, pooling="stratified")
    assert res.pooling_used == "stratified"
    assert set(res.model.keys()) == set(data.keys())


def test_shared_shape_pooling(data):
    res = _run(data, pooling="shared_shape")
    assert res.pooling_used == "shared_shape"


def test_invalid_pooling_raises(data):
    with pytest.raises(ValueError):
        _run(data, pooling="bogus")


def test_adaptive_alpha_requires_k_permute(data):
    with pytest.raises(ValueError):
        _run(data, adaptive_alpha=True, k_permute=0)


def test_adaptive_alpha_false_skips_alpha_curves(data):
    res = _run(data, adaptive_alpha=False)
    assert res.alpha_curves is None
