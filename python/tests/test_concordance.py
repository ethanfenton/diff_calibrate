import pandas as pd

from diff_calibrate.concordance import concordance, concordance_curve
from diff_calibrate.downsample import run_downsampling


def _de_df(genes, pval, logfc):
    return pd.DataFrame(
        {"gene": genes, "pval": pval, "padj": pval, "logFC": logfc}
    )


def test_concordance_identical_results_is_perfect():
    genes = [f"g{i}" for i in range(10)]
    pval = [0.001] * 5 + [0.9] * 5
    logfc = list(range(10))
    de = _de_df(genes, pval, logfc)
    full = _de_df(genes, pval, logfc)

    out = concordance(de, full, alpha=0.05).iloc[0]
    assert out["jaccard"] == 1.0
    assert out["f1"] == 1.0
    assert out["rank_cor"] == 1.0
    assert out["cat_overlap"] == 1.0


def test_concordance_disjoint_sig_sets():
    genes = [f"g{i}" for i in range(10)]
    de = _de_df(genes, [0.001] * 5 + [0.9] * 5, list(range(10)))
    full = _de_df(genes, [0.9] * 5 + [0.001] * 5, list(range(10)))

    out = concordance(de, full, alpha=0.05).iloc[0]
    assert out["jaccard"] == 0.0
    assert out["f1"] == 0.0


def test_concordance_both_empty_sig_sets_is_nan():
    genes = [f"g{i}" for i in range(6)]
    de = _de_df(genes, [0.9] * 6, list(range(6)))
    full = _de_df(genes, [0.9] * 6, list(range(6)))
    out = concordance(de, full, alpha=0.05).iloc[0]
    assert pd.isna(out["jaccard"])
    assert pd.isna(out["f1"])


def test_concordance_curve_shape(sim):
    d = sim(n_genes=60, n_per_group=40, seed=4)
    res = run_downsampling(
        d["counts"], d["group"], "wilcoxon", grid_args={"n_points": 3}, k=2, k_permute=0, seed=1
    )
    curve = concordance_curve(res, alpha=0.05)
    assert len(curve) == len(res["downsampled"])
    assert {"n", "rep", "jaccard", "f1", "rank_cor", "cat_overlap"} <= set(curve.columns)
