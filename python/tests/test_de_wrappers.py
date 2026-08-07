import numpy as np
import pytest

from diff_calibrate.de_wrappers import de_wilcoxon, resolve_de_fun


def test_de_wilcoxon_finds_de_genes(sim):
    d = sim(n_genes=150, n_per_group=80, n_de=15, effect=3.0, seed=1)
    res = de_wilcoxon(d["counts"], d["group"])
    assert set(res.columns) == {"gene", "pval", "padj", "logFC"}
    assert len(res) == 150
    sig = set(res.loc[res["padj"] < 0.05, "gene"])
    de_genes = {f"gene{i + 1}" for i in d["de_genes"]}
    # most true DE genes should be recovered at this effect size/sample size
    assert len(sig & de_genes) >= len(de_genes) * 0.5


def test_de_wilcoxon_requires_two_groups():
    counts = np.ones((5, 6))
    with pytest.raises(ValueError):
        de_wilcoxon(counts, ["A"] * 6)


def test_de_wilcoxon_constant_gene_gets_pval_one():
    counts = np.ones((2, 6)) * 5
    group = ["A", "A", "A", "B", "B", "B"]
    res = de_wilcoxon(counts, group, normalize=False)
    assert np.allclose(res["pval"], 1.0)


def test_resolve_de_fun_builtin_and_custom():
    assert resolve_de_fun("wilcoxon") is de_wilcoxon

    def custom(counts, group, **kwargs):
        return None

    assert resolve_de_fun(custom) is custom

    with pytest.raises(ValueError):
        resolve_de_fun("not_a_backend")
