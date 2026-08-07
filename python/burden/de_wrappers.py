"""Built-in and pluggable differential expression backends.

A DE backend is any callable ``f(counts, group, **kwargs) -> pd.DataFrame``
returning columns ``gene``, ``pval``, ``padj``, ``logFC``. ``counts`` is a
genes x cells numeric matrix (raw counts, unless the backend expects
otherwise), ``group`` is a length-``n_cells`` array-like giving the two-level
grouping to test. Any function following this contract can be passed as
``de_fun`` to :func:`burden.calculate_burden` or
:func:`burden.downsample.run_downsampling` directly; the wrappers below are
convenience implementations for common methods.
"""

from __future__ import annotations

import warnings
from typing import Callable, Union

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _check_two_groups(group) -> pd.Categorical:
    group = pd.Categorical(group)
    group = group.remove_unused_categories()
    if len(group.categories) != 2:
        raise ValueError(
            f"DE backends require exactly two groups, got {len(group.categories)}: "
            f"{', '.join(map(str, group.categories))}"
        )
    return group


def _gene_names(counts: np.ndarray, index=None) -> list:
    if index is not None:
        return list(index)
    return [f"gene{i + 1}" for i in range(counts.shape[0])]


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape, np.nan)
    mask = ~np.isnan(pvals)
    if mask.sum() > 0:
        out[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
    return out


def de_wilcoxon(counts, group, normalize: bool = True) -> pd.DataFrame:
    """Wilcoxon rank-sum DE.

    Gene-by-gene Wilcoxon rank-sum test on (by default) log1p-CPM normalized
    counts, with BH-adjusted p-values and log2 fold change of group means
    (pseudocount 1) as effect size. Dependency-free fallback DE method; fast
    enough to be practical inside the downsampling loop.

    Parameters
    ----------
    counts : array-like, genes x cells
    group : array-like, length n_cells, exactly two levels
    normalize : bool, CPM-normalize + log1p before testing (default True)

    Returns
    -------
    pd.DataFrame(gene, pval, padj, logFC)
    """
    counts = np.asarray(counts, dtype=float)
    index = getattr(counts, "index", None)
    group = _check_two_groups(group)

    if normalize:
        lib_size = counts.sum(axis=0)
        lib_size[lib_size == 0] = 1
        mat = np.log1p(counts / lib_size * 1e6)
    else:
        mat = counts

    lv = list(group.categories)
    idx1 = np.where(np.asarray(group) == lv[0])[0]
    idx2 = np.where(np.asarray(group) == lv[1])[0]

    n_genes = mat.shape[0]
    pval = np.empty(n_genes)
    for i in range(n_genes):
        x = mat[i, idx1]
        y = mat[i, idx2]
        if np.std(x) == 0 and np.std(y) == 0:
            pval[i] = 1.0
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        pval[i] = p

    logfc = np.log2((mat[:, idx2].mean(axis=1) + 1) / (mat[:, idx1].mean(axis=1) + 1))

    return pd.DataFrame(
        {
            "gene": _gene_names(mat, index),
            "pval": pval,
            "padj": _bh_adjust(pval),
            "logFC": logfc,
        }
    )


def de_edger(counts, group, **kwargs) -> pd.DataFrame:
    """edgeR quasi-likelihood F-test DE (via rpy2).

    Thin wrapper calling out to R's edgeR (``estimateDisp`` + ``glmQLFit`` +
    ``glmQLFTest``) for two-group comparisons. Requires ``rpy2`` and an R
    installation with the ``edgeR`` Bioconductor package (extra: ``burden[edger]``).
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.packages import importr
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "de_edger() requires rpy2 with a working R + edgeR installation. "
            "Install via `pip install burden[edger]` and install edgeR in R "
            "via BiocManager::install('edgeR')."
        ) from e

    counts = np.asarray(counts, dtype=float)
    index = getattr(counts, "index", None)
    group = _check_two_groups(group)
    gene_names = _gene_names(counts, index)

    edger = importr("edgeR")
    with (numpy2ri.converter + pandas2ri.converter).context():
        keep = counts.sum(axis=1) > 0
        r_counts = ro.r.matrix(
            counts[keep, :], nrow=int(keep.sum()), ncol=counts.shape[1]
        )
        r_group = ro.FactorVector(np.asarray(group).astype(str))
        y = edger.DGEList(counts=r_counts, group=r_group)
        y = edger.calcNormFactors(y)
        design = ro.r["model.matrix"](ro.Formula("~group"), data=ro.r("list")(group=r_group))
        y = edger.estimateDisp(y, design)
        fit = edger.glmQLFit(y, design)
        res = edger.glmQLFTest(fit, coef=2)
        tab = ro.r["as.data.frame"](res.rx2("table"))
        tab = pandas2ri.rpy2py(tab)

    kept_names = [g for g, k in zip(gene_names, keep) if k]
    out = pd.DataFrame({"gene": gene_names, "pval": 1.0, "padj": 1.0, "logFC": 0.0})
    idx = out.set_index("gene").index.get_indexer(kept_names)
    out.loc[idx, "pval"] = tab["PValue"].to_numpy()
    out.loc[idx, "logFC"] = tab["logFC"].to_numpy()
    out["padj"] = _bh_adjust(out["pval"].to_numpy())
    return out


def de_deseq2(counts, group, **kwargs) -> pd.DataFrame:
    """DESeq2-style Wald test DE (via pydeseq2).

    Thin wrapper around ``pydeseq2`` for two-group comparisons. Requires
    ``pydeseq2`` (extra: ``burden[deseq2]``). Substantially slower than
    :func:`de_wilcoxon`; expect longer runtimes inside the downsampling loop.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "de_deseq2() requires pydeseq2. Install via `pip install burden[deseq2]`."
        ) from e

    counts = np.asarray(counts, dtype=float)
    index = getattr(counts, "index", None)
    group = _check_two_groups(group)
    gene_names = _gene_names(counts, index)

    keep = counts.sum(axis=1) > 0
    cts = np.round(counts[keep, :]).T  # cells x genes, as pydeseq2 expects
    kept_names = [g for g, k in zip(gene_names, keep) if k]
    counts_df = pd.DataFrame(cts, columns=kept_names)
    metadata = pd.DataFrame({"group": np.asarray(group).astype(str)})

    dds = DeseqDataSet(counts=counts_df, metadata=metadata, design_factors="group", quiet=True)
    dds.deseq2()
    lv = list(group.categories)
    stat_res = DeseqStats(dds, contrast=["group", lv[1], lv[0]], quiet=True, **kwargs)
    stat_res.summary()
    res = stat_res.results_df

    out = pd.DataFrame({"gene": gene_names, "pval": 1.0, "padj": 1.0, "logFC": 0.0})
    idx = out.set_index("gene").index.get_indexer(kept_names)
    out.loc[idx, "pval"] = res["pvalue"].fillna(1.0).to_numpy()
    out.loc[idx, "logFC"] = res["log2FoldChange"].fillna(0.0).to_numpy()
    out["padj"] = _bh_adjust(out["pval"].to_numpy())
    return out


_DE_BACKENDS: dict[str, Callable] = {
    "wilcoxon": de_wilcoxon,
    "edger": de_edger,
    "deseq2": de_deseq2,
}


def resolve_de_fun(de_fun: Union[str, Callable]) -> Callable:
    """Resolve a DE backend from a name or function.

    Parameters
    ----------
    de_fun : str or callable
        Either the name of a built-in backend (``"wilcoxon"``, ``"edger"``,
        ``"deseq2"``) or a user-supplied function following the module-level
        DE backend contract.
    """
    if callable(de_fun):
        return de_fun
    if isinstance(de_fun, str) and de_fun in _DE_BACKENDS:
        return _DE_BACKENDS[de_fun]
    raise ValueError(
        f"de_fun must be a callable, or one of: {', '.join(_DE_BACKENDS)}"
    )
