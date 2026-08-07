"""Concordance between downsampled and full-data DE results."""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats


def concordance(
    de_res: pd.DataFrame, full_res: pd.DataFrame, alpha: float = 0.05, cat_k: Optional[int] = None
) -> pd.DataFrame:
    """Concordance measures between a downsampled DE result and the full-data result.

    Restricted to genes present in both tables:
    - jaccard: Jaccard index of DEG sets at alpha (NaN if both sets empty).
    - f1: F1 score of the downsampled DEG set against the full set as truth.
    - rank_cor: Spearman correlation of -log10(pval) between the two runs.
    - cat_overlap: concordance-at-the-top, fraction overlap of the top-cat_k
      genes by abs(logFC) in each run (cat_k defaults to min(50, n_genes // 4)).

    Returns a one-row DataFrame(jaccard, f1, rank_cor, cat_overlap).
    """
    common = np.intersect1d(de_res["gene"].to_numpy(), full_res["gene"].to_numpy())
    ds = de_res.set_index("gene").loc[common]
    fl = full_res.set_index("gene").loc[common]

    ds_sig = set(ds.index[(ds["padj"].notna()) & (ds["padj"] < alpha)])
    fl_sig = set(fl.index[(fl["padj"].notna()) & (fl["padj"] < alpha)])

    inter = len(ds_sig & fl_sig)
    union_n = len(ds_sig | fl_sig)
    jaccard = np.nan if union_n == 0 else inter / union_n

    tp = inter
    fp = len(ds_sig - fl_sig)
    fn = len(fl_sig - ds_sig)
    f1 = np.nan if (2 * tp + fp + fn) == 0 else (2 * tp) / (2 * tp + fp + fn)

    if len(common) >= 3:
        ds_p = -np.log10(np.maximum(ds["pval"].to_numpy(), np.finfo(float).tiny))
        fl_p = -np.log10(np.maximum(fl["pval"].to_numpy(), np.finfo(float).tiny))
        mask = ~(np.isnan(ds_p) | np.isnan(fl_p))
        if mask.sum() >= 3:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=stats.ConstantInputWarning)
                rank_cor = stats.spearmanr(ds_p[mask], fl_p[mask]).statistic
        else:
            rank_cor = np.nan
    else:
        rank_cor = np.nan

    if cat_k is None:
        cat_k = max(1, min(50, len(common) // 4))
    top_ds = set(ds["logFC"].abs().sort_values(ascending=False).index[: min(cat_k, len(ds))])
    top_fl = set(fl["logFC"].abs().sort_values(ascending=False).index[: min(cat_k, len(fl))])
    cat_overlap = len(top_ds & top_fl) / cat_k

    return pd.DataFrame(
        [{"jaccard": jaccard, "f1": f1, "rank_cor": rank_cor, "cat_overlap": cat_overlap}]
    )


def concordance_curve(
    downsampling_result: dict, alpha: float = 0.05, cat_k: Optional[int] = None
) -> pd.DataFrame:
    """Apply concordance() to every replicate in a run_downsampling() result.

    Returns DataFrame(n, rep, jaccard, f1, rank_cor, cat_overlap).
    """
    ds = downsampling_result["downsampled"]
    full_res = downsampling_result["full"]
    rows = []
    for _, row in ds.iterrows():
        conc = concordance(row["de_res"], full_res, alpha=alpha, cat_k=cat_k).iloc[0]
        rows.append({"n": row["n"], "rep": row["rep"], **conc.to_dict()})
    return pd.DataFrame(rows)
