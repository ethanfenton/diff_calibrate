"""Empirical FDR-calibrated significance threshold as a function of n.

At small cell counts, DE p-values can be miscalibrated (typically
anti-conservative), so a fixed padj < 0.05 cutoff does not carry the same
false-positive rate at every sample size. This calibrates a per-n threshold
alpha(n) from the label-permuted arm of run_downsampling (k_permute > 0):
permuted labels give a DE run with no true group difference, so its padj
distribution is an empirical null. alpha(n) is set to the empirical quantile
of that null padj distribution at target_fpr.

See docs/DESIGN.md (R package) for alternative calibration strategies
(q-values, power-matching) not yet implemented.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .downsample import n_degs


def calibrate_alpha(
    downsampling_result: dict, target_fpr: float = 0.05, cap_at_target: bool = False
) -> pd.DataFrame:
    """Calibrate alpha(n) from the permuted (null) arm of a downsampling result.

    Returns DataFrame(n, alpha, n_null_tests).
    """
    perm = downsampling_result["permuted"]
    if perm is None:
        raise ValueError(
            "downsampling_result has no permuted arm; rerun run_downsampling() "
            "with k_permute > 0."
        )

    rows = []
    for n, group in perm.groupby("n"):
        pooled_padj = np.concatenate([r["padj"].to_numpy() for r in group["de_res"]])
        pooled_padj = pooled_padj[~np.isnan(pooled_padj)]
        if len(pooled_padj) == 0:
            a = target_fpr
        else:
            a = float(np.quantile(pooled_padj, target_fpr, method="linear"))
        if cap_at_target:
            a = min(a, target_fpr)
        rows.append({"n": int(n), "alpha": a, "n_null_tests": len(pooled_padj)})
    return pd.DataFrame(rows).sort_values("n").reset_index(drop=True)


def predict_alpha(alpha_curve: pd.DataFrame, n) -> np.ndarray:
    """Linearly interpolate a calibrated alpha(n) curve to arbitrary cell counts.

    Endpoints are held constant outside the observed range.
    """
    n = np.atleast_1d(np.asarray(n, dtype=float))
    return np.interp(n, alpha_curve["n"].to_numpy(), alpha_curve["alpha"].to_numpy())


def apply_adaptive_alpha(downsampling_result: dict, alpha_curve: pd.DataFrame) -> pd.DataFrame:
    """Recompute DEG counts using a calibrated alpha(n) instead of a fixed alpha.

    Reuses the DE result tables already stored in downsampling_result, so no
    DE refitting is needed. Returns the `downsampled` DataFrame with n_degs
    replaced and a new `alpha_used` column.
    """
    ds = downsampling_result["downsampled"].copy()
    alpha_used = predict_alpha(alpha_curve, ds["n"].to_numpy())
    ds["n_degs"] = [
        n_degs(res, alpha=a) for res, a in zip(ds["de_res"], alpha_used)
    ]
    ds["alpha_used"] = alpha_used
    return ds
