"""DEGs ~ n modeling: nonlinearity test, model fit/predict, pooling test.

Mirrors the R package's approach (see docs/DESIGN.md in the repo root):
quasi-Poisson response family on the *identity* link (so "linear" means
DEGs literally linear in n, matching the natural-scale intuition and giving
safer extrapolation behavior than a log link), fit via `statsmodels` GLM
with `scale="X2"` (Pearson-dispersion quasi-likelihood scaling), and nested
models compared with a quasi-F test analogous to R's
`anova.glm(..., test="F")` / `anova.gam(..., test="F")`.

One deliberate difference from the R implementation: splines here are
unpenalized natural cubic regression-spline bases (`patsy`'s `cr()`, the
same basis family `mgcv`'s `bs="cr"` uses) rather than `mgcv`'s penalized
smooths. This package does not shrink a spline's effective df below its
nominal df the way a penalized GAM can, so the "spline collapses back to
linear" degenerate-df failure mode from the R implementation is less likely
here, but the `min_df` guard is kept for parity and safety.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy.mgcv_cubic_splines import cr  # noqa: F401  (used inside formula strings)
from scipy import stats


def _spline_k(n_unique: int, k_max: int = 8) -> int:
    return max(3, min(k_max, n_unique - 1))


def _family():
    """Quasi-Poisson response family: Poisson deviance/variance, identity link."""
    return sm.families.Poisson(link=sm.families.links.Identity())


def _fit_glm(formula: str, data: pd.DataFrame, weights: np.ndarray):
    """Fit a GLM, suppressing the harmless identity-link domain warning.

    The identity-link Poisson deviance can transiently evaluate at a
    negative mean during IRLS iteration; the optimizer recovers, but the
    warning is noisy for every caller (mirrors `.gam_quiet()` in the R
    implementation).
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=sm.tools.sm_exceptions.DomainWarning)
        warnings.filterwarnings("ignore", message="invalid value encountered")
        model = smf.glm(formula, data=data, family=_family(), var_weights=weights)
        return model.fit(scale="X2")


def quasi_f_test(result_null, result_alt, min_df: float = 0.4) -> dict:
    """Quasi-F nested-model test between two GLM fits (result_null nested in result_alt).

    F = ((D0 - D1) / (df0 - df1)) / dispersion_alt, dispersion_alt from the
    (larger) alternative model's Pearson-chi2-based scale estimate --
    mirrors R's `anova.glm(..., test="F")` convention. Returns
    dict(p_value, f_stat, df).
    """
    d0, d1 = result_null.deviance, result_alt.deviance
    df0, df1 = result_null.df_resid, result_alt.df_resid
    df_num = df0 - df1
    df_den = df1
    if df_den <= 0:
        return {"p_value": np.nan, "f_stat": np.nan, "df": df_num}
    if df_num < min_df:
        return {"p_value": 1.0, "f_stat": np.nan, "df": df_num}
    dispersion = result_alt.scale
    f_stat = ((d0 - d1) / df_num) / dispersion
    p_value = float(stats.f.sf(f_stat, df_num, df_den)) if f_stat >= 0 else 1.0
    return {"p_value": p_value, "f_stat": float(f_stat), "df": float(df_num)}


def test_nonlinearity(
    n: Sequence[float],
    degs: Sequence[float],
    weights: Optional[Sequence[float]] = None,
    cell_type: Optional[Sequence] = None,
    sig_level: float = 0.05,
) -> dict:
    """Test whether a spline term is needed for DEGs ~ n.

    Fits a linear quasi-Poisson GLM (`degs ~ n`) and a spline quasi-Poisson
    GLM (`degs ~ cr(n, df=k)`), both via `statsmodels`, and compares them
    with :func:`quasi_f_test`. A significant result (by default p < 0.05)
    indicates the DEGs~n relationship isn't well captured by a straight
    line, and a spline should be used.

    If `cell_type` is supplied, both models include `+ C(cell_type)` as an
    additive term (shared curve shape, per-cell-type intercept) -- use this
    to run the nonlinearity test on a `shared_shape` pooled fit rather than
    a single cell type.

    Returns dict(p_value, decision, model_linear, model_spline), decision
    one of "spline" / "linear".
    """
    n = np.asarray(n, dtype=float)
    degs = np.asarray(degs, dtype=float)
    w = np.ones_like(n) if weights is None else np.asarray(weights, dtype=float)
    k = _spline_k(len(np.unique(n)))

    df = pd.DataFrame({"n": n, "degs": degs})
    if cell_type is None:
        f_linear = "degs ~ n"
        f_spline = f"degs ~ cr(n, df={k})"
    else:
        df["cell_type"] = pd.Categorical(cell_type)
        f_linear = "degs ~ n + C(cell_type)"
        f_spline = f"degs ~ cr(n, df={k}) + C(cell_type)"

    model_linear = _fit_glm(f_linear, df, w)
    model_spline = _fit_glm(f_spline, df, w)

    lrt = quasi_f_test(model_linear, model_spline)
    decision = "spline" if (not np.isnan(lrt["p_value"]) and lrt["p_value"] < sig_level) else "linear"

    return {
        "p_value": lrt["p_value"],
        "decision": decision,
        "model_linear": model_linear,
        "model_spline": model_spline,
    }


@dataclass
class BurdenModel:
    model: object
    model_type: str
    decision: str
    nonlinearity_test: dict


def fit_deg_model(
    n: Sequence[float],
    degs: Sequence[float],
    weights: Optional[Sequence[float]] = None,
    cell_type: Optional[Sequence] = None,
    model_type: str = "auto",
    sig_level: float = 0.05,
) -> BurdenModel:
    """Fit the DEGs ~ n model for a single cell type (or a shared-shape pooled fit).

    model_type: "auto" (default; run test_nonlinearity() and pick), "spline"
    (force spline), "linear" (force linear).

    Returns a BurdenModel(model, model_type, decision, nonlinearity_test),
    nonlinearity_test always computed (even when model_type is forced) so
    the decision that *would* have been made is visible.
    """
    if model_type not in ("auto", "spline", "linear"):
        raise ValueError('model_type must be one of "auto", "spline", "linear"')

    nl_test = test_nonlinearity(n, degs, weights=weights, cell_type=cell_type, sig_level=sig_level)

    decision = {"auto": nl_test["decision"], "spline": "spline", "linear": "linear"}[model_type]
    model = nl_test["model_spline"] if decision == "spline" else nl_test["model_linear"]

    return BurdenModel(model=model, model_type=model_type, decision=decision, nonlinearity_test=nl_test)


def predict_degs(
    fit: Union[BurdenModel, object],
    n_new: Union[float, Sequence[float]],
    cell_type_new: Optional[Union[str, Sequence]] = None,
    se: bool = False,
) -> Union[np.ndarray, pd.DataFrame]:
    """Predict expected DEG count at a given cell count.

    Parameters
    ----------
    fit : BurdenModel (from fit_deg_model) or a raw fitted GLMResults.
    n_new : cell count(s) to predict at.
    cell_type_new : optional cell type label(s), required (recycled to
        len(n_new)) if fit's model includes a cell_type term.
    se : if True, also return the standard error of the prediction.

    Returns a numpy array of predicted DEG counts, or if se=True, a
    DataFrame(n, fit, se).
    """
    model_result = fit.model if isinstance(fit, BurdenModel) else fit
    n_new = np.atleast_1d(np.asarray(n_new, dtype=float))
    new_df = pd.DataFrame({"n": n_new})
    if cell_type_new is not None:
        ct = np.broadcast_to(np.asarray(cell_type_new, dtype=object), n_new.shape)
        new_df["cell_type"] = ct

    pred = model_result.get_prediction(new_df)
    if not se:
        return np.asarray(pred.predicted_mean)
    summary = pred.summary_frame()
    return pd.DataFrame(
        {"n": n_new, "fit": summary["mean"].to_numpy(), "se": summary["mean_se"].to_numpy()}
    )


def test_pooling(
    n: Sequence[float],
    degs: Sequence[float],
    cell_type: Sequence,
    weights: Optional[Sequence[float]] = None,
    sig_level: float = 0.05,
) -> dict:
    """Test whether a pooled (across cell type) model is justified.

    Compares three nested GLM fits of DEGs ~ n across all cell types:
    - pooled: `degs ~ cr(n, df=k)` -- single shared curve, no cell-type effect.
    - shared_shape: `degs ~ cr(n, df=k) + C(cell_type)` -- shared curve
      shape, per-cell-type intercept.
    - stratified: `degs ~ cr(n, df=k) * C(cell_type)` -- fully separate curve
      per cell type.

    Two quasi-F tests: pooled vs shared_shape (offset needed?), shared_shape
    vs stratified (curve shape needed, beyond an offset?).

    Returns dict(recommendation, offset_test, shape_test, model_pooled,
    model_shared_shape, model_stratified). recommendation is one of
    "pooled", "shared_shape", "stratified".
    """
    n = np.asarray(n, dtype=float)
    degs = np.asarray(degs, dtype=float)
    w = np.ones_like(n) if weights is None else np.asarray(weights, dtype=float)
    k = _spline_k(len(np.unique(n)))

    df = pd.DataFrame({"n": n, "degs": degs, "cell_type": pd.Categorical(cell_type)})

    model_pooled = _fit_glm(f"degs ~ cr(n, df={k})", df, w)
    model_shared_shape = _fit_glm(f"degs ~ cr(n, df={k}) + C(cell_type)", df, w)
    model_stratified = _fit_glm(f"degs ~ cr(n, df={k}) * C(cell_type)", df, w)

    offset_test = quasi_f_test(model_pooled, model_shared_shape)
    shape_test = quasi_f_test(model_shared_shape, model_stratified)

    if np.isnan(offset_test["p_value"]) or offset_test["p_value"] >= sig_level:
        recommendation = "pooled"
    elif np.isnan(shape_test["p_value"]) or shape_test["p_value"] >= sig_level:
        recommendation = "shared_shape"
    else:
        recommendation = "stratified"

    return {
        "recommendation": recommendation,
        "offset_test": offset_test,
        "shape_test": shape_test,
        "model_pooled": model_pooled,
        "model_shared_shape": model_shared_shape,
        "model_stratified": model_stratified,
    }
