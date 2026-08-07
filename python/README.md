# diff_calibrate (Python)

Python port of the `diffcalibrate` R package: sample-size-adjusted
differential expression (DE) "burden" for single-cell data. Cell types with
more cells tend to show more significant DE genes for reasons of statistical
power alone, which biases naive DEG counts toward larger cell types.
`diff_calibrate` downsamples each cell type across a range of cell counts,
refits DE at each, and models DEG count as a function of cell count — so you
can ask "how many DEGs would this cell type show if it had as many cells as
a typical cell type in this dataset," rather than comparing raw DEG counts
across cell types of very different size.

See [`../docs/DESIGN.md`](../docs/DESIGN.md) for the full methodology,
rationale, and known caveats (in particular: extrapolation risk when a cell
type is much smaller than the reference cell count). The R and Python
implementations share the same design and defaults; see "Differences from
the R package" below for the few places they diverge for
language/library reasons.

## Install

```bash
# from a local checkout
pip install -e "path/to/diff_calibrate/python"

# optional DE backends
pip install -e "path/to/diff_calibrate/python[edger]"    # requires rpy2 + R + edgeR
pip install -e "path/to/diff_calibrate/python[deseq2]"   # requires pydeseq2
```

## Usage

```python
import diff_calibrate

# data: dict, one entry per cell type, each
#   {"counts": <genes x cells matrix>, "group": <two-level treatment vector>}
data = {
    "neurons":    {"counts": neuron_counts,    "group": neuron_group},
    "astrocytes": {"counts": astrocyte_counts, "group": astrocyte_group},
    "microglia":  {"counts": microglia_counts, "group": microglia_group},
}

res = diff_calibrate.calculate_burden(
    data,
    de_fun="wilcoxon",       # or "edger", "deseq2", or your own function
    model_type="auto",       # "auto" | "spline" | "linear"
    pooling="pooled",        # "pooled" | "auto" | "shared_shape" | "stratified"
    adaptive_alpha=True,     # empirical FDR-calibrated alpha(n) per cell type
    alpha=0.05,
)

res.burden          # cell_type, n_full, reference_n, burden, burden_se, extrapolated
res.pooling_test     # tests behind the pooled/shared_shape/stratified recommendation
res.alpha_curves      # per-cell-type alpha(n) from empirical FDR calibration
res.concordance        # per-cell-type concordance of downsampled vs full DE calls
```

`de_fun` can be any callable following the contract in
`diff_calibrate.de_wrappers`: `f(counts, group, **kwargs) ->
pd.DataFrame(gene, pval, padj, logFC)`, so any DE method (not just the three
built-ins) can be plugged in.

## Differences from the R package

- **Spline basis.** The R package uses `mgcv`'s penalized thin-plate/cubic
  regression smooths. This port uses `patsy`'s unpenalized natural cubic
  regression-spline basis (`cr()`, same underlying basis family as `mgcv`'s
  `bs="cr"`, just without the smoothing penalty). Practically: the Python
  spline won't shrink toward a straight line the way a penalized GAM can, so
  `test_nonlinearity()`'s degenerate-df edge case (a spline that collapses
  onto its nested linear model) is rarer here, though the `min_df` guard is
  kept for parity.
- **Quasi-Poisson fitting.** R fits via `mgcv::gam(family = quasipoisson(link
  = "identity"))`; Python fits a `statsmodels` `GLM` with
  `family=Poisson(link=Identity())` and `scale="X2"` (Pearson-dispersion
  scaling), which is the standard way to get quasi-Poisson behavior in
  `statsmodels`. Nested-model comparisons use a quasi-F test
  (`diff_calibrate.model_fit.quasi_f_test`), the same statistic family as
  R's `anova(..., test="F")`.
- **edgeR/DESeq2 backends.** `de_edger()` calls out to R via `rpy2` (needs a
  local R + Bioconductor edgeR install); `de_deseq2()` uses the `pydeseq2`
  port rather than calling R's DESeq2 directly. Both are optional extras;
  `de_wilcoxon()` has no non-Python dependencies.

## Development

```bash
pip install -e ".[test]"
pytest
```
