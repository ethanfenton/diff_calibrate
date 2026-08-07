# `burden` — design notes

Package for quantifying sample-size bias in single-cell differential expression
(DE): larger cell types get more DEGs partly because they have more cells, not
because they are biologically more affected. `burden` estimates, per cell type,
the expected number of DEGs *after adjusting for cell count*, by downsampling,
modeling DEGs as a function of n, and predicting at a common reference n.

This file is a living planning document. Revisit the "future work / alternatives"
sections periodically — several design choices below were picked as reasonable
starting points, not because alternatives were ruled out.

## Pipeline

1. **Downsample.** For each cell type x treatment arm, downsample to a grid of
   target cell counts (log-spaced, from a configurable minimum up to the full
   count), with `k` replicate draws per grid point (without replacement).
   Also keep the full-data DE run per cell type as ground truth.
2. **Run DE.** Pluggable DE backend: user supplies (or selects a built-in
   wrapper for) a function `(counts, group_labels, ...) -> data.frame(gene,
   pval, padj, logFC)`. Built-ins: Wilcoxon (rank-sum, `presto`-style or base
   `wilcox.test`), edgeR (QLF/LRT), DESeq2. Keeps downsampling/modeling logic
   method-agnostic.
3. **Record DEG counts** per (cell type, n, replicate) at nominal alpha (and,
   for the adaptive-alpha path, at the calibrated alpha(n) — see below).
4. **Concordance.** At each downsampled run, compare its DE result to the
   full-data DE result for the same cell type:
   - Set overlap: Jaccard / F1 of DEG sets at matched alpha.
   - Rank correlation of -log10(p) or logFC between downsampled and full runs.
   - Concordance-at-the-top (CAT): overlap of top-K genes by effect size across
     a range of K — alpha-free, robust choice.
   Produces a concordance(n) curve per cell type, used to identify sample
   sizes too small to trust (either hard-excluded from the model or used as a
   weight, so unstable/low-concordance points get down-weighted rather than
   deleted outright — lets the spline "see" the instability rather than being
   blind to it).
5. **Model DEGs ~ n.** Count model (quasi-Poisson or negative binomial; beta-
   binomial if bounding by total genes tested is wanted) regressing DEG count
   on cell count, per cell type or pooled (see below). `model_type` argument:
   - `"auto"`: fit linear and spline/GAM models, LRT between nested fits
     (spline edf vs linear), pick spline if p < threshold (configurable).
   - `"spline"`: force GAM (`mgcv::gam(degs ~ s(n), family = nb())`).
   - `"linear"`: force parametric GLM.
   Decision + test statistic are always retained in output for transparency,
   even when the choice is forced.
6. **Pooled vs cell-type-specific model.** Same nested-model-comparison
   pattern as step 5, applied to whether the n->DEG relationship needs to vary
   by cell type: compare `degs ~ s(n, by = cell_type)` (or fully separate
   fits) vs `degs ~ s(n) + cell_type` (shared shape, different intercept) vs
   fully pooled `degs ~ s(n)`, via LRT/AIC. A middle ground worth using when
   cell types are sparse: cell type as a random effect (`s(cell_type,
   bs="re")` in mgcv) — partial pooling instead of an all-or-nothing choice.

   **Default is `pooled`, not `auto`.** It is a reasonable working
   assumption that cell count affects detected DEG count in a broadly
   similar way regardless of cell type (same DE method, same underlying
   statistical-power mechanics), so pooling all cell types' downsampling
   replicates into one DEGs~n curve by default is both defensible and gives
   the model far more data per fit than any single cell type could supply
   alone. Letting `pooling = "auto"` pick per-dataset via the LRT sounds
   more rigorous but is not the safer default: with the modest replicate
   counts (`k`) typical of this pipeline, a single noisy LRT can
   false-positive into recommending `shared_shape`/`stratified` even when
   the true relationship is shared across cell types — this was observed
   directly during development (a `test_pooling()` unit test using
   identically-generated data for two "cell types" intermittently
   recommended `shared_shape` depending on the random draw). `stratified`
   also compounds with the extrapolation caveat above: a per-cell-type
   model fit to few points is exactly the situation where predicting at a
   reference n outside that cell type's own range is least trustworthy.
   `pooling = "auto"`/`"stratified"` remain available for users who want
   per-cell-type flexibility and are prepared for the added noise;
   `pooling_test` is always computed and returned regardless of which mode
   is used, so `"auto"`'s recommendation is visible even under the
   `"pooled"` default.
7. **Burden metric.** Burden for a cell type = model-predicted DEG count at a
   *reference* cell count (default: median cells-per-treatment-arm across all
   cell types in the dataset; user-overridable), evaluated on the cell type's
   own fitted curve. I.e., "how many DEGs would this cell type have shown if
   it had as many cells as a typical cell type in this dataset," using the
   model fit (pooled or cell-type-specific, per step 6) to extrapolate/
   interpolate from its actual (n, DEG) observations. This is the number
   reported to the user, together with a CI/SE from the model.
8. **Adaptive alpha (starting approach: empirical FDR calibration).** Add a
   permutation arm to the downsampling loop: at each n, also run DE on
   label-shuffled data to get a null distribution of DEG counts (or p-values)
   at that sample size. Choose alpha(n) such that the expected false-positive
   count (or expected FDR) implied by that null stays constant across n,
   rather than using a fixed alpha = 0.05 everywhere. Reuses the same
   downsampling infrastructure (just adds `permute = TRUE` runs). Output is an
   alpha(n) curve per cell type/method, plus DEG counts computed under it.

## Alternatives considered for adaptive alpha (not yet implemented)

Kept here so we don't forget them when revisiting:
- **q-value (Storey) instead of BH**, which adapts to the estimated null
  proportion at each n automatically — may reduce the need to hand-tune alpha
  at all; cheaper than permutation but less directly interpretable as "held
  quantity."
- **Power-matching**: choose alpha(n) so that power to detect a fixed effect
  size is held constant across n (via analytical power formulas per DE method,
  or empirically via spike-in/simulated-effect benchmarks), rather than
  holding false-positive rate constant. Closer to the actual goal (comparable
  *sensitivity*) but needs either method-specific power formulas or its own
  simulation study — bigger lift, good stretch goal once empirical-FDR
  calibration is validated.

## Known caveat: extrapolation beyond a cell type's observed range

If `reference_n` (e.g. the dataset median) exceeds a given cell type's own
full cell count, `stratified`/per-cell-type models must extrapolate the
fitted curve outside its training range to compute burden — splines in
particular can diverge sharply there. `calculate_burden()` flags this via
`burden$extrapolated` and emits a warning; treat `burden`/`burden_se` for
flagged cell types with real skepticism (large `burden_se` is usually the
tell), and consider forcing `model_type = "linear"` (safer extrapolation
behavior) or excluding very small cell types from `reference_n` entirely for
those rows. This is a fundamental limitation of predicting "what if this
cell type had more cells than it actually has," not a bug to be fully
engineered away — a cell type with 30 cells simply has no data informing
what would happen at 150.

## Other open questions / future work

- Exact reference-n rule: median cells-per-treatment across cell types is the
  default; consider also supporting "min across cell types" (conservative,
  never extrapolates) as a preset.
- Whether burden should be reported as a raw predicted DEG count, or
  normalized (e.g., by number of genes tested in that cell type, since gene
  detection also scales with cell count).
- ~~Python port~~: done, see `python/`. Mirrors the R API and defaults
  (quasi-Poisson identity-link count model, pooled-by-default, empirical FDR
  calibration). Two intentional implementation differences, documented in
  `python/README.md`: splines are `patsy`'s unpenalized natural cubic
  regression-spline basis (`cr()`) rather than `mgcv`'s penalized smooths,
  and nested-model comparisons use `statsmodels` GLM with `scale="X2"` and a
  quasi-F test rather than `mgcv::gam`/`anova.gam`. edgeR is called via
  `rpy2` (optional extra); DESeq2 uses the `pydeseq2` port rather than
  calling R's DESeq2. Direct AnnData support (rather than requiring the
  caller to pass raw matrices) is still open future work.
- Concordance-based exclusion vs weighting: currently planned as a weighting
  scheme by default; revisit whether a hard minimum-n cutoff should be the
  default instead (simpler to explain, more conservative).
