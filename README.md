# burden

Sample-size-adjusted differential expression (DE) "burden" for single-cell
data. Cell types with more cells tend to show more significant DE genes for
reasons of statistical power alone, which biases naive DEG counts toward
larger cell types. `burden` downsamples each cell type across a range of
cell counts, refits DE at each, and models DEG count as a function of cell
count — so you can ask "how many DEGs would this cell type show if it had as
many cells as a typical cell type in this dataset," rather than comparing
raw DEG counts across cell types of very different size.

See [`docs/DESIGN.md`](docs/DESIGN.md) for the full methodology, rationale,
and known caveats (in particular: extrapolation risk when a cell type is
much smaller than the reference cell count).

## Install

```r
# from a local checkout
install.packages("devtools")
devtools::install("path/to/burden")
```

## Usage

```r
library(burden)

# data: named list, one element per cell type, each
#   list(counts = <genes x cells matrix>, group = <two-level treatment vector>)
data <- list(
  neurons     = list(counts = neuron_counts,     group = neuron_group),
  astrocytes  = list(counts = astrocyte_counts,  group = astrocyte_group),
  microglia   = list(counts = microglia_counts,  group = microglia_group)
)

res <- calculate_burden(
  data,
  de_fun = "wilcoxon",       # or "edger", "deseq2", or your own function
  model_type = "auto",       # "auto" | "spline" | "linear"
  pooling = "auto",          # "auto" | "pooled" | "shared_shape" | "stratified"
  adaptive_alpha = TRUE,     # empirical FDR-calibrated alpha(n) per cell type
  alpha = 0.05
)

res$burden          # cell_type, n_full, reference_n, burden, burden_se, extrapolated
res$pooling_test     # LRTs behind the pooled/shared_shape/stratified recommendation
res$alpha_curves      # per-cell-type alpha(n) from empirical FDR calibration
res$concordance        # per-cell-type concordance of downsampled vs full DE calls
```

`de_fun` can be any function following the contract in `?de_backends`:
`function(counts, group, ...)` returning `data.frame(gene, pval, padj,
logFC)`, so any DE method (not just the three built-ins) can be plugged in.

## Development

```r
devtools::document()   # regenerate NAMESPACE / man pages after editing R/
devtools::test()       # run the testthat suite
devtools::check()      # full R CMD check
```
