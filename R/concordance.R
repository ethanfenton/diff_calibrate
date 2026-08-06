#' Concordance between a downsampled DE result and the full-data DE result
#'
#' Computes three complementary concordance measures between a DE result
#' table from a downsampled run and the reference (full-data) DE result for
#' the same cell type, restricted to genes present in both tables:
#' \describe{
#'   \item{jaccard}{Jaccard index of DEG sets at `alpha` (set overlap /
#'     union). `NaN` if both sets are empty.}
#'   \item{f1}{F1 score of the downsampled DEG set against the full DEG set
#'     treated as truth. `NaN` if both sets are empty.}
#'   \item{rank_cor}{Spearman correlation of `-log10(pval)` between the two
#'     runs — alpha-free, sensitive to overall rank agreement.}
#'   \item{cat_overlap}{Concordance-at-the-top: fraction overlap of the
#'     top-`cat_k` genes by `abs(logFC)` in each run (alpha-free; `cat_k`
#'     defaults to the smaller of 50 or a quarter of genes tested).}
#' }
#'
#' @param de_res `data.frame` DE result for the downsampled run.
#' @param full_res `data.frame` DE result for the full-data run (reference).
#' @param alpha Significance threshold for the set-based measures.
#' @param cat_k Number of top genes to use for concordance-at-the-top.
#'   Defaults to `min(50, floor(n_genes / 4))`.
#' @return A one-row `data.frame(jaccard, f1, rank_cor, cat_overlap)`.
#' @export
concordance <- function(de_res, full_res, alpha = 0.05, cat_k = NULL) {
  common <- intersect(de_res$gene, full_res$gene)
  ds <- de_res[match(common, de_res$gene), , drop = FALSE]
  fl <- full_res[match(common, full_res$gene), , drop = FALSE]

  ds_sig <- ds$gene[!is.na(ds$padj) & ds$padj < alpha]
  fl_sig <- fl$gene[!is.na(fl$padj) & fl$padj < alpha]

  inter <- length(intersect(ds_sig, fl_sig))
  union_n <- length(union(ds_sig, fl_sig))
  jaccard <- if (union_n == 0) NaN else inter / union_n

  tp <- inter
  fp <- length(setdiff(ds_sig, fl_sig))
  fn <- length(setdiff(fl_sig, ds_sig))
  f1 <- if ((2 * tp + fp + fn) == 0) NaN else (2 * tp) / (2 * tp + fp + fn)

  rank_cor <- if (length(common) >= 3) {
    suppressWarnings(stats::cor(-log10(pmax(ds$pval, .Machine$double.xmin)),
                                 -log10(pmax(fl$pval, .Machine$double.xmin)),
                                 method = "spearman", use = "pairwise.complete.obs"))
  } else {
    NA_real_
  }

  if (is.null(cat_k)) cat_k <- max(1, min(50, floor(length(common) / 4)))
  top_ds <- ds$gene[order(-abs(ds$logFC))][seq_len(min(cat_k, nrow(ds)))]
  top_fl <- fl$gene[order(-abs(fl$logFC))][seq_len(min(cat_k, nrow(fl)))]
  cat_overlap <- length(intersect(top_ds, top_fl)) / cat_k

  data.frame(jaccard = jaccard, f1 = f1, rank_cor = rank_cor,
             cat_overlap = cat_overlap)
}

#' Concordance curve across all downsampled replicates for a cell type
#'
#' Applies [concordance()] to every replicate in a [run_downsampling()]
#' result, against that cell type's full-data result, and returns one row
#' per (n, replicate) with `n`, `rep`, and the concordance measures. Intended
#' to be inspected/plotted as a function of `n`, and optionally used to
#' derive weights or an exclusion cutoff before fitting the DEG~n model (see
#' `min_concordance`/`concordance_weight` in [fit_deg_model()]).
#'
#' @param downsampling_result A single cell type's result from
#'   [run_downsampling()] (i.e. `run_downsampling_dataset(...)[[cell_type]]`).
#' @param alpha,cat_k Passed to [concordance()].
#' @return `data.frame(n, rep, jaccard, f1, rank_cor, cat_overlap)`.
#' @export
concordance_curve <- function(downsampling_result, alpha = 0.05, cat_k = NULL) {
  ds <- downsampling_result$downsampled
  full_res <- downsampling_result$full
  conc <- do.call(rbind, lapply(ds$de_res, concordance, full_res = full_res,
                                 alpha = alpha, cat_k = cat_k))
  cbind(n = ds$n, rep = ds$rep, conc)
}
