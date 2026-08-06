#' Empirical FDR-calibrated significance threshold as a function of n
#'
#' At small cell counts, DE p-values can be miscalibrated (typically
#' anti-conservative), so a fixed `padj < 0.05` cutoff does not carry the
#' same false-positive rate at every sample size. This calibrates a
#' per-`n` threshold `alpha(n)` from the label-permuted arm of
#' [run_downsampling()] (`k_permute > 0`): permuted labels give a DE run
#' with no true group difference, so its `padj` distribution is an empirical
#' null. `alpha(n)` is set to the empirical quantile of that null `padj`
#' distribution at `target_fpr`, i.e. the threshold that would call a
#' `target_fpr` fraction of genes "significant" under the null at that `n`.
#' Where the method is anti-conservative at small n, this yields
#' `alpha(n) < target_fpr` (a stricter cutoff, correcting for excess false
#' positives); where it is conservative, `alpha(n) > target_fpr` is allowed
#' (recovering power) unless `cap_at_target = TRUE`.
#'
#' See `docs/DESIGN.md` for alternative calibration strategies (q-values,
#' power-matching) not yet implemented.
#'
#' @param downsampling_result A single cell type's result from
#'   [run_downsampling()], run with `k_permute > 0`.
#' @param target_fpr Target false-positive rate to hold constant across n.
#'   Default 0.05 (matches the conventional fixed-alpha default, so
#'   `alpha(n)` should be read as "what alpha achieves the same effective
#'   FPR as a nominal 0.05 would at the full sample size").
#' @param cap_at_target If `TRUE`, never allow `alpha(n) > target_fpr`
#'   (i.e. only correct for anti-conservative bias, never spend extra power
#'   on conservative bias). Default `FALSE`.
#' @return `data.frame(n, alpha, n_null_tests)`, one row per grid point,
#'   `n_null_tests` giving the pooled number of null (gene x permutation
#'   replicate) tests the estimate at that `n` is based on (useful for
#'   judging how noisy `alpha(n)` is at sparse grid points).
#' @export
calibrate_alpha <- function(downsampling_result, target_fpr = 0.05,
                             cap_at_target = FALSE) {
  perm <- downsampling_result$permuted
  if (is.null(perm)) {
    stop("downsampling_result has no permuted arm; rerun run_downsampling() ",
         "with k_permute > 0.", call. = FALSE)
  }

  ns <- sort(unique(perm$n))
  out <- do.call(rbind, lapply(ns, function(n) {
    rows <- perm[perm$n == n, , drop = FALSE]
    pooled_padj <- unlist(lapply(rows$de_res, function(r) r$padj))
    pooled_padj <- pooled_padj[!is.na(pooled_padj)]
    a <- if (length(pooled_padj) == 0) {
      target_fpr
    } else {
      as.numeric(stats::quantile(pooled_padj, probs = target_fpr, type = 7))
    }
    if (cap_at_target) a <- min(a, target_fpr)
    data.frame(n = n, alpha = a, n_null_tests = length(pooled_padj))
  }))
  out
}

#' Interpolate a calibrated alpha(n) curve to arbitrary cell counts
#'
#' Linear interpolation on the [calibrate_alpha()] grid, with endpoints held
#' constant outside the observed range.
#'
#' @param alpha_curve Output of [calibrate_alpha()].
#' @param n Vector of cell counts to evaluate at.
#' @return Numeric vector of `alpha` values, same length as `n`.
#' @export
predict_alpha <- function(alpha_curve, n) {
  stats::approx(alpha_curve$n, alpha_curve$alpha, xout = n, rule = 2)$y
}

#' Recompute DEG counts using a calibrated alpha(n) instead of a fixed alpha
#'
#' Reuses the full DE result tables already stored in a [run_downsampling()]
#' result (`downsampled$de_res`), so no DE refitting is needed — just
#' recounts `padj < alpha(n)` with a per-row, n-specific alpha from
#' [calibrate_alpha()]/[predict_alpha()].
#'
#' @param downsampling_result A single cell type's result from
#'   [run_downsampling()].
#' @param alpha_curve Output of [calibrate_alpha()] for the same cell type.
#' @return The `downsampled` `data.frame` from `downsampling_result`, with
#'   `n_degs` replaced by counts at the calibrated alpha, and a new column
#'   `alpha_used`.
#' @export
apply_adaptive_alpha <- function(downsampling_result, alpha_curve) {
  ds <- downsampling_result$downsampled
  alpha_used <- predict_alpha(alpha_curve, ds$n)
  ds$n_degs <- vapply(seq_len(nrow(ds)), function(i) {
    n_degs(ds$de_res[[i]], alpha = alpha_used[i])
  }, numeric(1))
  ds$alpha_used <- alpha_used
  ds
}
