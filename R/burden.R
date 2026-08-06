#' Calculate sample-size-adjusted DE burden across cell types
#'
#' End-to-end pipeline: downsamples each cell type across a grid of cell
#' counts, runs DE at each, computes concordance against the full-data DE
#' result, optionally calibrates a per-n adaptive alpha from label-permuted
#' null runs, tests whether a spline is needed and whether cell types share
#' a common DEGs~n curve, fits the resulting model(s), and predicts the
#' expected DEG count ("burden") for every cell type at a common reference
#' cell count. See `docs/DESIGN.md` for the full methodology and rationale.
#'
#' @param data Named list, one element per cell type, each
#'   `list(counts = <genes x cells matrix>, group = <two-level vector>)`.
#'   Names give cell type labels.
#' @param de_fun DE backend: `"wilcoxon"`, `"edger"`, `"deseq2"`, or a
#'   user function following the [de_backends] contract.
#' @param model_type One of `"auto"`, `"spline"`, `"linear"` — see
#'   [fit_deg_model()].
#' @param pooling One of `"pooled"` (default: a single DEGs~n curve fit
#'   across all cell types, on the assumption that cell count affects
#'   detected DEG count similarly regardless of cell type), `"auto"` (run
#'   [test_pooling()] and use its LRT-based recommendation instead), or
#'   force `"shared_shape"`/`"stratified"` (per-cell-type curves). `"auto"`
#'   and `"stratified"` trade the pooled default's stability for
#'   flexibility: each cell type's own model is fit to far fewer
#'   downsampling replicates than the pooled fit, so it is noisier, and
#'   `"auto"`'s per-dataset LRT can itself false-positive into
#'   `shared_shape`/`stratified` on a single noisy comparison (observed
#'   during development, see `docs/DESIGN.md`) — which is why pooling is
#'   not auto-selected by default. [test_pooling()]'s result is always
#'   computed and returned as `pooling_test` regardless of this argument, so
#'   you can inspect what `"auto"` *would* have chosen without switching to
#'   it.
#' @param adaptive_alpha If `TRUE` (default), calibrate alpha(n) per cell
#'   type via [calibrate_alpha()] (empirical FDR calibration against
#'   label-permuted null runs) and use it in place of `alpha` when counting
#'   DEGs for the model. Requires `k_permute > 0`.
#' @param alpha Nominal significance threshold. Used directly if
#'   `adaptive_alpha = FALSE`; used as `target_fpr` for [calibrate_alpha()]
#'   otherwise.
#' @param reference_n Cell count at which burden is predicted. Default
#'   `NULL` uses the median, across cell types, of each cell type's full
#'   per-group cell count (`min(table(group))`) — i.e. "a typical cell
#'   type's cell count in this dataset." Pass a number to override.
#' @param min_concordance Downsampled replicates with a `concordance_metric`
#'   value below this are excluded from model fitting entirely (`NULL`
#'   disables hard exclusion; see `concordance_metric` and `weight_floor`
#'   for soft down-weighting instead/in addition).
#' @param concordance_metric Which [concordance()] measure to use for
#'   `min_concordance` and for weighting. Default `"cat_overlap"` (alpha-free).
#' @param weight_floor Minimum weight assigned to any retained replicate
#'   (prevents near-zero concordance from making a point contribute
#'   ~nothing while not being hard-excluded). Default 0.1.
#' @param sig_level Significance threshold for the nonlinearity and pooling
#'   LRTs.
#' @param grid_args,k,k_permute,seed Passed to [run_downsampling()].
#' @param ... Passed to `de_fun`.
#' @return A `burden_result` object (list) with:
#'   \item{burden}{`data.frame(cell_type, n_full, reference_n, burden,
#'     burden_se, n_degs_observed)` — the main output.}
#'   \item{downsampling}{Per-cell-type [run_downsampling()] results.}
#'   \item{concordance}{Per-cell-type [concordance_curve()] results.}
#'   \item{alpha_curves}{Per-cell-type [calibrate_alpha()] results, or
#'     `NULL` if `adaptive_alpha = FALSE`.}
#'   \item{pooling_test}{[test_pooling()] output.}
#'   \item{model}{The fitted model object(s) actually used for prediction.}
#' @export
calculate_burden <- function(data, de_fun,
                              model_type = c("auto", "spline", "linear"),
                              pooling = c("pooled", "auto", "shared_shape", "stratified"),
                              adaptive_alpha = TRUE,
                              alpha = 0.05,
                              reference_n = NULL,
                              min_concordance = NULL,
                              concordance_metric = "cat_overlap",
                              weight_floor = 0.1,
                              sig_level = 0.05,
                              grid_args = list(), k = 5, k_permute = 5, seed = 1,
                              ...) {
  model_type <- match.arg(model_type)
  pooling <- match.arg(pooling)
  if (adaptive_alpha && k_permute <= 0) {
    stop("adaptive_alpha = TRUE requires k_permute > 0.", call. = FALSE)
  }

  cell_types <- names(data)

  downsampling <- run_downsampling_dataset(
    data, de_fun, grid_args = grid_args, k = k, k_permute = k_permute,
    alpha = alpha, seed = seed, ...
  )

  concordance <- stats::setNames(
    lapply(downsampling, concordance_curve, alpha = alpha), cell_types
  )

  alpha_curves <- NULL
  if (adaptive_alpha) {
    alpha_curves <- stats::setNames(
      lapply(downsampling, calibrate_alpha, target_fpr = alpha), cell_types
    )
    downsampling <- stats::setNames(lapply(cell_types, function(ct) {
      ds <- downsampling[[ct]]
      ds$downsampled <- apply_adaptive_alpha(ds, alpha_curves[[ct]])
      ds
    }), cell_types)
  }

  combined <- do.call(rbind, lapply(cell_types, function(ct) {
    ds <- downsampling[[ct]]$downsampled
    conc <- concordance[[ct]]
    w <- pmax(conc[[concordance_metric]], weight_floor, na.rm = FALSE)
    w[is.na(w)] <- weight_floor
    data.frame(cell_type = ct, n = ds$n, degs = ds$n_degs, weight = w)
  }))

  if (!is.null(min_concordance)) {
    conc_all <- do.call(rbind, lapply(cell_types, function(ct) {
      cbind(cell_type = ct, concordance[[ct]])
    }))
    keep <- conc_all[[concordance_metric]] >= min_concordance
    combined <- combined[keep, , drop = FALSE]
  }

  pooling_test <- test_pooling(combined$n, combined$degs, combined$cell_type,
                                weights = combined$weight, sig_level = sig_level)
  pooling_used <- if (pooling == "auto") pooling_test$recommendation else pooling

  if (is.null(reference_n)) {
    n_full_per_ct <- vapply(downsampling, function(d) max(d$downsampled$n), numeric(1))
    reference_n <- stats::median(n_full_per_ct)
  }

  model <- NULL
  burden_rows <- NULL

  if (pooling_used == "pooled") {
    model <- fit_deg_model(combined$n, combined$degs, weights = combined$weight,
                            model_type = model_type, sig_level = sig_level)
    pred <- predict_degs(model, rep(reference_n, length(cell_types)), se = TRUE)
    burden_rows <- data.frame(cell_type = cell_types, burden = pred$fit, burden_se = pred$se)
  } else if (pooling_used == "shared_shape") {
    model <- fit_deg_model(combined$n, combined$degs, weights = combined$weight,
                            cell_type = combined$cell_type, model_type = model_type,
                            sig_level = sig_level)
    pred <- predict_degs(model, rep(reference_n, length(cell_types)),
                          cell_type_new = cell_types, se = TRUE)
    burden_rows <- data.frame(cell_type = cell_types, burden = pred$fit, burden_se = pred$se)
  } else {
    model <- stats::setNames(lapply(cell_types, function(ct) {
      sub <- combined[combined$cell_type == ct, , drop = FALSE]
      fit_deg_model(sub$n, sub$degs, weights = sub$weight, model_type = model_type,
                    sig_level = sig_level)
    }), cell_types)
    burden_rows <- do.call(rbind, lapply(cell_types, function(ct) {
      pred <- predict_degs(model[[ct]], reference_n, se = TRUE)
      data.frame(cell_type = ct, burden = pred$fit, burden_se = pred$se)
    }))
  }

  n_full <- vapply(downsampling, function(d) max(d$downsampled$n), numeric(1))
  n_degs_observed <- vapply(downsampling, function(d) n_degs(d$full, alpha), numeric(1))

  burden <- data.frame(
    cell_type = cell_types,
    n_full = n_full[cell_types],
    reference_n = reference_n,
    n_degs_observed = n_degs_observed[cell_types],
    stringsAsFactors = FALSE
  )
  burden <- merge(burden, burden_rows, by = "cell_type", sort = FALSE)
  burden <- burden[match(cell_types, burden$cell_type), , drop = FALSE]
  rownames(burden) <- NULL
  burden$extrapolated <- reference_n > burden$n_full
  if (any(burden$extrapolated)) {
    warning(
      "reference_n (", reference_n, ") exceeds the full cell count for: ",
      paste(burden$cell_type[burden$extrapolated], collapse = ", "),
      ". Burden for these cell types is extrapolated beyond observed data ",
      "(see the `extrapolated` column and treat burden_se accordingly) ",
      "-- this is most severe for stratified/per-cell-type models with a ",
      "spline, which can diverge sharply outside the training range; ",
      "consider model_type = \"linear\" or excluding these cell types.",
      call. = FALSE
    )
  }

  structure(list(
    burden = burden,
    downsampling = downsampling,
    concordance = concordance,
    alpha_curves = alpha_curves,
    pooling_test = pooling_test,
    pooling_used = pooling_used,
    model = model
  ), class = "burden_result")
}

#' @export
print.burden_result <- function(x, ...) {
  cat("<burden_result>", nrow(x$burden), "cell types, pooling:", x$pooling_used, "\n")
  print(x$burden, row.names = FALSE)
  if (any(x$burden$extrapolated)) {
    cat("Note: burden is extrapolated beyond observed cell counts for: ",
        paste(x$burden$cell_type[x$burden$extrapolated], collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
