#' Nested-model likelihood ratio test between two `mgcv::gam` fits
#'
#' Thin wrapper around `anova.gam(..., test = "Chisq")` used for both the
#' spline-need test ([test_nonlinearity()]) and the pooling-need test
#' ([test_pooling()]) — both reduce to "is a more flexible nested model
#' justified over a simpler one," so they share this helper.
#'
#' @param model_null,model_alt `mgcv::gam` objects, `model_null` nested
#'   within `model_alt`.
#' @param min_df Minimum extra (effective) degrees of freedom `model_alt`
#'   must use beyond `model_null` for the test to be considered meaningful.
#'   When a penalized smooth collapses back onto (or very near) its nested
#'   null — e.g. a spline that shrinks to an effectively straight line — the
#'   nominal F/chi-square test can become numerically degenerate (near-zero
#'   denominator df) and report spuriously tiny p-values despite the two
#'   models being practically identical. Below this threshold the test is
#'   short-circuited to `p_value = 1` rather than trusting that degenerate
#'   statistic. Default 0.4.
#' @return A one-row `data.frame(p_value, chisq, df)`.
#' @keywords internal
lrt_test <- function(model_null, model_alt, min_df = 0.4) {
  a <- stats::anova(model_null, model_alt, test = "F")
  # anova.gam names differ across mgcv versions; grab the last row's stats
  # generically by column-name pattern.
  pcol <- grep("^Pr\\(", colnames(a), value = TRUE)
  chisqcol <- grep("Chisq|Deviance", colnames(a), value = TRUE)
  dfcol <- grep("^Df$", colnames(a), value = TRUE)
  df <- if (length(dfcol)) a[[dfcol]][2] else NA_real_
  p_value <- if (length(pcol)) a[[pcol[1]]][2] else NA_real_
  if (!is.na(df) && df < min_df) p_value <- 1
  data.frame(
    p_value = p_value,
    chisq = if (length(chisqcol)) a[[chisqcol[length(chisqcol)]]][2] else NA_real_,
    df = df
  )
}

.spline_k <- function(n_unique, k_max = 8) {
  max(3, min(k_max, n_unique - 1))
}

#' Default count-response family for DEGs ~ n models
#'
#' Quasi-Poisson on the identity link: identity because "linear" is meant in
#' the natural sense of "DEG count grows linearly with cell count" (the
#' quantity a user reasons about and that [predict_degs()] extrapolates),
#' not "log(DEG count) grows linearly with cell count" as a log-link GLM
#' would test; quasi-Poisson (one dispersion parameter, comparable across
#' nested fits) rather than negative binomial because comparing two
#' negative-binomial fits that each estimate their own `theta` via nested
#' LRT is anti-conservative (the dispersion mismatch itself inflates the
#' apparent deviance reduction) — quasi-likelihood's F-test is the standard
#' fix. See `docs/DESIGN.md`.
#' @keywords internal
.default_family <- function() stats::quasipoisson(link = "identity")

#' Fit an `mgcv::gam` suppressing the transient "NaNs produced" warning
#'
#' The identity-link quasi-Poisson deviance can briefly evaluate at a
#' negative mean during IRLS iteration, which is harmless (the optimizer
#' recovers) but noisy; suppressed here rather than left for every caller of
#' [test_nonlinearity()]/[test_pooling()] to see repeatedly.
#' @keywords internal
.gam_quiet <- function(formula, family, data, weights) {
  withCallingHandlers(
    mgcv::gam(formula = formula, family = family, data = data, weights = weights),
    warning = function(w) {
      if (grepl("NaNs produced", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

#' Test whether a spline term is needed for DEGs ~ n
#'
#' Fits a linear negative-binomial GAM (`degs ~ n`) and a spline
#' negative-binomial GAM (`degs ~ s(n)`), both via `mgcv::gam`, and compares
#' them with a likelihood ratio test. A significant result (by default
#' p < 0.05) indicates the relationship between cell count and DEG count is
#' not well captured by a straight line on the model's link scale, and a
#' spline should be used.
#'
#' @param n Numeric vector of cell counts.
#' @param degs Integer vector of DEG counts, same length as `n`.
#' @param weights Optional prior weights (e.g. concordance-derived), same
#'   length as `n`.
#' @param cell_type Optional factor/character vector, same length as `n`.
#'   When supplied, both the linear and spline models include `+ cell_type`
#'   as an additive term (shared curve shape, per-cell-type intercept) — use
#'   this to run the nonlinearity test on a `shared_shape` pooled fit rather
#'   than a single cell type.
#' @param family Response family for both `mgcv::gam` fits. Default
#'   [`.default_family()`], quasi-Poisson on the identity link — see that
#'   function's docs for why.
#' @param sig_level Significance threshold for the nonlinearity call.
#' @return `list(p_value, decision, model_linear, model_spline)`, `decision`
#'   one of `"spline"` / `"linear"`.
#' @export
test_nonlinearity <- function(n, degs, weights = NULL, cell_type = NULL,
                               family = .default_family(), sig_level = 0.05) {
  df <- data.frame(n = n, degs = degs)
  w <- weights %||% rep(1, length(n))
  k <- .spline_k(length(unique(n)))

  if (is.null(cell_type)) {
    f_linear <- degs ~ n
    f_spline <- degs ~ s(n, k = k)
  } else {
    df$cell_type <- as.factor(cell_type)
    f_linear <- degs ~ n + cell_type
    f_spline <- stats::as.formula(paste0("degs ~ s(n, k = ", k, ") + cell_type"))
  }

  model_linear <- .gam_quiet(f_linear, family = family, data = df, weights = w)
  model_spline <- .gam_quiet(f_spline, family = family, data = df, weights = w)

  lrt <- lrt_test(model_linear, model_spline)
  decision <- if (!is.na(lrt$p_value) && lrt$p_value < sig_level) "spline" else "linear"

  list(p_value = lrt$p_value, decision = decision,
       model_linear = model_linear, model_spline = model_spline)
}

#' Fit the DEGs ~ n model for a single cell type
#'
#' @param n Numeric vector of cell counts (one per downsampling replicate).
#' @param degs Integer vector of DEG counts, same length as `n`.
#' @param weights Optional prior weights, e.g. derived from
#'   [concordance_curve()] (see `min_concordance`/`weight_from`).
#' @param cell_type Optional cell-type vector, see [test_nonlinearity()];
#'   used for `shared_shape` pooled fits.
#' @param model_type One of `"auto"` (default; run [test_nonlinearity()] and
#'   pick), `"spline"` (force spline), `"linear"` (force linear).
#' @param family Passed to [test_nonlinearity()].
#' @param sig_level Significance threshold used when `model_type = "auto"`.
#' @return A `burden_model` object: `list(model, model_type, decision,
#'   nonlinearity_test)`, where `nonlinearity_test` is the full
#'   [test_nonlinearity()] output (always computed, even when `model_type`
#'   is forced, so the decision that *would* have been made is visible).
#' @export
fit_deg_model <- function(n, degs, weights = NULL, cell_type = NULL,
                           model_type = c("auto", "spline", "linear"),
                           family = .default_family(), sig_level = 0.05) {
  model_type <- match.arg(model_type)
  nl_test <- test_nonlinearity(n, degs, weights = weights, cell_type = cell_type,
                                family = family, sig_level = sig_level)

  decision <- switch(model_type,
                      auto = nl_test$decision,
                      spline = "spline",
                      linear = "linear")

  model <- if (decision == "spline") nl_test$model_spline else nl_test$model_linear

  structure(list(model = model, model_type = model_type, decision = decision,
                  nonlinearity_test = nl_test),
            class = "burden_model")
}

#' Predict expected DEG count at a given cell count
#'
#' @param fit A `burden_model` from [fit_deg_model()] (or a fitted
#'   `mgcv::gam` directly).
#' @param n_new Cell count(s) to predict at.
#' @param cell_type_new Optional cell type label(s) to predict at (required,
#'   recycled to `length(n_new)`, if `fit`'s model includes a `cell_type`
#'   term, e.g. a `shared_shape` pooled fit).
#' @param se If `TRUE`, also return the standard error of the prediction
#'   (on the response scale, via the delta method).
#' @return Numeric vector of predicted DEG counts, or if `se = TRUE`, a
#'   `data.frame(n, fit, se)`.
#' @export
predict_degs <- function(fit, n_new, cell_type_new = NULL, se = FALSE) {
  model <- if (inherits(fit, "burden_model")) fit$model else fit
  newdata <- data.frame(n = n_new)
  if (!is.null(cell_type_new)) {
    newdata$cell_type <- factor(cell_type_new, levels = levels(model$model$cell_type))
  }
  if (!se) {
    return(as.numeric(stats::predict(model, newdata = newdata, type = "response")))
  }
  pr <- stats::predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  linkinv <- model$family$linkinv
  mu <- linkinv(pr$fit)
  # delta method: se on response scale ~ |d(linkinv)/d(eta)| * se(eta)
  d <- (linkinv(pr$fit + 1e-6) - linkinv(pr$fit - 1e-6)) / 2e-6
  data.frame(n = n_new, fit = mu, se = abs(d) * pr$se.fit)
}

#' Test whether a pooled (across cell type) model is justified
#'
#' Compares three nested `mgcv::gam` fits of DEGs ~ n across all cell types
#' in a combined dataset:
#' \itemize{
#'   \item `pooled`: `degs ~ s(n)` — single shared curve, no cell-type effect.
#'   \item `shared_shape`: `degs ~ s(n) + cell_type` — shared curve shape,
#'     per-cell-type intercept (random-intercept-like via a fixed factor).
#'   \item `stratified`: `degs ~ s(n, by = cell_type) + cell_type` — fully
#'     separate curve per cell type.
#' }
#' Two LRTs are run: `pooled` vs `shared_shape` (is a per-cell-type offset
#' needed at all?) and `shared_shape` vs `stratified` (is a per-cell-type
#' curve *shape* needed, beyond an offset?). Use the result to decide which
#' of the three to actually use for prediction (see `recommendation`).
#'
#' @param n,degs,cell_type Equal-length vectors pooling downsampling results
#'   across cell types.
#' @param weights Optional prior weights.
#' @param family Response family for all three `mgcv::gam` fits. Default
#'   [`.default_family()`].
#' @param sig_level Significance threshold for both LRTs.
#' @return `list(recommendation, offset_test, shape_test, model_pooled,
#'   model_shared_shape, model_stratified)`. `recommendation` is one of
#'   `"pooled"`, `"shared_shape"`, `"stratified"`.
#' @export
test_pooling <- function(n, degs, cell_type, weights = NULL,
                          family = .default_family(), sig_level = 0.05) {
  df <- data.frame(n = n, degs = degs, cell_type = as.factor(cell_type))
  w <- weights %||% rep(1, length(n))
  k <- .spline_k(length(unique(n)))

  model_pooled <- .gam_quiet(degs ~ s(n, k = k), family = family, data = df, weights = w)
  model_shared_shape <- .gam_quiet(degs ~ s(n, k = k) + cell_type, family = family,
                                    data = df, weights = w)
  model_stratified <- .gam_quiet(degs ~ s(n, k = k, by = cell_type) + cell_type,
                                  family = family, data = df, weights = w)

  offset_test <- lrt_test(model_pooled, model_shared_shape)
  shape_test <- lrt_test(model_shared_shape, model_stratified)

  recommendation <- if (is.na(offset_test$p_value) || offset_test$p_value >= sig_level) {
    "pooled"
  } else if (is.na(shape_test$p_value) || shape_test$p_value >= sig_level) {
    "shared_shape"
  } else {
    "stratified"
  }

  list(recommendation = recommendation, offset_test = offset_test, shape_test = shape_test,
       model_pooled = model_pooled, model_shared_shape = model_shared_shape,
       model_stratified = model_stratified)
}
