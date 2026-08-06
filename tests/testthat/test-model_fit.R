test_that("test_nonlinearity picks linear for a truly linear relationship", {
  set.seed(20)
  n <- rep(seq(20, 200, by = 20), each = 20)
  degs <- rpois(length(n), lambda = 0.5 * n)
  out <- test_nonlinearity(n, degs)
  expect_equal(out$decision, "linear")
  expect_true(out$p_value > 0.05)
})

test_that("test_nonlinearity picks spline for a saturating relationship", {
  set.seed(21)
  n <- rep(seq(10, 500, length.out = 12), each = 6)
  mu <- 100 * (1 - exp(-n / 50))
  degs <- rpois(length(n), lambda = mu)
  out <- test_nonlinearity(n, degs)
  expect_equal(out$decision, "spline")
})

test_that("fit_deg_model respects forced model_type", {
  set.seed(22)
  n <- rep(seq(20, 200, by = 20), each = 5)
  degs <- rpois(length(n), lambda = 0.5 * n)
  fit_lin <- fit_deg_model(n, degs, model_type = "linear")
  fit_spl <- fit_deg_model(n, degs, model_type = "spline")
  expect_equal(fit_lin$decision, "linear")
  expect_equal(fit_spl$decision, "spline")
})

test_that("predict_degs returns sensible predictions with se", {
  set.seed(23)
  n <- rep(seq(20, 200, by = 20), each = 5)
  degs <- rpois(length(n), lambda = 0.5 * n)
  fit <- fit_deg_model(n, degs, model_type = "linear")
  pred <- predict_degs(fit, 100, se = TRUE)
  expect_true(pred$fit > 0)
  expect_true(pred$se > 0)
})

test_that("test_pooling recommends pooled when cell types share a curve", {
  set.seed(5)
  n <- rep(seq(20, 200, by = 20), each = 10)
  degs1 <- rpois(length(n), lambda = 0.5 * n)
  degs2 <- rpois(length(n), lambda = 0.5 * n)
  out <- test_pooling(c(n, n), c(degs1, degs2),
                       c(rep("a", length(n)), rep("b", length(n))))
  expect_equal(out$recommendation, "pooled")
})

test_that("test_pooling recommends stratified when curves clearly differ", {
  set.seed(25)
  n <- rep(seq(20, 200, by = 20), each = 5)
  degs1 <- rpois(length(n), lambda = 0.1 * n)
  degs2 <- rpois(length(n), lambda = 2 * n)
  out <- test_pooling(c(n, n), c(degs1, degs2),
                       c(rep("a", length(n)), rep("b", length(n))))
  expect_true(out$recommendation %in% c("shared_shape", "stratified"))
})
