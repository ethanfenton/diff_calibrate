test_that("calibrate_alpha returns one row per grid point and errors without permuted arm", {
  d <- sim_cell_type(40, n_genes = 40, seed = 8)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 3, seed = 12)
  ac <- calibrate_alpha(out, target_fpr = 0.05)
  expect_equal(nrow(ac), length(unique(out$downsampled$n)))
  expect_true(all(ac$alpha >= 0 & ac$alpha <= 1))

  out_noperm <- run_downsampling(d$counts, d$group, "wilcoxon",
                                  grid_args = list(n_min = 10, n_points = 3),
                                  k = 2, k_permute = 0, seed = 13)
  expect_error(calibrate_alpha(out_noperm), "no permuted arm")
})

test_that("cap_at_target enforces alpha(n) <= target_fpr", {
  d <- sim_cell_type(40, n_genes = 40, seed = 9)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 3, seed = 14)
  ac <- calibrate_alpha(out, target_fpr = 0.05, cap_at_target = TRUE)
  expect_true(all(ac$alpha <= 0.05))
})

test_that("predict_alpha interpolates and holds endpoints constant", {
  ac <- data.frame(n = c(10, 50, 100), alpha = c(0.01, 0.05, 0.1))
  expect_equal(predict_alpha(ac, 50), 0.05)
  expect_equal(predict_alpha(ac, 5), 0.01)
  expect_equal(predict_alpha(ac, 200), 0.1)
})

test_that("apply_adaptive_alpha recomputes n_degs without rerunning DE", {
  d <- sim_cell_type(40, n_genes = 40, seed = 10)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 3, seed = 15)
  ac <- calibrate_alpha(out)
  ds2 <- apply_adaptive_alpha(out, ac)
  expect_true("alpha_used" %in% names(ds2))
  expect_equal(nrow(ds2), nrow(out$downsampled))
})
