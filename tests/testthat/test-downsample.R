test_that("downsample_grid always includes n_max and respects n_min", {
  g <- downsample_grid(200, n_min = 20, n_points = 5)
  expect_true(max(g) == 200)
  expect_true(min(g) >= 20)
  expect_true(!is.unsorted(g))

  g_small <- downsample_grid(10, n_min = 20)
  expect_equal(g_small, 10)
})

test_that("downsample_counts samples the requested number per group", {
  d <- sim_cell_type(50, seed = 3)
  ds <- downsample_counts(d$counts, d$group, n_per_group = 10, seed = 42)
  expect_equal(as.integer(table(ds$group)), c(10L, 10L))
  expect_equal(ncol(ds$counts), 20)
})

test_that("downsample_counts warns and keeps all cells if n_per_group too large", {
  d <- sim_cell_type(5, seed = 4)
  expect_warning(ds <- downsample_counts(d$counts, d$group, n_per_group = 100),
                  "has only")
  expect_equal(as.integer(table(ds$group)), c(5L, 5L))
})

test_that("run_downsampling produces the documented structure", {
  d <- sim_cell_type(30, n_genes = 40, seed = 5)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 2, seed = 10)
  expect_true(all(c("full", "downsampled", "permuted") %in% names(out)))
  expect_true(all(c("n", "rep", "n_degs", "de_res") %in% names(out$downsampled)))
  expect_true(all(out$downsampled$n <= 30))
  expect_equal(nrow(out$permuted), nrow(out$downsampled))
})

test_that("n_degs counts padj below alpha", {
  res <- data.frame(padj = c(0.01, 0.2, NA, 0.049, 0.05))
  expect_equal(n_degs(res, alpha = 0.05), 2)
})
