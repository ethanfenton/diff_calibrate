test_that("calculate_burden runs end-to-end and returns one row per cell type", {
  data <- list(
    small = sim_cell_type(25, n_genes = 60, seed = 30),
    large = sim_cell_type(150, n_genes = 60, seed = 31)
  )
  res <- calculate_burden(
    data, de_fun = "wilcoxon",
    grid_args = list(n_min = 12, n_points = 3),
    k = 2, k_permute = 2, seed = 100
  )
  expect_s3_class(res, "burden_result")
  expect_equal(nrow(res$burden), 2)
  expect_true(all(c("cell_type", "n_full", "reference_n", "burden", "burden_se") %in%
                    names(res$burden)))
  expect_true(all(res$burden$burden >= 0))
})

test_that("calculate_burden works with adaptive_alpha = FALSE and forced pooling/model_type", {
  data <- list(
    a = sim_cell_type(20, n_genes = 50, seed = 32),
    b = sim_cell_type(80, n_genes = 50, seed = 33)
  )
  res <- calculate_burden(
    data, de_fun = "wilcoxon", model_type = "linear", pooling = "stratified",
    adaptive_alpha = FALSE,
    grid_args = list(n_min = 10, n_points = 3),
    k = 2, k_permute = 0, seed = 101
  )
  expect_equal(res$pooling_used, "stratified")
  expect_true(is.list(res$model) && !inherits(res$model, "burden_model"))
  expect_null(res$alpha_curves)
})

test_that("calculate_burden errors if adaptive_alpha requested without permuted reps", {
  data <- list(a = sim_cell_type(20, n_genes = 30, seed = 34))
  expect_error(
    calculate_burden(data, de_fun = "wilcoxon", adaptive_alpha = TRUE, k_permute = 0),
    "k_permute > 0"
  )
})

test_that("reference_n defaults to the median full cell count across cell types", {
  data <- list(
    a = sim_cell_type(20, n_genes = 30, seed = 35),
    b = sim_cell_type(60, n_genes = 30, seed = 36),
    c = sim_cell_type(100, n_genes = 30, seed = 37)
  )
  res <- calculate_burden(
    data, de_fun = "wilcoxon", adaptive_alpha = FALSE,
    grid_args = list(n_min = 10, n_points = 3),
    k = 2, k_permute = 0, seed = 102
  )
  expect_equal(unique(res$burden$reference_n), 60)
})
