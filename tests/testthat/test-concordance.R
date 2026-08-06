test_that("concordance is 1 when a result is compared to itself", {
  d <- sim_cell_type(40, n_genes = 60, seed = 6)
  res <- de_wilcoxon(d$counts, d$group)
  conc <- concordance(res, res, alpha = 0.05)
  expect_equal(conc$jaccard, 1)
  expect_equal(conc$f1, 1)
  expect_equal(conc$cat_overlap, 1)
  expect_equal(conc$rank_cor, 1, tolerance = 1e-8)
})

test_that("concordance handles disjoint DEG sets", {
  full <- data.frame(gene = paste0("g", 1:10), pval = c(rep(0.001, 5), rep(0.9, 5)),
                      padj = c(rep(0.01, 5), rep(0.9, 5)), logFC = c(rep(2, 5), rep(0, 5)))
  ds <- data.frame(gene = paste0("g", 1:10), pval = c(rep(0.9, 5), rep(0.001, 5)),
                    padj = c(rep(0.9, 5), rep(0.01, 5)), logFC = c(rep(0, 5), rep(2, 5)))
  conc <- concordance(ds, full, alpha = 0.05)
  expect_equal(conc$jaccard, 0)
  expect_equal(conc$f1, 0)
})

test_that("concordance_curve returns one row per downsampled replicate", {
  d <- sim_cell_type(30, n_genes = 40, seed = 7)
  out <- run_downsampling(d$counts, d$group, "wilcoxon",
                           grid_args = list(n_min = 10, n_points = 3),
                           k = 2, k_permute = 0, seed = 11)
  cc <- concordance_curve(out)
  expect_equal(nrow(cc), nrow(out$downsampled))
  expect_true(all(cc$cat_overlap >= 0 & cc$cat_overlap <= 1))
})
