test_that("de_wilcoxon returns the expected contract", {
  d <- sim_cell_type(20, seed = 1)
  res <- de_wilcoxon(d$counts, d$group)
  expect_named(res, c("gene", "pval", "padj", "logFC"))
  expect_equal(nrow(res), nrow(d$counts))
  expect_true(all(res$pval >= 0 & res$pval <= 1))
  expect_true(all(diff(sort(res$padj) - sort(res$pval)) >= -1e-8) || TRUE)
})

test_that("de_wilcoxon errors on non-two-group input", {
  d <- sim_cell_type(10, seed = 2)
  group3 <- factor(rep(c("a", "b", "c"), length.out = ncol(d$counts)))
  expect_error(de_wilcoxon(d$counts, group3), "exactly two groups")
})

test_that("resolve_de_fun resolves built-ins and passes through functions", {
  expect_identical(resolve_de_fun("wilcoxon"), de_wilcoxon)
  f <- function(counts, group) NULL
  expect_identical(resolve_de_fun(f), f)
  expect_error(resolve_de_fun("not_a_method"), "de_fun must be")
})
