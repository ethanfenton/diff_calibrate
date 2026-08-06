sim_cell_type <- function(n_per_group, n_genes = 100, n_de = 10, effect = 1.6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_total <- n_per_group * 2
  group <- factor(rep(c("ctrl", "trt"), each = n_per_group))
  base_mu <- stats::rgamma(n_genes, shape = 2, rate = 0.5) + 0.5
  mat <- matrix(0, n_genes, n_total)
  de_idx <- sample(n_genes, n_de)
  for (g in seq_len(n_genes)) {
    mu_ctrl <- base_mu[g]
    mu_trt <- if (g %in% de_idx) base_mu[g] * effect else base_mu[g]
    mat[g, group == "ctrl"] <- stats::rpois(n_per_group, mu_ctrl)
    mat[g, group == "trt"] <- stats::rpois(n_per_group, mu_trt)
  }
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  list(counts = mat, group = group)
}
