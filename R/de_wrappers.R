#' Built-in and pluggable differential expression backends
#'
#' A DE backend is any function with signature `function(counts, group, ...)`
#' returning a `data.frame` with columns `gene`, `pval`, `padj`, `logFC`.
#' `counts` is a genes x cells numeric matrix (raw counts, not normalized,
#' unless the backend expects otherwise), `group` is a factor/character vector
#' of length `ncol(counts)` giving the two-level grouping to test. Any
#' function following this contract can be passed as `de_fun` to
#' [calculate_burden()] or [run_downsampling()] directly; the wrappers below
#' are convenience implementations for common methods.
#'
#' @name de_backends
NULL

.check_two_groups <- function(group) {
  group <- droplevels(as.factor(group))
  lv <- levels(group)
  if (length(lv) != 2) {
    stop("DE backends require exactly two groups, got ", length(lv), ": ",
         paste(lv, collapse = ", "), call. = FALSE)
  }
  group
}

#' Wilcoxon rank-sum DE
#'
#' Gene-by-gene Wilcoxon rank-sum test on (by default) log1p-CPM normalized
#' counts, with BH-adjusted p-values and log2 fold change of group means
#' (pseudocount 1) as effect size. Dependency-free fallback DE method; fast
#' enough to be practical inside the downsampling loop.
#'
#' @param counts Genes x cells numeric/count matrix.
#' @param group Length-`ncol(counts)` vector with exactly two levels.
#' @param normalize If `TRUE` (default), CPM-normalize and log1p-transform
#'   `counts` before testing. Set `FALSE` if `counts` is already normalized.
#' @return `data.frame(gene, pval, padj, logFC)`.
#' @export
de_wilcoxon <- function(counts, group, normalize = TRUE) {
  group <- .check_two_groups(group)
  if (normalize) {
    lib_size <- colSums(counts)
    lib_size[lib_size == 0] <- 1
    cpm <- t(t(counts) / lib_size) * 1e6
    mat <- log1p(cpm)
  } else {
    mat <- counts
  }

  lv <- levels(group)
  idx1 <- which(group == lv[1])
  idx2 <- which(group == lv[2])

  pval <- vapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, idx1]
    y <- mat[i, idx2]
    if (stats::sd(x) == 0 && stats::sd(y) == 0) return(1)
    suppressWarnings(stats::wilcox.test(x, y)$p.value)
  }, numeric(1))

  logFC <- log2((rowMeans(mat[, idx2, drop = FALSE]) + 1) /
                (rowMeans(mat[, idx1, drop = FALSE]) + 1))

  data.frame(
    gene = rownames(counts) %||% paste0("gene", seq_len(nrow(counts))),
    pval = pval,
    padj = stats::p.adjust(pval, method = "BH"),
    logFC = logFC,
    stringsAsFactors = FALSE
  )
}

#' edgeR quasi-likelihood F-test DE
#'
#' Thin wrapper around `edgeR`'s QLF pipeline (`estimateDisp` +
#' `glmQLFit` + `glmQLFTest`) for two-group comparisons. Requires the
#' `edgeR` package (Suggests).
#'
#' @inheritParams de_wilcoxon
#' @param ... Passed to `edgeR::glmQLFTest`.
#' @return `data.frame(gene, pval, padj, logFC)`.
#' @export
de_edger <- function(counts, group, ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for de_edger(). Install it via BiocManager.",
         call. = FALSE)
  }
  group <- .check_two_groups(group)
  keep_genes <- rowSums(counts) > 0
  y <- edgeR::DGEList(counts = counts[keep_genes, , drop = FALSE], group = group)
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  res <- edgeR::glmQLFTest(fit, coef = 2, ...)
  tab <- res$table

  out <- data.frame(
    gene = rownames(counts),
    pval = 1,
    padj = 1,
    logFC = 0,
    stringsAsFactors = FALSE
  )
  m <- match(rownames(tab), out$gene)
  out$pval[m] <- tab$PValue
  out$padj[m] <- stats::p.adjust(out$pval, method = "BH")
  out$logFC[m] <- tab$logFC
  out
}

#' DESeq2 Wald test DE
#'
#' Thin wrapper around `DESeq2` for two-group comparisons (`DESeq()` with
#' default Wald test, `lfcShrink`-free `results()` extraction). Requires the
#' `DESeq2` package (Suggests). Substantially slower than
#' [de_wilcoxon()]/[de_edger()]; expect longer runtimes inside the
#' downsampling loop.
#'
#' @inheritParams de_wilcoxon
#' @param ... Passed to `DESeq2::DESeq`.
#' @return `data.frame(gene, pval, padj, logFC)`.
#' @export
de_deseq2 <- function(counts, group, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for de_deseq2(). Install it via BiocManager.",
         call. = FALSE)
  }
  group <- .check_two_groups(group)
  keep_genes <- rowSums(counts) > 0
  cts <- round(counts[keep_genes, , drop = FALSE])
  coldata <- data.frame(group = group)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                         design = ~group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)
  res <- as.data.frame(DESeq2::results(dds))

  out <- data.frame(
    gene = rownames(counts),
    pval = 1,
    padj = 1,
    logFC = 0,
    stringsAsFactors = FALSE
  )
  m <- match(rownames(res), out$gene)
  out$pval[m] <- ifelse(is.na(res$pvalue), 1, res$pvalue)
  out$padj[m] <- stats::p.adjust(out$pval, method = "BH")
  out$logFC[m] <- ifelse(is.na(res$log2FoldChange), 0, res$log2FoldChange)
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.de_backends <- list(
  wilcoxon = de_wilcoxon,
  edger = de_edger,
  deseq2 = de_deseq2
)

#' Resolve a DE backend from a name or function
#'
#' @param de_fun Either a string naming a built-in backend ("wilcoxon",
#'   "edger", "deseq2") or a user-supplied function following the
#'   [de_backends] contract.
#' @return A function.
#' @keywords internal
resolve_de_fun <- function(de_fun) {
  if (is.function(de_fun)) return(de_fun)
  if (is.character(de_fun) && length(de_fun) == 1 && de_fun %in% names(.de_backends)) {
    return(.de_backends[[de_fun]])
  }
  stop("de_fun must be a function, or one of: ",
       paste(names(.de_backends), collapse = ", "), call. = FALSE)
}
