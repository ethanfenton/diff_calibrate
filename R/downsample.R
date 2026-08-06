#' Build a downsampling grid of per-group cell counts
#'
#' Log-spaced grid of target cell counts from `n_min` up to `n_max`
#' (inclusive), de-duplicated and sorted. `n_max` is always included so the
#' full-data point anchors the model.
#'
#' @param n_max Full available cell count for the group being downsampled.
#' @param n_min Smallest cell count to include. Default 20.
#' @param n_points Number of grid points requested (actual count may be lower
#'   after de-duplication/rounding, or capped if `n_max < n_min`).
#' @return Integer vector, ascending, always ending in `n_max`.
#' @export
downsample_grid <- function(n_max, n_min = 20, n_points = 8) {
  if (n_max < n_min) return(n_max)
  grid <- round(exp(seq(log(n_min), log(n_max), length.out = n_points)))
  grid <- sort(unique(pmin(grid, n_max)))
  if (grid[length(grid)] != n_max) grid <- c(grid, n_max)
  grid
}

#' Downsample a counts matrix to a target number of cells per group
#'
#' Samples `n_per_group` cells (columns) without replacement from each level
#' of `group` independently. If a group has fewer than `n_per_group` cells,
#' all of its cells are kept (with a warning) rather than erroring, since
#' upstream grid construction should generally prevent this.
#'
#' @param counts Genes x cells matrix.
#' @param group Length-`ncol(counts)` grouping vector.
#' @param n_per_group Target cell count per group level.
#' @param seed Optional integer seed for reproducible sampling.
#' @return `list(counts = <subsetted matrix>, group = <subsetted group>)`.
#' @export
downsample_counts <- function(counts, group, n_per_group, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  group <- as.factor(group)
  keep_idx <- unlist(lapply(levels(group), function(lv) {
    idx <- which(group == lv)
    n <- min(n_per_group, length(idx))
    if (n < n_per_group) {
      warning("Group '", lv, "' has only ", length(idx),
               " cells (< requested ", n_per_group, "); using all of them.",
               call. = FALSE)
    }
    sample(idx, n, replace = FALSE)
  }))
  list(counts = counts[, keep_idx, drop = FALSE], group = droplevels(group[keep_idx]))
}

#' Count significant genes in a DE result table
#'
#' @param de_res `data.frame` with a `padj` column (as returned by a
#'   [de_backends] function).
#' @param alpha Significance threshold.
#' @return Integer count of genes with `padj < alpha` (NAs excluded).
#' @export
n_degs <- function(de_res, alpha = 0.05) {
  sum(de_res$padj < alpha, na.rm = TRUE)
}

#' Run the downsampling + DE loop for a single cell type
#'
#' For a grid of per-group cell counts, draws `k` replicate downsamples at
#' each grid point and runs `de_fun` on each, plus (optionally) `k_permute`
#' label-permuted replicates per grid point for empirical null calibration
#' (see [calibrate_alpha()]). Always also runs DE once on the full data as
#' the reference/ground-truth result.
#'
#' @param counts Genes x cells matrix for this cell type.
#' @param group Two-level grouping vector (e.g. treatment) for this cell type.
#' @param de_fun DE backend: string name or function, see [resolve_de_fun()].
#' @param grid Integer vector of per-group cell counts; if `NULL`, built via
#'   [downsample_grid()] using `grid_args`.
#' @param grid_args List of extra arguments to [downsample_grid()] (`n_min`,
#'   `n_points`) when `grid` is not supplied directly.
#' @param k Number of downsampling replicates per grid point.
#' @param k_permute Number of label-permuted replicates per grid point, for
#'   empirical alpha calibration. Set to `0` to skip (saves compute if
#'   adaptive alpha is not needed).
#' @param alpha Nominal significance threshold used to compute `n_degs` in
#'   the returned table (in addition to the full `de_res` table, which lets
#'   downstream code recompute at other thresholds without rerunning DE).
#' @param seed Optional integer seed; downstream draws are seeded
#'   deterministically off of this (`seed + rep_index`) for reproducibility.
#' @param ... Passed to `de_fun`.
#' @return A list with:
#'   \item{full}{`data.frame` DE result at full cell count (the reference).}
#'   \item{downsampled}{`data.frame` with one row per (n, replicate):
#'     `n`, `rep`, `n_degs`, plus a list-column `de_res` holding the full DE
#'     table for that replicate (used by [concordance()]).}
#'   \item{permuted}{Same shape as `downsampled`, from label-shuffled runs,
#'     or `NULL` if `k_permute == 0`.}
#' @export
run_downsampling <- function(counts, group, de_fun,
                              grid = NULL, grid_args = list(),
                              k = 5, k_permute = 5, alpha = 0.05,
                              seed = 1, ...) {
  de_fun <- resolve_de_fun(de_fun)
  group <- as.factor(group)
  n_max <- min(table(group))

  if (is.null(grid)) {
    grid <- do.call(downsample_grid, c(list(n_max = n_max), grid_args))
  }

  full_res <- de_fun(counts, group, ...)

  rep_idx <- 0L
  downsampled <- do.call(rbind, lapply(grid, function(n) {
    do.call(rbind, lapply(seq_len(k), function(r) {
      rep_idx <<- rep_idx + 1L
      ds <- downsample_counts(counts, group, n, seed = seed + rep_idx)
      res <- de_fun(ds$counts, ds$group, ...)
      data.frame(n = n, rep = r, n_degs = n_degs(res, alpha),
                 de_res = I(list(res)))
    }))
  }))

  permuted <- NULL
  if (k_permute > 0) {
    perm_idx <- 0L
    permuted <- do.call(rbind, lapply(grid, function(n) {
      do.call(rbind, lapply(seq_len(k_permute), function(r) {
        perm_idx <<- perm_idx + 1L
        ds <- downsample_counts(counts, group, n, seed = seed + 100000L + perm_idx)
        shuffled_group <- sample(ds$group)
        res <- de_fun(ds$counts, shuffled_group, ...)
        data.frame(n = n, rep = r, n_degs = n_degs(res, alpha),
                   de_res = I(list(res)))
      }))
    }))
  }

  list(full = full_res, downsampled = downsampled, permuted = permuted)
}

#' Run the downsampling loop across all cell types in a dataset
#'
#' @param data Named list, one element per cell type, each an unnamed
#'   `list(counts = <genes x cells matrix>, group = <two-level vector>)`.
#'   Names of `data` are used as cell type labels throughout the package.
#' @param de_fun DE backend, see [resolve_de_fun()].
#' @param ... Passed to [run_downsampling()] (applied identically to every
#'   cell type; pass per-cell-type grids via `grid` if needed by calling
#'   [run_downsampling()] directly instead).
#' @return Named list (by cell type) of [run_downsampling()] results.
#' @export
run_downsampling_dataset <- function(data, de_fun, ...) {
  if (is.null(names(data)) || any(names(data) == "")) {
    stop("`data` must be a named list, with names giving cell type labels.",
         call. = FALSE)
  }
  stats::setNames(lapply(names(data), function(ct) {
    d <- data[[ct]]
    run_downsampling(d$counts, d$group, de_fun, ...)
  }), names(data))
}
