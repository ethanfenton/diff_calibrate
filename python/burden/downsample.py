"""Downsampling + DE loop, per cell type and across a dataset."""

from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from .de_wrappers import resolve_de_fun


def downsample_grid(n_max: int, n_min: int = 20, n_points: int = 8) -> np.ndarray:
    """Log-spaced grid of target per-group cell counts from n_min to n_max.

    De-duplicated and sorted; ``n_max`` is always included so the full-data
    point anchors the model.
    """
    if n_max < n_min:
        return np.array([n_max])
    grid = np.round(np.exp(np.linspace(np.log(n_min), np.log(n_max), n_points)))
    grid = np.unique(np.minimum(grid, n_max)).astype(int)
    if grid[-1] != n_max:
        grid = np.append(grid, n_max)
    return grid


def downsample_counts(
    counts: np.ndarray,
    group,
    n_per_group: int,
    rng: Optional[np.random.Generator] = None,
) -> dict:
    """Sample n_per_group cells without replacement from each group level.

    If a group has fewer than n_per_group cells, all of its cells are kept
    (with a warning) rather than erroring.
    """
    import warnings

    rng = rng or np.random.default_rng()
    group = pd.Categorical(group)
    keep_idx = []
    for lv in group.categories:
        idx = np.where(np.asarray(group) == lv)[0]
        n = min(n_per_group, len(idx))
        if n < n_per_group:
            warnings.warn(
                f"Group '{lv}' has only {len(idx)} cells (< requested "
                f"{n_per_group}); using all of them."
            )
        keep_idx.append(rng.choice(idx, size=n, replace=False))
    keep_idx = np.concatenate(keep_idx)
    sub_group = pd.Categorical(np.asarray(group)[keep_idx]).remove_unused_categories()
    return {"counts": counts[:, keep_idx], "group": sub_group}


def n_degs(de_res: pd.DataFrame, alpha: float = 0.05) -> int:
    """Count genes with padj < alpha (NAs excluded)."""
    return int((de_res["padj"] < alpha).sum())


def run_downsampling(
    counts: np.ndarray,
    group,
    de_fun: Union[str, callable],
    grid: Optional[Sequence[int]] = None,
    grid_args: Optional[dict] = None,
    k: int = 5,
    k_permute: int = 5,
    alpha: float = 0.05,
    seed: Optional[int] = 1,
    **kwargs,
) -> dict:
    """Downsampling + DE loop for a single cell type.

    For a grid of per-group cell counts, draws ``k`` replicate downsamples at
    each grid point and runs ``de_fun`` on each, plus (optionally)
    ``k_permute`` label-permuted replicates per grid point for empirical null
    calibration (see :func:`burden.adaptive_alpha.calibrate_alpha`). Always
    also runs DE once on the full data as the reference/ground-truth result.

    Returns
    -------
    dict with keys:
        full : pd.DataFrame, DE result at full cell count.
        downsampled : pd.DataFrame with columns n, rep, n_degs, de_res (list
            of per-replicate DE result DataFrames).
        permuted : same shape as downsampled, from label-shuffled runs, or
            None if k_permute == 0.
    """
    de_fun = resolve_de_fun(de_fun)
    group = pd.Categorical(group)
    n_max = int(pd.Series(group).value_counts().min())

    if grid is None:
        grid = downsample_grid(n_max=n_max, **(grid_args or {}))

    full_res = de_fun(counts, group, **kwargs)

    master_rng = np.random.default_rng(seed)

    rows = []
    for n in grid:
        for r in range(k):
            rng = np.random.default_rng(master_rng.integers(1 << 32))
            ds = downsample_counts(counts, group, int(n), rng=rng)
            res = de_fun(ds["counts"], ds["group"], **kwargs)
            rows.append({"n": int(n), "rep": r + 1, "n_degs": n_degs(res, alpha), "de_res": res})
    downsampled = pd.DataFrame(rows)

    permuted = None
    if k_permute > 0:
        prows = []
        for n in grid:
            for r in range(k_permute):
                rng = np.random.default_rng(master_rng.integers(1 << 32))
                ds = downsample_counts(counts, group, int(n), rng=rng)
                shuffled_group = pd.Categorical(rng.permutation(np.asarray(ds["group"])))
                res = de_fun(ds["counts"], shuffled_group, **kwargs)
                prows.append(
                    {"n": int(n), "rep": r + 1, "n_degs": n_degs(res, alpha), "de_res": res}
                )
        permuted = pd.DataFrame(prows)

    return {"full": full_res, "downsampled": downsampled, "permuted": permuted}


def run_downsampling_dataset(data: dict, de_fun: Union[str, callable], **kwargs) -> dict:
    """Run the downsampling loop across all cell types in a dataset.

    Parameters
    ----------
    data : dict
        Mapping cell_type -> {"counts": genes x cells matrix, "group": ...}.
    de_fun : DE backend, see resolve_de_fun.
    **kwargs : passed to run_downsampling identically for every cell type.

    Returns
    -------
    dict, cell_type -> run_downsampling() result.
    """
    if not data:
        raise ValueError("`data` must be a non-empty dict of cell_type -> {counts, group}.")
    return {
        ct: run_downsampling(d["counts"], d["group"], de_fun, **kwargs)
        for ct, d in data.items()
    }
