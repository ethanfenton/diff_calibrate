import numpy as np
import pytest


def sim_cell_type(n_genes=200, n_per_group=60, n_de=20, effect=2.5, seed=0):
    """Simulate a genes x cells count matrix with a known set of DE genes."""
    rng = np.random.default_rng(seed)
    n_cells = n_per_group * 2
    base = rng.gamma(2, 2, size=n_genes)
    counts = rng.poisson(base[:, None], size=(n_genes, n_cells)).astype(float)
    group = np.array(["A"] * n_per_group + ["B"] * n_per_group)
    de_idx = rng.choice(n_genes, n_de, replace=False)
    b_idx = np.where(group == "B")[0]
    counts[np.ix_(de_idx, b_idx)] = rng.poisson(
        (base[de_idx] * effect)[:, None], size=(n_de, len(b_idx))
    )
    return {"counts": counts, "group": group, "de_genes": de_idx}


@pytest.fixture
def sim():
    return sim_cell_type
