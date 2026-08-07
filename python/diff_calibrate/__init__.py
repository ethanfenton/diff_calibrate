"""diff_calibrate: sample-size-adjusted differential expression burden for single-cell data."""

from .adaptive_alpha import apply_adaptive_alpha, calibrate_alpha, predict_alpha
from .burden import BurdenResult, calculate_burden
from .concordance import concordance, concordance_curve
from .de_wrappers import de_deseq2, de_edger, de_wilcoxon, resolve_de_fun
from .downsample import (
    downsample_counts,
    downsample_grid,
    n_degs,
    run_downsampling,
    run_downsampling_dataset,
)
from .model_fit import BurdenModel, fit_deg_model, predict_degs, quasi_f_test, test_nonlinearity, test_pooling

__version__ = "0.1.0"

__all__ = [
    "calculate_burden",
    "BurdenResult",
    "BurdenModel",
    "de_wilcoxon",
    "de_edger",
    "de_deseq2",
    "resolve_de_fun",
    "downsample_grid",
    "downsample_counts",
    "n_degs",
    "run_downsampling",
    "run_downsampling_dataset",
    "concordance",
    "concordance_curve",
    "calibrate_alpha",
    "predict_alpha",
    "apply_adaptive_alpha",
    "test_nonlinearity",
    "fit_deg_model",
    "predict_degs",
    "test_pooling",
    "quasi_f_test",
]
