"""TGA signal denoising, event detection, mass-loss quantification, and
Frazer-Suzuki peak deconvolution for overlapping events."""
from .io import load_tga, read_batch_list, TGARun
from .smooth import smooth_savgol, SmoothingResult
from .derivative import dtg_savgol
from .events import detect_events, Event
from .quantify import quantify, MassLoss
from .plot import plot_run
from .deconvolve import (
    DeconvolvedPeak,
    baseline_dtg,
    deconvolve,
    fit_region_baseline,
    frazer_suzuki,
    frazer_suzuki_sum,
)

__all__ = [
    "load_tga", "read_batch_list", "TGARun",
    "smooth_savgol", "SmoothingResult",
    "dtg_savgol",
    "detect_events", "Event",
    "quantify", "MassLoss",
    "plot_run",
    "DeconvolvedPeak", "baseline_dtg", "deconvolve", "fit_region_baseline",
    "frazer_suzuki", "frazer_suzuki_sum",
]
