"""TGA signal denoising and mass-loss event quantification."""
from .io import load_tga, read_batch_list, TGARun
from .smooth import smooth_savgol, SmoothingResult
from .derivative import dtg_savgol
from .events import detect_events, Event
from .quantify import quantify, MassLoss
from .plot import plot_run

__all__ = [
    "load_tga", "read_batch_list", "TGARun",
    "smooth_savgol", "SmoothingResult",
    "dtg_savgol",
    "detect_events", "Event",
    "quantify", "MassLoss",
    "plot_run",
]
