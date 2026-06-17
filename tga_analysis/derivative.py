"""DTG (dW/dT) via Savitzky-Golay derivative.

Differentiating noisy data amplifies noise; SG sidesteps that by fitting
a local polynomial and analytically differentiating the polynomial,
which is much better behaved than a finite difference on smoothed data.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter


def dtg_savgol(
    w: np.ndarray,
    dT: float,
    window_pts: int,
    poly_order: int = 3,
) -> np.ndarray:
    """First derivative dW/dT (% per °C) on the uniform T grid.

    Uses the same SG window/order chosen for smoothing the mass signal
    so the smoothed mass and DTG are mutually consistent.
    """
    return savgol_filter(w, window_pts, poly_order, deriv=1, delta=dT)
