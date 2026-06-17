"""Optimal smoothing of TGA mass signal.

Primary smoother: Savitzky-Golay (SG) with polynomial order fixed at 3.
SG is the standard in thermal analysis because the same local polynomial
fit gives both the smoothed signal and its derivative without an extra
differentiation step that would re-amplify noise.

Window selection is data-driven via residual whiteness. The idea: when
the smoother is correctly tuned, residuals (raw - smoothed) should look
like white noise. If the window is too wide, the smoother removes part
of the signal, and the residuals retain low-frequency structure — they
show positive lag-1 autocorrelation. If the window is too narrow, the
residuals are essentially noise and the autocorrelation is near zero.

We sweep candidate windows in °C-equivalent units (so the criterion is
consistent across runs with different sampling densities), and choose
the *largest* window whose residual |lag-1 autocorrelation| stays below
a small threshold. That is the most aggressive smoothing that has not
yet started to remove real signal. The Durbin-Watson statistic on the
residuals is reported alongside for cross-checking.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import savgol_filter


@dataclass
class SmoothingResult:
    w_smooth: np.ndarray
    window_pts: int
    window_C: float
    poly_order: int
    rho1: float
    durbin_watson: float
    residual_std: float
    candidates: list[dict] = field(default_factory=list)


def _lag1_autocorr(r: np.ndarray) -> float:
    r = r - r.mean()
    denom = float(np.sum(r * r))
    if denom == 0.0:
        return 0.0
    return float(np.sum(r[:-1] * r[1:]) / denom)


def _durbin_watson(r: np.ndarray) -> float:
    denom = float(np.sum(r * r))
    if denom == 0.0:
        return 2.0
    d = np.diff(r)
    return float(np.sum(d * d) / denom)


def smooth_savgol(
    w: np.ndarray,
    dT: float,
    poly_order: int = 3,
    window_C_range: tuple[float, float] = (0.3, 5.0),
    n_candidates: int = 20,
    rho1_threshold: float = 0.05,
) -> SmoothingResult:
    """Smooth w (mass vs. T on a uniform grid) with an auto-selected SG window.

    Parameters
    ----------
    w : array
        Weight signal on a uniform T grid (% of initial mass).
    dT : float
        Grid spacing in °C, matching the spacing used in io.load_tga.
    poly_order : int
        SG polynomial order. Order 3 is the thermal-analysis default and
        works well for the first derivative (DTG) used downstream.
    window_C_range : (float, float)
        Lower and upper bounds for the SG window length in °C. The
        defaults span the range used in practice for 10 K/min TGA runs.
    n_candidates : int
        Number of candidate windows sampled across the range.
    rho1_threshold : float
        Maximum acceptable magnitude of the residual lag-1 autocorrelation.
        Below this, residuals are considered close to white.

    Returns
    -------
    SmoothingResult
        Includes the chosen window plus per-candidate diagnostics so the
        user can verify the choice was sensible.
    """
    n = len(w)
    if n < 7:
        raise ValueError(f"Signal too short for SG smoothing: n={n}")

    w_min = max(poly_order + 2, int(round(window_C_range[0] / dT)))
    w_max = min(n - 1, int(round(window_C_range[1] / dT)))
    w_min |= 1  # SG needs odd window
    w_max |= 1
    if w_max <= w_min:
        raise ValueError(
            f"Window search range too narrow: [{w_min}, {w_max}] pts"
        )

    cand_pts = np.unique(
        np.clip(
            (np.round(np.linspace(w_min, w_max, n_candidates)).astype(int) | 1),
            w_min,
            w_max,
        )
    )

    candidates: list[dict] = []
    for wp in cand_pts:
        w_s = savgol_filter(w, int(wp), poly_order)
        r = w - w_s
        candidates.append(
            {
                "window_pts": int(wp),
                "window_C": float(wp * dT),
                "rho1": _lag1_autocorr(r),
                "durbin_watson": _durbin_watson(r),
                "residual_std": float(r.std()),
            }
        )

    under_threshold = [c for c in candidates if abs(c["rho1"]) <= rho1_threshold]
    if under_threshold:
        chosen = max(under_threshold, key=lambda c: c["window_pts"])
    else:
        chosen = min(candidates, key=lambda c: abs(c["rho1"]))

    w_smooth = savgol_filter(w, chosen["window_pts"], poly_order)
    return SmoothingResult(
        w_smooth=w_smooth,
        window_pts=chosen["window_pts"],
        window_C=chosen["window_C"],
        poly_order=poly_order,
        rho1=chosen["rho1"],
        durbin_watson=chosen["durbin_watson"],
        residual_std=chosen["residual_std"],
        candidates=candidates,
    )
