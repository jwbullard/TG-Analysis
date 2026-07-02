"""Frazer-Suzuki peak-sum deconvolution for overlapping DTG events.

The Frazer-Suzuki (FS) function is the standard asymmetric peak form
used in thermal-analysis kinetics. Its asymmetry parameter b smoothly
interpolates between Gaussian (b → 0) and asymmetric-tail shapes
characteristic of real thermal-decomposition events:

    y(T) = A · exp[ -(ln 2 / b²) · ln²(1 + 2b(T - T₀)/W) ]

Parameters
----------
A   — peak amplitude (DTG units, % per °C)
T₀  — peak temperature (°C)
W   — full width at half maximum (°C); exact in the Gaussian limit b → 0
b   — asymmetry (b > 0: tail to the right; b < 0: tail to the left).
      The function is zero outside the domain 1 + 2b(T - T₀)/W > 0;
      with |b| ≲ 0.95 this domain comfortably exceeds typical event widths.

For a region containing N overlapping peaks, fit
    y_total(T) = Σ_i FS(T; A_i, T₀_i, W_i, b_i)
to the baseline-corrected positive-going DTG (i.e. -dW/dT with the
chosen drift baseline subtracted). Each fitted peak's area is its mass
loss; the sum of all peak areas equals the integral of the baseline-
corrected DTG over the region.

This module pairs with quantify.py: the same three baseline modes
(constant/linear/alpha) are applied region-wide, the FS sum is fit
separately for each, and the deconvolved per-peak mass losses can be
compared across modes. For an isolated peak in a region with quiet
baselines on both sides, the per-peak result is invariant to mode (to
the precision of the fit); strong divergence flags either an ill-chosen
region or a deconvolution that has not properly separated the peaks.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

_LN2 = float(np.log(2.0))


@dataclass
class DeconvolvedPeak:
    A: float                     # amplitude (% per °C)
    T0: float                    # peak temperature (°C)
    W: float                     # FWHM (°C)
    b: float                     # asymmetry
    area: float                  # numerical integral over the fit T grid (% mass loss)


def frazer_suzuki(T: np.ndarray, A: float, T0: float, W: float, b: float) -> np.ndarray:
    """Single Frazer-Suzuki peak evaluated on T. Zero outside the valid domain."""
    if abs(b) < 1e-8:
        return A * np.exp(-_LN2 * (2.0 * (T - T0) / W) ** 2)
    arg = 1.0 + 2.0 * b * (T - T0) / W
    out = np.zeros_like(T, dtype=float)
    valid = arg > 0.0
    out[valid] = A * np.exp(-(_LN2 / b**2) * np.log(arg[valid]) ** 2)
    return out


def frazer_suzuki_sum(T: np.ndarray, *params: float) -> np.ndarray:
    if len(params) % 4 != 0:
        raise ValueError(f"FS sum expects 4 params per peak, got {len(params)}")
    y = np.zeros_like(T, dtype=float)
    for i in range(0, len(params), 4):
        y = y + frazer_suzuki(T, *params[i : i + 4])
    return y


def _initial_peak_positions(T: np.ndarray, y: np.ndarray, n: int) -> list[float]:
    """Seed peak T₀ guesses from scipy.find_peaks; pad evenly if too few found."""
    mag = np.where(y > 0, y, 0.0)
    peaks, props = find_peaks(mag, prominence=0.0)
    if len(peaks) >= n:
        order = np.argsort(props["prominences"])[::-1]
        chosen = sorted(peaks[order[:n]].tolist())
        return [float(T[i]) for i in chosen]
    positions = [float(T[i]) for i in sorted(peaks.tolist())]
    T_min, T_max = float(T[0]), float(T[-1])
    while len(positions) < n:
        positions.append(T_min + (T_max - T_min) * (len(positions) + 1) / (n + 1))
    return sorted(positions[:n])


def deconvolve(
    T: np.ndarray,
    y: np.ndarray,
    n_peaks: int,
    W_init: float | None = None,
    b_init: float = 0.0,
) -> list[DeconvolvedPeak]:
    """Fit a sum of n_peaks Frazer-Suzuki functions to y(T).

    y should be the positive-going baseline-corrected DTG signal over a
    contiguous T-range (i.e. -(dW/dT - baseline) with all values nominally
    ≥ 0 inside the peaks). T is in °C on the run's uniform grid.

    Initial peak positions are taken from scipy.find_peaks; if fewer
    peaks are found than requested the remaining seeds are evenly spaced.
    Bounds: A ≥ 0, T₀ within [T_min, T_max], W ∈ [1 °C, span], |b| < 0.95.

    Returns a list of DeconvolvedPeak ordered by increasing T₀.
    """
    if len(T) < 4 * n_peaks:
        raise ValueError(
            f"Region has {len(T)} grid points; need ≥ {4*n_peaks} for {n_peaks}-peak fit"
        )

    T_min, T_max = float(T[0]), float(T[-1])
    span = T_max - T_min
    if W_init is None:
        W_init = max(span / (4.0 * n_peaks), 10.0)

    T0_init = _initial_peak_positions(T, y, n_peaks)
    A_init = max(float(y.max()), 1e-6)

    p0: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    for T0 in T0_init:
        p0.extend([A_init, T0, W_init, b_init])
        lower.extend([0.0, T_min, 1.0, -0.95])
        upper.extend([10.0 * A_init, T_max, span, 0.95])

    try:
        popt, _ = curve_fit(
            frazer_suzuki_sum,
            T, y,
            p0=p0,
            bounds=(lower, upper),
            maxfev=20000,
        )
    except (RuntimeError, ValueError) as exc:
        raise RuntimeError(f"Deconvolution failed to converge: {exc}") from exc

    peaks: list[DeconvolvedPeak] = []
    for i in range(0, len(popt), 4):
        A, T0, W, b = popt[i : i + 4]
        area = float(np.trapezoid(frazer_suzuki(T, A, T0, W, b), T))
        peaks.append(DeconvolvedPeak(A=float(A), T0=float(T0), W=float(W), b=float(b), area=area))
    peaks.sort(key=lambda p: p.T0)
    return peaks


def fit_region_baseline(
    T: np.ndarray,
    w: np.ndarray,
    T_lo: float,
    T_hi: float,
    pad_C: float = 20.0,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Linear pre- and post-baselines for a deconvolution region.

    Fits w(T) ≈ a*T + b in two 20 °C-wide windows immediately outside
    [T_lo, T_hi]. Returns ((a_pre, b_pre), (a_post, b_post)).

    Raises ValueError if either window falls off the data array.
    """
    dT = float(np.median(np.diff(T)))
    pad_pts = max(3, int(round(pad_C / dT)))

    i_lo = int(np.searchsorted(T, T_lo))
    i_hi = int(np.searchsorted(T, T_hi))
    pre_end = i_lo
    pre_start = pre_end - pad_pts
    post_start = i_hi
    post_end = post_start + pad_pts

    if pre_start < 0 or post_end > len(T):
        raise ValueError(
            f"Region [{T_lo}, {T_hi}] too close to data edge; need {pad_C} °C of pad"
        )
    a_pre, b_pre = np.polyfit(T[pre_start:pre_end], w[pre_start:pre_end], 1)
    a_post, b_post = np.polyfit(T[post_start:post_end], w[post_start:post_end], 1)
    return (float(a_pre), float(b_pre)), (float(a_post), float(b_post))


def baseline_dtg(
    T_evt: np.ndarray,
    w_evt: np.ndarray,
    a_pre: float,
    b_pre: float,
    a_post: float,
    b_post: float,
    mode: str,
) -> np.ndarray:
    """Baseline dB/dT in the DTG domain for a deconvolution region.

    mode ∈ {"constant", "linear", "alpha"} — see quantify.py for the
    physical justification of each. Same definitions used here.
    """
    if mode == "constant":
        return np.full_like(T_evt, a_pre)
    if mode == "linear":
        span = float(T_evt[-1] - T_evt[0])
        alpha_T = (T_evt - T_evt[0]) / span
        return (1.0 - alpha_T) * a_pre + alpha_T * a_post
    if mode == "alpha":
        dm = float(w_evt[0] - w_evt[-1])
        if abs(dm) < 1e-12:
            return np.full_like(T_evt, 0.5 * (a_pre + a_post))
        alpha = (w_evt[0] - w_evt) / dm
        return (1.0 - alpha) * a_pre + alpha * a_post
    raise ValueError(f"Unknown baseline mode: {mode!r}")
