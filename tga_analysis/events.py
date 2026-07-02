"""Mass-loss event detection and sloped-tangent construction.

For each event detected as a peak in -DTG, build the ICTAC-style sloped
tangent construction:

  1. Pre-event baseline: a linear fit to the mass-vs-T data in a quiet
     region just before the event (|DTG| below a small fraction of the
     event's peak amplitude).
  2. Post-event baseline: same idea, just after the event.
  3. Inflection tangent: the line through (T_peak, w_smooth(T_peak))
     with slope dW/dT(T_peak), i.e., the tangent to the mass curve at
     the steepest-slope point of the event.
  4. Onset T  = intersection of pre-baseline with the inflection tangent.
     Endset T = intersection of inflection tangent with post-baseline.

This is the standard ICTAC tangent construction and is the form most
referees expect for reported mass losses.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks


@dataclass
class Event:
    peak_idx: int
    peak_T: float
    peak_dwdT: float
    onset_T: float                           # ICTAC tangent intersection (pre-baseline ∩ inflection)
    endset_T: float                          # ICTAC tangent intersection (inflection ∩ post-baseline)
    event_start_T: float                     # where DTG first leaves baseline (threshold crossing)
    event_end_T: float                       # where DTG first returns to baseline
    pre_baseline: tuple[float, float]        # (slope, intercept) in w = a*T + b
    post_baseline: tuple[float, float]
    inflection_tangent: tuple[float, float]
    pre_window: tuple[int, int]              # index range [start, end) of pre-baseline fit
    post_window: tuple[int, int]             # index range [start, end) of post-baseline fit


def _line_fit(T: np.ndarray, w: np.ndarray, i0: int, i1: int) -> tuple[float, float]:
    a, b = np.polyfit(T[i0:i1], w[i0:i1], 1)
    return float(a), float(b)


def _intersect(l1: tuple[float, float], l2: tuple[float, float]) -> float:
    a1, b1 = l1
    a2, b2 = l2
    if abs(a1 - a2) < 1e-15:
        return float("nan")
    return (b2 - b1) / (a1 - a2)


def detect_events(
    T: np.ndarray,
    w: np.ndarray,
    dwdT: np.ndarray,
    prominence_frac: float = 0.02,
    min_separation_C: float = 20.0,
    baseline_exit_threshold_frac: float = 0.02,
    baseline_fit_width_C: float = 20.0,
) -> list[Event]:
    """Detect mass-loss events and construct sloped tangents around each.

    Parameters
    ----------
    T, w, dwdT : arrays
        Temperature, smoothed weight, and smoothed dW/dT on a uniform grid.
    prominence_frac : float
        Required prominence of a peak in |DTG|, expressed as a fraction of
        the run's maximum |DTG|. Default 2% suppresses noise peaks.
    min_separation_C : float
        Minimum °C separation between detected peaks. Peaks closer than
        this are merged into the higher one. For events that genuinely
        overlap within this window, peak deconvolution is required and is
        not done here.
    baseline_exit_threshold_frac : float
        Walking outward from a peak, the location where |DTG| first drops
        below this fraction of the peak amplitude defines the edge of the
        event. The baseline fit window is anchored at that edge.
    baseline_fit_width_C : float
        Width (°C) of the linear baseline fit window adjacent to each
        event edge. Wide enough to suppress noise, narrow enough to keep
        the slope local. Clipped if a neighboring peak or array boundary
        is closer than this width.
    """
    dT = float(np.median(np.diff(T)))
    if dT <= 0:
        raise ValueError("Non-uniform or descending T grid passed to detect_events.")
    mag = np.where(-dwdT > 0, -dwdT, 0.0)
    max_mag = float(mag.max())
    if max_mag <= 0:
        return []

    peaks, _ = find_peaks(
        mag,
        prominence=prominence_frac * max_mag,
        distance=max(1, int(round(min_separation_C / dT))),
    )

    fit_width_pts = max(3, int(round(baseline_fit_width_C / dT)))
    events: list[Event] = []

    for k, p_idx in enumerate(peaks):
        peak_mag = mag[p_idx]
        thr = baseline_exit_threshold_frac * peak_mag

        # Pre-baseline: walk left from the peak until |DTG| drops below thr,
        # then fit a fixed-width window ending at that crossing point.
        # If no quiet point exists before the adjacent peak (crowded events,
        # e.g. cement paste below 200 °C), fall back to the inter-peak DTG
        # minimum — the closest thing to a baseline that exists locally.
        left_bound = 0 if k == 0 else int(peaks[k - 1])
        i = int(p_idx) - 1
        while i > left_bound and mag[i] > thr:
            i -= 1
        if i <= left_bound and mag[max(i, left_bound) + 1] > thr:
            # No quiet point found; use inter-peak DTG minimum as the anchor.
            i = int(left_bound + 1 + np.argmin(mag[left_bound + 1 : int(p_idx)]))
        end_pre = i + 1
        start_pre = max(left_bound, end_pre - fit_width_pts)

        # Post-baseline: walk right until |DTG| drops below thr; fit a
        # fixed-width window starting at the crossing point. Same fallback
        # as above when no quiet point exists before the next peak.
        right_bound = len(T) - 1 if k == len(peaks) - 1 else int(peaks[k + 1])
        j = int(p_idx) + 1
        while j < right_bound and mag[j] > thr:
            j += 1
        if j >= right_bound and mag[min(j, right_bound) - 1] > thr:
            j = int(p_idx + 1 + np.argmin(mag[int(p_idx) + 1 : right_bound]))
        start_post = j
        end_post = min(right_bound + 1, start_post + fit_width_pts)

        if end_pre - start_pre < 2 or end_post - start_post < 2:
            continue

        pre_line = _line_fit(T, w, start_pre, end_pre)
        post_line = _line_fit(T, w, start_post, end_post)

        slope_inflect = float(dwdT[p_idx])
        intercept_inflect = float(w[p_idx] - slope_inflect * T[p_idx])
        inflect_line = (slope_inflect, intercept_inflect)

        T_onset = _intersect(pre_line, inflect_line)
        T_endset = _intersect(inflect_line, post_line)

        # Sanity: if intersections fall outside the data range, fall back
        # to threshold-crossing T's (less precise but always well-defined).
        if not np.isfinite(T_onset) or T_onset < T[start_pre] or T_onset > T[p_idx]:
            T_onset = float(T[end_pre - 1])
        if not np.isfinite(T_endset) or T_endset < T[p_idx] or T_endset > T[end_post - 1]:
            T_endset = float(T[start_post])

        events.append(
            Event(
                peak_idx=int(p_idx),
                peak_T=float(T[p_idx]),
                peak_dwdT=float(dwdT[p_idx]),
                onset_T=float(T_onset),
                endset_T=float(T_endset),
                event_start_T=float(T[end_pre]),
                event_end_T=float(T[min(start_post, len(T) - 1)]),
                pre_baseline=pre_line,
                post_baseline=post_line,
                inflection_tangent=inflect_line,
                pre_window=(int(start_pre), int(end_pre)),
                post_window=(int(start_post), int(end_post)),
            )
        )

    return events
