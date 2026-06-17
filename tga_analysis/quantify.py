"""Three mass-loss quantification methods.

Each method has its own conventional range of integration. For a clean
event with flat baselines the three should agree to within a few percent;
disagreements flag baseline drift, overlap with a neighboring event, or
poor baseline-region selection.

Method 1, "naive":
    Smoothed-mass difference between the event-start and event-end
    temperatures (where |DTG| crosses the baseline-exit threshold). No
    baseline correction. Sensitive to ongoing slow drift under the event.

Method 2, "tangent" (ICTAC sloped-tangent construction):
    Each baseline is a linear fit (sloped, not horizontal) to a quiet
    region adjacent to the event. The mass loss is the difference between
    the pre-baseline value at the onset T and the post-baseline value at
    the endset T, where onset and endset are the intersections of the
    baseline lines with the tangent to the mass curve at the inflection
    point (the DTG peak). This corrects for slow drifts.

Method 3, "dtg_area":
    Integrate the baseline-subtracted DTG over the full event extent
    [event_start_T, event_end_T]. The DTG baseline is linearly
    interpolated from a_pre (pre-baseline slope) at the start to a_post
    (post-baseline slope) at the end. With baseline subtraction this is
    the same physical quantity as "tangent", computed in the derivative
    domain — useful as a self-consistency check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .events import Event


@dataclass
class MassLoss:
    naive: float
    tangent: float
    dtg_area: float


def _interp_at(T: np.ndarray, y: np.ndarray, T_target: float) -> float:
    return float(np.interp(T_target, T, y))


def quantify(
    T: np.ndarray,
    w: np.ndarray,
    dwdT: np.ndarray,
    event: Event,
) -> MassLoss:
    a_pre, b_pre = event.pre_baseline
    a_post, b_post = event.post_baseline

    # Naive: smoothed mass at the event-extent boundaries.
    naive = _interp_at(T, w, event.event_start_T) - _interp_at(T, w, event.event_end_T)

    # Tangent (baseline-corrected) at ICTAC onset/endset.
    tangent = (a_pre * event.onset_T + b_pre) - (a_post * event.endset_T + b_post)

    # DTG area (baseline-corrected) over the full event extent.
    mask = (T >= event.event_start_T) & (T <= event.event_end_T)
    T_evt = T[mask]
    dwdT_evt = dwdT[mask]
    if len(T_evt) < 2 or event.event_end_T <= event.event_start_T:
        dtg_area = float("nan")
    else:
        span = event.event_end_T - event.event_start_T
        alpha = (T_evt - event.event_start_T) / span
        dtg_baseline = (1.0 - alpha) * a_pre + alpha * a_post
        integrand = -(dwdT_evt - dtg_baseline)   # mass loss as positive integrand
        dtg_area = float(np.trapezoid(integrand, T_evt))

    return MassLoss(naive=naive, tangent=tangent, dtg_area=dtg_area)
