"""Mass-loss quantification: naive, ICTAC tangent, and DTG-area (3 baselines).

For the DTG-area method the under-event baseline is hidden by the event
signal and must be inferred from the quiet pre/post regions. Three
constructions are computed in parallel so the user can compare; on a
clean isolated event they agree to within a fraction of a percent, and
disagreements diagnose how drift behaves under the event.

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
    point (the DTG peak). Independent of the under-event baseline choice.

Method 3, "dtg_area" — three baseline variants, three distinct answers:

    "constant":  dB/dT = a_pre, the pre-event slope held across the event.
        Drift is extrapolated forward from the pre-event region only;
        any change in drift slope between pre and post is ignored.
        Most common "single linear baseline" choice in commercial
        software. Justified when post-event drift is poorly defined
        (short or absent post-event quiet region) or when one trusts
        only the pre-event slope estimate.

    "linear":    dB/dT = (1-α_T)·a_pre + α_T·a_post,
                 α_T = (T - T_start)/(T_end - T_start).
        Linear interpolation of slopes in T. Justified when the drift
        mechanism varies smoothly with T and is uncoupled from the
        reaction itself.
        (Note: for the integrated mass loss, this is equivalent to
        holding dB/dT at the midpoint slope (a_pre + a_post)/2 across
        the event — the integral of a linear function equals its
        midpoint times its width. The pointwise baseline shape differs,
        but the area does not.)

    "alpha":     extent-of-reaction blend (Borchardt-Daniels / ICTAC
                 kinetics convention):
        dB/dT(T) = (1 - α(T))·a_pre + α(T)·a_post
        where α(T) = (w[T_start] - w[T]) / (w[T_start] - w[T_end])
        is the cumulative fractional mass loss. Same structure as
        "linear" but α(T) tracks reaction extent rather than T; for a
        peaked DTG, α rises fastest near the inflection, so the slope
        transition concentrates there rather than spreading linearly
        in T. Self-consistent when drift behavior is tied to the state
        of the sample (large mass change altering buoyancy, evolved gas,
        sample geometry). Computed non-iteratively using the already-
        smoothed mass.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .events import Event


@dataclass
class MassLoss:
    naive: float
    tangent: float
    dtg_area_constant: float
    dtg_area_linear: float
    dtg_area_alpha: float


def _interp_at(T: np.ndarray, y: np.ndarray, T_target: float) -> float:
    return float(np.interp(T_target, T, y))


def _dtg_area(
    T_evt: np.ndarray, dwdT_evt: np.ndarray, dtg_baseline: np.ndarray
) -> float:
    integrand = -(dwdT_evt - dtg_baseline)
    return float(np.trapezoid(integrand, T_evt))


def quantify(
    T: np.ndarray,
    w: np.ndarray,
    dwdT: np.ndarray,
    event: Event,
) -> MassLoss:
    a_pre, b_pre = event.pre_baseline
    a_post, b_post = event.post_baseline

    naive = _interp_at(T, w, event.event_start_T) - _interp_at(T, w, event.event_end_T)
    tangent = (a_pre * event.onset_T + b_pre) - (a_post * event.endset_T + b_post)

    mask = (T >= event.event_start_T) & (T <= event.event_end_T)
    T_evt = T[mask]
    dwdT_evt = dwdT[mask]
    w_evt = w[mask]
    if len(T_evt) < 2 or event.event_end_T <= event.event_start_T:
        nan = float("nan")
        return MassLoss(
            naive=naive, tangent=tangent,
            dtg_area_constant=nan, dtg_area_linear=nan, dtg_area_alpha=nan,
        )

    span = event.event_end_T - event.event_start_T
    alpha_T = (T_evt - event.event_start_T) / span

    base_const = np.full_like(T_evt, a_pre)
    dtg_area_constant = _dtg_area(T_evt, dwdT_evt, base_const)

    base_linear = (1.0 - alpha_T) * a_pre + alpha_T * a_post
    dtg_area_linear = _dtg_area(T_evt, dwdT_evt, base_linear)

    dm_total = float(w_evt[0] - w_evt[-1])
    if abs(dm_total) < 1e-12:
        dtg_area_alpha = float("nan")
    else:
        alpha = (w_evt[0] - w_evt) / dm_total           # 0 → 1 across event
        base_alpha = (1.0 - alpha) * a_pre + alpha * a_post
        dtg_area_alpha = _dtg_area(T_evt, dwdT_evt, base_alpha)

    return MassLoss(
        naive=naive,
        tangent=tangent,
        dtg_area_constant=dtg_area_constant,
        dtg_area_linear=dtg_area_linear,
        dtg_area_alpha=dtg_area_alpha,
    )
