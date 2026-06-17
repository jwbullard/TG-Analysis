"""Diagnostic plot: raw + smoothed mass, DTG, and per-event tangent overlays."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .events import Event


def plot_run(
    T: np.ndarray,
    w_raw: np.ndarray,
    w_smooth: np.ndarray,
    dwdT: np.ndarray,
    events: list[Event],
    title: str,
    out_path: Path,
) -> None:
    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, figsize=(9.5, 7.0),
        gridspec_kw={"hspace": 0.05}, constrained_layout=True,
    )

    ax1.plot(T, w_raw, color="lightgrey", lw=0.6, label="raw")
    ax1.plot(T, w_smooth, color="C0", lw=1.3, label="smoothed")
    ax1.set_ylabel("Weight (%)")
    ax1.set_title(title, fontsize=10)
    ax1.legend(loc="best", fontsize=8)

    ax2.plot(T, dwdT, color="C3", lw=1.0)
    ax2.axhline(0, color="k", lw=0.4)
    ax2.set_xlabel("Temperature (°C)")
    ax2.set_ylabel("dW/dT (% per °C)")

    for ev in events:
        # Shade the baseline fit regions.
        for (i0, i1) in (ev.pre_window, ev.post_window):
            i1 = min(i1, len(T))
            if i1 > i0:
                ax1.axvspan(T[i0], T[i1 - 1], color="C2", alpha=0.10)
                ax2.axvspan(T[i0], T[i1 - 1], color="C2", alpha=0.10)
        # Pre/post baselines extended into the event for visibility.
        Tp_pre = np.array([T[ev.pre_window[0]], ev.onset_T])
        ax1.plot(
            Tp_pre, ev.pre_baseline[0] * Tp_pre + ev.pre_baseline[1],
            "--", color="C2", lw=1.0,
        )
        Tp_post = np.array(
            [ev.endset_T, T[min(ev.post_window[1] - 1, len(T) - 1)]]
        )
        ax1.plot(
            Tp_post, ev.post_baseline[0] * Tp_post + ev.post_baseline[1],
            "--", color="C2", lw=1.0,
        )
        # Inflection tangent.
        Tp_inf = np.array([ev.onset_T, ev.endset_T])
        ax1.plot(
            Tp_inf,
            ev.inflection_tangent[0] * Tp_inf + ev.inflection_tangent[1],
            "--", color="C1", lw=1.0,
        )
        # Markers: onset/peak/endset (tangent construction).
        for Tm in (ev.onset_T, ev.peak_T, ev.endset_T):
            for ax in (ax1, ax2):
                ax.axvline(Tm, color="k", lw=0.4, ls=":")
        # Event-extent markers (where DTG departs/returns to baseline).
        for Tm in (ev.event_start_T, ev.event_end_T):
            for ax in (ax1, ax2):
                ax.axvline(Tm, color="C2", lw=0.4, ls="-.")
        ax2.plot([ev.peak_T], [ev.peak_dwdT], "o", color="C3", ms=4)

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
