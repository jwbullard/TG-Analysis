"""Batch CLI: process a list of TGA files, write summary + per-file outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .deconvolve import (
    baseline_dtg,
    deconvolve,
    fit_region_baseline,
    frazer_suzuki,
)
from .derivative import dtg_savgol
from .events import detect_events
from .io import load_tga, read_batch_list
from .plot import plot_run
from .quantify import quantify
from .smooth import smooth_savgol


_BASELINE_MODES = ("constant", "linear", "alpha")
_BASELINE_COLS = {
    "constant": "mass_loss_dtg_area_constant_pct",
    "linear":   "mass_loss_dtg_area_linear_pct",
    "alpha":    "mass_loss_dtg_area_alpha_pct",
}


def _empty_event_row(file_name: str, label: str) -> dict:
    """Skeleton row with all columns present and NaN — keeps the events CSV
    schema uniform across detected and deconvolved entries."""
    return {
        "file": file_name,
        "label": label,
        "event": "",
        "event_kind": "",
        "peak_T_C": float("nan"),
        "onset_T_C": float("nan"),
        "endset_T_C": float("nan"),
        "event_start_T_C": float("nan"),
        "event_end_T_C": float("nan"),
        "peak_dwdT_pct_per_C": float("nan"),
        "mass_loss_naive_pct": float("nan"),
        "mass_loss_tangent_pct": float("nan"),
        "mass_loss_dtg_area_constant_pct": float("nan"),
        "mass_loss_dtg_area_linear_pct": float("nan"),
        "mass_loss_dtg_area_alpha_pct": float("nan"),
        "fwhm_C": float("nan"),
        "asymmetry_b": float("nan"),
    }


def _plot_deconvolution(
    T_evt: np.ndarray,
    y_by_mode: dict[str, np.ndarray],
    peaks_by_mode: dict[str, list],
    label: str,
    region_idx: int,
    out_path: Path,
) -> None:
    """Stack one panel per baseline mode showing the corrected DTG, each
    fitted Frazer-Suzuki peak, and the sum-of-peaks reconstruction."""
    fig, axes = plt.subplots(
        len(_BASELINE_MODES), 1,
        figsize=(8.5, 3.0 * len(_BASELINE_MODES)),
        sharex=True, constrained_layout=True,
    )
    if len(_BASELINE_MODES) == 1:
        axes = [axes]
    for ax, mode in zip(axes, _BASELINE_MODES):
        y = y_by_mode[mode]
        peaks = peaks_by_mode[mode]
        ax.plot(T_evt, y, color="k", lw=0.8, label="corrected DTG")
        total = np.zeros_like(T_evt)
        for k, p in enumerate(peaks):
            y_p = frazer_suzuki(T_evt, p.A, p.T0, p.W, p.b)
            ax.fill_between(T_evt, 0, y_p, alpha=0.25,
                            label=f"P{k+1} @ {p.T0:.0f} °C  area = {p.area:.2f} %")
            total = total + y_p
        ax.plot(T_evt, total, color="C3", lw=1.0, ls="--", label="Σ FS peaks")
        ax.axhline(0, color="k", lw=0.4)
        ax.set_ylabel("-dW/dT (% / °C)")
        ax.set_title(f"baseline mode: {mode}", fontsize=10, loc="left")
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Temperature (°C)")
    fig.suptitle(f"{label}    deconvolution region {region_idx}", fontsize=11)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def process_one(
    file_path: Path,
    label: str | None,
    out_dir: Path,
    *,
    grid_dT: float = 0.05,
    sg_poly: int = 3,
    sg_window_C_range: tuple[float, float] = (0.3, 5.0),
    sg_n_candidates: int = 20,
    sg_rho1_threshold: float = 0.05,
    peak_prominence_frac: float = 0.02,
    peak_min_separation_C: float = 20.0,
    baseline_exit_threshold_frac: float = 0.02,
    baseline_fit_width_C: float = 20.0,
    decon_regions: list[tuple[float, float, int]] | None = None,
) -> dict:
    run = load_tga(file_path, label=label, dT=grid_dT)
    sm = smooth_savgol(
        run.w,
        run.dT,
        poly_order=sg_poly,
        window_C_range=sg_window_C_range,
        n_candidates=sg_n_candidates,
        rho1_threshold=sg_rho1_threshold,
    )
    dwdT = dtg_savgol(run.w, run.dT, sm.window_pts, sm.poly_order)
    events = detect_events(
        run.T,
        sm.w_smooth,
        dwdT,
        prominence_frac=peak_prominence_frac,
        min_separation_C=peak_min_separation_C,
        baseline_exit_threshold_frac=baseline_exit_threshold_frac,
        baseline_fit_width_C=baseline_fit_width_C,
    )

    event_rows: list[dict] = []
    for i, ev in enumerate(events, start=1):
        ml = quantify(run.T, sm.w_smooth, dwdT, ev)
        row = _empty_event_row(run.source.name, run.label)
        row.update({
            "event": str(i),
            "event_kind": "detected",
            "peak_T_C": ev.peak_T,
            "onset_T_C": ev.onset_T,
            "endset_T_C": ev.endset_T,
            "event_start_T_C": ev.event_start_T,
            "event_end_T_C": ev.event_end_T,
            "peak_dwdT_pct_per_C": ev.peak_dwdT,
            "mass_loss_naive_pct": ml.naive,
            "mass_loss_tangent_pct": ml.tangent,
            "mass_loss_dtg_area_constant_pct": ml.dtg_area_constant,
            "mass_loss_dtg_area_linear_pct": ml.dtg_area_linear,
            "mass_loss_dtg_area_alpha_pct": ml.dtg_area_alpha,
        })
        event_rows.append(row)

    decon_rows: list[dict] = []
    for r_idx, (T_lo, T_hi, n_peaks) in enumerate(decon_regions or [], start=1):
        try:
            (a_pre, b_pre), (a_post, b_post) = fit_region_baseline(
                run.T, sm.w_smooth, T_lo, T_hi
            )
        except ValueError as exc:
            print(f"[warn] {run.label}: decon region {r_idx} skipped: {exc}",
                  file=sys.stderr)
            continue

        mask = (run.T >= T_lo) & (run.T <= T_hi)
        T_evt = run.T[mask]
        w_evt = sm.w_smooth[mask]
        dwdT_evt = dwdT[mask]
        if len(T_evt) < 4 * n_peaks:
            print(f"[warn] {run.label}: decon region {r_idx} too narrow for "
                  f"{n_peaks} peaks", file=sys.stderr)
            continue

        peaks_by_mode: dict[str, list] = {}
        y_by_mode: dict[str, np.ndarray] = {}
        for mode in _BASELINE_MODES:
            B = baseline_dtg(T_evt, w_evt, a_pre, b_pre, a_post, b_post, mode)
            y = -(dwdT_evt - B)
            y_by_mode[mode] = y
            try:
                peaks_by_mode[mode] = deconvolve(T_evt, y, n_peaks)
            except RuntimeError as exc:
                print(f"[warn] {run.label}: decon region {r_idx} mode {mode}: {exc}",
                      file=sys.stderr)
                peaks_by_mode[mode] = []

        # Pair deconvolved peaks across modes by ascending T0 (already sorted).
        max_peaks = max((len(p) for p in peaks_by_mode.values()), default=0)
        for k in range(max_peaks):
            row = _empty_event_row(run.source.name, run.label)
            row.update({
                "event": f"decon-R{r_idx}-P{k+1}",
                "event_kind": "deconvolved",
                "event_start_T_C": T_lo,
                "event_end_T_C": T_hi,
            })
            # Take canonical T0/W/b/A from the linear mode (or first non-empty).
            canonical = (
                peaks_by_mode["linear"]
                or peaks_by_mode["constant"]
                or peaks_by_mode["alpha"]
            )
            if k < len(canonical):
                p = canonical[k]
                row["peak_T_C"] = p.T0
                row["peak_dwdT_pct_per_C"] = -p.A
                row["fwhm_C"] = p.W
                row["asymmetry_b"] = p.b
            for mode in _BASELINE_MODES:
                if k < len(peaks_by_mode[mode]):
                    row[_BASELINE_COLS[mode]] = peaks_by_mode[mode][k].area
            decon_rows.append(row)

        _plot_deconvolution(
            T_evt, y_by_mode, peaks_by_mode,
            run.label, r_idx,
            out_dir / f"{run.label}_decon_R{r_idx}.png",
        )

    pd.DataFrame(event_rows + decon_rows).to_csv(
        out_dir / f"{run.label}_events.csv", index=False
    )

    plot_run(
        run.T,
        run.w,
        sm.w_smooth,
        dwdT,
        events,
        title=(
            f"{run.label}    "
            f"SG window = {sm.window_C:.2f} °C, order = {sm.poly_order}, "
            f"|ρ₁| = {abs(sm.rho1):.3f}, DW = {sm.durbin_watson:.2f}"
        ),
        out_path=out_dir / f"{run.label}.png",
    )

    return {
        "file": run.source.name,
        "label": run.label,
        "n_points": len(run.T),
        "T_min_C": float(run.T.min()),
        "T_max_C": float(run.T.max()),
        "sg_window_C": sm.window_C,
        "sg_window_pts": sm.window_pts,
        "sg_poly_order": sm.poly_order,
        "residual_lag1_autocorr": sm.rho1,
        "durbin_watson": sm.durbin_watson,
        "residual_std": sm.residual_std,
        "n_events": len(events),
        "n_decon_peaks": len(decon_rows),
        "_event_rows": event_rows + decon_rows,
    }


def _parse_decon_region(s: str) -> tuple[float, float, int]:
    parts = [tok.strip() for tok in s.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"--deconvolve-region expects 'Tlo,Thi,Npeaks', got {s!r}"
        )
    try:
        Tlo, Thi = float(parts[0]), float(parts[1])
        N = int(parts[2])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--deconvolve-region parse error: {exc}")
    if Thi <= Tlo:
        raise argparse.ArgumentTypeError(f"Thi must exceed Tlo in {s!r}")
    if N < 1:
        raise argparse.ArgumentTypeError(f"Npeaks must be ≥ 1 in {s!r}")
    return Tlo, Thi, N


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="TGA batch analysis: smooth, detect events, quantify mass loss.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("batch_list", type=Path, help="Text file listing TGA files.")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: <batch-list-dir>/Results/<batch-list-stem>/).",
    )

    g_grid = p.add_argument_group("grid")
    g_grid.add_argument(
        "--grid-dT",
        type=float,
        default=0.05,
        help="Uniform T-grid spacing (°C) for resampling.",
    )

    g_sg = p.add_argument_group("Savitzky-Golay smoothing")
    g_sg.add_argument("--sg-poly", type=int, default=3, help="SG polynomial order.")
    g_sg.add_argument(
        "--sg-window-min",
        type=float,
        default=0.3,
        help="Lower bound (°C) of the SG window search range.",
    )
    g_sg.add_argument(
        "--sg-window-max",
        type=float,
        default=5.0,
        help="Upper bound (°C) of the SG window search range.",
    )
    g_sg.add_argument(
        "--sg-n-candidates",
        type=int,
        default=20,
        help="Number of candidate windows sampled across the range.",
    )
    g_sg.add_argument(
        "--sg-rho1-threshold",
        type=float,
        default=0.05,
        help="Max |lag-1 residual autocorrelation| for the chosen window.",
    )

    g_peak = p.add_argument_group("peak detection / baseline")
    g_peak.add_argument(
        "--peak-prominence",
        type=float,
        default=0.02,
        help="Required |DTG| peak prominence as a fraction of run max |DTG|.",
    )
    g_peak.add_argument(
        "--peak-min-sep",
        type=float,
        default=20.0,
        help="Minimum °C separation between detected peaks.",
    )
    g_peak.add_argument(
        "--baseline-threshold",
        type=float,
        default=0.02,
        help="|DTG| fraction of peak below which an event edge is declared.",
    )
    g_peak.add_argument(
        "--baseline-width",
        type=float,
        default=20.0,
        help="Width (°C) of the linear baseline fit window at each event edge.",
    )

    g_dec = p.add_argument_group("deconvolution")
    g_dec.add_argument(
        "--deconvolve-region",
        action="append",
        type=_parse_decon_region,
        default=[],
        metavar="Tlo,Thi,N",
        help="Fit a sum of N Frazer-Suzuki peaks to the baseline-corrected DTG "
             "over [Tlo, Thi] °C. May be repeated for multiple regions. "
             "Results land in <label>_events.csv as event_kind='deconvolved' "
             "rows, and a diagnostic plot is written per region per file.",
    )

    args = p.parse_args(argv)

    if args.sg_window_max <= args.sg_window_min:
        print(
            "--sg-window-max must exceed --sg-window-min",
            file=sys.stderr,
        )
        return 2

    batch = read_batch_list(args.batch_list)
    if not batch:
        print(f"No entries in batch list {args.batch_list}", file=sys.stderr)
        return 1

    out_dir = args.out or (
        args.batch_list.parent / "Results" / args.batch_list.stem
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    tunables = dict(
        grid_dT=args.grid_dT,
        sg_poly=args.sg_poly,
        sg_window_C_range=(args.sg_window_min, args.sg_window_max),
        sg_n_candidates=args.sg_n_candidates,
        sg_rho1_threshold=args.sg_rho1_threshold,
        peak_prominence_frac=args.peak_prominence,
        peak_min_separation_C=args.peak_min_sep,
        baseline_exit_threshold_frac=args.baseline_threshold,
        baseline_fit_width_C=args.baseline_width,
        decon_regions=args.deconvolve_region,
    )

    summary_rows: list[dict] = []
    all_event_rows: list[dict] = []
    for file_path, label in batch:
        try:
            info = process_one(file_path, label, out_dir, **tunables)
            all_event_rows.extend(info.pop("_event_rows"))
            summary_rows.append(info)
            decon_note = (
                f", {info.get('n_decon_peaks', 0)} decon peak(s)"
                if args.deconvolve_region else ""
            )
            print(
                f"[ok]   {file_path.name}: {info['n_events']} event(s){decon_note}, "
                f"SG window = {info['sg_window_C']:.2f} °C"
            )
        except Exception as exc:
            print(f"[fail] {file_path.name}: {exc}", file=sys.stderr)
            summary_rows.append(
                {"file": file_path.name, "label": label, "error": str(exc)}
            )

    pd.DataFrame(summary_rows).to_csv(out_dir / "summary.csv", index=False)
    pd.DataFrame(all_event_rows).to_csv(out_dir / "events_all.csv", index=False)
    print(f"\nResults written to: {out_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
