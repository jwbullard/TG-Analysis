"""Batch CLI: process a list of TGA files, write summary + per-file outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .derivative import dtg_savgol
from .events import detect_events
from .io import load_tga, read_batch_list
from .plot import plot_run
from .quantify import quantify
from .smooth import smooth_savgol


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
        event_rows.append(
            {
                "file": run.source.name,
                "label": run.label,
                "event": i,
                "peak_T_C": ev.peak_T,
                "onset_T_C": ev.onset_T,
                "endset_T_C": ev.endset_T,
                "event_start_T_C": ev.event_start_T,
                "event_end_T_C": ev.event_end_T,
                "peak_dwdT_pct_per_C": ev.peak_dwdT,
                "mass_loss_naive_pct": ml.naive,
                "mass_loss_tangent_pct": ml.tangent,
                "mass_loss_dtg_area_pct": ml.dtg_area,
            }
        )

    pd.DataFrame(event_rows).to_csv(
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
        "_event_rows": event_rows,
    }


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
    )

    summary_rows: list[dict] = []
    all_event_rows: list[dict] = []
    for file_path, label in batch:
        try:
            info = process_one(file_path, label, out_dir, **tunables)
            all_event_rows.extend(info.pop("_event_rows"))
            summary_rows.append(info)
            print(
                f"[ok]   {file_path.name}: {info['n_events']} event(s), "
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
