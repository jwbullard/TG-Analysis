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


def process_one(file_path: Path, label: str | None, out_dir: Path) -> dict:
    run = load_tga(file_path, label=label)
    sm = smooth_savgol(run.w, run.dT)
    dwdT = dtg_savgol(run.w, run.dT, sm.window_pts, sm.poly_order)
    events = detect_events(run.T, sm.w_smooth, dwdT)

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
        description="TGA batch analysis: smooth, detect events, quantify mass loss."
    )
    p.add_argument("batch_list", type=Path, help="Text file listing TGA files.")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: ./Results/<batch-list-stem>/).",
    )
    args = p.parse_args(argv)

    batch = read_batch_list(args.batch_list)
    if not batch:
        print(f"No entries in batch list {args.batch_list}", file=sys.stderr)
        return 1

    out_dir = args.out or (
        args.batch_list.parent / "Results" / args.batch_list.stem
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    all_event_rows: list[dict] = []
    for file_path, label in batch:
        try:
            info = process_one(file_path, label, out_dir)
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
