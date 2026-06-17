"""I/O for TGA data files and batch lists.

Reader accepts Excel (.xlsx, .xls) or CSV. The expected columns are
temperature (°C), weight (% of initial mass), and optionally an
instrument-supplied derivative dW/dT. Column names are matched by keyword
so vendor variations ("Sample Temperature", "Weight %", "Deriv. Weight",
"DTG", etc.) are picked up automatically.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class TGARun:
    """A single TGA ramp, after cleanup and uniform-T resampling."""

    label: str
    source: Path
    T: np.ndarray                          # uniform-grid temperature (°C)
    w: np.ndarray                          # weight (% of initial mass)
    w_instr_dwdT: np.ndarray | None        # instrument DTG on the same grid, or None
    dT: float                              # grid spacing (°C)


def _detect_columns(df: pd.DataFrame) -> tuple[str, str, str | None]:
    cols = {str(c).lower().strip(): c for c in df.columns}
    T_col = next((cols[k] for k in cols if "temp" in k), None)
    W_col = next(
        (cols[k] for k in cols if "weight" in k and "deriv" not in k and "dtg" not in k),
        None,
    )
    D_col = next(
        (cols[k] for k in cols if "deriv" in k or "dtg" in k or "dw/dt" in k),
        None,
    )
    if T_col is None or W_col is None:
        raise ValueError(
            f"Could not identify temperature and weight columns in {list(df.columns)}"
        )
    return T_col, W_col, D_col


def _read_raw(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported file extension: {path.suffix}")


def load_tga(
    path: str | Path,
    label: str | None = None,
    dT: float = 0.05,
) -> TGARun:
    """Load a TGA file and resample to a uniform T grid.

    Steps:
      1. Read the file and identify T, weight, and (optional) DTG columns.
      2. Drop rows with NaN in T or weight.
      3. Truncate at the temperature maximum (drop the post-ramp isothermal
         hold or cooling tail that shows up as flat or decreasing T).
      4. Average duplicate-T rows, sort, and interpolate onto a uniform grid
         of spacing dT (°C). The default 0.05 °C ≈ 20 samples per °C is fine
         for SG smoothing with the windows used downstream.
    """
    path = Path(path)
    raw = _read_raw(path)
    T_col, W_col, D_col = _detect_columns(raw)

    keep_cols = [T_col, W_col] + ([D_col] if D_col else [])
    df = raw[keep_cols].copy()
    df = df.rename(
        columns={T_col: "T", W_col: "w", **({D_col: "d"} if D_col else {})}
    )
    df = df.dropna(subset=["T", "w"]).astype({"T": float, "w": float})
    if df.empty:
        raise ValueError(f"No valid (T, weight) rows in {path}")

    i_max = int(df["T"].values.argmax())
    df = df.iloc[: i_max + 1]
    df = df.groupby("T", as_index=False).mean(numeric_only=True)
    df = df.sort_values("T").reset_index(drop=True)
    if len(df) < 2:
        raise ValueError(f"Not enough monotonic ramp data in {path}")

    T_uniform = np.arange(df["T"].iloc[0], df["T"].iloc[-1], dT)
    w_uniform = np.interp(T_uniform, df["T"].values, df["w"].values)
    d_uniform = (
        np.interp(T_uniform, df["T"].values, df["d"].values)
        if "d" in df.columns
        else None
    )

    return TGARun(
        label=label or path.stem,
        source=path,
        T=T_uniform,
        w=w_uniform,
        w_instr_dwdT=d_uniform,
        dT=dT,
    )


def read_batch_list(path: str | Path) -> list[tuple[Path, str | None]]:
    """Parse a batch list. Format: one file per line, with an optional label.

        Fortera_Calcite_1.xlsx
        Fortera_Aragonite_1.xlsx  Aragonite-rep1
        # Fortera_Vaterite_4.xlsx   ← skipped
        ../OtherDir/Run.xlsx      Other-run

    Paths are resolved relative to the batch-list file's parent directory.
    A label, if given, is whatever follows the path on the same line
    (whitespace-separated, first split only — labels with spaces are kept).
    """
    path = Path(path)
    base = path.parent
    entries: list[tuple[Path, str | None]] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        file_path = (base / parts[0]).resolve()
        label = parts[1].strip() if len(parts) > 1 else None
        entries.append((file_path, label))
    return entries
