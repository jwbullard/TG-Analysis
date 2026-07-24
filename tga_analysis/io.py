"""I/O for TGA data files and batch lists.

Two CSV flavors are handled transparently:

  * Plain table — header on row 1, data from row 2.
  * TA Instruments TRIOS export — ~130 lines of [Section] metadata,
    then a column-header row containing both 'Temperature' and 'Weight',
    then a units row ('min,°C,mg,%,...'), then data. Two columns named
    'Weight' (mg and %) are disambiguated by the units row.

Excel files (.xlsx, .xls) are passed through pandas' Excel reader.
TA Instruments SDT 650 full exports carry two 'Weight' columns (mg and
%). Some exports include a units row directly under the header
('min,°C,mg,%,...') and some do not; both are handled. When the units
row is present it is used to rename the weight columns; otherwise the
first non-NaN value picks the % (larger) vs the mg (smaller).

Column names are matched by keyword so vendor variations ('Sample
Temperature', 'Weight %', 'Deriv. Weight', 'DTG', etc.) are picked up
automatically. When both 'Weight (mg)' and 'Weight (%)' are present, the
percentage column is preferred.
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
    # Prefer Weight (%) over Weight (mg) when both are present (TRIOS export
    # carries both; the % column is what downstream code assumes).
    pct_keys = [
        k for k in cols
        if "weight" in k and "%" in k and "deriv" not in k and "dtg" not in k
    ]
    if pct_keys:
        W_col = cols[pct_keys[0]]
    else:
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


def _looks_numeric(row: pd.Series) -> bool:
    try:
        pd.to_numeric(row)
        return True
    except (ValueError, TypeError):
        return False


def _read_csv(path: Path) -> pd.DataFrame:
    """Read a CSV that may have a TRIOS-style metadata preamble.

    Scans the first ~500 lines for a header row that contains both a
    'Temperature' and a 'Weight' token; if found, reads from there and
    drops a units row if present. Falls back to a plain pd.read_csv if
    no preamble is detected.
    """
    header_row: int | None = None
    with path.open("r", encoding="utf-8-sig", errors="replace") as f:
        for i, line in enumerate(f):
            tokens = [t.strip().lower() for t in line.split(",")]
            has_T = any("temperature" in t for t in tokens)
            has_W = any("weight" in t for t in tokens)
            if has_T and has_W and len(tokens) >= 2:
                header_row = i
                break
            if i > 500:
                break

    if header_row is None or header_row == 0:
        df = pd.read_csv(path, encoding="utf-8-sig")
    else:
        df = pd.read_csv(path, skiprows=header_row, header=0, encoding="utf-8-sig")

    if len(df) > 0 and not _looks_numeric(df.iloc[0]):
        units = df.iloc[0].astype(str).tolist()
        df = df.iloc[1:].reset_index(drop=True)
        new_cols: list[str] = []
        for col, unit in zip(df.columns, units):
            base = str(col).strip()
            unit = unit.strip()
            if base.lower().startswith("weight") and not (
                "deriv" in base.lower() or "dtg" in base.lower()
            ):
                if "%" in unit:
                    new_cols.append("Weight (%)")
                elif "mg" in unit.lower():
                    new_cols.append("Weight (mg)")
                else:
                    new_cols.append(f"{base} ({unit})")
            else:
                new_cols.append(base)
        df.columns = new_cols

    return df.apply(pd.to_numeric, errors="coerce")


def _disambiguate_weight_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename bare 'Weight' columns to 'Weight (%)' or 'Weight (mg)'.

    TA Instruments SDT 650 xlsx exports emit sample mass twice: once in
    mg and once as a percentage of the initial mass. Both columns are
    literally headed 'Weight' with no units row, so pandas renames the
    second to 'Weight.1'. Neither carries a token that the downstream
    column detector can key on. Given exactly two bare 'Weight' columns,
    the one with the larger first non-NaN value is the % (bounded at
    100 and typically starting there, though hydrated pastes may have
    lost a few percent during the pre-ramp equilibration hold); the
    smaller one is mg (typical sample loading 20-70 mg on the SDT 650).
    Derivative and heat-flow columns are ignored.
    """
    mass_cols = [
        c for c in df.columns
        if str(c).lower().strip().startswith("weight")
        and "deriv" not in str(c).lower()
        and "dtg" not in str(c).lower()
        and "heat" not in str(c).lower()
        and "%" not in str(c)
        and "mg" not in str(c).lower()
    ]
    if len(mass_cols) != 2:
        return df
    first_vals: list[tuple[str, float]] = []
    for col in mass_cols:
        series = df[col].dropna()
        if series.empty:
            return df
        first_vals.append((col, float(series.iloc[0])))
    (col_a, val_a), (col_b, val_b) = first_vals
    pct_col, mg_col = (col_a, col_b) if val_a > val_b else (col_b, col_a)
    return df.rename(columns={pct_col: "Weight (%)", mg_col: "Weight (mg)"})


def _strip_units_row(df: pd.DataFrame) -> pd.DataFrame:
    """If row 0 is a units row, use it to rename Weight columns and drop it.

    Some TA Instruments SDT 650 xlsx exports emit a units row directly
    under the column-header row ('min, °C, mg, %, % / °C, mW, ...').
    Others emit only the header. When present, the units row lets us
    disambiguate the two identically-named 'Weight' columns
    unambiguously (rather than relying on the value-range heuristic),
    but pandas otherwise treats it as string data and downstream numeric
    coercion fails on values like 'mg'.
    """
    if len(df) == 0 or _looks_numeric(df.iloc[0]):
        return df
    units = df.iloc[0].astype(str).tolist()
    df = df.iloc[1:].reset_index(drop=True)
    new_cols: list[str] = []
    for col, unit in zip(df.columns, units):
        base = str(col).strip()
        unit_s = unit.strip()
        if base.lower().startswith("weight") and not (
            "deriv" in base.lower() or "dtg" in base.lower()
        ):
            if "%" in unit_s:
                new_cols.append("Weight (%)")
            elif "mg" in unit_s.lower():
                new_cols.append("Weight (mg)")
            else:
                new_cols.append(f"{base} ({unit_s})")
        else:
            new_cols.append(base)
    df.columns = new_cols
    return df.apply(pd.to_numeric, errors="coerce")


def _read_raw(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xls"):
        return _disambiguate_weight_columns(_strip_units_row(pd.read_excel(path)))
    if path.suffix.lower() == ".csv":
        return _read_csv(path)
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
