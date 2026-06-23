# TG-Analysis

Thermogravimetric (TGA) signal denoising and mass-loss event quantification.
Reads vendor Excel/CSV exports, picks an optimal Savitzky-Golay smoother
from residual whiteness, detects mass-loss events in the DTG, and reports
mass loss by three methods (naive, ICTAC sloped tangent, DTG-area).

Author: Jeffrey W. Bullard. Code lives at `~/Code/TG-Analysis/`.

## Requirements

- Python ≥ 3.10
- numpy, scipy, pandas, matplotlib, openpyxl (pulled in by `pyproject.toml`)

The package is installed editable into the master venv at
`~/Code/Python/Envs/Default/`. Activate it with the `py` alias from
`~/.zshrc`, or directly:

```sh
source ~/Code/Python/Envs/Default/bin/activate
```

You'll know the venv is live when `which tga-analyze` resolves to a path
under `~/Code/Python/Envs/Default/bin/`.

## Installation

If the master venv ever loses the install (fresh machine, etc.):

```sh
source ~/Code/Python/Envs/Default/bin/activate
cd ~/Code/TG-Analysis
pip install -e .
```

This registers the `tga-analyze` console script and keeps the source
editable in place.

## Quick start

```sh
source ~/Code/Python/Envs/Default/bin/activate
cd ~/Research/TG-CaCO3-Meihe/Data/Raw
tga-analyze batch_calcite.txt
```

Outputs land at `<data-dir>/Results/<batch-list-stem>/`. Override with
`--out <dir>`.

## Input data

Vendor Excel (`.xlsx`, `.xls`) or CSV with at minimum:

- a temperature column (anything containing "temp"),
- a weight column in % of initial mass (anything containing "weight",
  but not "deriv" or "dtg"),
- optionally an instrument-supplied DTG column (matched on "deriv",
  "dtg", or "dw/dt").

Column names are picked up case-insensitively by keyword, so Meihe's
columns (`Temperature/°C`, `Weight/%`, `Deriv. Weight/(%/°C)`) work
without configuration.

The loader drops the post-ramp isothermal hold (truncates at the
temperature max), averages duplicate-T rows, and resamples onto a
uniform °C grid (default 0.05 °C; tune with `--grid-dT`).

## Batch list format

One file per line. Optional whitespace-separated label. Lines starting
with `#` are skipped. Paths resolve relative to the batch-list file's
parent.

```
Fortera_Calcite_1.xlsx        Calcite-1
Fortera_Aragonite_1.xlsx      Aragonite-rep1
# Fortera_Vaterite_4.xlsx     ← excluded
../OtherDir/Run.xlsx          Other-run
```

## Output files

Per run (one set per file in the batch list):

- `<label>.png` — plot of mass and DTG vs. T with detected events,
  baselines, and onset/endset markers
- `<label>_events.csv` — per-event onset/endset/peak temperatures plus
  the three mass-loss values

Batch-level:

- `summary.csv` — one row per input file (smoothing diagnostics, event
  count, error if any)
- `events_all.csv` — all event rows across all files, concatenated

## CLI flags

`tga-analyze --help` is the source of truth; the table below is a
snapshot.

**Grid**

| Flag | Default | What it controls |
|---|---|---|
| `--grid-dT` | 0.05 | Uniform T-grid spacing (°C) for resampling |

**Savitzky-Golay smoothing**

| Flag | Default | What it controls |
|---|---|---|
| `--sg-poly` | 3 | SG polynomial order |
| `--sg-window-min` | 0.3 | Lower bound of SG window search (°C) |
| `--sg-window-max` | 5.0 | Upper bound of SG window search (°C) |
| `--sg-n-candidates` | 20 | Number of candidate windows sampled |
| `--sg-rho1-threshold` | 0.05 | Max acceptable \|residual lag-1 autocorrelation\| |

Window selection picks the *largest* candidate whose residual lag-1
autocorrelation stays below the threshold — i.e., the most aggressive
smoothing that has not yet removed real signal.

**Peak detection / baseline**

| Flag | Default | What it controls |
|---|---|---|
| `--peak-prominence` | 0.02 | Required \|DTG\| peak prominence (fraction of run max) |
| `--peak-min-sep` | 20.0 | Minimum °C separation between peaks |
| `--baseline-threshold` | 0.02 | \|DTG\| fraction of peak below which an event edge is declared |
| `--baseline-width` | 20.0 | Width (°C) of linear baseline fit window at each edge |

**Output location**

| Flag | Default | What it controls |
|---|---|---|
| `--out` | `<batch-dir>/Results/<batch-stem>/` | Output directory |

## Mass-loss quantification

Each event reports three values:

1. **naive** — smoothed-mass difference between event-start and
   event-end T (no baseline correction). Sensitive to ongoing slow
   drift under the event.
2. **tangent** — ICTAC sloped-tangent construction. Sloped baselines
   (not horizontal) are fit on each side of the event; mass loss is
   the difference between pre-baseline at the onset T and post-baseline
   at the endset T, where onset/endset are intersections with the
   inflection tangent. Corrects for slow drift; this is the value to
   report.
3. **dtg_area** — baseline-subtracted DTG integrated over the full
   event extent. Same physical quantity as "tangent", computed in the
   derivative domain. Useful as a self-consistency check.

On clean Fortera samples the three agree within 1–2 %. Larger spread
flags baseline drift, overlap with a neighboring event, or a
poorly-chosen baseline window.

Horizontal-tangent extrapolation (what most instrument software does)
is mathematically equivalent to *naive* subtraction under drift —
report the sloped tangent.

## Limitations

- **Partially overlapping events** (e.g., C-S-H + AFt + AFm below
  200 °C in hydrated cement paste) are detected as a single peak if the
  shoulder is shallower than `--peak-prominence` or closer than
  `--peak-min-sep`. Deconvolution via asymmetric peak fits is a planned
  extension.
- **Heavily overlapping events** with no shoulder cannot be separated
  from single-ramp TGA in principle. Use multi-rate (isoconversional)
  experiments or coupled techniques (TGA-MS, TGA-FTIR).
- The loader expects a single monotonic ramp; segments past the
  temperature max are dropped.

## Programmatic use

The CLI is a thin wrapper. To use from a notebook:

```python
from tga_analysis.io import load_tga
from tga_analysis.smooth import smooth_savgol
from tga_analysis.derivative import dtg_savgol
from tga_analysis.events import detect_events
from tga_analysis.quantify import quantify

run = load_tga("Fortera_Calcite_1.xlsx")
sm = smooth_savgol(run.w, run.dT)
dwdT = dtg_savgol(run.w, run.dT, sm.window_pts, sm.poly_order)
events = detect_events(run.T, sm.w_smooth, dwdT)
for ev in events:
    print(quantify(run.T, sm.w_smooth, dwdT, ev))
```

All module-level functions accept the same tunables as the CLI flags
(keyword arguments, not hyphenated).
