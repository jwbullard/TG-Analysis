# TG-Analysis

Thermogravimetric (TGA) signal denoising, mass-loss event quantification,
and Frazer-Suzuki peak deconvolution for overlapping events. Reads vendor
Excel/CSV exports (including TA Instruments TRIOS CSV with metadata
preamble), picks an optimal Savitzky-Golay smoother from residual
whiteness, detects mass-loss events in the DTG, and reports mass loss
under five constructions: naive (no correction), ICTAC sloped tangent,
and DTG-area under three under-event baseline choices (constant, linear,
alpha-weighted). Optional `--deconvolve-region` fits a sum of
Frazer-Suzuki peaks to a user-specified T range for partially
overlapping events.

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

**TA Instruments TRIOS CSV exports** (~130 lines of `[Section]` metadata
followed by a data header, then a units row, then data) are detected
automatically: the loader skips the preamble, drops the units row, and
disambiguates the duplicate `Weight` columns by units (the percentage
column is preferred over milligrams).

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
- `<label>_events.csv` — per-event row with peak/onset/endset
  temperatures and five mass-loss values (`naive`, `tangent`, and three
  DTG-area variants: `constant`, `linear`, `alpha`). Each row carries
  `event_kind` ∈ {`detected`, `deconvolved`}; deconvolved rows also
  populate `fwhm_C` and `asymmetry_b`.
- `<label>_decon_R<i>.png` — only when `--deconvolve-region` is used.
  One panel per baseline mode showing the corrected DTG, each fitted
  Frazer-Suzuki peak, and the sum-of-peaks reconstruction. One file per
  region.

Batch-level:

- `summary.csv` — one row per input file (smoothing diagnostics, event
  count, deconvolved-peak count, error if any)
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

**Deconvolution**

| Flag | Default | What it controls |
|---|---|---|
| `--deconvolve-region Tlo,Thi,N` | none | Fit a sum of `N` Frazer-Suzuki peaks to the baseline-corrected DTG over `[Tlo, Thi]` °C. Repeatable for multiple regions. |

The region must include 20 °C of **truly quiet** baseline pad on each
side: the pre- and post-baselines are linear fits to those pads, so if
the pad falls inside a tail of an event the fit is contaminated and the
three baseline modes will scatter wildly. Picking the region wider than
the event extent (with ~20-50 °C of quiet pad) is usually sufficient.

**Output location**

| Flag | Default | What it controls |
|---|---|---|
| `--out` | `<batch-dir>/Results/<batch-stem>/` | Output directory |

## Mass-loss quantification

Each detected event reports five values:

1. **naive** — smoothed-mass difference between event-start and
   event-end T (no baseline correction). Sensitive to ongoing slow
   drift under the event.
2. **tangent** — ICTAC sloped-tangent construction. Sloped baselines
   (not horizontal) are fit on each side of the event; mass loss is
   the difference between pre-baseline at the onset T and post-baseline
   at the endset T, where onset/endset are intersections with the
   inflection tangent. Independent of the under-event baseline choice.
3. **dtg_area_constant** — baseline-subtracted DTG area with the
   under-event baseline held constant at the pre-event slope
   (`dB/dT = a_pre`). "Drift continues forward from what we saw
   before the event."
4. **dtg_area_linear** — baseline-subtracted DTG area with the
   under-event baseline interpolating linearly in T between the pre-
   and post-event slopes (`dB/dT = (1-α_T)·a_pre + α_T·a_post`). For
   the *integrated* mass loss this equals holding the mean slope
   `(a_pre + a_post)/2` constant — only the pointwise baseline shape
   differs.
5. **dtg_area_alpha** — baseline-subtracted DTG area with the slope
   blend weighted by extent of reaction
   (`dB/dT = (1-α(T))·a_pre + α(T)·a_post`,
   `α(T) = (w_start - w(T)) / (w_start - w_end)`). Borchardt-Daniels /
   ICTAC-kinetics convention. Self-consistent when drift behaviour is
   tied to the state of the sample (large mass change altering
   buoyancy, evolved gas, sample geometry).

Horizontal-tangent extrapolation (what most instrument software does)
is mathematically equivalent to *naive* subtraction under drift —
report the sloped tangent instead.

### Method agreement as a diagnostic

On a clean isolated event with quiet baseline on both sides, the three
DTG-area variants agree to within ~1 %. Substantial divergence between
them — or large differences relative to `tangent` — is a *diagnostic*
that the event extent is being contaminated by overlapping chemistry.
In that case no choice of under-event baseline gives the right answer;
the integration limits themselves are wrong. Use `--deconvolve-region`
to fit asymmetric peaks within the region and report per-peak masses.

A demonstrated validation case: in a PC + zeolite ternary blend, the
detected event at 715 °C (well-crystallized calcite decarbonation)
reported 0.6/2.3/2.5 % across the three baseline modes — a factor-of-4
spread. After deconvolving the 540–870 °C region into four
Frazer-Suzuki peaks, the calcite peak gave 4.15/4.13/4.12 %, a 0.7 %
spread invariant to the baseline choice.

## Deconvolution

`--deconvolve-region Tlo,Thi,N` fits a sum of `N` Frazer-Suzuki peaks
to the baseline-corrected DTG over `[Tlo, Thi]`:

    y(T) = A · exp[ -(ln 2 / b²) · ln²(1 + 2b(T - T₀)/W) ]

Parameters per peak: amplitude `A` (% per °C), peak temperature `T₀`
(°C), full width at half maximum `W` (°C; exact in the Gaussian limit
`b → 0`), asymmetry `b` (positive = right tail, negative = left tail —
real thermal-decomposition events typically have left tails).

The fit is repeated under each of the three baseline modes so the
per-peak mass loss can be compared. For a properly-isolated peak the
three modes agree to <1 % per peak.

**Practical guidance for region selection:**

- The region must include 20 °C of **truly quiet** baseline pad on
  each side. If the pad falls inside an event tail, `a_pre` or
  `a_post` is contaminated and the three modes scatter.
- Pick the region wider than the event extent. For a peak with a long
  left tail, extend the lower bound well below the apparent rise.
- Validate by checking that the three baseline modes agree to <1 % per
  peak. If they don't, widen the region.
- The diagnostic PNG (`<label>_decon_R<i>.png`) shows the fit per
  mode — visual inspection is fast and catches most issues
  (e.g., a fit that over-resolved one feature into two narrow peaks).

## Limitations

- **Partially overlapping events** are detected as a single peak if
  the shoulder is shallower than `--peak-prominence` or closer than
  `--peak-min-sep`. Use `--deconvolve-region Tlo,Thi,N` to fit the
  region as a sum of `N` asymmetric peaks.
- **Heavily overlapping events** with no shoulder cannot be separated
  from single-ramp TGA in principle. Use multi-rate (isoconversional)
  experiments or coupled techniques (TGA-MS, TGA-FTIR).
- **Automatic deconvolution region selection** is not implemented;
  `Tlo`, `Thi`, and `N` must be specified explicitly.
- **Number-of-peaks selection** (AIC/BIC) is not implemented; `N` is
  user-fixed. Compare fits at successive `N` to judge whether an
  additional peak is justified.
- The loader expects a single monotonic ramp; segments past the
  temperature max are dropped.

## Programmatic use

The CLI is a thin wrapper. To use from a notebook:

```python
from tga_analysis import (
    load_tga, smooth_savgol, dtg_savgol,
    detect_events, quantify,
    fit_region_baseline, baseline_dtg, deconvolve,
)

run = load_tga("Fortera_Calcite_1.xlsx")
sm = smooth_savgol(run.w, run.dT)
dwdT = dtg_savgol(run.w, run.dT, sm.window_pts, sm.poly_order)

# Five-method mass loss per detected event
events = detect_events(run.T, sm.w_smooth, dwdT)
for ev in events:
    print(quantify(run.T, sm.w_smooth, dwdT, ev))

# Deconvolve a user-specified region
T_lo, T_hi, N = 540.0, 830.0, 1
(a_pre, b_pre), (a_post, b_post) = fit_region_baseline(
    run.T, sm.w_smooth, T_lo, T_hi
)
mask = (run.T >= T_lo) & (run.T <= T_hi)
T_evt, w_evt, dwdT_evt = run.T[mask], sm.w_smooth[mask], dwdT[mask]
for mode in ("constant", "linear", "alpha"):
    B = baseline_dtg(T_evt, w_evt, a_pre, b_pre, a_post, b_post, mode)
    peaks = deconvolve(T_evt, -(dwdT_evt - B), N)
    print(mode, [(p.T0, p.area) for p in peaks])
```

All module-level functions accept the same tunables as the CLI flags
(keyword arguments, not hyphenated).
