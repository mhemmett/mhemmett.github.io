---
type: project-hub
project: Axial Seamount Shear-Wave Splitting
repo: ~/Seismology/axial-splitting-ml
created: 2026-04-27
updated: 2026-04-27
tags: [axial, shear-wave-splitting, anisotropy, ooi-rca, project-hub]
---

# Axial Splitting — Project Hub

Central note for the shear-wave splitting (SWS) project at Axial Seamount. Source code lives in `~/Seismology/axial-splitting-ml`. This note summarizes scope, methods, code structure, current status, and links to topic notes in this folder.

## Scientific Goal

Use shear-wave splitting on local seismicity to track the temporal evolution of stress and crack-induced anisotropy in Axial Seamount's caldera, with a focus on the **April 2015 eruption** and the broader **2014–2021** inflation/deflation cycles. Splitting parameters — fast polarization direction (φ) and delay time (δt) — encode the orientation and density of aligned cracks in the shear-wave window. Tracking how φ and δt change in space and time should constrain how the volcanic stress field reorganizes around eruptive episodes and feed into volcanic hazard / forecasting work on the Juan de Fuca Ridge.

Key questions:
- Does φ rotate before/during/after the 2015 eruption?
- Do δt magnitudes track inflation cycles (more cracks open → larger δt)?
- Are the anisotropy patterns spatially coherent across the caldera or station-dependent?
- How sensitive are these signals to methodological choices (windowing, clustering, QC)?

## Data

- **Network**: OOI Regional Cabled Array (RCA), University of Washington
- **Period**: 2014–2021
- **Stations** (caldera + flanks): AXAS1, AXAS2, AXCC1, AXEC1, AXEC2, AXEC3, AXID1
- **Instrumentation**: short-period and broadband OBS
- **Catalogs in repo**:
	- `data/AXIAL.PHASE.FINAL_3D_V2.nlloc` — NonLinLoc 3D hypocenter relocations with P/S picks (primary catalog)
	- `data/Axial.MLDD.v202112.2` and `data/Axial.DD.pha.v20221129.1` — machine-learning + double-difference relocated catalogs
	- `data/stations_axial.llz`, `data/AXIAL_stations.xml`
	- Various filtered CSVs (e.g. `final_catalog.csv`, `mldd_catalog_2015_2021.csv`, `extended_catalog_for_snr.csv`)
- **Per-station pickled QC'd catalogs**: `data/AXAS1.clean.cat.pickle`, `AXAS2.clean.cat.pickle`, ... `AXID1.clean.cat.pickle`

## Methods Implemented

The repo carries **three** SWS implementations side-by-side so results can be cross-validated.

### 1. Modified SWSPy (primary)
Silver & Chan (1991) eigenvalue minimization with Teanby et al. (2004) cluster analysis, vendored in `swspy/` with corrections relative to upstream:
- Replaced DBSCAN with **hierarchical Ward clustering** (matches MFAST).
- **Circular-safe coordinate transform** `(δt·cos 2φ, δt·sin 2φ)` for clustering.
- Implements **Teanby representative variance** (Eq. 13–14) for cluster selection.
- Picks best observation inside the chosen cluster by minimum observation variance.
- MFAST-compliant Ncmin=1, Mmax=15.
- Core: `swspy/swspy/splitting/split.py::_sws_win_clustering` (~lines 1029–1150).

### 2. Baillard eigenvalue method (reference / validation)
From Christian Baillard's unpublished Python implementation, retained for cross-checking:
- **Single adaptive window** sized to dominant period × [0.5, 2.0].
- 100 lags × 100 angles exhaustive grid search per event.
- Multiple-minima detection on the eigenvalue surface (morphological filtering).
- RMS-based ranking instead of F-statistic confidence regions.
- P-coda contamination handling on the noise window.
- Lives in `scripts/splitting_functions.py` (~lines 1350–2094) and `scripts/run_sws.py` / `scripts/shearwavesplit.py` (Baillard's originals).

### 3. Standard SWSPy (kept only for reference)
Original upstream behavior — DBSCAN, no Teanby weighting, largest-cluster selection. Not used for production runs; documented in README for traceability of why we forked.

### Key methodological contrasts
| | Modified SWSPy | Baillard |
|---|---|---|
| Window strategy | multiple fixed windows, cluster in window space | single adaptive window, cluster in (φ,δt) space |
| Window length | fixed (≈0.5–2.0 s) | T_dom × [0.5, 2.0] |
| Max lag | fixed (e.g. 20 samples) | T_dom × sampling rate |
| Quality metric | F-statistic 95% region | RMS difference + grid bounds |

## Quality Control Pipeline

Applied uniformly before either splitting method:
1. **SNR** (horizontal) > 2.0; signal window S → S+2.0 s; noise window adjusted to dodge P-coda.
2. **P-wave rectilinearity** > 0.7 over P ± 0.12 s (Jurkevics 1988 covariance).
3. **Incidence angle** < 30° from vertical (from P polarization) — keeps events inside the shear-wave window.
4. **Magnitude** filter (default M > 0).

QC is computed and attached per-event in `organized_waveforms` dicts before splitting is attempted.

## Dynamic Parameter Pipeline

Recent refactor (see `docs/splitting_workflow_update.md` and `docs/splitting_integration_update.md`) replaces hardcoded windowing with per-event, data-driven parameters:

```
estimate_dominant_period(trace)        # PPSD or multitaper → T_dom
calculate_dynamic_parameters(event)    # window = 3–5 × T_dom; band = 1/T_dom × [0.67, 1.5]
create_splitting_analysis(event)       # builds configured SWSPy splitting object
perform_splitting_analysis(event)      # .split() and extract φ, δt, errors, quality
```

Each result row carries the dynamic parameters used (T_dom, window length, filter band), so post-hoc QC and sensitivity analysis are tractable.

Note: deprecated `mtspec` was replaced with the maintained `multitaper` (`MTSpec` class). A `general_stub.py` covers the missing `general` utility module from Baillard's original codebase (`cov_eig`, `smooth_curve`, `ll2xy`).

## Repo Layout (what to look at where)

```
axial-splitting-ml/
├── README.md                                # overview (source for much of this note)
├── data/                                     # catalogs, station metadata, per-station QC pickles
├── docs/
│   ├── splitting_workflow_update.md         # dynamic-parameter integration
│   ├── splitting_integration_update.md      # create_/perform_splitting_analysis split
│   └── package_updates.md                   # mtspec → multitaper, general_stub
├── notebooks/                                # earlier versions of the workflows
├── scripts/                                  # active code — see below
├── swspy/                                    # vendored, modified SWSPy (Teanby clustering fix)
├── results/                                  # CSV outputs, rose plots, batched results
├── pickles/                                  # cached intermediate objects
└── references/                               # papers
```

### Active scripts
- `scripts/splitting_functions.py` — core analysis library (~3500 lines). Both Baillard and SWSPy entry points; QC; dynamic-parameter helpers.
- `scripts/seismic_geometry.py` — back-azimuth, epicentral distance, incidence helpers.
- `scripts/get_all_traces.py` — automated waveform retrieval.
- `scripts/run_sws.py`, `scripts/shearwavesplit.py` — Baillard's reference implementation.

### Per-station notebooks (current production runs)
SWSPy-MFAST + MLDD catalog, one per station, with `_batched` variants for chunked processing:
- `axial_splitting_mldd_AXAS1.ipynb`, `…_AXAS2.ipynb`, `…_AXCC1.ipynb`, `…_AXEC1.ipynb`, `…_AXEC2.ipynb` (+ `_batched`, `_full_catalog`), `…_AXEC3.ipynb`, `…_AXID1.ipynb`
- NonLinLoc-catalog runs: `axial_splitting_nonlinloc_swspy_AX*.ipynb` and `…_swspy-mfast*.ipynb`
- Parameter-sensitivity notebooks: `axial_splitting_nonlinloc_AXAS2_parameter_comparison.ipynb`, `…_AXEC2_parameter_comparison*`

### Plot/figure generators (for poster/paper)
`scripts/ARTICLE_*.py` and `scripts/POSTER_*.py` — rose diagrams, density polar plots, lag vs distance, time variations, earthquake maps, multi-rose, deformation overlays.

## Current Status (April 2026)

- **Active task**: full-catalog SWS run on **AXEC2** using the 2015–2021 MLDD catalog. As of 2026-04-27 the batched output has progressed through `splitting_results_mldd_2015_2021_axec2_batch_41.csv` (~15-minute cadence per batch in `results/`).
- **Recent commits**:
	- `5709d29` Additional parameter test on AXAS2, AXEC3 vs. Baillard rose plots and 24-hr density plots from poster
	- `25696b6` Parameter comparison with adaptation of Baillard plotting interpolation + 24 hrs around eruption
	- `1b9b393` Workflow split into separate scripts for AXAS2, AXEC2, plotting
	- `2a626f9` First version with output across 10k earthquakes (NonLinLoc + AXAS2)
	- `ed13aad` Initial commit
- **Validation status**:
	- Baillard method matches Christian Baillard's reference: SNR < 5% diff, geometry within 1°, φ/δt consistent.
	- Cross-method (modified SWSPy ↔ Baillard) comparison ongoing; main open question is sensitivity to windowing strategy (multi-window cluster vs single adaptive).
- **Known artefacts in repo root** (working figures, not authoritative):
	- `AXAS2-2015-01-22-*.png` — pre-eruption single-event diagnostics
	- `baillard_24hr_scatter*.png`, `baillard_april.png`, `baillard_rose_plot_april.png` — Baillard pipeline outputs
	- `swspy_april.png`, `swspy_rose_plot_april.png` — modified SWSPy comparison
	- `results_swspy_sigma_1_5_0_5_tmid_1_8_2_2_rose_plot.png` — parameter sweep
	- `two_clusters.png` — cluster-selection sanity check

## Open Threads / Next Steps

- Finish AXEC2 MLDD full-catalog batched run; merge batches into a single `axial-mldd-2015-2021-axec2-full.csv`.
- Re-run AXAS1, AXCC1, AXEC1, AXEC3, AXID1 with the dynamic-parameter pipeline so all stations use the same QC and windowing logic.
- Quantify systematic offsets between modified SWSPy and Baillard on the same events (φ, δt, error) to bound method dependence in any temporal-rotation claim.
- Build temporal stacks (rose + density polar) per station around April 2015 with consistent event selection.
- Spatial inversion of φ/δt for stress orientation (handoff to the [[Axial-Stress-Inversion/|stress-inversion]] thread once results are stable).
- Convert results dictionaries to a tidy DataFrame schema across all stations for joint plotting.

## Related Notes

- [[Axial-Stress-Inversion/]] — downstream stress-field inversion using φ
- (add per-station analysis notes here as they're created)

## Citation Pointers

- Hudson, Asplet, Walker (2023), *Seismica* — SWSPy.
- Silver & Chan (1991), *JGR* — eigenvalue method.
- Teanby, Kendall, Jones, Barkved (2004), *GJI* — clustering.
- Wessel, Savage, Teanby (2017) — MFAST reference.
