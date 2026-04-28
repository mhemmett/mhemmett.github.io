
# Axial Splitting — Project Hub

Central note for the shear-wave splitting (SWS) project at Axial Seamount. Source code lives in `~/Seismology/axial-splitting-ml`. This note tracks scope, methods, code structure, current status, and links out to topic notes and daily logs.

## Scientific Goal

Track the temporal evolution of stress and crack-induced anisotropy in Axial Seamount's caldera through shear-wave splitting on local seismicity. The long-term arc starts **just before the April 2015 eruption** and runs through to the present, currently re-inflated state of stress. Splitting parameters — fast polarization direction (φ) and delay time (δt) — encode the orientation and density of aligned cracks in the shear-wave window. Tracking how these change in space and time should constrain how the volcanic stress field reorganizes around eruptive episodes and feed into forecasting work on the Juan de Fuca Ridge.

Key questions:
- Does φ rotate before/during/after the 2015 eruption?
- Do δt magnitudes track the inflation cycle (more cracks open → larger δt)?
- Are the anisotropy patterns spatially coherent across the caldera or station-dependent?
- How sensitive are these signals to methodological choices (windowing, clustering, QC)?

## Data

- **Network**: OOI Regional Cabled Array (RCA), University of Washington
- **Period**: 2015–2021 (current); long-term goal is pre-2015 → present
- **Instrumentation**: short-period and broadband ocean-bottom seismometers (OBS)

### Stations
Five caldera-floor stations carry the eruption-window analysis:
- **AXAS1, AXAS2** — south caldera flank
- **AXEC1, AXEC2, AXEC3** — eastern caldera floor

**AXCC1** (central caldera) tipped over in **March 2015** due to caldera volumetric inflation; it'll be folded back in once orientation/tilt corrections are settled. **AXID1** sits on the cabled array but isn't part of the eruption-window five.

Future additions: the **15-station 2022–2024 temporary OBS deployment along the North Rift Zone**.

### Catalogs

The **production catalog** is the **MLdd (machine-learning + double-difference) catalog from Kaiwen Wang's** work at Axial Seamount (`data/Axial.MLDD.v202112.2`, parsed into `data/mldd_catalog_2015_2021.csv`). **S-pick uncertainties from Wang (2024)** drive the per-event windowing.

Method development and ongoing **cross-validation** use a **NonLinLoc 3D relocated catalog with kurtosis-based picks from Wilcock & Baillard** (`data/AXIAL.PHASE.FINAL_3D_V2.nlloc`). Results between the two pipelines are consistent — the NonLinLoc catalog now serves mostly as a regression check.

| Catalog | Role | File |
|---|---|---|
| Wang MLdd, 2015–2021 | **primary** | `data/Axial.MLDD.v202112.2` → `mldd_catalog_2015_2021.csv` |
| Wilcock/Baillard NonLinLoc | method dev + cross-validation | `data/AXIAL.PHASE.FINAL_3D_V2.nlloc` |

Station metadata in `data/stations_axial.llz` and `data/AXIAL_stations.xml`. Per-station QC'd catalogs are pickled to `data/AX*.clean.cat.pickle`.

## Method

Production analysis uses a **modified SWSPy** (Silver & Chan 1991 eigenvalue minimization with Teanby et al. 2004 cluster analysis) with several MFAST-style adjustments and per-event dynamic parameters. **Christian Baillard's** original single-window Python implementation is the historical foundation of this work — many of his supporting scripts are still active dependencies — and provides the cross-validation reference.

### Modified SWSPy (production)

Vendored in `swspy/` with the following corrections relative to upstream:
- **Hierarchical Ward clustering** (replacing upstream DBSCAN) in a circular-safe `(δt cos 2φ, δt sin 2φ)` coordinate space — matches MFAST.
- **Teanby (2004) representative variance** (Eq. 13–14) for cluster selection; best observation within cluster is the one with minimum observation variance.
- **Adapted windowing and clustering parameters** relative to upstream.
- **Dynamic windowing on `T_dom`** of each event (similar to MFAST), rather than fixed window lengths.
- **Wang (2024) S-pick uncertainty** drives the analysis-window placement around the S arrival.

Core: `swspy/swspy/splitting/split.py::_sws_win_clustering`, plus the `splitting_functions.py` workflow in `scripts/`.

### Christian Baillard's single-window method (foundation + validation)

Much of the early method work, and a large fraction of the supporting scripts in this repo, draws on Baillard's unpublished Python implementation:
- **Single adaptive window**, sized to `[0.02, T_dom × 2.0]`.
- 100 lags × 100 angles exhaustive grid search per event.
- Multiple-minima detection on the eigenvalue surface (morphological filtering).
- RMS-based ranking of solutions; P-coda contamination handling.

His scripts (`shearwavesplit.py`, `sws_methods.py`, `plotwaveform.py`, `GMT.py`, `projection.py`, `util.py`, `sws_vertices.py`) are still active dependencies — `splitting_functions.py` imports them transitively. Cross-validation against his reference results: SNR within ~5%, geometry within 1°, φ/δt consistent.

### Method comparison

| Feature             | Modified SWSPy (production)               | Baillard (validation)  |
| ------------------- | ----------------------------------------- | ---------------------- |
| Window strategy     | multiple windows, cluster in window space | single adaptive window |
| Window length       | dynamic, based on `T_dom` (MFAST-style)   | `[0.02, T_dom × 2.0]`  |
| S-arrival placement | Wang (2024) S-pick uncertainty            | Baillard adaptive      |
| Clustering          | Ward, circular-safe `(δt, φ)`             | parameter-space minima |
| Quality metric      | F-statistic 95% region                    | RMS difference         |

### Quality control

Applied uniformly before splitting:
1. **SNR** (horizontal) > 2.0; signal window S → S+2.0 s; noise window dodges P-coda.
2. **P-wave rectilinearity** > 0.7 over P ± 0.12 s (Jurkevics 1988).
3. **Incidence angle** < 30° from vertical (from P polarization).
4. **Magnitude** filter (default M > 0).

QC is computed and attached per-event in `organized_waveforms` dicts before splitting.

## Repo Layout (post-cleanup)

```
axial-splitting-ml/
├── README.md                      # main doc — kept in sync with this note
├── env.sh                         # one-shot venv + deps + Jupyter kernel
├── requirements.txt
├── .gitignore                     # ignores data/, results/, pickles/, *.mseed, figures
├── cleanup.sh                     # 2026-04-27 repo-cleanup helper
├── data/                          # catalogs, station metadata, per-station QC pickles  (untracked)
├── pickles/                       # cached intermediate objects                          (untracked)
├── results/                       # CSV outputs, rose plots, batched results             (untracked)
├── docs/
│   ├── splitting_workflow_update.md
│   ├── splitting_integration_update.md
│   └── package_updates.md
├── swspy/                         # vendored, modified SWSPy submodule
└── scripts/                       # active code — see below
```

### Active scripts
- `splitting_functions.py` — main workflow: QC, dynamic-parameter pipeline, both Baillard and SWSPy entry points.
- `teanby_clustering.py` — clustering helpers used by `splitting_functions.py`.
- `get_all_traces.py` — automated waveform retrieval.
- Baillard-derived modules still in active use: `sws_methods.py`, `shearwavesplit.py`, `plotwaveform.py`, `GMT.py`, `projection.py`, `util.py`, `sws_vertices.py`.

### Active notebooks

Per-station MLdd processing (production):
- `axial_splitting_mldd_AXAS1.ipynb`, `…_AXAS2.ipynb`, `…_AXCC1.ipynb`, `…_AXEC1.ipynb`, `…_AXEC2.ipynb` (+ `_batched`, `_full_catalog`), `…_AXEC3.ipynb`, `…_AXID1.ipynb`.

Per-station NonLinLoc processing (cross-validation):
- `axial_splitting_nonlinloc_swspy_AX*.ipynb` (+ `_batched`).
- `axial_splitting_nonlinloc_swspy-mfast*.ipynb` for MFAST-style parameter exploration.
- `axial_splitting_nonlinloc_AXAS2_parameter_comparison.ipynb`, `…_AXEC2_parameter_comparison.ipynb` for sensitivity sweeps.

Plotting and catalog prep:
- `nonlinloc_apr_14_jun_01_plots.ipynb` — current plotter for combined cleaned per-station results.
- `mldd_apr_20_28_plots.ipynb` — MLdd 24-hr-window plots around the eruption.
- `parse_mldd_catalog.ipynb` — converts Wang's MLdd file into the workflow's catalog format.
- `tmid_analysis.ipynb` — sensitivity to analysis-window center.

### Files kept locally but not tracked
- All of Baillard's plotting / figure framework: `scripts/ARTICLE_*.py`, `scripts/POSTER_*.py`, `scripts/test_*.py`, `scripts/untitled*.py`. Reference only.
- Orphan utilities and superseded notebooks (`master_shear_wave_splitting_workflow.ipynb`, the `notebooks/` folder, etc.) — see `[[04-27-26 Notes]]` for the full list.
- All data, results, pickles, waveforms (`*.mseed`), and figures (`*.png/pdf/jpg`) — see `.gitignore`.

## Current Status (April 2026)

- **Active production run**: AXEC2 across the full 2015–2021 MLdd catalog, batched, currently in the high-30s of ~100-event batches (`results/splitting_results_mldd_2015_2021_axec2_batch_*.csv`).
- **Validation status**: modified SWSPy ↔ Baillard agreement holds — SNR within ~5%, geometry within 1°, φ/δt consistent across the NonLinLoc cross-check catalog.
- **Recent commits**:
	- `5709d29` Additional parameter test on AXAS2, AXEC3 vs. Baillard rose plots and 24-hr density plots from poster
	- `25696b6` Parameter comparison with adaptation of Baillard plotting interpolation + 24 hrs around eruption
	- `1b9b393` Workflow split into separate scripts for AXAS2, AXEC2, plotting
	- `2a626f9` First version with output across 10k earthquakes (NonLinLoc + AXAS2)
	- `ed13aad` Initial commit

## Open Threads / Next Steps

- Finish AXEC2 MLDD full-catalog batched run; merge batches into a single `axial-mldd-2015-2021-axec2-full.csv`.
- Re-run AXAS1, AXAS2, AXEC1, AXEC3 with the same dynamic-parameter / Wang-S-pick pipeline.
- Fold **AXCC1** in once tilt-correction (post-March-2015 tip-over) is settled.
- Bring in the 2022–2024 **North Rift Zone** OBS deployment (15 stations) once the cabled-array workflow is stable across the existing five.
- Quantify systematic offsets between modified SWSPy and Baillard on the same events (φ, δt, error) to bound method dependence in any temporal-rotation claim.
- Build temporal stacks (rose + density polar) per station around April 2015 with consistent event selection.
- Spatial inversion of φ/δt for stress orientation — handoff to the [[Axial-Stress-Inversion/|stress-inversion]] thread once results are stable.
- Convert results dictionaries to a tidy DataFrame schema across all stations for joint plotting.

## Daily Notes

- [[04-27-26 Notes]] — repo cleanup: `.gitignore`, untracked legacy/orphan files, tracked active MLDD + NonLinLoc + plotting notebooks.

## Related Notes

- [[Axial-Stress-Inversion/]] — downstream stress-field inversion using φ.
- (add per-station analysis notes here as they're created)

## Citation Pointers

- **Wang, K. (2024)** — MLdd catalog and S-pick uncertainty model (production input).
- Wilcock, W.S.D. & Baillard, C. — NonLinLoc / kurtosis-pick catalog (method dev + cross-validation).
- Hudson, T.S., Asplet, J., Walker, A.M. (2023). *Seismica* — SWSPy.
- Silver, P.G., Chan, W.W. (1991). *JGR*, 96(B10), 16429–16454 — eigenvalue method.
- Teanby, N.A., Kendall, J.M., Jones, R.H., Barkved, O. (2004). *GJI*, 156(3), 459–466 — clustering.
- Wessel, A., Savage, M., Teanby, N. (2017). MFAST — DOI: 10.5281/zenodo.1042760.
- Jurkevics, A. (1988). *BSSA*, 78(5), 1725–1743 — polarization analysis.

## Acknowledgements

This work builds substantially on **Christian Baillard's** earlier shear-wave-splitting code at Axial Seamount; many supporting scripts in `scripts/` are his and are still in use. The **MLdd catalog and S-pick uncertainty model** are from **Kaiwen Wang's** work on Axial seismicity. NonLinLoc relocations and kurtosis picks are from **William Wilcock and Christian Baillard**.
