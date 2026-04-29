# Joint Stress Inversion at Axial Seamount: Integrated Project Plan

**Michael Hemmett — Working Document, April 2026**

---

## 1. Project Overview

This project develops a 4D model of the stress field and seismic anisotropy at Axial Seamount across a complete eruption cycle, using a joint inversion of three seismological observables — shear-wave splitting (SWS), relative seismic velocity variations (dv/v), and P-wave focal mechanisms — constrained by geodetic strain measurements. The inversion is unified through the crack density tensor α_ij, which serves as the common physical quantity linking all observables to the elastic stiffness tensor C_ijkl and ultimately to the stress tensor σ_ij.

### 1.1 Core Theoretical Insight

All three observables are projections of the elastic stiffness tensor C_ijkl, which is itself a function of the crack state of the medium. The crack state is controlled by the stress tensor σ_ij and pore pressure P_p through the mechanics of stress-dependent crack opening and closure (Nur, 1971; Hudson, 1981; Sayers and Kachanov, 1995). The forward model is:

```
σ_ij(x,t), P_p(x,t)  →  α_ij(x,t)  →  C_ijkl(x,t)  →  {SWS, dv/v, focal mechanisms}
```

The inverse problem recovers α_ij(x,t) jointly from all three data types, then converts to C_ijkl via Hudson's perturbation theory, and finally to σ_ij using geodetic strain and orientation constraints from focal mechanisms and SWS.

### 1.2 What is Novel

No published work jointly inverts all three observables simultaneously for a crack density field. Existing studies treat these measurements independently (e.g., Gerst and Savage, 2004 for SWS; Brenguier et al., 2008 for dv/v; Martínez-Garzón et al., 2014 for focal mechanisms) or compare pairs qualitatively (Boness and Zoback, 2006 for SWS vs. focal mechanisms). The specific innovations here are: (1) using the crack density tensor α_ij as the primary inversion target that both dv/v and splitting intensity constrain simultaneously; (2) converting the recovered C_ijkl to absolute stress using geodetic strain measurements; and (3) applying this framework to 10 years of continuous data at the best-instrumented submarine volcano in the world.

---

## 2. Data

All data come from the Ocean Observatories Initiative (OOI) Regional Cabled Array at Axial Seamount (Wilcock et al., 2018) and associated instruments.

**Seismological data:**

The earthquake catalog is a machine learning-based detection catalog relocated with the double-difference algorithm (Wang et al., 2024), spanning 2015–2025. This catalog provides the earthquake locations and waveforms needed for all three seismological analyses. The key derived products are: shear-wave splitting parameters (φ, δt) and splitting intensity from S-wave arrivals (Kendall et al., 2025); P-wave focal mechanisms from first-motion polarities and/or waveform inversion (Zhang et al., in prep); and dv/v time series from ambient noise cross-correlations between station pairs (Lee et al., 2024; Wang et al., in prep).

**Geodetic data:**

Seafloor geodetic measurements include vertical displacement from bottom pressure recorders (BPRs), which record the caldera inflation/deflation with sub-centimeter precision (Nooner and Chadwick, 2016), and horizontal displacement from acoustic ranging between seafloor benchmarks. Together these constrain the surface deformation field, from which strain at depth can be modeled.

**Velocity model:**

A 3D P- and S-wave velocity model is required for ray tracing (SWS forward modeling), sensitivity kernel computation (dv/v), and earthquake location. Existing models from OBS tomography at Axial (e.g., West et al., 2001; Baillard et al.) provide the starting structure, supplemented by any updated models from the OOI data (Hefner et al., 2020; Slead et al., 2024).

---

## 3. Theoretical Framework

### 3.1 The Crack Density Tensor α_ij as the Unifying Variable

The central inversion target is the second-rank crack density tensor α_ij(x,t) (Sayers and Kachanov, 1995), defined as:

```
α_ij = (1/V) Σ_r B_T^(r) a_r³ n_i^(r) n_j^(r)
```

where the sum is over all open cracks in volume V, B_T is the tangential compliance of each crack, a_r is the crack radius, and n^(r) is the crack normal. This symmetric tensor encodes both the scalar crack density (through its trace, tr(α) ∝ ε) and the crack orientation distribution (through its eigenvectors and eigenvalue ratios).

Additionally, the fourth-rank crack density tensor β_ijkl captures deviations from elliptical anisotropy:

```
β_ijkl = (1/V) Σ_r (B_N^(r) − B_T^(r)) a_r³ n_i^(r) n_j^(r) n_k^(r) n_l^(r)
```

For fluid-saturated cracks at Axial Seamount, B_N → 0 (fluid stiffens the normal compliance), so β_ijkl → 0 and the anisotropy is controlled primarily by α_ij. This simplification reduces the problem significantly while remaining physically appropriate for the submarine volcanic setting.

**Why α_ij rather than scalar ε or σ_ij:** Scalar crack density ε discards orientation information that splitting intensity contains. The full stress tensor σ_ij is further from the observables and requires additional modeling assumptions to predict splitting and velocity. The tensor α_ij sits at the natural intermediate level — close enough to the observables to enable a nearly linear inversion, while retaining the directional information that makes the joint approach powerful.

### 3.2 From α_ij to C_ijkl: Hudson's Perturbation Theory

The elastic stiffness tensor is constructed from the isotropic background plus crack-induced perturbations (Hudson, 1981):

```
C_ijkl = C⁰_ijkl + C¹_ijkl(α) + C²_ijkl(α) + ...
```

where C⁰ is the isotropic (crack-free) background stiffness from the reference velocity model, and C¹ and C² are the first- and second-order corrections. In the Sayers-Kachanov formulation, these corrections enter through the excess compliance:

```
ΔS_ijkl = (1/4)(α_ik δ_jl + α_il δ_jk + α_jk δ_il + α_jl δ_ik) + β_ijkl
```

so that C = (S⁰ + ΔS)⁻¹. For weak anisotropy, C ≈ C⁰ − C⁰ ΔS C⁰ to first order. At crack densities above ε ~ 0.05–0.1 (which the data suggest are reached co-eruptively at Axial), the second-order term matters and should be retained.

### 3.3 Observable 1: dv/v → Isotropic Projection of α_ij

The dv/v measured from ambient noise interferometry (Shapiro and Campillo, 2004; Sens-Schönfelder and Wegler, 2006) responds primarily to the isotropic average of the stiffness perturbation, which is proportional to the trace of the crack density tensor:

```
dv/v ∝ −f(ν) × tr(α)
```

where f(ν) is a function of Poisson's ratio that depends on whether P- or S-wave velocities dominate the measurement (for Rayleigh-wave-dominated noise cross-correlations, this is a combination of both). The proportionality follows from the O'Connell and Budiansky (1974) effective medium expressions:

```
K_eff/K₀ = 1 − (16/9)[(1 − ν²)/(1 − 2ν)] tr(α)/B_T
μ_eff/μ₀ = 1 − (32/45)[(1 − ν)(5 − ν)/(2 − ν)] tr(α)/B_T
```

so that changes in tr(α) map linearly onto changes in elastic moduli and hence velocities.

The observed dv/v for a station pair is a spatial integral weighted by the sensitivity kernel K:

```
(dv/v)_obs = ∫∫∫ K(x; freq, station pair) × [−f(ν) Δtr(α(x,t))] d³x
```

The sensitivity kernel K depends on the frequency band (higher frequency → shallower sensitivity) and the station pair geometry. Frequency-dependent dv/v measurements provide crude depth resolution (Brenguier et al., 2014).

**What dv/v constrains:** tr(α), the total open crack density summed over all orientations. This is the isotropic part of the crack density tensor. It does not constrain the anisotropic part (eigenvalue differences and eigenvector orientations of α_ij).

### 3.4 Observable 2: Splitting Intensity → Anisotropic Projection of α_ij

The splitting intensity SI (Chevrot, 2000) for a shear wave propagating through a weakly anisotropic medium is a linear path integral of the local anisotropy:

```
SI_k = ∫_ray K_SI(ℓ; ray geometry) × δC_ijkl(ℓ) dℓ
```

where δC_ijkl = C¹_ijkl(α) is the first-order Hudson correction, which is a linear function of α_ij. Since SI ≈ δt sin(2(θ − φ)), the splitting intensity is sensitive to the anisotropic part of the stiffness perturbation — the eigenvalue difference and orientation of α_ij in the plane perpendicular to the ray.

For a 3D discretized model, this becomes:

```
SI_k = Σ_j G_kj^(SI) m_j
```

where m_j parameterizes the anisotropy in each grid cell (e.g., as γ sin(2ψ) and γ cos(2ψ), where γ is the Thomsen shear-wave anisotropy parameter and ψ is the fast direction azimuth). The sensitivity matrix G^(SI) is determined by the ray path geometry through each cell.

**What SI constrains:** The anisotropic part of α_ij — specifically, the eigenvalue difference (related to percent anisotropy / differential crack density) and the eigenvector orientation (related to stress orientation). This is complementary to dv/v, which constrains the isotropic part.

### 3.5 Observable 3: Focal Mechanisms → Stress Orientation and R

Focal mechanism stress inversion (Michael, 1984, 1987; Hardebeck and Michael, 2006) constrains the deviatoric stress tensor orientation and the stress ratio R = (σ₁ − σ₂)/(σ₁ − σ₃) via the Wallace-Bott hypothesis (Wallace, 1951; Bott, 1959): slip occurs in the direction of maximum resolved shear stress on the fault plane.

The stress tensor recovered from focal mechanisms provides:

- Three Euler angles specifying the orientations of σ₁, σ₂, σ₃
- The stress ratio R (shape of the deviatoric tensor)
- Temporal variations in these quantities (with ≥15–20 mechanisms per time window)

These are **constraints on the eigenvectors and eigenvalue ratios of σ_ij**, which, through the stress–crack relationship, are constraints on the eigenvectors and eigenvalue ratios of α_ij. The physical link is that the principal directions of α_ij are determined by the principal directions of σ_ij (cracks open normal to σ₃, close normal to σ₁; Crampin, 1981), and the eigenvalue ratios of α_ij reflect R.

**What focal mechanisms constrain:** The orientation of α_ij eigenvectors and the ratio of its eigenvalues. They do NOT constrain the absolute magnitude of α_ij or its trace.

### 3.6 Complementarity Summary

```
                        dv/v            Splitting Intensity      Focal Mechanisms
─────────────────────────────────────────────────────────────────────────────────
tr(α)  [total density]    ●●●              ○                       ○
Δα eigenvalues            ○                ●●●                     ●●
  [anisotropy magnitude]
α eigenvectors            ○                ●●●                     ●●●
  [orientation]
Temporal resolution       ●●● (daily)      ●● (per event)          ● (per window)
Depth resolution          ●● (freq-dep.)   ●● (ray geometry)       ● (event depth)
Spatial coverage          ●● (interstation) ●●● (ray paths)        ●● (seismogenic)
```

Together, the three observables constrain all components of the symmetric tensor α_ij: its trace (total crack density), its deviatoric part (anisotropy), and its orientation (stress direction).

---

## 4. Inversion Methodology

### 4.1 Model Parameterization

The model domain is the Axial Seamount caldera and its immediate surroundings, discretized into a 3D grid of cells P(x,y,z). At each cell and each time step t, the model parameters are:

**Primary unknowns (what the seismological inversion recovers):**

At each cell, α_ij is a symmetric 3×3 tensor with 6 independent components. For practical parameterization, decompose it into:

- tr(α) = α₁₁ + α₂₂ + α₃₃ (scalar crack density, 1 parameter — constrained primarily by dv/v)
- Anisotropic part: eigenvalue differences and orientation (up to 5 additional parameters, but for a medium with predominantly sub-vertical cracks in the shallow crust, this can often be reduced to 2–3: horizontal anisotropy magnitude, fast azimuth ψ, and tilt — constrained primarily by SI and focal mechanisms)

**Structural parameters (determined separately or fixed):**

- C⁰_ijkl(x): isotropic background stiffness from the 3D velocity model
- B_T, B_N: crack compliances (determined by assuming fluid-saturated penny-shaped cracks in basalt, following Hudson, 1981)
- K(x; freq, pair): dv/v sensitivity kernels (computed from the velocity model)
- G^(SI): splitting intensity sensitivity matrix (from ray tracing through the velocity model)

### 4.2 Joint Inversion: Stacked Linear System

Rather than solving each data type independently and then combining (weighted stacking), formulate a single linear system that simultaneously fits all data:

```
┌              ┐       ┌           ┐
│  G_dv/v      │       │  d_dv/v   │
│              │       │           │
│  G_SI        │ × m = │  d_SI     │
│              │       │           │
│  G_FM        │       │  d_FM     │
│              │       │           │
│  λ_s L_s     │       │  0        │
│  λ_t L_t     │       │  0        │
└              ┘       └           ┘
```

where:

**d_dv/v** is the vector of observed dv/v values (one per station pair per time window per frequency band), and **G_dv/v** maps the model parameters m (the α_ij field) to predicted dv/v via the sensitivity kernels and the O'Connell-Budiansky relations. G_dv/v primarily senses tr(α).

**d_SI** is the vector of observed splitting intensities (one per earthquake per station), and **G_SI** maps m to predicted SI via ray-traced sensitivity kernels through the anisotropic model. G_SI primarily senses the anisotropic components of α_ij.

**d_FM** encodes the focal mechanism constraints. For each time-space window where a stress inversion has been performed (Michael, 1984; Arnold and Townend, 2007), the recovered stress orientations and R provide constraints on the eigenvectors and eigenvalue ratios of α_ij in those cells. These enter as penalty terms: G_FM × m penalizes models where the α_ij eigenvectors deviate from the focal-mechanism-derived stress directions, and where the eigenvalue ratios are inconsistent with R. This is the coupling between the orientation information from focal mechanisms and the α_ij field.

**L_s and L_t** are spatial and temporal regularization operators (e.g., Laplacian smoothing in space, first-difference smoothing in time), with damping parameters λ_s and λ_t chosen by L-curve analysis or cross-validation. These enforce the physical expectation that the crack density field varies smoothly except across known structural boundaries (ring faults, caldera walls).

**m** is the model vector containing the α_ij parameters at every grid cell and time step.

This combined system is solved by damped least squares (LSQR or conjugate gradient), yielding the α_ij(x,t) field that jointly explains all three data types.

This formulation follows the standard approach for multi-data-type geophysical inversions (cf. Julià et al., 2000 for joint surface wave + receiver function inversion; Hardebeck and Michael, 2006 for spatially damped stress inversion).

### 4.3 From α_ij to C_ijkl

At each grid cell, the recovered α_ij is inserted into the Sayers-Kachanov compliance perturbation:

```
ΔS_ijkl = (1/4)(α_ik δ_jl + α_il δ_jk + α_jk δ_il + α_jl δ_ik)
```

(with β_ijkl ≈ 0 for fluid-saturated cracks), and the full stiffness tensor is:

```
C_ijkl(x,t) = [S⁰_ijkl(x) + ΔS_ijkl(x,t)]⁻¹
```

retaining second-order terms where crack densities are large (ε > 0.05). This gives C_ijkl at every point in the model at every time step.

### 4.4 From C_ijkl to σ_ij: The Geodetic Strain Bridge

This is the step that converts the seismologically-determined C_ijkl into absolute stress in physical units. From Hooke's law:

```
σ_ij(x,t) = C_ijkl(x,t) × ε_kl(x,t)
```

where ε_kl is the strain tensor. The strain field is derived from geodetic observations:

**Surface strain:** The BPR-measured vertical displacement and acoustic-ranging horizontal displacements directly constrain the surface strain tensor at the seafloor. At the instrument locations, the strain is measured, not modeled.

**Strain at depth:** Extending strain to depth requires a source model for the deformation. Candidate models include: a Mogi point source (simplest, spherically symmetric pressure source; Mogi, 1958), a horizontal penny-shaped sill (more realistic for Axial's shallow magma reservoir; Fialko et al., 2001), or a finite element model incorporating the caldera geometry, ring faults, and heterogeneous elastic properties from the velocity model (cf. Masterlark, 2007).

The source model dependence is a systematic uncertainty. To quantify it, run the σ_ij calculation with multiple plausible source geometries and report the spread as a modeling uncertainty on the absolute stress magnitudes. The stress *orientations* (from focal mechanisms and SWS) and the relative stress *changes* (from dv/v and δt temporal variations) are much less sensitive to the geodetic source model than the absolute magnitudes.

**Empirical calibration as a cross-check:** The co-eruption deflation event (2015) provides an independent calibration opportunity. The BPR-measured 4 m of subsidence, combined with the observed Δ(δt) ≈ 0.05–0.10 s, gives an empirically calibrated splitting-stress sensitivity of S ≈ 0.002–0.004 s/MPa. This can be used to check whether the model-predicted stress changes at the surface are consistent with the observed splitting changes, independent of the source model.

### 4.5 Pore Pressure Separation

The dv/v responds to effective stress σ_ij^eff = σ_ij − ηP_p δ_ij (where η is the Biot coefficient), so a change in dv/v could reflect either a stress change or a pore pressure change. To separate these:

**From the joint inversion itself:** If dv/v changes (which sense tr(α), the isotropic part) occur without corresponding changes in SI (which senses the anisotropic part) or focal mechanism orientations, the most parsimonious explanation is an isotropic effective stress change — i.e., a pore pressure change. Conversely, if dv/v, SI, and focal mechanisms all change consistently, the dominant cause is deviatoric stress evolution.

**As an additional model parameter:** Include P_p(x,t) as an inversion parameter at each grid cell. The dv/v is sensitive to the combination (tr(α) + pressure effect), while the SI is sensitive only to the anisotropic part of α_ij. The pore pressure enters only through the isotropic part, so the anisotropic observables (SI, focal mechanisms) break the trade-off.

**From hydrothermal modeling:** Independent constraints on P_p may come from diffuse flow measurements, hydrothermal vent temperature data, or coupled thermo-hydro-mechanical models of the Axial hydrothermal system.

---

## 5. Implementation Plan

### Phase 1: Individual Observable Pipelines (Prerequisite)

These analyses are underway or completed by the research group:

- **SWS measurement:** Measure (φ, δt) for all events in the Wang et al. (2024) catalog. Compute splitting intensity SI for each event-station pair. Assess temporal changes at each station. (Hemmett; Kendall et al., 2025)
- **dv/v monitoring:** Compute noise cross-correlations for all station pairs. Extract dv/v time series at multiple frequency bands (0.5–1 Hz, 1–2 Hz, 2–4 Hz) for depth resolution. (Lee et al., 2024; Wang et al., in prep)
- **Focal mechanism catalog:** Determine P-wave first-motion focal mechanisms. Perform stress inversion in spatial and temporal bins. (Zhang et al., in prep)
- **Velocity model:** Assemble or update the 3D V_p, V_s model for Axial from existing tomographic studies.

### Phase 2: Sensitivity Kernel Computation

- **dv/v kernels:** Compute 2D (and ideally 3D) sensitivity kernels for each station pair at each frequency band, using the velocity model. Methods: Born approximation for surface waves (Tromp et al., 2010) or adjoint methods if feasible.
- **SI kernels:** Ray-trace S-wave paths from each event to each station through the 3D velocity model. Compute the splitting intensity sensitivity kernel along each ray (Chevrot, 2000; Favier and Chevrot, 2003). For finite-frequency effects, use the banana-doughnut kernels of Favier and Chevrot (2003).
- **Resolution analysis:** At each grid cell, assess the resolution from each data type (dv/v, SI, focal mechanisms). Cells with poor resolution from all data types will be heavily damped in the inversion.

### Phase 3: Joint Inversion for α_ij(x,t)

- **Grid design:** Define a 3D grid covering the caldera with cell sizes appropriate to the resolution (~500 m horizontal, ~500 m vertical, based on ray coverage). Time stepping at ~monthly intervals for the inter-eruption period, sub-daily during eruption sequences.
- **Build the combined G matrix:** Assemble G_dv/v, G_SI, and G_FM into the stacked system (Section 4.2).
- **Solve:** Use LSQR with spatial and temporal damping. Test sensitivity to damping parameters via L-curve analysis.
- **Output:** α_ij(x,y,z,t) at every grid cell and time step.

### Phase 4: Forward Modeling of C_ijkl and σ_ij

- **C_ijkl computation:** At each cell, compute C_ijkl from α_ij using the Hudson/Sayers-Kachanov perturbation (Section 4.3).
- **Geodetic strain modeling:** Invert the BPR and acoustic ranging data for a source model at each time step. Compute the 3D strain field ε_kl(x,t) from the source model. Test multiple source geometries (Mogi, sill, FEM) and propagate uncertainty.
- **Stress computation:** Apply σ_ij = C_ijkl × ε_kl at each cell. Compare the resulting stress orientations with the independently determined focal mechanism orientations as a consistency check.
- **Empirical calibration check:** Compare the model-predicted co-eruption stress change with the empirically calibrated value from the δt–BPR comparison (~20–30 MPa).

### Phase 5: Interpretation and Validation

- **Eruption cycle dynamics:** Map the evolution of σ_ij(x,t) across the full inflation–eruption–reinflation cycle. Identify spatial patterns of stress concentration (ring faults, reservoir roof, rift zones).
- **Pre-eruption precursors:** Analyze the hours-to-days before the 2015 eruption onset for precursory signals in the α_ij field, particularly at the eastern caldera stations where splitting changes are observed before mixed-frequency seismicity.
- **Validation against independent data:** Check consistency with: (a) bottom pressure recorder time series, (b) earthquake rate changes, (c) known lava flow locations and volumes, (d) hydrothermal vent observations.
- **Sensitivity analysis:** Quantify how results depend on: the choice of velocity model, the geodetic source geometry, the assumed crack aspect ratio distribution, the damping parameters, and the frequency bands used for dv/v.

---

## 6. Key Assumptions and Limitations

**Assumptions that are well-justified at Axial:**

- Crack-induced anisotropy dominates over intrinsic mineral anisotropy in the shallow crust (supported by the temporal variability of splitting, which rules out a static structural origin).
- Cracks are fluid-saturated (the submarine setting ensures pervasive seawater saturation), so β_ijkl ≈ 0 and α_ij controls the anisotropy.
- The Hudson perturbation theory is valid (crack densities ε < 0.15, aspect ratios ≪ 1).
- The velocity model is sufficiently accurate for ray tracing and kernel computation (verifiable against travel time residuals).

**Assumptions that introduce uncertainty:**

- The geodetic source model for strain at depth (Mogi vs. sill vs. FEM) — mitigated by testing multiple geometries.
- Uniform crack compliance B_T (in reality, crack compliance varies with depth, temperature, and rock type) — mitigated by depth-dependent parameterization if needed.
- The crack closure relationship is approximately linear over the relevant stress range — justified for Δσ ≪ πEα_max/4(1−ν²) ≈ 100–500 MPa, which is satisfied at Axial.

**Known limitations:**

- Temporal resolution varies by data type: dv/v is daily, SI is per-earthquake (irregular in time), focal mechanisms require ≥15 events per window. The joint inversion naturally handles this through the time-stepping scheme, but the effective temporal resolution at any cell is limited by the sparsest data type with meaningful sensitivity there.
- Spatial resolution is limited by station geometry (7 OOI broadband seismometers) and earthquake distribution. The caldera interior is well-sampled; the distal rift zones are not.
- The inversion recovers effective stress (σ − ηP_p), not total stress and pore pressure independently, unless the pore pressure separation scheme (Section 4.5) provides sufficient constraint.

---

## 7. References

Angelier, J. (1979). Determination of the mean principal directions of stresses for a given fault population. *Tectonophysics*, 56(3-4), T17-T26.

Arnold, R. and Townend, J. (2007). A Bayesian approach to estimating tectonic stress from seismological data. *Geophys. J. Int.*, 170, 1336-1356.

Boness, N.L. and Zoback, M.D. (2006). Mapping stress and structurally controlled crustal shear velocity anisotropy in California. *Geology*, 34(10), 825-828.

Bott, M.H.P. (1959). The mechanics of oblique slip faulting. *Geological Magazine*, 96(2), 109-117.

Brenguier, F. et al. (2008). Towards forecasting volcanic eruptions using seismic noise. *Nature Geoscience*, 1, 126-130.

Brenguier, F. et al. (2014). Mapping pressurized volcanic fluids from induced crustal seismic velocity drops. *Science*, 345, 80-82.

Chevrot, S. (2000). Multichannel analysis of shear wave splitting. *J. Geophys. Res.*, 105(B9), 21579-21590. DOI: 10.1029/2000JB900199.

Crampin, S. (1981). A review of wave motion in anisotropic and cracked elastic-media. *Wave Motion*, 3, 343-391.

Donaldson, C. et al. (2017). Relative seismic velocity variations correlate with deformation at Kīlauea volcano. *Sci. Adv.*, 3, e1700219. DOI: 10.1126/sciadv.1700219.

Favier, N. and Chevrot, S. (2003). Sensitivity kernels for shear wave splitting in transverse isotropic media. *Geophys. J. Int.*, 153, 213-228.

Fialko, Y. et al. (2001). Deformation due to a pressurized horizontal circular crack in an elastic half-space. *Geophys. J. Int.*, 146, 181-190.

Gerst, A. and Savage, M.K. (2004). Seismic anisotropy beneath Ruapehu volcano: A possible eruption forecasting tool. *Science*, 306, 1543-1547.

Hardebeck, J.L. and Michael, A.J. (2006). Damped regional-scale stress inversions. *J. Geophys. Res.*, 111, B11310.

Hefner, W.L. et al. (2020). *J. Geophys. Res.*, 125, e2020JB019356. DOI: 10.1029/2020JB019356.

Hudson, J.A. (1981). Wave speeds and attenuation of elastic waves in material containing cracks. *Geophys. J. R. astr. Soc.*, 64, 133-150. DOI: 10.1111/j.1365-246X.1981.tb02662.x.

Julià, J. et al. (2000). Joint inversion of receiver function and surface wave dispersion observations. *Geophys. J. Int.*, 143, 99-112.

Kendall, M.J. et al. (2025). *Seismica*, 4(1). DOI: 10.26443/seismica.v4i1.1101.

Lee, M.K. et al. (2024). *Geophys. Res. Lett.*, 51, e2024GL108883. DOI: 10.1029/2024GL108883.

Martínez-Garzón, P. et al. (2014). Stress tensor changes related to fluid injection at The Geysers geothermal field. *Geophys. Res. Lett.*, 41, 8441-8449.

Massa, B. et al. (2016). *Front. Earth Sci.*, 4:103. DOI: 10.3389/feart.2016.00103.

Mavko, G., Mukerji, T., and Dvorkin, J. (2009). *The Rock Physics Handbook*, 2nd ed. Cambridge University Press.

Michael, A.J. (1984). Determination of stress from slip data: Faults and folds. *J. Geophys. Res.*, 89(B13), 11517-11526. DOI: 10.1029/JB089iB13p11517.

Mogi, K. (1958). Relations between the eruptions of various volcanoes and the deformations of the ground surfaces around them. *Bull. Earthq. Res. Inst.*, 36, 99-134.

Nooner, S.L. and Chadwick, W.W. (2016). Inflation-predictable behavior and co-eruption deformation at Axial Seamount. *Science*, 354, 1399-1403. DOI: 10.1126/science.aah4666.

Nur, A. (1971). Effects of stress on velocity anisotropy in rocks with cracks. *J. Geophys. Res.*, 76(8), 2022-2034. DOI: 10.1029/JB076i008p02022.

O'Connell, R.J. and Budiansky, B. (1974). Seismic velocities in dry and saturated cracked solids. *J. Geophys. Res.*, 79(35), 5412-5426.

Sayers, C.M. and Kachanov, M. (1995). Microcrack-induced elastic wave anisotropy of brittle rocks. *J. Geophys. Res.*, 100(B3), 4149-4156.

Sens-Schönfelder, C. and Wegler, U. (2006). Passive image interferometry and seasonal variations of seismic velocities at Merapi Volcano. *Geophys. Res. Lett.*, 33, L21302.

Shapiro, N.M. and Campillo, M. (2004). Emergence of broadband Rayleigh waves from correlations of the ambient seismic noise. *Geophys. Res. Lett.*, 31, L07614.

Silver, P.G. and Chan, W.W. (1991). Shear wave splitting and subcontinental mantle deformation. *J. Geophys. Res.*, 96(B10), 16429-16454. DOI: 10.1029/91JB00899.

Slead, S.R. et al. (2024). *J. Geophys. Res.*, 129, e2023JB028414. DOI: 10.1029/2023JB028414.

Thomsen, L. (1986). Weak elastic anisotropy. *Geophysics*, 51(10), 1954-1966. DOI: 10.1190/1.1442051.

Wallace, R.E. (1951). Geometry of shearing stress and relation to faulting. *Journal of Geology*, 59, 118-130.

Walsh, J.B. (1965). The effect of cracks on the compressibility of rock. *J. Geophys. Res.*, 70(2), 381-389.

Wang, K. et al. (2024). *Geophys. Res. Lett.*, 51, e2024GL108631. DOI: 10.1029/2024GL108631.

Wilcock, W.S.D. et al. (2018). The cabled observatory at Axial Seamount. *Oceanography*, 31(1), 114-123. DOI: 10.5670/oceanog.2018.117.

Zatsepin, S.V. and Crampin, S. (1997). Modelling the compliance of crustal rock: I. Response of shear-wave splitting to differential stress. *Geophys. J. Int.*, 129, 477-494.

Zhu, J. et al. (2026). *Geology*. DOI: 10.1130/G54254.1.
