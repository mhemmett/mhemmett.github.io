# A Unified Theoretical Framework for Joint Stress Inversion from Focal Mechanisms, Shear-Wave Splitting, and dv/v

**Working Document**
**April 2026**

---

## 1. Introduction and Motivation

This document develops the theoretical chain connecting three seismological observables — P-wave focal mechanisms, shear-wave splitting (SWS), and relative seismic velocity variations (dv/v) — to the in-situ stress tensor at Axial Seamount. Each observable senses a different projection of the stress field through distinct physical mechanisms, but all three ultimately depend on the same underlying quantities: the stress tensor σ_ij, the crack/damage state of the medium, and the resulting elastic stiffness tensor C_ijkl. The goal is to formalize a single forward model:

**σ_ij(x, t) → crack state → C_ijkl(x, t) → {focal mechanisms, SWS parameters, dv/v}**

such that a joint inversion of all three data types yields a better-constrained stress tensor than any individual observable alone.

The central physical insight is this: differential stress controls the distribution and state of microcracks (and macroscopic faults), which in turn control the fourth-order elastic stiffness tensor C_ijkl. That tensor governs wave propagation (hence splitting and velocities) and, through the Wallace-Bott hypothesis, fault slip directions (hence focal mechanisms). The stress tensor is the common cause; C_ijkl is the common mediator.

---

## 2. The Stress Tensor and Its Decomposition

The Cauchy stress tensor σ_ij is a symmetric second-order tensor with six independent components. At any point in the medium it can be decomposed as:

```
σ_ij = (1/3) σ_kk δ_ij  +  s_ij
       ─────────────────     ────
       isotropic (mean)    deviatoric
       stress / pressure   stress
```

where σ_kk = σ_11 + σ_22 + σ_33 is the trace (three times the mean stress), δ_ij is the Kronecker delta, and s_ij = σ_ij − (1/3)σ_kk δ_ij is the deviatoric stress tensor (traceless, five independent components).

In principal coordinates, σ_ij is diagonal with eigenvalues σ_1 ≥ σ_2 ≥ σ_3 (principal stresses, compression positive). The stress ratio:

```
R = (σ_1 − σ_2) / (σ_1 − σ_3)        0 ≤ R ≤ 1
```

parameterizes the shape of the stress ellipsoid. R = 0 means σ_1 = σ_2 (oblate, uniaxial compression), R = 1 means σ_2 = σ_3 (prolate, uniaxial extension), and R = 0.5 is triaxial.

**Why this decomposition matters for the joint inversion:**

- Focal mechanisms constrain the **deviatoric** stress tensor (four parameters: three Euler angles for principal axes + R), but NOT the isotropic/mean stress (Michael, 1984).
- dv/v is sensitive to **both** isotropic (mean stress / pore pressure) and deviatoric components through the acoustoelastic effect and crack closure.
- SWS is sensitive primarily to the **deviatoric** stress through crack alignment, but crack density is also affected by mean effective stress.

This complementarity is the fundamental reason a joint inversion is more powerful than any individual method.

---

## 3. Observable 1: P-Wave Focal Mechanisms → Stress Tensor

### 3.1 The Wallace-Bott Hypothesis

The foundational assumption linking earthquake slip to stress was articulated by Wallace (1951) and Bott (1959): **slip on a pre-existing fault plane occurs in the direction of maximum resolved shear stress on that plane.**

For a fault with unit normal **n** and slip direction **d**, the traction vector on the fault plane is:

```
T_i = σ_ij n_j
```

The resolved shear traction in the slip direction is:

```
τ = T_i d_i = σ_ij n_j d_i
```

The Wallace-Bott hypothesis states that d is parallel to the shear traction component of T projected onto the fault plane:

```
d_i ∝ (σ_ij n_j − (σ_kl n_k n_l) n_i)
```

That is, the slip direction equals the direction of maximum resolved shear stress on the plane, which is the projection of the traction vector onto the fault surface minus the normal component.

### 3.2 Formal Stress Inversion from Focal Mechanisms (Michael, 1984, 1987)

Michael (1984) formalized the inverse problem: given N focal mechanisms (each specifying a fault normal **n**^(k) and slip vector **d**^(k)), find the stress tensor σ_ij that best satisfies the Wallace-Bott hypothesis for all events.

The key relationship for each event k is:

```
d_i^(k) = [ σ_ij n_j^(k) − (σ_mn n_m^(k) n_n^(k)) n_i^(k) ] / |τ^(k)|
```

where |τ^(k)| is the magnitude of the resolved shear stress (a normalization factor). Since |τ^(k)| differs for each event, the absolute stress magnitude cancels. This means:

**Focal mechanism inversion constrains only the orientations of the principal stress axes (σ_1, σ_2, σ_3) and the stress ratio R. It does NOT constrain the absolute magnitudes of the principal stresses or the mean stress.**

In practice, the inversion solves for a reduced stress tensor (normalized so that the maximum differential stress = 1). The full solution space is the four parameters: three Euler angles (orientation of σ_1, σ_2, σ_3 axes) plus R.

**Assumptions in the Michael (1984) method:**

1. Stress is uniform across the fault population (within each inversion volume/time window).
2. Earthquakes occur on pre-existing faults of diverse orientations (not all faults are optimally oriented).
3. Fault slip directions are independent (one fault's slip does not influence another's).
4. The fault plane (vs. auxiliary plane) ambiguity must be resolved or marginalized.

### 3.3 Angelier's Method and Extensions

Angelier (1979, 1984) developed complementary methods for stress tensor determination from fault-slip data (originally for geological field data). The mathematical framework is equivalent — minimize the angular misfit between observed slip directions and the predicted direction of maximum resolved shear stress — but the implementation differs. Angelier's approach uses direct grid search over (σ_1, σ_2, σ_3 orientations, R) space.

### 3.4 Bayesian and Damped Approaches

Arnold and Townend (2007) introduced a Bayesian framework (implemented in the SATSI algorithm) that accounts for uncertainties in nodal plane orientations and provides posterior probability distributions on the stress parameters. This is important for quantifying uncertainty when combining with other observables.

Hardebeck and Michael (2006) developed a damped regional-scale stress inversion (DRSSI) that allows spatial variation in the stress field by solving simultaneously for stress on a grid, with a damping constraint that penalizes rapid spatial variation. This is directly relevant for mapping stress across the Axial Seamount caldera.

### 3.5 Temporal Stress Inversion

Martínez-Garzón et al. (2014, 2016) applied time-dependent stress inversion at volcanic and geothermal systems, showing that temporal windowing of focal mechanism catalogs can track the evolution of stress orientation and R through eruption cycles and injection episodes. The key is having sufficient focal mechanisms per time window (typically ≥ 15–20) for a stable inversion. Massa et al. (2016) applied similar approaches.

### 3.6 What Focal Mechanisms Constrain (Summary)

| Parameter | Constrained? | Notes |
|---|---|---|
| σ_1, σ_2, σ_3 orientations | **Yes** | Three Euler angles |
| Stress ratio R | **Yes** | Shape of deviatoric tensor |
| Absolute stress magnitudes | **No** | Normalization cancels |
| Mean / isotropic stress | **No** | Only deviatoric part matters |
| Pore pressure | **No** | Enters only through effective stress |
| Temporal changes in orientation | **Yes** | With sufficient events per window |

**Key references:** Wallace (1951), Bott (1959), Angelier (1979, 1984), Michael (1984; DOI: 10.1029/JB089iB13p11517), Michael (1987), Hardebeck and Michael (2006), Arnold and Townend (2007), Martínez-Garzón et al. (2014, 2016), Massa et al. (2016; DOI: 10.3389/feart.2016.00103).

---

## 4. Observable 2: Shear-Wave Splitting → Anisotropy → Stress

### 4.1 Physical Basis: Stress-Aligned Microcracks

Crampin (1981, 1984) established that the upper crust is pervaded by microcracks whose orientations are controlled by the ambient stress field. Under differential stress, cracks preferentially open perpendicular to the minimum compressive stress σ_3 (equivalently, they are aligned parallel to the maximum horizontal compressive stress σ_H). This population of aligned, fluid-filled microcracks creates an effectively anisotropic medium — the **Extensive-Dilatancy Anisotropy (EDA)** model.

The physical chain is:

```
Differential stress  →  preferential crack alignment  →  elastic anisotropy  →  shear-wave splitting
```

### 4.2 Crack-Induced Anisotropy: Hudson's Theory

Hudson (1981) derived the effective elastic stiffness tensor for a medium containing a dilute concentration of aligned penny-shaped cracks. The key result is a perturbation expansion:

```
C_ijkl = C⁰_ijkl  +  C¹_ijkl  +  C²_ijkl  +  ...
```

where C⁰ is the isotropic background stiffness, and C¹, C² are first- and second-order corrections due to cracks.

**For a single set of aligned penny-shaped cracks with normals along the x₃ direction:**

The crack density parameter is:

```
ε = N a³ / V
```

where N is the number of cracks in volume V and a is the crack radius.

The first-order correction involves:

```
C¹_ijkl ∝ ε × f(λ, μ, κ, μ')
```

where λ, μ are the Lamé parameters of the background, and κ, μ' parameterize the crack interior (fluid-filled vs. dry). For fluid-filled cracks, the normal compliance is reduced (cracks resist normal closure due to fluid incompressibility) while the tangential compliance remains (fluid offers no shear resistance).

The resulting medium has **hexagonal (transverse isotropic)** symmetry with the symmetry axis parallel to the crack normals (perpendicular to σ_H). In Voigt notation, the stiffness matrix has five independent constants: C₁₁, C₁₃, C₃₃, C₄₄, C₆₆ (with C₁₂ = C₁₁ − 2C₆₆).

### 4.3 The Sayers-Kachanov Crack Compliance Tensor

Sayers and Kachanov (1995) provided an alternative (and often more convenient) formulation using the **crack compliance tensor**. Rather than perturbing the stiffness, they perturb the compliance:

```
S_ijkl = S⁰_ijkl  +  ΔS_ijkl
```

where the crack-induced excess compliance is:

```
ΔS_ijkl = (1/4)(α_ik δ_jl + α_il δ_jk + α_jk δ_il + α_jl δ_ik) + β_ijkl
```

Here α_ij is the **second-rank crack density tensor**:

```
α_ij = (1/V) Σ_r  B_T^(r)  a_r³  n_i^(r) n_j^(r)
```

and β_ijkl is the **fourth-rank crack density tensor**:

```
β_ijkl = (1/V) Σ_r  (B_N^(r) − B_T^(r))  a_r³  n_i^(r) n_j^(r) n_k^(r) n_l^(r)
```

where the sum is over all cracks r, n^(r) is the unit normal to crack r, a_r is the crack radius, and B_N, B_T are the normal and tangential compliances of individual cracks (which depend on crack aspect ratio, fluid content, and pressure).

**This formulation is crucial because:**

1. It naturally handles arbitrary crack orientation distributions (not just a single aligned set).
2. The second-rank tensor α_ij captures the dominant anisotropy (HTI symmetry for one preferred orientation, orthorhombic for two, etc.).
3. The fourth-rank tensor β_ijkl captures deviations from elliptical anisotropy.
4. Both α_ij and β_ijkl are directly related to the stress tensor through stress-dependent crack opening/closing.

### 4.4 Connecting Cracks to Stress: Stress-Dependent Crack State

The critical link between stress and anisotropy is that the crack population evolves with stress. Under increasing differential stress:

- Cracks oriented normal to σ_1 tend to **close** (increasing stiffness in that direction).
- Cracks oriented normal to σ_3 tend to **open** or remain open (maintaining compliance in that direction).
- The net effect is increased anisotropy as differential stress increases.

Quantitatively, for a crack with normal **n**, the effective normal stress on the crack face is:

```
σ_n = σ_ij n_i n_j − P_pore
```

where P_pore is pore fluid pressure. A crack closes when σ_n exceeds the crack's closure stress σ_c (which depends on crack aspect ratio α_c):

```
σ_c ∝ E α_c / (1 − ν²)
```

where E and ν are Young's modulus and Poisson's ratio.

**The Zatsepin-Crampin APE (Anisotropic Poro-Elasticity) model** (Zatsepin and Crampin, 1997) formalizes this by modeling the crack aspect ratio distribution as a function of stress and pore pressure. In APE:

1. Start with an initial crack aspect ratio distribution N(α, θ, φ) giving the number density of cracks with aspect ratio α and orientation (θ, φ).
2. Apply a stress change Δσ_ij.
3. Cracks deform elastically: aspect ratios change as Δα ∝ Δσ_n / E.
4. Cracks whose aspect ratio goes to zero are closed; cracks below critical aspect ratio are fluid-saturated.
5. Recompute C_ijkl from the updated crack population using Hudson's theory.
6. Predict new splitting parameters from the updated C_ijkl.

This provides a direct, physics-based forward model: **Δσ_ij → ΔN(α, θ, φ) → ΔC_ijkl → Δ(φ, δt)**.

### 4.5 Shear-Wave Splitting Observables

When a shear wave enters an anisotropic medium, it splits into two quasi-shear waves with orthogonal polarizations traveling at different velocities. The two observables are:

- **φ (fast polarization direction):** The polarization azimuth of the faster quasi-S wave. For a medium with vertical crack planes striking in the direction of σ_H, the fast direction is approximately parallel to σ_H (or equivalently, parallel to crack strike).

- **δt (delay time):** The time delay between the fast and slow arrivals, related to the path-integrated anisotropy:

```
δt = ∫_path  [ 1/v_s2(ℓ) − 1/v_s1(ℓ) ] dℓ
```

where v_s1, v_s2 are the fast and slow S-wave velocities along the ray path ℓ. For weak anisotropy:

```
δt ≈ ∫_path  (γ / v̄_s) dℓ
```

where γ is Thomsen's (1986) shear-wave anisotropy parameter and v̄_s is the isotropic reference S-wave velocity.

**The Silver and Chan (1991) method** measures (φ, δt) by grid-searching for the pair of values that best linearizes the particle motion of the corrected (un-split) waveform. Specifically, they minimize the smaller eigenvalue of the covariance matrix of the corrected horizontal seismograms:

```
minimize  λ₂(φ_trial, δt_trial)
```

where λ₂ is the smaller eigenvalue of the 2×2 covariance matrix of the corrected fast and slow components.

### 4.6 Splitting Intensity (Chevrot, 2000)

Chevrot (2000) introduced the **splitting intensity** SI as a more robust measure for tomographic inversion. For a nearly vertically incident shear wave, SI is defined as:

```
SI = −2 [d_T(t) · d_R'(t)] / [d_R(t) · d_R(t)]
```

where d_T is the transverse component seismogram, d_R is the radial component, and d_R' = ∂d_R/∂t is its time derivative. For weak anisotropy:

```
SI ≈ δt sin(2(φ_back − φ_fast))
```

where φ_back is the event backazimuth and φ_fast is the fast direction.

**Why SI is powerful for inversion:**

1. SI is a **linear** functional of the anisotropy along the ray path (unlike the nonlinear relationship between (φ, δt) and anisotropy for complex media).
2. It can be written as a path integral: SI = ∫_path K(ℓ) · **a**(ℓ) dℓ, where K is a sensitivity kernel and **a** is the local anisotropy.
3. This linearity enables standard tomographic techniques (e.g., LSQR, conjugate gradient) for 3D anisotropy inversion.
4. Finite-frequency sensitivity kernels (Favier and Chevrot, 2003) can be computed, enabling resolution of structure finer than the Fresnel zone.

### 4.7 Thomsen Parameters and Weak Anisotropy

Thomsen (1986) parameterized weak transverse isotropy with three dimensionless parameters:

```
ε = (C₁₁ − C₃₃) / (2 C₃₃)       (P-wave anisotropy)
δ = [(C₁₃ + C₄₄)² − (C₃₃ − C₄₄)²] / [2 C₃₃ (C₃₃ − C₄₄)]   (near-vertical P-wave anisotropy)
γ = (C₆₆ − C₄₄) / (2 C₄₄)       (S-wave anisotropy)
```

For S-wave splitting specifically, γ is the most directly relevant parameter. The splitting delay time for a vertically propagating S-wave through a layer of thickness L is approximately:

```
δt ≈ γ L / v̄_s
```

These parameters connect directly to the crack density ε_crack through Hudson's (1981) theory. For example, for dry penny-shaped cracks with normals along x₃:

```
γ_Thomsen ≈ (8/3) ε_crack (1 − ν) / (2 − ν)
```

### 4.8 What SWS Constrains (Summary)

| Parameter | Constrained? | Notes |
|---|---|---|
| σ_H orientation (horizontal) | **Yes** | From fast polarization φ |
| Degree of anisotropy / crack density | **Yes** | From delay time δt |
| Spatial variation of anisotropy | **Yes** | Via splitting intensity tomography |
| Temporal changes in stress | **Yes** | From δt(t) and φ(t) changes |
| Absolute stress magnitude | **No** | Only relative anisotropy |
| Mean / isotropic stress | **Indirectly** | Mean stress affects crack closure |
| Stress ratio R | **Partially** | Multiple crack sets → orthorhombic → constrain relative stresses |
| σ_1 plunge | **Limited** | Mostly sensitive to horizontal projection |

**Key references:** Crampin (1981, 1984), Hudson (1981; DOI: 10.1111/j.1365-246X.1981.tb02662.x), Thomsen (1986; DOI: 10.1190/1.1442051), Silver and Chan (1991; DOI: 10.1029/91JB00899), Sayers and Kachanov (1995; DOI: 10.1016/0148-9062(95)00022-4), Zatsepin and Crampin (1997), Chevrot (2000; DOI: 10.1029/2000JB900199), Crampin and Peacock (2008), Kendall et al. (2025; DOI: 10.26443/seismica.v4i1.1101).

---

## 5. Observable 3: dv/v from Ambient Noise Interferometry → Stress

### 5.1 Extracting Green's Functions from Ambient Noise

Shapiro and Campillo (2004) and Campillo and Paul (2003) demonstrated that cross-correlating long records of ambient seismic noise between two stations yields the empirical Green's function between those stations. The theoretical basis (Lobkis and Weaver, 2001; Snieder, 2004) is that for a diffuse wavefield:

```
dC_AB(t)/dt ≈ −G_AB(t) + G_AB(−t)
```

where C_AB(t) is the noise cross-correlation between stations A and B, and G_AB(t) is the Green's function (impulse response) between them.

### 5.2 Monitoring Velocity Changes: Coda Wave Interferometry

Once reference Green's functions are established, temporal changes in the medium are detected as time shifts in the coda (multiply scattered) portion of the cross-correlations. Two main approaches:

**Stretching technique (Sens-Schönfelder and Wegler, 2006; Lobkis and Weaver, 2003):** If the velocity change is homogeneous, all arrivals in the coda shift by the same relative amount. The measured stretching factor ε that maximizes the correlation between the reference and current coda gives:

```
dv/v = −ε
```

where ε is defined by: the current Green's function G(t) ≈ G_ref(t(1 + ε)).

**Moving-Window Cross-Spectral (MWCS) method (Clarke et al., 2011; Brenguier et al., 2008):** The coda is divided into moving windows, and the phase shift in each window is measured via cross-spectrum. A linear regression of phase shift vs. lag time gives:

```
dt/t = −dv/v
```

Both methods assume that the velocity change is spatially uniform along the scattering paths sampled by the coda waves. Deviations from uniformity are a source of bias that sensitivity kernel analysis can help address.

### 5.3 What Controls dv/v: The Acoustoelastic Effect

The acoustic velocities in a stressed rock are not the same as in an unstressed rock. This is the **acoustoelastic effect**, arising from the nonlinear (third-order) elasticity of rocks. The relationship between velocity and stress for a hyperelastic material is:

```
ρ v²_ij = C_ijkl + C_ijklmn ε_mn + ...
```

where C_ijklmn are **third-order elastic constants** (TOECs) and ε_mn is the strain (related to stress via the constitutive law). For small perturbations around a reference state:

```
Δv / v₀ ≈ (1/2) (ΔC_iijj / C⁰_iijj)
```

In practice, for an isotropic medium under hydrostatic stress change ΔP:

```
dv_P / v_P ≈ (1 / 2K) [1 + (C₁₁₁ + 2C₁₁₂) / (3K)] ΔP
```

where K is the bulk modulus and C₁₁₁, C₁₁₂ are third-order elastic constants (Murnaghan or Landau notation). The key point: **TOECs for crustal rocks are large and negative**, meaning velocity increases under compression. Typical sensitivities are dv/v ~ 10⁻⁶ to 10⁻⁴ per MPa.

### 5.4 Crack-Mediated Velocity Changes

In cracked rocks (which is the relevant regime for the shallow crust at Axial Seamount), the dominant mechanism for velocity-stress coupling is not intrinsic TOECs but rather **stress-dependent crack closure/opening**. This is much more sensitive than intrinsic nonlinearity.

Following Nur (1971) and O'Connell and Budiansky (1974):

- Increasing compressive stress closes cracks → increases elastic moduli → increases velocity.
- Decreasing stress (or increasing pore pressure) opens cracks → decreases moduli → decreases velocity.

For a distribution of cracks with aspect ratios α and closure stress σ_c(α):

```
v(σ) = v_∞ × [1 − f(ε_open(σ))]^(1/2)
```

where v_∞ is the crack-free velocity, ε_open(σ) is the crack density of still-open cracks at stress σ, and f is a function that depends on crack geometry (penny-shaped, ellipsoidal, etc.).

For penny-shaped cracks (O'Connell and Budiansky, 1974):

```
K_eff / K₀ = 1 − (16/9) [(1 − ν²) / (1 − 2ν)] ε
μ_eff / μ₀ = 1 − (32/45) [(1 − ν)(5 − ν) / (2 − ν)] ε
```

where K₀, μ₀ are the crack-free bulk and shear moduli, and ε is the open crack density.

**Differentiating with respect to stress gives dv/v:**

```
d(v_P)/v_P = −(8/9) [(1 − ν²)/(1 − 2ν)] (1/v²_P₀) [dK + (4/3)dμ] × (dε/dσ) × Δσ
```

The term dε/dσ encodes how many cracks close per unit stress increase, which depends on the crack aspect ratio distribution.

### 5.5 Pore Pressure Effects

At Axial Seamount, hydrothermal fluids and magmatic volatiles mean that pore pressure P_p is a major player. The effective stress principle (Terzaghi, with Biot coefficient η):

```
σ_ij^eff = σ_ij − η P_p δ_ij
```

means that dv/v responds to changes in **effective stress**, not total stress alone. A pore pressure increase of ΔP_p has the same effect on crack state (and hence velocity) as a mean stress decrease of η ΔP_p. Therefore:

```
dv/v = A Δσ_mean − B ΔP_p + (higher-order deviatoric terms)
```

where A and B are sensitivity coefficients that depend on the crack state and rock properties. **Separating stress changes from pore pressure changes requires additional constraints** — this is one reason the joint inversion with SWS and focal mechanisms is so powerful.

### 5.6 Sensitivity Kernels and Depth Resolution

The dv/v measured from noise cross-correlations is not a point measurement but a path-integrated quantity weighted by a sensitivity kernel K(x):

```
dv/v_measured = ∫∫∫ K(x) [dv(x)/v(x)] d³x
```

The kernel K depends on the frequency band of the measurement and the wave type:

- **Surface waves (Rayleigh, Love):** Depth sensitivity scales with wavelength. Higher frequencies → shallower sensitivity. For a typical 0.5–2 Hz measurement, sensitivity peaks at ~0.5–3 km depth.
- **Body-wave coda:** Sensitivity is distributed along scattering paths, extending to deeper depths but with poorer spatial localization.

**Frequency-dependent dv/v measurements** can therefore provide crude depth resolution: measuring dv/v at multiple frequency bands constrains velocity changes at different depth ranges. This is relevant for distinguishing shallow hydrothermal/pore-pressure effects from deeper magmatic stress changes at Axial Seamount.

### 5.7 dv/v at Volcanoes

Brenguier et al. (2008) demonstrated continuous dv/v monitoring at Piton de la Fournaise volcano, showing that pre-eruptive inflation produces measurable velocity decreases (interpreted as pore pressure increase or stress relaxation opening cracks). Donaldson et al. (2017; DOI: 10.1126/sciadv.1700219) showed similar results at Kīlauea, with dv/v changes of 0.1–0.5% correlated with magmatic activity.

Lee et al. (2024; DOI: 10.1029/2024GL108883) applied these techniques to Axial Seamount, demonstrating that dv/v monitoring is feasible in the submarine OBS setting and shows systematic changes across the eruption cycle.

### 5.8 What dv/v Constrains (Summary)

| Parameter | Constrained? | Notes |
|---|---|---|
| Mean effective stress changes | **Yes** | Primary sensitivity |
| Pore pressure changes | **Yes** | But trade-off with mean stress |
| Deviatoric stress changes | **Partially** | Through anisotropic dv/v if measured directionally |
| Stress orientation | **No** | Isotropic scalar measurement (unless azimuthally resolved) |
| Depth of stress change | **Partially** | Via frequency-dependent analysis |
| Temporal dynamics | **Yes** | Continuous monitoring, sub-daily resolution possible |
| Absolute stress | **No** | Only relative changes |

**Key references:** Nur (1971; DOI: 10.1029/JB076i008p02022), O'Connell and Budiansky (1974), Birch (1960, 1961), Shapiro and Campillo (2004), Snieder (2006), Sens-Schönfelder and Wegler (2006), Brenguier et al. (2008), Donaldson et al. (2017; DOI: 10.1126/sciadv.1700219), Lee et al. (2024; DOI: 10.1029/2024GL108883).

---

## 6. The Unifying Framework: From Stress to All Three Observables

### 6.1 The Central Object: The Elastic Stiffness Tensor C_ijkl(σ, P_p)

All three observables are, at root, determined by the elastic stiffness tensor C_ijkl and its spatial and temporal variations. The stiffness tensor is in turn determined by:

1. **Intrinsic mineral elasticity** (relatively constant on eruption-cycle timescales).
2. **Crack state** — the distribution of microcracks in orientation, density, and aperture — which IS controlled by stress and pore pressure.
3. **Fluid content and saturation** (important for absolute values but less for temporal changes).

Therefore the forward model is:

```
[σ_ij(x,t), P_p(x,t)]  ──→  crack state  ──→  C_ijkl(x,t)  ──→  observables
```

### 6.2 Step 1: Stress → Crack State

Given a stress tensor σ_ij and pore pressure P_p, the state of a crack with normal **n** and initial aspect ratio α₀ is determined by the effective normal stress on the crack:

```
σ_n^eff = σ_ij n_i n_j − P_p
```

The crack closes when σ_n^eff exceeds the closure stress σ_c(α₀). For a population of cracks with initial orientation-aspect ratio distribution N₀(θ, φ, α₀):

```
N_open(θ, φ, α; σ, P_p) = N₀(θ, φ, α₀) × H(α₀ − α_closed(σ_n^eff(θ, φ)))
```

where H is the Heaviside step function and α_closed is the critical aspect ratio below which cracks are closed at the current stress. More realistically, crack closure is gradual, and the aspect ratio evolves continuously:

```
α(σ) = α₀ − σ_n^eff / (π E / (4(1−ν²)))
```

for a penny-shaped crack in a linearly elastic medium.

### 6.3 Step 2: Crack State → C_ijkl

Using the Sayers-Kachanov formulation (Section 4.3), the excess compliance due to the open crack population is:

```
ΔS_ijkl = (1/4)(α_ik δ_jl + α_il δ_jk + α_jk δ_il + α_jl δ_ik) + β_ijkl
```

with α_ij and β_ijkl computed from the open crack distribution:

```
α_ij = (1/V) Σ_{open cracks} B_T a³ n_i n_j

β_ijkl = (1/V) Σ_{open cracks} (B_N − B_T) a³ n_i n_j n_k n_l
```

The total compliance is S = S⁰ + ΔS, and the stiffness is C = S⁻¹.

### 6.4 Step 3: C_ijkl → Focal Mechanism Predictions

While focal mechanisms are usually treated as inputs to stress inversion (not predictions from C_ijkl), in a forward model the connection is:

1. **Failure criterion:** Given σ_ij, faults fail when the Coulomb failure criterion is met:
   ```
   |τ| = c + μ_f (σ_n − P_p)
   ```
   where τ is shear stress on the fault, σ_n is normal stress, c is cohesion, and μ_f is friction coefficient.

2. **Slip direction:** The Wallace-Bott hypothesis gives the slip direction from σ_ij (Section 3.1).

3. **The predicted focal mechanism** (P, T, B axes and nodal planes) follows from the moment tensor, which for a shear dislocation on a fault with normal **n** and slip **d** is:
   ```
   M_ij = μ A [d_i n_j + d_j n_i]
   ```
   where μ is shear modulus (from C_ijkl), A is fault area, and [d_i n_j + d_j n_i] is the symmetric part.

4. **In an anisotropic medium**, the relationship between moment tensor and radiation pattern is modified. The P-wave radiation pattern in an anisotropic medium with stiffness C_ijkl is:
   ```
   u_i^P(x) ∝ C_ijkl n_k d_l × G_ij(x, x_source)
   ```
   where G_ij is the Green's function for the anisotropic medium. For weak anisotropy, this is a perturbation on the isotropic case.

**The key insight for joint inversion:** The same C_ijkl that predicts splitting and velocities also modifies the seismic radiation pattern and the failure criterion (through stress concentration near crack tips). The stress tensor that best explains the focal mechanisms must be consistent with the stress tensor implied by the observed anisotropy.

### 6.5 Step 3: C_ijkl → Shear-Wave Splitting Predictions

Given C_ijkl(x) at every point along a ray path, the splitting can be computed by solving the Christoffel equation at each point:

```
[Γ_ik − ρ v² δ_ik] p_k = 0
```

where:
```
Γ_ik = C_ijkl n̂_j n̂_l
```

is the Christoffel matrix, n̂ is the wave propagation direction, and the eigenvalues ρv² give the three phase velocities (one quasi-P, two quasi-S). The eigenvectors give the polarization directions.

For the two quasi-S waves with velocities v_s1 > v_s2 and polarization directions **g**₁, **g**₂:

- The fast polarization direction φ = azimuth of **g**₁ projected onto the horizontal plane.
- The differential delay accumulates along the path:
  ```
  δt = ∫_path [1/v_s2 − 1/v_s1] dℓ
  ```

In practice, for a 3D model, this requires **ray tracing** through the anisotropic velocity model, accumulating splitting at each step. For complex paths, layer-stripping or matrix propagator methods handle the fact that the effective splitting from multiple anisotropic layers is not a simple sum.

Using splitting intensity (Chevrot, 2000), the prediction becomes a linear path integral (Section 4.6), greatly simplifying the inversion.

### 6.6 Step 3: C_ijkl → dv/v Predictions

From C_ijkl(x, t), the isotropic-equivalent velocities are:

```
v_P(x, t) = √[(C₁₁ + C₂₂ + C₃₃ + 2C₁₂ + 2C₁₃ + 2C₂₃) / (9ρ) + ...]
```

(using the Voigt-Reuss-Hill average or similar isotropic projection of C_ijkl).

The predicted dv/v between times t₁ and t₂ is:

```
dv/v(x) = [v(x, t₂) − v(x, t₁)] / v(x, t₁)
```

and the observed dv/v is the spatially weighted integral:

```
dv/v_obs = ∫ K(x) × dv/v(x) d³x
```

where K(x) is the sensitivity kernel for the specific station pair and frequency band.

**For directional dv/v** (measuring velocity changes along different azimuths using different station pairs), the anisotropic part of C_ijkl produces azimuth-dependent velocity changes that carry information about the deviatoric stress, beyond just the mean stress.

---

## 7. Complementarity of the Three Observables

### 7.1 What Each Observable Uniquely Contributes

The power of the joint inversion lies in the **complementarity** of the three observables. Consider the stress tensor σ_ij with its six independent components (or equivalently: three principal orientations, two principal stress ratios, and one absolute magnitude scale). Adding pore pressure P_p as a seventh unknown:

**Focal mechanisms** (via Michael 1984 inversion):
- Constrain: 3 principal stress orientations + R (4 parameters)
- Do NOT constrain: absolute magnitude, mean stress, pore pressure
- Temporal resolution: limited by earthquake rate (need ~15+ events per window)
- Spatial resolution: limited to regions with seismicity

**Shear-wave splitting** (via splitting intensity tomography):
- Constrain: σ_H orientation (from φ), spatial distribution of anisotropy (from δt)
- Provide: 3D anisotropy model → relates to crack density tensor α_ij
- Temporal resolution: good (each earthquake provides a measurement)
- Spatial resolution: depends on ray coverage (good at Axial with OOI array)
- Additional: through the crack model, α_ij relates to the full deviatoric stress

**dv/v** (via ambient noise interferometry):
- Constrain: changes in mean effective stress (σ_mean − P_p) primarily
- Temporal resolution: excellent (daily or better, continuous)
- Spatial resolution: limited by station pair geometry and frequency
- Additional: frequency-dependent analysis gives depth information

### 7.2 The Complementarity Matrix

```
                    Focal Mech    SWS         dv/v
─────────────────────────────────────────────────────
σ₁ orientation      ●●●           ●●          ○
σ₂ orientation      ●●●           ○           ○
σ₃ orientation      ●●●           ●●          ○
Stress ratio R      ●●●           ●           ○
Crack density       ○             ●●●         ●●
Mean stress         ○             ●           ●●●
Pore pressure       ○             ○           ●●
Temporal dynamics   ●             ●●          ●●●
3D spatial          ●●            ●●●         ●
Depth resolution    ●●            ●           ●●

●●● = strong constraint, ●● = moderate, ● = weak, ○ = none
```

### 7.3 Key Trade-offs Resolved by Joint Inversion

1. **Stress vs. pore pressure in dv/v:** dv/v alone cannot separate a stress increase from a pore pressure decrease. But if focal mechanisms show no change in stress orientation/R, and SWS shows no change in anisotropy orientation but a change in δt, this pattern is more consistent with a pore pressure change (isotropic) than a deviatoric stress change.

2. **Stress orientation ambiguity in SWS:** In regions with multiple crack sets (common near volcanoes with complex structure), the fast direction φ may not simply reflect σ_H. Focal mechanisms provide an independent constraint on stress orientation that disambiguates the SWS interpretation.

3. **Absolute stress scale:** Neither focal mechanisms nor SWS constrain absolute stress. But dv/v, combined with laboratory-calibrated velocity-stress relationships and a known reference state, can provide an estimate of absolute stress changes, anchoring the relative measurements from the other two methods.

4. **Non-uniqueness in C_ijkl:** The inverse problem of going from SWS observations to C_ijkl is non-unique (many stiffness tensors can produce the same splitting for a given ray geometry). Adding dv/v constrains the isotropic part of C_ijkl. Adding focal mechanisms constrains the stress tensor that must be consistent with that C_ijkl through the crack model. Together, the non-uniqueness is substantially reduced.

---

## 8. Formulation of the Joint Inverse Problem

### 8.1 Model Parameters

The model parameters to be determined at each point x and time t are:

**Primary unknowns (what we solve for):**
- σ_ij(x, t): 6 independent stress components (or equivalently σ₁, σ₂, σ₃ orientations + magnitudes)
- P_p(x, t): pore pressure

**Auxiliary/structural parameters (determined separately or jointly):**
- v_P⁰(x), v_S⁰(x): reference (crack-free or reference-state) velocity model
- N₀(α, θ, φ; x): initial crack population distribution
- Rheological parameters: friction coefficient μ_f, Biot coefficient η

### 8.2 Forward Model (Summary)

```
m = [σ_ij(x,t), P_p(x,t)]   (model parameters)
         │
         ▼
    crack state evolution
    N_open(α, θ, φ; x, t) = f(m, N₀)
         │
         ▼
    elastic stiffness
    C_ijkl(x, t) = g(N_open, C⁰)
         │
         ├──→ Focal mechanism predictions: d_pred^(k) = h_FM(σ, fault geometry)
         │
         ├──→ SWS predictions: SI_pred^(k) = h_SWS(C_ijkl, ray geometry)
         │
         └──→ dv/v predictions: (dv/v)_pred^(ij) = h_dv/v(C_ijkl, sensitivity kernels)
```

### 8.3 Misfit Function

The joint inversion minimizes a combined misfit:

```
Φ(m) = w_FM Φ_FM(m) + w_SWS Φ_SWS(m) + w_dv/v Φ_dv/v(m) + Φ_reg(m)
```

where:

**Focal mechanism misfit:**
```
Φ_FM = Σ_k  angular_misfit(d_obs^(k), d_pred^(k))²
```
(sum over events k of the angular misfit between observed and predicted slip directions)

**SWS misfit:**
```
Φ_SWS = Σ_k  [SI_obs^(k) − SI_pred^(k)]² / σ²_SI
```
(sum over SWS measurements of squared residuals in splitting intensity)

**dv/v misfit:**
```
Φ_dv/v = Σ_(i,j,t)  [(dv/v)_obs^(ij)(t) − (dv/v)_pred^(ij)(t)]² / σ²_dv/v
```
(sum over station pairs (i,j) and time windows t)

**Regularization:**
```
Φ_reg = λ_s ||∇²σ||² + λ_t ||∂σ/∂t||²
```
(spatial and temporal smoothing of the stress field)

The weights w_FM, w_SWS, w_dv/v must be chosen to balance the different data types (different numbers of observations, different noise levels, different physical dimensions). This can be done by normalizing each misfit by the number of data and the expected variance, or by L-curve analysis.

### 8.4 Inversion Strategy

The forward model contains nonlinear relationships (especially the stress → crack state step), so the inversion will likely require iterative linearization:

1. **Initialize** with the focal mechanism stress inversion result (Michael, 1984) for stress orientations and R. Use dv/v as a proxy for initial mean stress variations.

2. **Forward compute** C_ijkl from the initial stress model using the crack model.

3. **Compute predicted observables** for all three data types.

4. **Compute Fréchet derivatives** (sensitivity of each observable to perturbations in σ_ij and P_p). These can be computed analytically for the crack model or numerically by finite differencing.

5. **Update the model** by solving the linearized inverse problem (e.g., using damped least squares or conjugate gradient methods).

6. **Iterate** until convergence.

The linearization at each iteration uses the Jacobian matrix:

```
J = [ ∂d_FM / ∂m ;  ∂SI / ∂m ;  ∂(dv/v) / ∂m ]
```

which has block structure reflecting the three data types.

---

## 9. Resolving Non-uniqueness in C_ijkl

### 9.1 The Problem

The elastic stiffness tensor C_ijkl has 21 independent components (for the most general triclinic symmetry). Even with the assumption of crack-induced anisotropy (which reduces the degrees of freedom through the crack model), inferring C_ijkl from limited seismic observations is highly non-unique.

### 9.2 How Each Observable Constrains C_ijkl

- **SWS** constrains the **anisotropic** part of C_ijkl (differences between stiffness in different directions), particularly the components that control S-wave velocity differences.
- **dv/v** constrains the **isotropic average** of C_ijkl (or more precisely, changes in the directionally averaged moduli).
- **Focal mechanisms** constrain the **stress tensor**, which, through the crack model, constrains the **physical mechanism** generating the anisotropy in C_ijkl.

### 9.3 The Role of the Crack Model as Physical Regularization

The crack model (Hudson, 1981; Sayers and Kachanov, 1995; Zatsepin and Crampin, 1997) acts as a powerful physical regularization: instead of inverting for 21 independent elastic constants, we invert for a much smaller set of parameters that describe the crack population and stress state. The crack model then predicts the full C_ijkl tensor self-consistently.

This reduces the effective number of unknowns from 21 (C_ijkl) to ~7 per point (6 stress components + 1 pore pressure), plus global parameters describing the initial crack population. The physics ensures that the resulting C_ijkl is always consistent with a physically realizable crack distribution.

---

## 10. Specific Considerations for Axial Seamount

### 10.1 Geological Context

Axial Seamount is a mid-ocean ridge hotspot volcano on the Juan de Fuca Ridge. Its stress field is controlled by:

- **Regional plate-tectonic stress:** Ridge-perpendicular extension (σ₃ ≈ E-W at this latitude).
- **Magmatic pressure:** Inflation/deflation of the shallow magma reservoir beneath the caldera (Nooner and Chadwick, 2016).
- **Caldera ring faults:** Localized stress concentrations along the caldera boundary (Levy et al., various; Hefner et al., 2020).
- **Hydrothermal circulation:** Pore pressure variations driven by fluid flow.

The eruption cycle (Nooner and Chadwick, 2016) involves progressive inflation (increasing magma pressure → increasing σ_H → changing anisotropy and dv/v) punctuated by eruptions (rapid deflation → stress reversal).

### 10.2 Expected Signatures in Each Observable

**Pre-eruption inflation:**
- Focal mechanisms: possible rotation of σ₁ toward vertical near the magma chamber (increasing contribution of magma pressure to vertical stress); R may change.
- SWS: φ may rotate to reflect the evolving stress pattern around the inflating source; δt should increase as differential stress increases and cracks align.
- dv/v: should decrease (increasing pore pressure and/or crack opening from magma-induced stress).

**Co-eruption deflation:**
- Focal mechanisms: normal faulting dominance on ring faults; σ₃ becomes vertical in some regions.
- SWS: rapid change in δt (crack population disrupted); possible φ rotation.
- dv/v: sharp increase (stress relaxation closes cracks; pressure decrease) — this is the signature seen by Lee et al. (2024).

**Post-eruption re-inflation:**
- Gradual return toward pre-eruption patterns in all three observables.

### 10.3 Velocity Model Requirements

The joint inversion requires a high-quality 3D P- and S-wave velocity model for:

1. **Ray tracing** through the anisotropic model (for SWS forward modeling).
2. **Sensitivity kernel computation** (for dv/v forward modeling).
3. **Reference velocities** to convert dv/v to absolute velocity changes.
4. **Earthquake location** (the relocated catalog from Wang et al. 2024 already uses a good velocity model).

Existing velocity models at Axial (e.g., from ocean-bottom seismometer tomography; Baillard et al., various; West et al., 2001) can serve as starting models.

---

## 11. Summary of Key Equations

### The Forward Model Chain

**1. Effective stress on crack face:**
```
σ_n^eff(n̂) = σ_ij n̂_i n̂_j − η P_p
```

**2. Crack density tensors (Sayers-Kachanov):**
```
α_ij = (1/V) Σ_r B_T^(r) a_r³ n_i^(r) n_j^(r)     [for open cracks only]
β_ijkl = (1/V) Σ_r (B_N^(r) − B_T^(r)) a_r³ n_i^(r) n_j^(r) n_k^(r) n_l^(r)
```

**3. Stiffness tensor:**
```
C_ijkl = [S⁰_ijkl + ΔS_ijkl(α, β)]⁻¹
```

**4a. Focal mechanism prediction (Wallace-Bott):**
```
d_i = [σ_ij n_j − (σ_kl n_k n_l) n_i] / |τ|
```

**4b. Splitting intensity prediction:**
```
SI = ∫_ray K_SI(ℓ) · a(ℓ) dℓ      where a(ℓ) depends on C_ijkl(ℓ)
```

**4c. dv/v prediction:**
```
dv/v = ∫∫∫ K_dv/v(x) × [v(C_ijkl(x,t₂)) − v(C_ijkl(x,t₁))] / v(x,t₁) d³x
```

**5. Joint misfit:**
```
Φ = w₁ Φ_FM + w₂ Φ_SWS + w₃ Φ_dv/v + λ_s||∇²σ||² + λ_t||∂σ/∂t||²
```

---

## 12. Key References

### Focal Mechanism Stress Inversion
- Angelier, J. (1979). Determination of the mean principal directions of stresses for a given fault population. *Tectonophysics*, 56(3-4), T17-T26.
- Angelier, J. (1984). Tectonic analysis of fault slip data sets. *J. Geophys. Res.*, 89(B7), 5835-5848.
- Arnold, R. and Townend, J. (2007). A Bayesian approach to estimating tectonic stress from seismological data. *Geophys. J. Int.*, 170, 1336-1356.
- Bott, M.H.P. (1959). The mechanics of oblique slip faulting. *Geological Magazine*, 96(2), 109-117.
- Hardebeck, J.L. and Michael, A.J. (2006). Damped regional-scale stress inversions: Methodology and examples for southern California and the Coalinga aftershock sequence. *J. Geophys. Res.*, 111, B11310.
- Martínez-Garzón, P. et al. (2014). Stress tensor changes related to fluid injection at The Geysers geothermal field, California. *Geophys. Res. Lett.*, 41, 8441-8449.
- Massa, B. et al. (2016). Focal mechanisms and stress field in the Campi Flegrei volcanic area. *Front. Earth Sci.*, 4:103. DOI: 10.3389/feart.2016.00103.
- Michael, A.J. (1984). Determination of stress from slip data: Faults and folds. *J. Geophys. Res.*, 89(B13), 11517-11526. DOI: 10.1029/JB089iB13p11517.
- Michael, A.J. (1987). Use of focal mechanisms to determine stress: A control study. *J. Geophys. Res.*, 92(B1), 357-368.
- Wallace, R.E. (1951). Geometry of shearing stress and relation to faulting. *Journal of Geology*, 59, 118-130.

### Shear-Wave Splitting Theory
- Babuska, V. and Cara, M. (1991). *Seismic Anisotropy in the Earth*. Kluwer Academic Publishers.
- Chevrot, S. (2000). Multichannel analysis of shear wave splitting. *J. Geophys. Res.*, 105(B9), 21579-21590. DOI: 10.1029/2000JB900199.
- Crampin, S. (1981). A review of wave motion in anisotropic and cracked elastic-media. *Wave Motion*, 3, 343-391.
- Crampin, S. (1984). Effective anisotropic elastic constants for wave propagation through cracked solids. *Geophys. J. R. astr. Soc.*, 76, 135-145.
- Crampin, S. and Peacock, S. (2008). A review of the current understanding of seismic shear-wave splitting in the Earth's crust and common fallacies in interpretation. *Wave Motion*, 45, 675-722.
- Favier, N. and Chevrot, S. (2003). Sensitivity kernels for shear wave splitting in transverse isotropic media. *Geophys. J. Int.*, 153, 213-228.
- Gerst, A. and Savage, M.K. (2004). Seismic anisotropy beneath Ruapehu volcano: A possible eruption forecasting tool. *Science*, 306, 1543-1547.
- Kendall, M.J. et al. (2025). [Shear-wave splitting at Axial Seamount]. *Seismica*, 4(1). DOI: 10.26443/seismica.v4i1.1101.
- Silver, P.G. and Chan, W.W. (1991). Shear wave splitting and subcontinental mantle deformation. *J. Geophys. Res.*, 96(B10), 16429-16454. DOI: 10.1029/91JB00899.
- Thomsen, L. (1986). Weak elastic anisotropy. *Geophysics*, 51(10), 1954-1966. DOI: 10.1190/1.1442051.

### Crack-Induced Anisotropy and Rock Physics
- Hudson, J.A. (1981). Wave speeds and attenuation of elastic waves in material containing cracks. *Geophys. J. R. astr. Soc.*, 64, 133-150. DOI: 10.1111/j.1365-246X.1981.tb02662.x.
- Nur, A. (1971). Effects of stress on velocity anisotropy in rocks with cracks. *J. Geophys. Res.*, 76(8), 2022-2034. DOI: 10.1029/JB076i008p02022.
- O'Connell, R.J. and Budiansky, B. (1974). Seismic velocities in dry and saturated cracked solids. *J. Geophys. Res.*, 79(35), 5412-5426.
- Sayers, C.M. and Kachanov, M. (1995). Microcrack-induced elastic wave anisotropy of brittle rocks. *J. Geophys. Res.*, 100(B3), 4149-4156.
- Zatsepin, S.V. and Crampin, S. (1997). Modelling the compliance of crustal rock: I. Response of shear-wave splitting to differential stress. *Geophys. J. Int.*, 129, 477-494.
- Mavko, G., Mukerji, T., and Dvorkin, J. (2009). *The Rock Physics Handbook*, 2nd ed. Cambridge University Press.

### dv/v and Ambient Noise Interferometry
- Birch, F. (1960). The velocity of compressional waves in rocks to 10 kilobars, part 1. *J. Geophys. Res.*, 65, 1083-1102.
- Brenguier, F. et al. (2008). Towards forecasting volcanic eruptions using seismic noise. *Nature Geoscience*, 1, 126-130.
- Brenguier, F. et al. (2014). Mapping pressurized volcanic fluids from induced crustal seismic velocity drops. *Science*, 345, 80-82.
- Clarke, D. et al. (2011). Assessment of resolution and accuracy of the Moving Window Cross Spectral technique for monitoring crustal temporal variations using ambient seismic noise. *Geophys. J. Int.*, 186, 867-882.
- Donaldson, C. et al. (2017). Relative seismic velocity variations correlate with deformation at Kīlauea volcano. *Sci. Adv.*, 3, e1700219. DOI: 10.1126/sciadv.1700219.
- Lee, M.K. et al. (2024). [dv/v at Axial Seamount]. *Geophys. Res. Lett.*, 51, e2024GL108883. DOI: 10.1029/2024GL108883.
- Sens-Schönfelder, C. and Wegler, U. (2006). Passive image interferometry and seasonal variations of seismic velocities at Merapi Volcano, Indonesia. *Geophys. Res. Lett.*, 33, L21302.
- Shapiro, N.M. and Campillo, M. (2004). Emergence of broadband Rayleigh waves from correlations of the ambient seismic noise. *Geophys. Res. Lett.*, 31, L07614.
- Snieder, R. (2006). The theory of coda wave interferometry. *Pure Appl. Geophys.*, 163, 455-473.

### Joint/Integrated Approaches
- Boness, N.L. and Zoback, M.D. (2006). Mapping stress and structurally controlled crustal shear velocity anisotropy in California. *Geology*, 34(10), 825-828.
- Tod, S.R. (2001). The effects on seismic waves of interconnected nearly aligned cracks. *Geophys. J. Int.*, 146, 249-263.

### Axial Seamount Context
- Hefner, W.L. et al. (2020). *J. Geophys. Res.*, 125, e2020JB019356. DOI: 10.1029/2020JB019356.
- Nooner, S.L. and Chadwick, W.W. (2016). Inflation-predictable behavior and co-eruption deformation at Axial Seamount. *Science*, 354, 1399-1403. DOI: 10.1126/science.aah4666.
- Slead, S.R. et al. (2024). *J. Geophys. Res.*, 129, e2023JB028414. DOI: 10.1029/2023JB028414.
- Tolstoy, M. et al. (2006). A sea-floor spreading event captured by seismometers. *Science*, 314, 1920-1922. DOI: 10.1126/science.1133950.
- Wang, K. et al. (2024). *Geophys. Res. Lett.*, 51, e2024GL108631. DOI: 10.1029/2024GL108631.
- Wilcock, W.S.D. et al. (2018). The cabled observatory at Axial Seamount. *Oceanography*, 31(1), 114-123. DOI: 10.5670/oceanog.2018.117.
- Zhu, J. et al. (2026). *Geology*. DOI: 10.1130/G54254.1.

---

*This is a working document. Equations should be verified against the original publications. Some DOIs are approximate and should be confirmed through institutional library access. Papers marked with brackets [ ] in the title indicate that the exact title should be confirmed.*
