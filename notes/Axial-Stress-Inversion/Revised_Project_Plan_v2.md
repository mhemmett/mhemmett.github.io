# Joint Stress Inversion at Axial Seamount: Revised Project Plan (v2)

**Michael Hemmett — Working Document, April 2026**

---

## 1. Overview: What Changed and Why

This revision makes three fundamental changes to the project plan:

**Change 1: The stress inversion route is now through crack closure physics, not geodetic strain.** The previous plan proposed computing σ_ij = C_ijkl × ε_kl using geodetically derived strain. As discussed, this fails because the geodetic data (a handful of BPRs and acoustic ranging pairs) cannot resolve a 3D strain tensor field at depth — the strain you compute is entirely controlled by the assumed source model (Mogi, sill, etc.), not by the data. The revised plan instead recovers Δσ from Δα_ij using the stress-dependent crack closure relationship, which is constrained by the seismological observations themselves. Geodetic data are retained as an independent validation, not as an input to the inversion.

**Change 2: The Hudson second-order term is replaced by the Cheng Padé approximation.** Mavko et al. (2009, §4.10) explicitly warn that Hudson's second-order expansion "is not a uniformly converging series and predicts increasing moduli with crack density beyond the formal limit" (Cheng, 1993). At the crack densities observed at Axial (ε ~ 0.05–0.13), the second-order correction can produce artifacts. The Cheng Padé approximation avoids this problem and should be used throughout.

**Change 3: The fluid effect on C_ijkl is treated properly using Gassmann's anisotropic equations.** The previous plan set β_ijkl ≈ 0 for fluid-saturated cracks using the high-frequency Hudson result (U₃ = 0 for infinitely thin fluid-filled cracks). But at seismic frequencies (1–5 Hz), pore pressures equilibrate within cracks and the correct treatment is the low-frequency (Gassmann) limit. Following Mavko et al. (2009, §4.10 and §6.5): compute the dry-rock C_ijkl using Hudson's theory with dry crack parameters, then saturate using Gassmann's anisotropic fluid substitution equations (Gassmann, 1951; Brown and Korringa, 1975). This gives the correct low-frequency saturated stiffness tensor.

---

## 2. Revised Theoretical Framework

### 2.1 The Forward Model Chain

```
σ_ij(x,t), P_p(x,t)
    │
    ▼  [crack closure physics, §2.2]
α_ij(x,t)   (crack density tensor)
    │
    ▼  [Hudson + Padé, dry frame, §2.3]
C^dry_ijkl(x,t)
    │
    ▼  [Gassmann anisotropic fluid substitution, §2.4]
C^sat_ijkl(x,t)
    │
    ├──▶ dv/v predictions  [via sensitivity kernels, §2.5]
    ├──▶ Splitting intensity predictions  [via ray tracing, §2.6]
    └──▶ Focal mechanism predictions  [via Wallace-Bott, §2.7]
```

The inverse problem recovers α_ij(x,t) jointly from all three data types, then maps Δα_ij → Δσ using the crack closure model.

### 2.2 Stress → Crack State: The Closure Model

The crack closure stress for a penny-shaped crack with aspect ratio α in a matrix with shear modulus μ₀ and Poisson's ratio ν is (Mavko et al., 2009, §2.7, eq. from Walsh, 1965):

```
σ_close = (π / 2(1 − ν)) α μ₀
```

Note: this uses shear modulus μ₀, not Young's modulus E. For basalt with μ₀ ~ 25 GPa and ν ~ 0.25:

```
σ_close = (π / 1.5) × α × 25 GPa = 52.4 α GPa
```

So a crack with α = 10⁻³ closes at σ_close ≈ 52 MPa, and a crack with α = 10⁻⁴ closes at ≈ 5.2 MPa. At Axial, the stress changes we care about (co-eruption ~20–30 MPa, inter-eruption inflation ~1–10 MPa) affect cracks in the aspect ratio range α ~ 10⁻⁴ to 5 × 10⁻⁴.

For a crack with normal **n̂**, the effective normal stress that closes it is:

```
σ_n^eff(n̂) = σ_ij n̂_i n̂_j − P_p
```

A crack is open when σ_n^eff < σ_close(α). For a population of cracks with initial aspect ratio distribution N₀(α, n̂), the open crack density tensor at stress state σ is:

```
α_ij(σ) = (1/V) Σ_{open} B_T a_r³ n_i n_j

where "open" = { cracks r : σ_n^eff(n̂_r) < σ_close(α_r) }
```

**The inverse relationship — Δα_ij → Δσ:** A change in α_ij between times t₁ and t₂ implies that cracks have opened or closed, which requires a change in the effective stress on those crack faces. For a uniform aspect ratio distribution between 0 and α_max:

```
Δα_ij ∝ −κ_ij Δσ_ij^eff
```

where the proportionality κ depends on the crack population parameters. This is the relationship that replaces the geodetic strain bridge. The key parameter is α_max (or equivalently, the aspect ratio distribution), which controls the stress sensitivity. At Axial, α_max can be calibrated empirically using the co-eruption event where both Δ(δt) and the pressure change are independently measured.

### 2.3 Crack State → Dry-Frame Stiffness: Hudson + Padé

For a single set of aligned cracks with normals along the 3-axis, the first-order Hudson corrections to the isotropic background are (Mavko et al., 2009, §4.10):

```
c¹₁₁ = −(λ²/μ) ε U₃
c¹₁₃ = −(λ(λ+2μ)/μ) ε U₃
c¹₃₃ = −((λ+2μ)²/μ) ε U₃
c¹₄₄ = −μ ε U₁
c¹₆₆ = 0
```

where ε = Na³/V is the crack density, and U₁, U₃ depend on the crack content.

**For dry cracks:**
```
U₁ = 16(λ+2μ) / 3(3λ+4μ)
U₃ = 4(λ+2μ) / 3(λ+μ)
```

**For fluid-filled cracks (high-frequency, isolated):**
```
U₁ = [16(λ+2μ) / 3(3λ+4μ)] × 1/(1+M)
U₃ = [4(λ+2μ) / 3(λ+μ)] × 1/(1+κ)
```

where M = 4μ'(λ+2μ) / (παμ(3λ+4μ)) and κ = [K'+4μ'/3](λ+2μ) / (παμ(λ+μ)), with K', μ' being the bulk and shear moduli of the inclusion material.

**For infinitely thin, fluid-filled cracks (μ' = 0):**
```
U₁ = 16(λ+2μ) / 3(3λ+4μ)     (same as dry — fluid has no shear resistance)
U₃ = 0                         (no change in compressional modulus!)
```

**Critical implication:** For thin fluid-filled cracks at high frequency, the P-wave modulus is unaffected by crack density. Only the shear modulus changes. This is NOT the correct result for low-frequency (seismic) conditions. At seismic frequencies, pore pressure equilibrates and we must use the **dry frame + Gassmann** approach instead of the high-frequency fluid-filled Hudson result.

**The correct procedure for seismic frequencies** (Mavko et al., 2009, §4.10, Caution box):
1. Compute C^dry_ijkl using Hudson's theory with **dry crack** parameters (U₁^dry, U₃^dry).
2. Saturate using Gassmann's anisotropic equations (§2.4 below).

**Padé approximation (Cheng, 1993):** Instead of the Hudson expansion c^eff = c⁰ + c¹ + c², use:

```
c^eff_Padé = c⁰ + c¹ / (1 − c²/c¹)
```

This avoids the divergence of the second-order term at moderate crack densities and is well-behaved up to ε ~ 0.1–0.15.

**Multiple crack sets:** For three sets of cracks with normals along the 1-, 2-, and 3-axes with crack densities ε₁, ε₂, ε₃, the first-order corrections superpose linearly (Mavko et al., 2009, §4.10, p.197). The general case with arbitrary crack orientations uses the orientation distribution function ε(θ,φ) and the normalized second- and fourth-rank crack density tensors ε̃_ij and ε̃_ijpq (Mavko et al., 2009, §4.10, p.198).

### 2.4 Fluid Substitution: Gassmann's Anisotropic Equations

For a dry anisotropic rock saturated with a fluid of bulk modulus K_fl, the saturated stiffness is (Gassmann, 1951; Mavko et al., 2009, §6.5):

```
c^sat_ijkl = c^dry_ijkl + [(K₀ δ_ij − c^dry_ijαα/3)(K₀ δ_kl − c^dry_bbkl/3)] / [(K₀/K_fl)φ(K₀ − K_fl) + (K₀ − c^dry_ccdd/9)]
```

where K₀ is the mineral bulk modulus, K_fl is the fluid bulk modulus, and φ is the porosity. The repeated index notation c_ijαα means Σ_α c_ijαα. This is the anisotropic generalization of the familiar isotropic Gassmann equation.

For the inverse (computing dry from saturated), the analogous expression applies with K_fl → 0 (or using the compliance-domain Brown-Korringa equations, Mavko et al., 2009, §6.4).

**Why this matters for the plan:** The dv/v measurements sense the saturated rock, and the splitting measurements sense the saturated rock. We need C^sat_ijkl to predict the observables. But the crack closure physics acts on the dry frame (cracks open and close mechanically, and the fluid response is a separate effect). So the correct forward model is: stress → dry crack state → dry C_ijkl → Gassmann → saturated C_ijkl → observables.

**Fluid properties at Axial:** K_fl for seawater at the conditions relevant to Axial (~1500 m water depth, ~2°C bottom temperature to ~350°C at hydrothermal vents, ~15 MPa hydrostatic pressure) ranges from ~2.2 GPa (cold seawater) to ~0.5–1.0 GPa (hot hydrothermal fluid). Near the magma reservoir, supercritical fluids and magmatic volatiles (CO₂, SO₂) have much lower K_fl. This variation in K_fl is itself a source of spatial heterogeneity in C^sat_ijkl that is not crack-related, and the inversion should account for it (or at least assess its magnitude).

### 2.5 The Mavko et al. (1995) Shortcut: Empirical Stress-Velocity Mapping

An alternative to the inclusion-model approach (Hudson → crack parameters → C_ijkl) is the empirical recipe of Mavko et al. (1995), described in §2.7 of the Handbook. The procedure estimates stress-induced velocity anisotropy directly from measured isotropic V_P and V_S versus hydrostatic pressure, without assuming crack shapes or aspect ratio distributions.

The idea: if V_P(P) and V_S(P) are known as functions of hydrostatic pressure P, then the stress-dependent compliance can be written as integrals over the measured pressure-dependent compliances W'₃₃₃₃(σ_n) and W'₂₃₂₃(σ_n), where σ_n = σ₀ cos²θ is the normal stress resolved on crack faces at angle θ from the applied stress axis (Mavko et al., 2009, §2.7, p.50).

**Why this is relevant to Axial:** If laboratory velocity-pressure data exist for Axial basalts (or for representative young oceanic basalts), this method provides a way to calibrate the Δα → Δσ relationship without assuming α_max or a crack geometry model. The measured V_P(P) and V_S(P) curves encode the full crack closure spectrum empirically. This would be a powerful complement to the theoretical approach.

### 2.6 Effective Stress: The n-Value Problem

Mavko et al. (2009, §2.6) emphasize that the effective stress coefficient n (in σ^eff = P_C − nP_P) is NOT necessarily 1, is NOT necessarily the Biot-Willis coefficient α, and is DIFFERENT for different rock properties. For velocities specifically:

```
δV/V = f(δP_C − θ δP_P)
```

where θ is the effective stress coefficient for velocity, and "it is generally observed that n_velocity is sometimes close to 1, and sometimes less than one" (Mavko et al., 2009, §2.6, p.46).

**Implication:** Our pore pressure separation scheme assumed that dv/v responds to (σ − P_p), i.e., n = 1. If n < 1, then a given pore pressure change produces a smaller velocity change than the same change in confining stress. This biases the pore pressure estimates if not accounted for. The value of n should be treated as an unknown (or calibrated from laboratory data) rather than assumed to be 1.

### 2.7 Symmetry Considerations for Axial Seamount

Table 2.7.1 in Mavko et al. (2009, §2.7, p.49) shows that the resulting anisotropy symmetry depends on both the initial crack distribution and the applied stress:

- Random initial cracks + uniaxial stress → axial (VTI/HTI) symmetry, 5 elastic constants
- Random initial cracks + triaxial stress → orthorhombic symmetry, 9 elastic constants
- Axial initial cracks + stress inclined to symmetry → monoclinic, 13 constants
- Orthorhombic initial + inclined triaxial → triclinic, 21 constants

At Axial, the initial crack distribution is likely NOT random — the sheeted dike complex imposes a pre-existing fabric (axial symmetry with a horizontal axis parallel to the ridge). When the volcanic stress field (approximately radial from the magma reservoir) is superimposed on this, the combined symmetry is at best orthorhombic and possibly monoclinic.

The plan should parameterize the anisotropy appropriately — not assume HTI everywhere. For the splitting intensity tomography, the Tsvankin extended Thomsen parameters for orthorhombic media (Mavko et al., 2009, §2.4) may be needed, requiring up to 9 independent elastic constants per cell rather than 5.

---

## 3. The Revised Inverse Problem

### 3.1 What We Invert For

**Primary model parameters at each grid cell and time step:**

The model vector m contains the crack density tensor α_ij(x,t), decomposed as:

- tr(α): total crack density (1 parameter) — constrained primarily by dv/v
- Deviatoric part of α: eigenvalue differences and eigenvector orientations (up to 5 parameters, reduced to 2–3 for approximately horizontal crack systems) — constrained primarily by SI and focal mechanisms
- A static structural anisotropy field α^(structural)_ij(x): the time-independent component from the sheeted dike fabric and other geological structure — identified as the component of α that does not vary over the 10-year observation window

**Secondary parameters recovered post-inversion:**

- C^dry_ijkl(x,t) from α_ij via Hudson + Padé (§2.3)
- C^sat_ijkl(x,t) from C^dry via Gassmann (§2.4)
- Δσ_ij(x,t) from Δα_ij via the crack closure model (§2.2)

### 3.2 The Joint Linear System

The inversion solves:

```
┌              ┐       ┌           ┐
│  G_dv/v      │       │  d_dv/v   │
│  G_SI        │ × m = │  d_SI     │
│  G_FM        │       │  d_FM     │
│  λ_s L_s     │       │  0        │
│  λ_t L_t     │       │  0        │
└              ┘       └           ┘
```

**G_dv/v:** Maps α_ij to predicted dv/v. The forward operator is: α → C^dry (Hudson/Padé) → C^sat (Gassmann) → isotropic-average velocity → dv/v via sensitivity kernel. For small perturbations, this is approximately linear in α through the O'Connell-Budiansky relations: dv/v ∝ −f(ν, K_fl) × Δtr(α). The sensitivity kernel K(x; freq, station pair) determines the spatial weighting.

**G_SI:** Maps α_ij to predicted splitting intensity. The forward operator is: α → C^sat → Christoffel equation along ray → v_s1, v_s2, polarization → SI. For weak anisotropy, SI is linear in the anisotropic part of α_ij (Chevrot, 2000). The sensitivity is determined by the ray geometry through each cell.

**G_FM:** Encodes the focal mechanism constraint on α_ij eigenvectors. The stress orientations from Michael (1984) inversion constrain the eigenvectors of σ_ij, which must be consistent with the eigenvectors of α_ij (the maximum eigenvector of α is parallel to σ₃, the minimum to σ₁). This enters as a penalty term on the angular misfit between the α_ij eigenvectors and the focal-mechanism stress axes.

**Regularization:** Spatial Laplacian L_s and temporal first-difference L_t with damping parameters λ_s, λ_t chosen by L-curve analysis. Known structural boundaries (ring faults, caldera walls) can be excluded from the smoothing.

**Multi-scale temporal approach:** Coarse time steps (months) for the full three-data-type inversion. Fine time steps (days) for dv/v-only or dv/v + SI sub-inversion, anchored to the coarser joint solution.

### 3.3 From Δα_ij to Δσ: The Crack Closure Inversion

After recovering α_ij(x,t), the stress change between times t₁ and t₂ is computed at each cell using the crack closure model:

**Step 1:** Compute Δα_ij = α_ij(t₂) − α_ij(t₁).

**Step 2:** The eigenvalues of Δα_ij give the change in crack density along each principal direction. An increase in α along direction n̂ means cracks with normals along n̂ have opened, requiring a decrease in σ_n^eff along n̂.

**Step 3:** Using the crack closure stress σ_close = (π/2(1−ν)) α μ₀, the stress change along each principal direction of Δα is:

```
Δσ_n^eff(n̂_k) = −Δα_k / κ_k
```

where κ_k = (2(1−ν))/(πμ₀ α_max) × ε₀,k is the stress sensitivity along direction k, and ε₀,k is the reference crack density in that direction.

**Step 4:** Reconstruct the full Δσ_ij from its principal values and the eigenvectors of Δα_ij.

**Calibration:** The co-eruption event provides an empirical calibration point. The BPR-measured deflation (Δz ≈ 4 m), combined with the observed Δ(δt) ≈ 0.05–0.10 s, gives:

```
S_empirical = Δ(δt) / Δσ ≈ 0.003 s/MPa
```

This constrains the effective α_max without needing to assume it. As a cross-check, the Mavko et al. (1995) recipe can provide an independent estimate of the stress-velocity sensitivity if laboratory V_P(P), V_S(P) data are available for Axial basalts.

### 3.4 Geodetic Comparison (Validation, Not Input)

The geodetically derived deformation field provides an independent check on the seismologically determined stress changes:

- Fit a source model (Mogi, sill, FEM) to the BPR and acoustic ranging data at each epoch.
- Compute the predicted stress change at the surface and at depth from the source model.
- Compare with the seismologically determined Δσ_ij at the same locations and times.
- Agreement validates the inversion. Disagreement diagnoses either a wrong source model or unaccounted physics in the seismological inversion (e.g., pore pressure effects).

This comparison is valuable precisely because it is independent — the geodetic and seismological stress estimates use completely different physics and data.

---

## 4. What the Rock Physics Handbook Revealed Was Missing

### 4.1 The Frequency Regime Problem

Hudson's theory with fluid-filled crack parameters gives the **high-frequency (isolated pore)** result. At this limit, for infinitely thin fluid-filled cracks, U₃ = 0 — meaning the P-wave modulus is unaffected by cracks. This is physically because at high frequency, the fluid in each crack is isolated and resists compression independently.

At **seismic frequencies** (1–5 Hz), pore pressure has time to equilibrate between cracks and with the surrounding pore space. The correct approach (Mavko et al., 2009, §4.10, Caution) is:

1. Compute C^dry using Hudson with dry crack parameters.
2. Saturate using Gassmann's anisotropic equations.

This gives a DIFFERENT (and correct) result for the P-wave sensitivity to crack density at seismic frequencies. Specifically, the saturated P-wave modulus IS sensitive to crack density (unlike the high-frequency Hudson result), because the Gassmann fluid term couples the volumetric deformation of cracks to the bulk stiffness.

**Impact on dv/v interpretation:** If you use the high-frequency Hudson result (U₃ = 0) to relate dv/v to crack density, you will systematically underestimate the P-wave sensitivity to cracks and overestimate the stress changes needed to explain the observed dv/v. The Gassmann correction is essential.

### 4.2 The Padé Approximation vs. Second-Order Hudson

The previous plan proposed retaining the Hudson second-order term c² for crack densities ε > 0.05. Mavko et al. (2009, §4.10, Caution box) explicitly warn against this: "Better results will be obtained by using just the first-order correction rather than inappropriately using the second-order correction."

The Cheng (1993) Padé approximation resolves this:
```
c^eff = c⁰ + c¹/(1 − c²/c¹)
```

This should be used throughout the plan wherever ε > 0.05.

### 4.3 The Effective Stress Coefficient for Velocity

The plan assumed n = 1 (Terzaghi effective stress) for the pore pressure separation. Mavko et al. (2009, §2.6) show that the effective stress coefficient for velocity, θ, is often less than 1 and varies with rock type, stress state, and even the specific elastic property being measured. For cracked rocks where the crack compliance is dominated by thin cracks (the relevant regime at Axial), θ is expected to be close to 1 but may deviate, especially at high pore pressures where fluid compressibility effects become important.

This should be treated as a free parameter in the inversion (or at minimum, a sensitivity test should be run for θ = 0.8, 0.9, 1.0).

### 4.4 Structural Anisotropy and Symmetry

The previous plan assumed HTI (transverse isotropy with a horizontal symmetry axis) throughout. Mavko et al. (2009, §2.7, Table 2.7.1) show that the actual symmetry depends on the combination of initial crack fabric and applied stress. At Axial:

- The sheeted dike complex creates a pre-existing axial fabric (cracks/contacts aligned parallel to the ridge axis).
- The volcanic stress field (radial from the magma reservoir) is NOT aligned with the ridge.
- The combination produces at least **orthorhombic** symmetry (9 independent elastic constants), and possibly **monoclinic** (13 constants) where the stress axes are inclined to the dike fabric.

The inversion should allow for orthorhombic symmetry near the caldera center (where the volcanic stress dominates) and at least acknowledge the possibility of lower symmetry near the rift zones (where tectonic and volcanic stresses compete). In practice, the data may not resolve more than 5 constants (VTI/HTI) at any given cell, but the plan should not theoretically restrict itself to HTI when the physics demands more.

### 4.5 Nonlinear Elasticity as an Alternative/Complement

Mavko et al. (2009, §2.5 and §2.7) present the third-order elasticity formulation as an alternative to the crack-based models. The effective stiffness under stress is:

```
c^eff_ijkl = c⁰_ijkl + c_ijklmn ε_mn
```

where c_ijklmn are third-order elastic constants (TOECs) and ε_mn is the applied strain. For an isotropic medium, three independent TOECs (c₁₁₁, c₁₁₂, c₁₂₃) describe the full stress dependence to first order.

The stress-induced Thomsen parameters become (Sarkar et al., 2003; Prioul et al., 2004):

```
ε^(1) = ε⁰ + (c₁₅₅ / c⁰₃₃ c⁰₅₅) (σ₂₂ − σ₃₃)
γ^(1) = γ⁰ + (c₄₅₆ / 2c⁰₅₅ c⁰₅₅) (σ₂₂ − σ₃₃)
```

where c₁₅₅ = (c₁₁₁ − c₁₁₂)/4 and c₄₅₆ = (c₁₁₁ − 3c₁₁₂ + 2c₁₂₃)/8.

This provides a DIRECT relationship between stress differences and Thomsen parameter changes, without needing to go through crack density as an intermediate. If TOECs can be measured or estimated for Axial basalts, this offers a complementary (and potentially more robust) route from observables to stress. The crack-based approach and the TOEC approach should agree in the weak-anisotropy, low-crack-density limit; where they diverge provides a diagnostic for when the assumptions of one or both are breaking down.

---

## 5. Implementation Plan

### Phase 1: Data Preparation and Individual Pipelines
- SWS parameters (φ, δt) and splitting intensity SI for all catalog events (Hemmett; Kendall et al., 2025)
- dv/v time series at multiple frequency bands from ambient noise (Lee et al., 2024; Wang et al., in prep)
- Focal mechanism catalog and stress inversions in spatial/temporal bins (Zhang et al., in prep)
- Assemble 3D V_P, V_S model; compute λ, μ, ν at each grid point

### Phase 2: Forward Model Setup
- Compute dry-frame Hudson corrections (with Padé) for each cell as functions of α_ij
- Implement Gassmann anisotropic fluid substitution (§6.5 of Mavko et al.) to get C^sat
- Compute dv/v sensitivity kernels (coda wave kernels if using late coda; surface wave kernels if using early arrivals — must determine which regime applies)
- Ray-trace S-wave paths and compute SI sensitivity kernels (Chevrot, 2000; Favier and Chevrot, 2003)

### Phase 3: Joint Inversion for α_ij(x,t)
- Build stacked G matrix incorporating all three data types
- Solve by LSQR with spatial/temporal damping
- Include static structural anisotropy term α^(structural)_ij
- Multi-scale temporal stepping (monthly for full joint; daily for dv/v + SI sub-inversion)
- Compute resolution matrix R and data-space model covariance

### Phase 4: Stress Recovery
- Compute Δα_ij(x,t) relative to a reference epoch
- Map to Δσ_ij using crack closure model, calibrated against co-eruption BPR data
- Cross-check with Mavko et al. (1995) empirical recipe if lab data available
- Cross-check with third-order elasticity (TOEC) approach if constants can be estimated

### Phase 5: Validation
- Compare Δσ at the surface with geodetically derived stress changes (from Mogi/sill models fit to BPR + acoustic ranging)
- Check internal consistency: do the α_ij eigenvectors agree with focal mechanism stress orientations at every cell and time step?
- Sensitivity analysis: vary θ (effective stress coefficient), K_fl (fluid modulus), α_max (crack aspect ratio limit), and the anisotropy symmetry assumption

---

## 6. Uncertainty Framework

### 6.1 Data Covariance

Each observable has a measurement uncertainty:
- dv/v: typically ±0.01–0.05% from noise cross-correlation stability
- SI: scatter from waveform quality, typically ±0.01–0.05 s
- Focal mechanism orientations: ±10–20° on nodal planes

These define the data covariance matrix C_d in the inverse problem.

### 6.2 Model Resolution

The resolution matrix R = G(G^T C_d⁻¹ G + λ²I)⁻¹ G^T C_d⁻¹ should be computed and reported at every cell. Features in the recovered α_ij field are only interpretable where the diagonal elements of R are close to 1. Where R is small, the result is controlled by the regularization, not the data.

### 6.3 The tr(α) − P_p Trade-off

The formal model covariance matrix (G^T C_d⁻¹ G + λ²I)⁻¹ contains off-diagonal elements between tr(α) and any P_p parameter, quantifying the trade-off. This should be computed and shown to be small enough that the two effects are resolvable given the data.

### 6.4 The α_max Calibration Uncertainty

The mapping Δα → Δσ depends on α_max (or equivalently, the aspect ratio distribution). The co-eruption calibration gives one constraint on α_max. The uncertainty in α_max propagates directly into uncertainty in absolute stress magnitudes (but NOT into stress orientations or relative temporal changes, which are independent of α_max).

---

## 7. Key References

Angelier, J. (1979). *Tectonophysics*, 56, T17-T26.
Arnold, R. and Townend, J. (2007). *Geophys. J. Int.*, 170, 1336-1356.
Bott, M.H.P. (1959). *Geological Magazine*, 96, 109-117.
Brenguier, F. et al. (2008). *Nature Geoscience*, 1, 126-130.
Brown, R.J.S. and Korringa, J. (1975). *Geophysics*, 40, 608-616.
Cheng, C.H. (1993). *J. Geophys. Res.*, 98, 4209-4219.
Chevrot, S. (2000). *J. Geophys. Res.*, 105, 21579-21590.
Crampin, S. (1981). *Wave Motion*, 3, 343-391.
Donaldson, C. et al. (2017). *Sci. Adv.*, 3, e1700219.
Favier, N. and Chevrot, S. (2003). *Geophys. J. Int.*, 153, 213-228.
Gassmann, F. (1951). *Vierteljahrsschrift der Naturforschenden Gesellschaft in Zürich*, 96, 1-23.
Hardebeck, J.L. and Michael, A.J. (2006). *J. Geophys. Res.*, 111, B11310.
Hudson, J.A. (1981). *Geophys. J. R. astr. Soc.*, 64, 133-150.
Kendall, M.J. et al. (2025). *Seismica*, 4(1).
Lee, M.K. et al. (2024). *Geophys. Res. Lett.*, 51, e2024GL108883.
Martínez-Garzón, P. et al. (2014). *Geophys. Res. Lett.*, 41, 8441-8449.
Mavko, G., Mukerji, T., and Dvorkin, J. (2009). *The Rock Physics Handbook*, 2nd ed. Cambridge University Press.
Mavko, G., Mukerji, T., and Dvorkin, J. (1995). Predicting stress-induced velocity anisotropy in rocks. *Geophysics*, 60, 1081-1087.
Michael, A.J. (1984). *J. Geophys. Res.*, 89, 11517-11526.
Nooner, S.L. and Chadwick, W.W. (2016). *Science*, 354, 1399-1403.
Nur, A. (1971). *J. Geophys. Res.*, 76, 2022-2034.
O'Connell, R.J. and Budiansky, B. (1974). *J. Geophys. Res.*, 79, 5412-5426.
Prioul, R. et al. (2004). *Geophysics*, 69, 415-425.
Sarkar, D. et al. (2003). *Geophysics*, 68, 1285-1296.
Sayers, C.M. and Kachanov, M. (1995). *J. Geophys. Res.*, 100, 4149-4156.
Sens-Schönfelder, C. and Wegler, U. (2006). *Geophys. Res. Lett.*, 33, L21302.
Silver, P.G. and Chan, W.W. (1991). *J. Geophys. Res.*, 96, 16429-16454.
Thomsen, L. (1986). *Geophysics*, 51, 1954-1966.
Wallace, R.E. (1951). *Journal of Geology*, 59, 118-130.
Walsh, J.B. (1965). *J. Geophys. Res.*, 70, 381-389.
Wang, K. et al. (2024). *Geophys. Res. Lett.*, 51, e2024GL108631.
Wilcock, W.S.D. et al. (2018). *Oceanography*, 31(1), 114-123.
Zatsepin, S.V. and Crampin, S. (1997). *Geophys. J. Int.*, 129, 477-494.
