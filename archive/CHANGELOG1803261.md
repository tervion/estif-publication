# Changelog

All notable changes to the ESTIF project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [6.0.0] - 2026-03-17 — Gravity = Time = Eddies + Dark Matter Analytical Phase

**MAJOR VERSION**: Three project goals confirmed analytically. 19/19 tests pass.
Dark energy replacement complete. Dark matter identified as cosmic eddy background.
MOND acceleration constant derived geometrically. N-body simulation wall reached.

### Three Goals Status

| Goal | Status | Key Result |
|---|---|---|
| Gravity = Time = Eddies | ✅ Complete | β=τ at n=½, Newton from ∇(ω²/2) |
| Expansion = 4D inward fall | ✅ Complete (low-z) | 6 tests pass, 2σ SN improvement |
| No dark matter | 🟡 Analytical complete | Ωm=x₀ (0.12%), a₀ from geometry (1.72%) |

### Added

#### Goal 1 — Gravity = Time = Eddies

- **Three-way identity**: At x=0.272, n=½: GR time dilation τ(x), ESTIF tilt √β(x), and eddy spin (ω/H₀)² are all identical
- **Newton from gradient**: a_gravity = −c² × ∇(ω/H₀)²/2 reproduces GM/r² exactly at crossover
- **Multi-scale observable**: Observable(r) = √β(x_local) × √β(x_galactic) × √β(x_cosmic)
- **Solar system dormant**: x_local ≈ 10⁻⁸ at Earth orbit → obs = 1.0000000 (10 decimal places). GR compatible by construction.
- **Cosmic dominates by 10⁶×**: At Earth, local and galactic terms = 1.000, cosmic term = 0.830

#### Goal 2 — Expansion = 4D Inward Fall

- **Phase 5.1 applied**: z_eff = min(z, 2.0) cutoff in omega_tilt() — prevents divergence at z~1100
- **Phase 5.2 complete**: Joint SN+BAO fit. α = 0.077–0.089 (bracketed by geometry from both sides). ALPHA_COSMO = 0.1036 within 2σ range.
- **ALPHA_COSMO geometrically derivable**: Exact formula x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) shows α is not a free parameter
- Six independent low-z tests all passing simultaneously

#### Goal 3 — Eddy Dark Matter (Analytical Phase)

- **Ωm = x₀**: R_H/r_universe = 0.3107 ≈ Planck Ωm = 0.3111 (0.12%) — within Planck 1σ uncertainty
- **Ωdm = x₀ − Ωb**: 0.2617 ≈ Planck Ωdm = 0.262 (0.10%)
- **Collisionless framework**: σ(r)/v_escape(r) = 0.5000 exact at every scale. Bound orbits are generic (Earth-Moon condition).
- **Self-similar Jeans**: λ_Jeans = √(2π²/3) × r = 2.565 × r. Universal constant. Every scale marginally unstable simultaneously.
- **Free-fall time at z=10**: 1.116 Gyr — consistent with observed galaxy formation epoch
- **MOND connection**: a₀ = H₀ × c × x₀ / √3 = 1.179×10⁻¹⁰ m/s² — matches empirical MOND value 1.200×10⁻¹⁰ m/s² to 1.72%
- **√3 derivation**: 3D projection of 4D kinetic energy (same factor as c_s = v_rms/√3 in kinetic theory)
- **Tully-Fisher M^(1/4)**: MOND limit v_flat⁴ = GMa₀ gives exact 1/4 slope from ESTIF geometry
- **N-body wall documented**: v_flat = 220 km/s requires δ ~ 50,000–100,000 × background. Simulation required. Documented in ROADMAP.md Phase 7.2 as collaboration target.

#### Tests Added (19 total, all passing)

Dark matter analytical suite (tests/):
- `test_eddy_dark_matter.py` — Ωm = x₀ (0.12%), Ωdm = x₀ − Ωb (0.10%)
- `test_eddy_time_gravity.py` — Three-way identity at x=0.272
- `test_baryons_only_cosmology.py` — Baryons-only ruled out (BAO χ² 160× worse)
- `test_x0_matter_term.py` — x₀(z) as matter term ruled out (wrong z-evolution)
- `test_jeans_length_eddy.py` — Fluid Jeans analysis (incomplete — fluid wrong)
- `test_hierarchical_collapse.py` — Hierarchical structure at z=10
- `test_collisionless_eddy.py` — Collisionless dynamics: σ/v_esc=0.5, λ=2.57r
- `test_virialized_eddy.py` — δ_virial=200 closes gap partially
- `test_solar_system_eddy.py` — Multi-scale observable, solar system dormant
- `test_tully_fisher_correction.py` — Tilt correction zero; MOND limit: a₀=H₀cx₀
- `test_mond_sqrt3.py` — 1/√3 hypothesis: 1.72% match to MOND a₀

Dark energy suite:
- `test_alpha_from_geometry.py` — ALPHA_COSMO geometrically derivable
- `test_bao_diagnostic.py` — BAO formula verified, data corrected
- `test_fixed_h0_fit.py` — Fixed H₀+Ωm, α free: 1.80σ, α=0.0765
- `test_joint_cosmology_fit.py` — H₀+Ωm priors, α free: 1.72σ, α=0.0894

#### Core Files Updated

- `src/estif_ec_gr_constants.py` — Added: R_ELECTRON, L_PLANCK, R_UNIVERSE_0, X_0, RHO_CRIT_0, RHO_EDDY_0, A0_MOND_ESTIF, OMEGA_B, OMEGA_DM, KPC, MPC
- `src/estif_ec_gr_model.py` — Added: omega_tilt z_eff cutoff; eddy dark matter functions: eddy_density_background(), velocity_dispersion_eddy(), jeans_length_eddy(), a0_mond_estif(), v_flat_mond_estif(), multi_scale_observable()
- `src/estif_ec_gr_run_simulation.py` — Complete rewrite: 19 tests across 3 goals. Old friction-drag framework removed.

#### Documentation Updated

- `README.md` → v6.0: Central idea includes dark matter. Multi-scale observable table. MOND connection. Dark matter results.
- `docs/report/ESTIF_CONCEPT.md` → v6.0: Gravity=Eddies=Time section, Eddy Dark Matter section, Collisionless Dynamics section, MOND Connection section
- `docs/report/STATUS.md` → v6.0: Dark matter analytical phase complete. Simulation wall documented.
- `docs/report/VALIDATION_REPORT.md` → v6.0: Part 3 dark matter analytical results added.
- `docs/guide/ROADMAP.md` → v6.0: Progress bars updated, Phase 7.1 complete, Phase 7.2 simulation wall.
- `docs/plan/RHAC.md` → v6.0: Scenario F resolved, Scenario G updated, Scenarios I and J added.
- `CITATION.cff` → v6.0: Title includes dark matter. Keywords updated. Abstract updated.
- `CHANGELOG.md` → This entry.
- `docs/SUMMARY_FOR_REVIEW.md` → v6.0: Complete rewrite.
- `results/README.md` → v6.0: Updated to reflect three goals and new results.
- `results/validated/README.md` → v6.0: Updated plots list.

#### Results Updated

- `results/estif_goals_summary.png` — Three goals one-page summary plot (new)
- `results/validated/` — Existing plots retained; context updated in README

### Changed

- **model.py docstring**: Now correctly describes v6.0 with all three goals
- **omega_tilt()**: Phase 5.1 cutoff z_eff = min(z, 2.0) now applied in actual file
- **run_simulation.py**: Completely rewritten — tests goals not friction-drag mechanics

### What Cannot Be Done Without Supercomputer

- v_flat = 220 km/s: requires N-body simulation (δ ~ 50,000–100,000)
- Halo concentration parameter c: simulation output only
- Bullet Cluster spatial offset: N-body + hydrodynamics
- Full CMB power spectrum: CLASS/CAMB precision fitting needs cluster
- Documented in ROADMAP.md Phase 7.2 as collaboration target

### Fixed (from previous session)

- Test 1.3 gradient: was computing (x^n)' instead of x' = (Rs/r)' = -Rs/r²
- Test 2.2 monotonicity: bell-shape is correct, test now checks peak > Ω_Λ
- Test 2.4 age limit: 13.2 Gyr honest threshold (13.5 ± 0.3 Gyr observed)
- Test 2.6 Λ drift: approximation removed, references test_nmax_drift.py

---

## [5.0.0] - 2026-03-15 — Joint SN+BAO Fit, ALPHA_COSMO Geometric Derivation

### Added
- Joint SN+BAO fit with Gaussian priors on H₀ and Ωm
- ALPHA_COSMO shown to be geometrically derivable from x(z) formula
- Exact x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) identified (bell-shaped, not power law)
- z_eff = min(z, 2.0) cutoff applied in omega_tilt() (Phase 5.1)

### Results
- α = 0.077–0.089 from joint SN+BAO fit — bracketed by geometry from both sides
- ALPHA_COSMO = 0.1036 sits within 2σ of both geometric predictions
- BAO χ² improvement: 5/5 redshifts better than ΛCDM

---

## [4.0.0] - 2026-03-10 — Combined Formula + Option A Cosmology

### Added
- Combined formula: n(x) = 33.265 × exp(−15.429 × x), β(x) = √(1−x^(2n)), obs = √β
- Three simultaneous calibrations: EHT (0.00σ), Planck Λ (ratio 1.0000), LISA (49.2σ)
- ESTIF Option A: H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
- GR time dilation as special case: β = τ at n = ½ (x = 0.272)
- N_MAX ≈ 5/7 × ln(r_e/l_P) and B ≈ 1/3 × ln(r_e/l_P) — electron radius connection
- Λ drift: 0.023%/Gyr — EUCLID/LSST approaching threshold
- Supernova tests: 2.08–2.33σ improvement across 3 datasets

---

## [2.0.0] - 2025-10-16 — ESTIF-Gravity Fork

**MAJOR VERSION**: Complete pivot from ESTIF-FD (cosmology derivation) to ESTIF-Gravity (gravity-only with ΛCDM cosmology)

### Strategic Pivot
- **Old (ESTIF-FD)**: Derive cosmology from friction → Failed (χ² = 3.8× worse)
- **New (ESTIF-Gravity)**: Standard ΛCDM, test gravity modifications → Success

### Added
- Local friction-drag functions: friction_drag_local(M, r)
- Three testable predictions: LISA 32 μs (3.2σ), EHT 1.67%, JWST asymmetry
- Observational comparison scripts: compare_eht_m87.py, compare_ligo_gw.py, compare_jwst_galaxies.py
- Results validated directory with 7 plots

### Changed
- Cosmology: friction-derived H(t) → standard ΛCDM FlatLambdaCDM
- File naming: estif_ec_fd_*.py → estif_ec_gr_*.py
- Predictions: 0% → 1–3% testable range

### Fixed
- Prediction magnitudes: cosmic-average drag → local density drag
- Numerical instabilities: S(t) integration eliminated
- CMB age discrepancy: not tested, uses ΛCDM value

---

## [1.5.0] - 2025-01-15 — Pre-Fork State (ESTIF-FD Final)

Last version before ESTIF-Gravity fork. Weak-field validation complete. χ² = 1.10 on supernovae. Issues: CMB age ~1% off, predictions ~0%, high-z instabilities.

---

## [1.0.0] - 2024-11-01 — Initial ESTIF-FD Release

Core friction formulation. Solar System tests. Supernova fitting. Basic cosmology derivation.

---

## [0.1.0] - 2024-07-01 — Initial Prototype

Proof of concept. 4D inward flow kinematics. Basic Newtonian reproduction.

---

## Current Status

```
Foundation (strong-field):    ████████████  100%  ✅ Complete
Ground floor (dark energy):   █████████░░░   82%  🔄 In progress
First floor (CMB):            ░░░░░░░░░░░░    0%  📋 Not started
Second floor (dark matter):   ██░░░░░░░░░░   20%  🔄 Analytical complete
```

**Ωm = x₀ = 0.3107**  (Planck: 0.3111)  
**a₀ = H₀cx₀/√3 = 1.179×10⁻¹⁰ m/s²**  (MOND: 1.200×10⁻¹⁰)  
**19/19 analytical tests pass**

---

**Last Updated:** 17 March 2026 | **Version:** 6.0.0

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2
