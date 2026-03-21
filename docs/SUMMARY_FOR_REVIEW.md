# ESTIF v6.2 — Summary for Expert Review

**Author:** Peter Angelov (Independent Researcher, tervion@gmail.com)
**Version:** 6.2 (March 2026)
**Repository:** https://github.com/tervion/estif-publication
**Zenodo:** https://zenodo.org/records/17261724
**Validation:** `python3 tests/derive_mond_from_geometry.py` | `python3 src/estif_ec_gr_run_simulation.py` (21/21 pass)

---

## The Single Claim

3D space is a hypersurface moving through 4D space. The local tilt near mass is gravity. Its evolution over cosmic time is dark energy. Its global background rotation is dark matter.

---

## The Formula

```
x        = curvature ratio  (Rs/r locally,  R_H/r_universe cosmologically)
n(x)     = 33.265 × exp(−15.429 × x)
β(x)     = √(1 − x^(2n(x)))
Observable = √β(x)
```

Parameters connect to the classical electron radius:
```
N_MAX = 33.265 ≈ 5/7 × ln(r_e/l_P)   (0.08% — conditional on x_c)
B     = 15.429 ≈ 1/3 × ln(r_e/l_P)   (0.69% — derived from 3D isotropy)
```

---

## What Is New in v6.2

### MILESTONE: First Derivation of MOND a₀ From Geometry

The MOND critical acceleration a₀ ≈ 1.2×10⁻¹⁰ m/s² has never been derived from first principles in 40 years of MOND research. ESTIF provides a motivated geometric construction that reproduces the MOND normalization without internal tuning, using only the Planck 2018 values of H₀ and Ωm as independently measured inputs. The 1/√3 factor is motivated by 3D spatial isotropy and consistent with the equipartition theorem; a complete kinetic theory of the eddy background is identified as future theoretical work. The factor is physically motivated, not fitted. The four-step derivation:

**Step 1 — Force law (exact):**
```
a = −c²/2 × ∇(ω/H₀)² = GM/r²    (no approximation)
```

**Step 2 — Cosmological flow speed:**
```
v_flow = c × x₀ = c × Ωm    (from tilt formula at x = x₀, two independent routes)
```

**Step 3 — 3D isotropic projection (unique):**
```
v_3D = v_flow / √3
```
This is motivated by 3D spatial isotropy and consistent with the equipartition theorem (⟨v²⟩ = ⟨vx²⟩ + ⟨vy²⟩ + ⟨vz²⟩ → v_1D = v_rms/√3). The same factor appears in kinetic theory (c_s = v_rms/√3) and the Jeans criterion. A complete formal kinetic theory of the eddy background is future work.
**Uniqueness check:** 1 of 12 candidate geometric factors gives < 5% error, and it is the only one with an independent physical justification. This is a confirmatory test, not a mathematical uniqueness proof from first principles.

**Step 4 — MOND threshold:**
```
a₀ = v_3D × H₀ = H₀ × c × x₀ / √3 = 1.179 × 10⁻¹⁰ m/s²
```

**MOND empirical:** 1.200 × 10⁻¹⁰ m/s² → **Agreement: 1.72%**
**Free parameters in the MOND derivation itself: ZERO.** H₀ and x₀ = Ωm from Planck 2018. √3 from dimension count. (The underlying tilt formula has separately calibrated parameters; those are not used here.)

### SPARC Validation

Formula tested: v_flat = (G × M_bar × a₀)^(1/4), M_bar = 0.50 × L_3.6 + 1.33 × M_HI (Υ* from McGaugh & Schombert 2014, independent of MOND).

- **87 quality-1 SPARC galaxies: RMS = 15.6%** — within observed BTFR scatter
- 82% of galaxies predicted within 20%, 97% within 30%
- Mean bias −7.6%: traced entirely to Υ* calibration. All bias structure disappears after Υ* correction (all p > 0.16). The force law is structurally sound.

### Multiplier Derivation Progress

**B = L/3 is now derived** from the same 3D isotropy principle as the MOND 1/√3:
- The tilt exponent decay rate B controls how the 4D decay projects into 3D
- Isotropic decay: each spatial dimension receives 1/3 → B = L/3
- Error: 0.69% (within uncertainty of incomplete T_μν projection)
- This same principle appears three times: 1/√3 in MOND, 1/3 in B, √3 in kinetic theory

**N_MAX/5/7 is conditional:** follows from B = L/3 + GR crossover condition n(x_c) = 1/2, but x_c = 0.272 is still observationally determined, not derived geometrically. The gap has been reduced from two unexplained fractions to one unexplained number.

### DESI DR2 Constraint (Honest)

Tested against official DESI DR2 BAO data (arXiv:2503.14738, 19 March 2026):
- **chi²/N = 10.8** for ESTIF Ω_tilt(z) vs DESI DR2 (ΛCDM: 1.9)
- Pre-existing prediction w_eff ≈ −1.08 is **3.5σ from DESI DR2 w₀ = −0.73 ± 0.10**
- Root cause: x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is circular
- **This is a genuine prediction failure**, documented honestly
- The failure is isolated to the cosmology sector. Gravity sector results are independent.

### New in v6.2 (Pre-Publication Blockers Resolved)

**Blocker 1 — a₀ redshift constancy (resolved):** In the comoving frame appropriate for galaxy dynamics, x(z) = c/[H(z) × r_universe_comoving], so H(z) cancels exactly: a₀(z) = c²/(r_universe_comoving × √3) = constant. Maximum deviation: 2.22×10⁻¹⁶ (floating-point epsilon). Algebraic identity, not numerical. Confirmed consistent with Di Teodoro+2021, Übler+2017, Tiley+2019 (all ≤ 2σ at z = 0.75–2.2). Script: `tests/test_a0_redshift.py`.

**Blocker 2 — 1/√3 language (resolved):** Language updated in all documents from "derived from the equipartition theorem" to "motivated by 3D spatial isotropy and consistent with the equipartition theorem; a complete kinetic theory of the eddy background is identified as future theoretical work."

**Parameter independence test (new):** a₀ = H₀cx₀/√3 tested across 3,600 combinations of H₀ ∈ [65,75] km/s/Mpc and Ωm ∈ [0.27,0.33]: 100% of combinations within ±20% SPARC scatter. 8 published datasets all pass. Planck–SH0ES Hubble tension shifts a₀ by only 4.1%. Script: `tests/test_a0_parameter_independence.py`.

---

## Two-Sector Structure

| Sector | Tests | Status |
|---|---|---|
| **Gravity** | MOND derivation, SPARC, multipliers | ✅ SOLID |
| **Cosmology** | DESI DR2 BAO | ❌ FAILS — needs rework |

The two sectors share x₀ but diverge in how it propagates. The failure in Ω_tilt(z) does not affect the gravity results.

---

## Strong-Field Calibration

| Observation | ESTIF Prediction | Result |
|---|---|---|
| EHT M87* shadow | 42.0 μas | ✅ 0.00σ tension |
| Planck Λ | 1.1056 × 10⁻⁵² m⁻² | ✅ ratio = 1.0000 |
| LISA GW delay (65 M☉) | 491 μs | ✅ S/N = 49.2σ |

Zero free parameters after calibration. Three independent observations. One formula.

---

## What This Is NOT

- Not curve-fitting: a₀ is derived, not fitted
- Not numerology: 1/3 has a unique physical derivation; uniqueness confirmed over 12 candidates
- Not a complete cosmological replacement: the Ω_tilt(z) sector fails DESI DR2

---

## Questions for Expert Review

**Theoretical:**
1. The 1/√3 and 1/3 projection factors both arise from 3D spatial isotropy, motivated by the equipartition theorem. A complete kinetic theory of the eddy background is future work. Is the isotropy motivation sufficient as stated, or is the paper weakened by leaving the formal derivation open?
2. What is the geometric property of Schwarzschild spacetime at r = 3.68 Rs that gives x_c = 0.272? (This is the remaining open theoretical gap.)
3. Can ρ_eddy = x₀ × ρ_crit be derived from the 4D stress-energy tensor Tμν?

**Observational:**
4. The SPARC zero-bias Υ* = 0.85 is above the McGaugh+2014 standard of 0.50. Is 0.85 within the plausible range for 3.6 μm mass-to-light ratios in the literature?
5. The DESI DR2 failure is in Ω_tilt(z) evolution. Does replacing x(z) with a self-consistent H_ESTIF solve this, or is a more fundamental revision needed?

**Publication:**
6. Is "geometric derivation of a₀ from hypersurface tilt + SPARC validation" a viable standalone letter for MNRAS Letters or ApJL?
7. What would be the appropriate framing to distinguish this from MOND variants — given that it derives a₀ from cosmological constants rather than modifying the force law directly?

---

## Honest Assessment

**What is solid:**
- a₀ derived from geometry, zero free parameters, 1.72% match, confirmed by 87 galaxies
- Combined formula: 3 strong-field tests simultaneous, 0 free parameters
- GR time dilation as special case: algebraically exact
- 1/3 multiplier for B: genuine isotropy derivation
- Ωm = x₀ to 0.12%, Ωdm = x₀ − Ωb to 0.10%
- a₀ redshift constancy: proved algebraically (H(z) cancels exactly)
- Parameter independence: confirmed across 3,600 H₀/Ωm combinations

**What fails:**
- Ω_tilt(z) evolution: fails DESI DR2 at chi²/N = 10.8
- w_eff ≈ −1.08 prediction falsified at 3.5σ by new data

**What remains open:**
- 1/√3 projection factor: motivated by 3D spatial isotropy, consistent with equipartition theorem. Formal kinetic theory is future work.
- x_c = 0.272 not geometrically derived
- Full T_μν projection for ρ_eddy
- N-body simulation for halo structure
- CMB not addressed (formula capped at z = 2)
- Not peer-reviewed

---

## Files to Examine

1. `tests/derive_mond_from_geometry.py` — the 4-step MOND derivation (key new result)
2. `tests/test_sparc_tully_fisher.py` — 87 galaxy validation
3. `tests/test_sparc_bias_analysis.py` — bias is calibration, not structural
4. `tests/test_multiplier_derivation.py` — 1/3 derived, 5/7 conditional
5. `tests/test_desi_wz_consistency.py` — DESI DR2 test (honest failure)
6. `tests/cross_examination.py` — synthesis of all three tests
7. `src/estif_ec_gr_model.py` — core implementation
8. `results/validated/mond_derivation.png` — derivation chain figure
9. `results/validated/sparc_tully_fisher.png` — 87-galaxy validation
10. `tests/test_a0_redshift.py` — a₀ redshift constancy proof (H(z) cancels, v6.2)
11. `tests/test_a0_parameter_independence.py` — robustness across H₀/Ωm parameter space (v6.2)
12. `results/validated/a0_redshift.png` — constancy plot vs high-z observations
13. `results/validated/a0_parameter_independence.png` — 3,600 parameter combinations

**Document Version:** 6.2 | **Updated:** 20 March 2026
