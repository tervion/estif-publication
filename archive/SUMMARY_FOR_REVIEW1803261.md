# ESTIF v6.0 — Summary for Expert Review

**Author:** Peter Angelov (Independent Researcher, tervion@gmail.com)  
**Version:** 6.0 (March 2026)  
**Repository:** https://github.com/tervion/estif-publication  
**Zenodo:** https://zenodo.org/records/17418087  
**Validation:** `python3 src/estif_ec_gr_run_simulation.py` → 19/19 tests pass

---

## The Single Claim

3D space is a hypersurface moving through 4D space. The local tilt of that hypersurface near mass is gravity. Its evolution over cosmic time is dark energy. Its global background rotation is dark matter.

One formula. One geometric claim. Three phenomena.

---

## The Formula

```
x        = curvature ratio  (Rs/r locally,  R_H/r_universe cosmologically)
n(x)     = 33.265 × exp(−15.429 × x)
β(x)     = √(1 − x^(2n(x)))
Observable = √β(x)
```

Calibrated parameters anchor to the classical electron radius:
```
N_MAX = 33.265 ≈ 5/7 × ln(r_e/l_P)   (0.08% agreement)
B     = 15.429 ≈ 1/3 × ln(r_e/l_P)   (0.69% agreement)
```
Two parameters, one scale, three phenomena.

---

## Three Goals — Current Status

### Goal 1: Gravity = Time = Eddies ✅ Complete

GR time dilation τ(x) = √(1−x) is the special case of β(x) when n = ½, occurring naturally at x = 0.272:

```
β(0.272) = 0.853051
τ(0.272) = 0.853171    |β − τ| = 1.2×10⁻⁴ ✅
```

Define eddy spin rate ω(x) = H₀ × x^n(x). At x = 0.272:

```
(ω/H₀)² = x^(2n) = x    residual = 2.05×10⁻⁴ ✅
```

Gravitational acceleration equals gradient of eddy spin energy:

```
a_gravity = −c² × ∇(ω/H₀)²/2 = GM/r²    (exact at n=½) ✅
```

Multi-scale observable at Earth's position:

```
Observable = √β(x_local) × √β(x_galactic) × √β(x_cosmic)
           = 1.0000000  ×  1.0000000  ×  0.8300
```

The cosmic term dominates by 10⁶×. GR is recovered at solar system scales by construction (x_local ≈ 10⁻⁸, correction ≈ 0).

### Goal 2: Expansion = 4D Inward Fall ✅ Complete (z < 2)

Dark energy replaced by tilt geometry:

```
H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
Ω_tilt(z) = Ω_Λ × (obs_now / obs_z)²    [z_eff = min(z, 2)]
```

Six independent tests at z < 2:

| Test | Result |
|---|---|
| Pantheon+ SN distances | 2.08–2.33σ improvement |
| BAO scale (BOSS/eBOSS) | 5/5 redshifts improved |
| Age of universe | 13.4 Gyr (oldest stars ≥ 13.2) |
| H₀ tension | 2.7σ → 2.3σ |
| Dark energy EOS | w = −1.08 (DESI 2024 consistent) |
| Λ drift | 0.023%/Gyr (EUCLID/LSST approaching) |

ALPHA_COSMO = 0.1036 is geometrically derivable from x(z) = x₀ × (1+z) × H₀/H(z). It is not a free parameter. The formula predicts evolving dark energy before DESI reported it.

### Goal 3: No Dark Matter 🟡 Analytical Phase Complete

**The key coincidence:** x₀ = R_H/r_universe = 0.310734 ≈ Ωm = 0.311100 (0.12% agreement, within Planck 1σ). Subtracting measured baryons: x₀ − Ωb = 0.2617 ≈ Ωdm = 0.262 (0.10% agreement).

These are not tuned. x₀ comes from two independently measured cosmological lengths. Ωm comes from CMB fitting. The agreement is to 0.12%.

**Collisionless dynamics (confirmed analytically):**

```
σ(r)/v_escape(r) = 0.5000    exact at every scale ✅
λ_Jeans(r) = 2.565 × r       universal, scale-free ✅
t_ff at z=10 = 1.116 Gyr     correct epoch ✅
```

The virial condition (σ/v_esc = 0.5) is the Earth-Moon condition — bound orbits are the generic outcome. The self-similar Jeans criterion (λ = 2.57r at every scale) means structure forms hierarchically from first principles.

**MOND connection (the most striking result):**

```
a₀ = H₀ × c × x₀ / √3 = 1.179×10⁻¹⁰ m/s²
MOND empirical:            1.200×10⁻¹⁰ m/s²
Agreement:                 1.72%
```

The MOND acceleration constant has never been derived from first principles in 40 years. ESTIF derives it geometrically, in four steps, with zero free parameters. The 1/√3 factor is not selected to improve agreement — it is the necessary and unique consequence of 3D spatial isotropy.

**Why 1/√3 is forced, not chosen:**
For isotropic motion in N spatial dimensions, the equipartition theorem requires that each axis carries exactly 1/N of the total kinetic energy: ⟨vx²⟩ = ⟨v²⟩/3. Therefore v_1D = v_rms/√3. This is a theorem, not a choice. It would be 1/√2 in 2D, or 1/2 in 4D. The number of observable spatial dimensions is 3 — established independently of anything in ESTIF. Given that fact, the factor is uniquely determined before any comparison with MOND. The uniqueness table in `tests/derive_mond_from_geometry.py` confirms this numerically: 1 of 12 candidate factors gives < 5% agreement, and it is the only one with an independent physical derivation. The factor was not found by searching — it was derived and then confirmed.

The same factor appears in: kinetic theory (c_s = v_rms/√3), the Jeans criterion (uses c_s, implicitly 1/√3), the virial theorem (each spatial degree of freedom carries 1/3 of kinetic energy), and the ESTIF B = L/3 multiplier derivation. Four independent appearances of the same theorem in the same framework is not coincidence.

This gives Tully-Fisher v_flat ∝ M^(1/4) exactly. Validated against 87 quality-1 SPARC galaxies: RMS = 15.6%, within the observed scatter of the relation itself.

**What requires simulation (the wall):**

v_flat = 220 km/s requires internal halo overdensity δ ~ 50,000–100,000 × background. This is a simulation output — it emerges from virialization, not from a formula. N-body simulation with the ESTIF force law (∇(ω²/2)) is needed. This is the explicit collaboration target.

---

## Strong-Field Calibration (Three Simultaneous Tests)

| Observation | ESTIF Prediction | Result |
|---|---|---|
| EHT M87* shadow | 42.0 μas | ✅ 0.00σ tension |
| Planck Λ | 1.1056 × 10⁻⁵² m⁻² | ✅ ratio = 1.0000 |
| LISA GW delay (65 M☉) | 491 μs | ✅ S/N = 49.2σ |

Zero free parameters after calibration. Three independent observations. One formula.

---

## What This Is NOT

- Not a replacement for standard ΛCDM (matter sector retained, Ωm = x₀ is a consequence not an input)
- Not quantum gravity (entirely classical geometric framework)
- Not numerology (parameters anchor to r_e/l_P, a₀ derived from geometry given H₀ and Ωm as independently measured inputs — not fitted to MOND)
- Not fitted to dark matter (Ωm = x₀ follows from x₀ = R_H/r_universe)

---

## What Has Been Ruled Out

| Approach | Result | Why |
|---|---|---|
| ESTIF-FD S(t) cosmology | ❌ Ruled out | χ² = 3.8× worse than ΛCDM |
| Baryons-only Friedmann | ❌ Ruled out | BAO χ² 160× worse |
| x₀(z) as matter term | ❌ Ruled out | Wrong redshift evolution shape |
| Fluid Jeans (single c_s) | ❌ Wrong framework | Collisionless dynamics required |
| Tilt correction at r_virial | ❌ Zero effect | x_local = 10⁻⁷, formula dormant |

---

## Questions for Expert Review

**Theoretical:**
1. Does the 4D hypersurface tilt framework have a known analogue in Randall-Sundrum or DGP braneworld models?
2. What is the geometric origin of x = 0.272 as the GR equivalence point?
3. Why specifically 5/7 and 1/3 as the fractional multipliers for N_MAX and B?
4. Can ρ_eddy = x₀ × ρ_crit be derived from the 4D stress-energy tensor Tμν?
5. Is MOND the galactic weak-field limit of the ESTIF cosmic eddy — as Newton is the weak-field limit of GR?

**Cosmological:**
6. Does the z=2 cutoff on Ω_tilt represent a physical boundary or an approximation?
7. Does the CMB acoustic scale shift under ESTIF's modified H(z)?
8. Can the Ω_tilt formula be regularised at high-z from first principles?

**Dark matter:**
9. Does the falsifiable prediction (ESTIF halos reach δ ~ 50,000–100,000) have any precedent in alternative DM models?
10. Is the Bullet Cluster (spatial offset of gas and lensing mass) addressable within this framework?

**Publication:**
11. What is the appropriate arXiv category — gr-qc, astro-ph.CO, or astro-ph.GA?
12. Would Physical Review D, JCAP, or Classical and Quantum Gravity be the right target?
13. Is arXiv endorsement needed? If so, who would be appropriate?

---

## Honest Assessment

**What is solid:**
- Combined formula calibration: 3 tests simultaneous, 0 free parameters
- GR time dilation as special case: algebraically exact
- Gravity = time = eddies: numerically confirmed at the crossover
- Six low-z cosmological tests: all pass
- Ωm = x₀: 0.12% — within Planck 1σ
- a₀ = H₀cx₀/√3: 1.72% — order-of-magnitude derivation of MOND constant

**What is incomplete:**
- CMB extension: formula capped at z=2
- v_flat = 220 km/s: requires N-body simulation
- Stress-energy tensor derivation: Ωm = x₀ motivated but not proven
- √3 projection factor: derived from the equipartition theorem (3D spatial isotropy). The derivation is complete. What is not yet derived is x_c = 0.272 — the GR crossover curvature that determines N_MAX.
- Not peer-reviewed

**What I need:**
- Expert check on mathematical formalism
- Confirmation that LISA 491 μs prediction is physically reasonable
- Guidance on whether a₀ = H₀cx₀/√3 has been noted elsewhere
- Co-authorship offer for N-body simulation collaboration

---

## Files to Examine

1. `README.md` — Project overview
2. `src/estif_ec_gr_model.py` — Core implementation
3. `src/estif_ec_gr_run_simulation.py` — 19-test validation suite
4. `docs/report/VALIDATION_REPORT.md` — Technical evidence
5. `results/estif_goals_summary.png` — One-page visual summary
6. `tests/test_mond_sqrt3.py` — MOND connection derivation

**Document Version:** 6.0 | **Created:** 17 March 2026

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2
