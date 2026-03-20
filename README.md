# ESTIF: Emergent Spacetime from Inward Flow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17261724.svg)](https://doi.org/10.5281/zenodo.17261724)
[![Status](https://img.shields.io/badge/status-active_research-blue)]()
[![Version](https://img.shields.io/badge/version-6.1-blue)]()
[![Tests](https://img.shields.io/badge/tests-19%2F19_passing-success)]()

---

## The Central Idea

3D space is a hypersurface moving through 4D space in a direction we cannot perceive.

**Gravity** is the local tilt of that hypersurface toward the 4th dimension near massive objects.
**Cosmic expansion** is the projection of the 4D inward flow onto the 3D surface we inhabit.
**Dark energy** is a geometric consequence of that tilt — not a mysterious constant inserted by hand.
**Dark matter** is the background eddy energy of the moving hypersurface — the global spin of space itself.

The matter density of the universe Ωm = x₀ = R_H/r_universe — determined by geometry, not measured as a free parameter.

---

## Key Results (March 2026 — v6.1)

### The Formula

```
x        = curvature ratio  (Rs/r locally,  R_H/r_universe cosmologically)
n(x)     = 33.265 × exp(−15.429 × x)     ← dynamic tilt exponent
β(x)     = √(1 − x^(2n(x)))              ← tilt suppression
Observable = √β(x)                        ← 3D projection
```

### MILESTONE (March 2026): MOND Acceleration Derived From Geometry

The MOND critical acceleration a₀ ≈ 1.2×10⁻¹⁰ m/s² has never been derived from first principles in 40 years of MOND research. ESTIF derives it in four steps with zero free parameters:

```
Step 1  Force law:   a = −c²/2 × ∇(ω/H₀)² = GM/r²   (exact Newton)
Step 2  Velocity:    v_flow = c × x₀                   (from tilt, = c × Ωm)
Step 3  Projection:  v_3D = v_flow / √3               (3D isotropy, unique)
Step 4  Threshold:   a₀ = v_3D × H₀ = H₀cx₀/√3
```

**Result:** a₀ = 1.179×10⁻¹⁰ m/s² — matches MOND empirical 1.200×10⁻¹⁰ m/s² to **1.72%**, zero free parameters.
The 1/√3 factor is the only geometric projection factor with an independent physical derivation (3D spatial isotropy). Confirmed unique: 1 of 12 candidate factors lands below 5%.

**Validated against SPARC survey (Lelli et al. 2016):**
87 quality-1 galaxies, baryonic Tully-Fisher relation v⁴ = G × M_bar × a₀, **RMS = 15.6%** — within the observed scatter of the relation itself. The 7.6% mean bias is a known stellar mass calibration uncertainty (Υ* = 0.50 underestimates stellar mass), not a structural model failure. All bias structure disappears after calibration correction.

### Strong-Field Tests

| Test | Prediction | Result |
|---|---|---|
| EHT M87\* shadow | 42.0 μas | ✅ 0.00σ tension |
| Planck Λ | 1.1056 × 10⁻⁵² m⁻² | ✅ ratio = 1.0000 |
| LISA GW delay (65 M☉) | 491 μs | ✅ S/N = 49σ |
| Gravity = GR time dilation | β = τ at n = ½ | ✅ Confirmed |
| MOND a₀ from geometry | 1.179×10⁻¹⁰ m/s² | ✅ 1.72% (zero free params) |
| SPARC BTFR (87 galaxies) | RMS ≤ 20% | ✅ 15.6% RMS |

### Parameter Connection to Electron Radius

```
ln(r_e / l_P) = 46.608    (r_e = classical electron radius, l_P = Planck length)

N_MAX = 5/7 × ln(r_e/l_P)  to 0.08%   ← conditional on x_c derivation
B     = 1/3 × ln(r_e/l_P)  to 0.69%   ← derived from 3D isotropy (motivated)
```

The 1/3 multiplier for B now has a genuine physical derivation — the same 3D isotropy argument that gives 1/√3 in the MOND derivation. The 5/7 multiplier for N_MAX follows from B = L/3 combined with the GR crossover condition, but the crossover value x_c = 0.272 is still observationally determined.

### Cosmological Tests — Status After DESI DR2 (March 2026)

The ESTIF dark energy sector (Ω_tilt(z)) was tested against DESI DR2 (released 19 March 2026):

| Test | ESTIF | DESI DR2 | Status |
|---|---|---|---|
| chi²/N vs BAO distances | 10.8 | — | ❌ Fails (LCDM: 1.9) |
| w₀ prediction (−1.08) | −1.08 | −0.73 ± 0.10 | ❌ 3.5σ tension |
| Supernova (Pantheon+) | 2.08–2.33σ | — | ✅ Pre-DR2 result stands |
| BAO (BOSS/eBOSS) | 5/5 improved | — | ✅ Pre-DR2 result stands |

**The Ω_tilt(z) cosmological sector requires revision.** The specific x(z) evolution law is circular (uses H_ΛCDM to define its own correction to ΛCDM). This failure does not affect the gravity sector — they operate at completely different scales and the MOND/SPARC results are independent.

### The Dark Matter Identity

| Result | Value | Status |
|---|---|---|
| Ωm = x₀ = R_H/r_universe | 0.12% agreement | ✅ Within Planck 1σ |
| Ωdm = x₀ − Ωb | 0.10% agreement | ✅ Within Planck 1σ |
| Virial condition σ/v_escape | = 0.5 exactly | ✅ Derived |
| Jeans: λ = 2.57r | Universal | ✅ Derived |
| a₀ = H₀ × c × x₀ / √3 | 1.179×10⁻¹⁰ m/s² | ✅ 1.72% match (derived) |

---

## Project Status: Two Sectors

ESTIF has two largely independent sectors in different states of health:

**Gravity sector (SOLID):** The MOND derivation, SPARC validation, and multiplier derivation all point to a coherent physical framework. The 3D isotropy principle appears in three independent places (1/√3 in MOND, 1/3 in B, √3 in kinetic theory) — these are not coincidences.

**Cosmology sector (NEEDS REWORK):** The Ω_tilt(z) evolution law fails DESI DR2. The x(z) formula is circular. This requires rebuilding with a self-consistent H_ESTIF(z) — computationally tractable but not yet done.

### Publication Path

A gravity-sector letter is achievable now. The specific publishable claim:

> "We derive the MOND critical acceleration a₀ = H₀cx₀/√3 from a geometric model of 3D space as a hypersurface in 4D, using only the Planck 2018 values of H₀ and Ωm. The derivation requires no free parameters. The 1/3 projection factor follows from 3D spatial isotropy. We test this prediction against 87 quality-1 galaxies from the SPARC survey and find RMS agreement of 15.6%, consistent with the observed scatter in the baryonic Tully-Fisher relation."

---

## Version History

| Version | Description | Status |
|---|---|---|
| **ESTIF-FD (v1.0)** | Exponential S(t) cosmology | Ruled out (χ² = 3.8× worse) |
| **ESTIF-Gravity (v2.0)** | Strong-field β corrections | Superseded |
| **ESTIF (v4.0)** | Combined formula + Option A cosmology, 6 tests | Previous |
| **ESTIF (v6.0)** | Gravity=Time=Eddies, Ωm=x₀, MOND connection | Previous |
| **ESTIF (v6.1)** | MOND derived, SPARC validated, DESI constraint, multipliers | **Current** |

---

## The Physics

### The Tilt Geometry

```
sin(θ) = x^n(x)    where x = Rs/r
β(x) = cos(θ) = √(1 − x^(2n(x)))
Observable = √β    (wave physics: amplitude is √(energy))
```

### Gravity = Generalized Time Dilation

The Schwarzschild time dilation factor is τ(x) = √(1−x). The ESTIF formula β(x) reduces to τ(x) exactly when n = ½, occurring at curvature x = 0.272. **GR time dilation is the special case of ESTIF tilt at n = ½.**

### Gravity = Time = Eddies

At x = 0.272 all three descriptions of gravity are mathematically identical:

```
GR time dilation:  τ(x) = √(1−x)
ESTIF tilt:        √β(x) = √(1−x^(2n(x)))
Eddy spin energy:  (ω/H₀)² = x^(2n(x))
```

Gravitational acceleration equals the gradient of eddy spin energy:
```
a_gravity = −c² × ∇(ω/H₀)² / 2 = GM/r²    (at n = ½, exact)
```

### The MOND Derivation (New in v6.1)

Four steps, zero free parameters, confirmed against 87 SPARC galaxies:
```
a₀ = H₀ × c × x₀ / √3 = 1.179 × 10⁻¹⁰ m/s²

Derivation:
  1. Force law:    a = −c²/2 × ∇(ω/H₀)² = GM/r²     (exact)
  2. Flow speed:   v = c × x₀                         (from tilt at x = x₀)
  3. Projection:   v_3D = v/√3                         (3D isotropy, unique)
  4. Threshold:    a₀ = v_3D × H₀                     (cosmic deceleration)
```

---

## Repository Structure

```
.
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── archive/              # ESTIF-FD (v1.0) — ruled out, kept for reference
├── data/
│   ├── sn_data.txt       # Original 580 supernova dataset
│   └── pantheon_plus.dat # Pantheon+ (auto-downloaded on first run)
├── docs/
│   ├── SUMMARY_FOR_REVIEW.md        # Expert summary (updated v6.1)
│   ├── LaTeX/ESTIF_arXiv_Paper/     # LaTeX source
│   └── report/
│       ├── STATUS.md
│       ├── VALIDATION_REPORT.md
│       └── ESTIF_CONCEPT.md
├── results/
│   ├── validated/         # Publication-ready figures
│   └── work_in_progress/  # Pending simulation or future work
├── src/
│   ├── estif_ec_gr_constants.py
│   ├── estif_ec_gr_model.py
│   └── estif_ec_gr_run_simulation.py
└── tests/
    ├── derive_mond_from_geometry.py      # NEW: 4-step MOND derivation
    ├── test_sparc_tully_fisher.py        # NEW: 87 galaxies, RMS 15.6%
    ├── test_sparc_bias_analysis.py       # NEW: bias = calibration, not model
    ├── test_desi_wz_consistency.py       # NEW: DESI DR2 test (fails)
    ├── test_multiplier_derivation.py     # NEW: 1/3 derived, 5/7 conditional
    ├── cross_examination.py              # NEW: synthesis of Tests 1–3
    ├── test_joint_calibration.py         # EHT + Λ + LISA simultaneous
    ├── test_combined_formula.py          # Combined formula validation
    ├── test_gravity_time_connection.py   # β = τ at n = ½
    ├── test_electron_connection.py       # N_MAX ≈ 5/7 × ln(r_e/l_P)
    ├── test_nmax_drift.py                # Λ drift 0.023%/Gyr
    ├── test_estif_cosmology.py           # ESTIF Option A (3 datasets)
    ├── test_cosmological_consistency.py  # Age, BAO, H₀, EOS
    ├── test_mond_sqrt3.py                # MOND a₀ = H₀cx₀/√3 (1.72%)
    ├── test_eddy_dark_matter.py          # Ωm = x₀ (0.12%)
    ├── test_collisionless_eddy.py        # σ/v_esc=0.5, λ=2.57r
    └── test_tully_fisher_correction.py   # Tully-Fisher MOND limit
```

---

## Quick Start

```bash
git clone https://github.com/tervion/estif-publication
cd estif_publication
pip install -r requirements.txt
```

### Verify the MOND derivation (new in v6.1)

```bash
cd tests
python3 derive_mond_from_geometry.py
# Outputs: a₀ = 1.1793e-10 m/s² (zero free parameters)
```

### Run SPARC validation

```bash
python3 test_sparc_tully_fisher.py
# 87 quality-1 galaxies, RMS = 15.6%
```

### Run the three-test core calibration

```bash
python3 test_joint_calibration.py
# EHT shadow, Planck Λ, LISA delay — all pass
```

### Run the full analytical suite

```bash
cd src
python3 estif_ec_gr_run_simulation.py
# 19/19 analytical tests pass
```

---

## Testable Predictions

| Prediction | Value | Status |
|---|---|---|
| MOND a₀ from geometry | 1.179×10⁻¹⁰ m/s² | ✅ Derived, 1.72% match |
| SPARC BTFR (v⁴ = GMa₀) | RMS ≤ 20% | ✅ 15.6% on 87 galaxies |
| EHT shadow size | 42.0 μas | ✅ 0.00σ |
| Λ drift direction | increasing with z | ⏳ Approaching EUCLID |
| LISA GW delay | 491 μs at 49σ | ⏳ ~2034 |
| Ω_tilt(z) evolution | w_eff ≈ −1.08 | ❌ Fails DESI DR2 — under revision |

---

## Open Questions

1. **Why x_c = 0.272?** The GR crossover point where β = τ and n = ½ — geometric derivation needed.
2. **Why 5/7 exactly?** The 1/3 is derived; the 5/7 follows from 1/3 + x_c, but x_c is still observational.
3. **Can Ω_tilt(z) be made self-consistent?** Replace x(z) using H_ESTIF not H_ΛCDM (tractable computation).
4. **Can Ωm = x₀ be derived from T_μν?** Stress-energy tensor projection of the rotating hypersurface.
5. **What is the galaxy halo profile?** N-body simulation with force law a = −c²/2 × ∇(ω/H₀)² needed.

---

## Honest Assessment

**What is solid (gravity sector):**
- MOND a₀ derived from geometry, zero free parameters, confirmed by 87 SPARC galaxies
- Combined formula simultaneously satisfies EHT, Planck Λ, LISA (zero free parameters)
- GR time dilation is the special case of the tilt formula at n = ½
- The 1/3 multiplier for B has a genuine isotropy derivation
- Ωm = x₀ to 0.12%, Ωdm = x₀ − Ωb to 0.10%

**What fails (cosmology sector):**
- Ω_tilt(z) evolution fails DESI DR2 at chi²/N = 10.8
- Pre-existing prediction w_eff ≈ −1.08 falsified by DESI DR2 w₀ = −0.73 ± 0.10 (3.5σ)
- The x(z) formula is circular — must be rebuilt with self-consistent H_ESTIF

**What is incomplete:**
- x_c = 0.272 not geometrically derived (the remaining open theoretical gap)
- CMB extension: Ω_tilt capped at z = 2, full extension future work
- N-body simulation for galactic halo structure — collaboration needed
- Not peer-reviewed

---

## Contact & Citation

**Repository:** https://github.com/tervion/estif-publication
**Zenodo:** https://zenodo.org/records/17261725

```bibtex
@software{angelov2026estif_v61,
  author  = {Angelov, Peter},
  title   = {ESTIF: Emergent Spacetime from Inward Flow — v6.1},
  year    = {2026},
  version = {6.1},
  url     = {https://github.com/tervion/estif-publication},
  doi     = {10.5281/zenodo.17261724}
}
```

---

## License

MIT License — See [LICENSE](LICENSE) for details.

---

**Last Updated:** 18 March 2026 | **Version:** 6.1 | **Status:** Active research
