# ESTIF Development Roadmap

**Version:** 6.2  
**Last Updated:** 20 March 2026  
**Status:** Gravity letter ready for submission. Cosmology sector under revision.

---

## Mission Statement

Develop a complete geometric model of gravity and cosmology from first principles —
deriving all observable consequences from a single physical claim: that 3D space is
a hypersurface moving through 4D space.

The goal is not to adjust ΛCDM. The goal is to replace it.

---

## Progress Overview

```
Foundation (strong-field gravity):    ████████████  100%  ✅ Complete
Ground floor (dark energy / Λ):       █████████░░░   75%  🔄 In progress
First floor (CMB / early universe):   ░░░░░░░░░░░░    0%  📋 Not started
Second floor (dark matter):           ██░░░░░░░░░░   20%  🔄 Analytical phase complete
```

---

## What Has Been Achieved (16 March 2026)

Before looking forward, the full record of what stands:

### Strong-Field Gravity — Complete

The combined formula:
```
n(x)     = 33.265 × exp(−15.429 × x)
β(x)     = √(1 − x^(2n(x)))
Observable = √β(x)
```

simultaneously satisfies three independent observations with zero free parameters
after calibration:

| Test | Result |
|---|---|
| EHT M87\* shadow | 0.00σ tension |
| Planck Λ | ratio = 1.0000 |
| LISA GW delay (65 M☉) | 49.2σ S/N |

Additional confirmed results:
- GR time dilation is the special case β = τ at n = ½ (x = 0.272)
- N_MAX ≈ 5/7 × ln(r_e/l_P) and B ≈ 1/3 × ln(r_e/l_P) to within 0.7%
- Λ drifts at 0.023%/Gyr — approaching EUCLID/LSST detection threshold

### Gravity = Time = Eddies — Confirmed Analytically

The three descriptions of gravity — GR time dilation, ESTIF tilt, and 4D eddy
spin — are the same phenomenon at different scales. At x = 0.272 all three
become identical. The gravitational acceleration equals the spatial gradient of
the eddy spin energy. The formula is correctly dormant at solar system scales
(GR compatible) and active at galactic and cosmic scales.

The multi-scale observable formula:
```
Observable(r) = √β(x_local) × √β(x_galactic) × √β(x_cosmic)
```
At Earth's position: local and galactic terms = 1.000000, cosmic term = 0.830.
Dark matter IS the cosmic eddy — not the solar or galactic eddies.

### Dark Energy — Partial

ESTIF Option A replaces ΩΛ with Ω_tilt(z):
```
H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
Ω_tilt(z) = Ω_Λ × (obs_now / obs_z)²
```

Six low-redshift tests pass simultaneously:

| Test | Result |
|---|---|
| Supernova distances (Pantheon+) | 2.08–2.33σ improvement |
| Age of universe | 13.63 Gyr ✅ (oldest stars ≥ 13.5 Gyr) |
| BAO scale | 5/5 redshifts improved |
| H₀ tension | 2.7σ → 2.3σ |
| Dark energy EOS | w = −1.08 (DESI 2024 consistent) |
| Λ drift | 0.023%/Gyr (EUCLID/LSST testable) |

**What remains on the ground floor:**
The Ω_tilt formula diverges at high redshift. A physically motivated cutoff
is applied at z = 2 (Phase 5.1 complete). CMB extension is Phase 6.

### Dark Matter — Analytical Phase Complete, Simulation Required

Five analytical results confirmed:

| Result | Status |
|---|---|
| Ωm = x₀ = R_H/r_universe (0.12% agreement) | ✅ Confirmed |
| Ωdm = x₀ − Ωb (0.10% agreement) | ✅ Confirmed |
| Virial condition σ/v_escape = 0.5 (exact) | ✅ Confirmed |
| Self-similar Jeans: λ_Jeans = 2.57 × r (universal) | ✅ Confirmed |
| Free-fall time at z=10: ~1 Gyr (correct epoch) | ✅ Confirmed |
| Formula dormant at solar system scales | ✅ Confirmed (GR compatible) |

**The wall hit:** v_flat = 220 km/s requires internal halo overdensity
δ ~ 50,000–100,000 × ρ_eddy. This is a simulation output. The Tully-Fisher
exponent is 1/3 (ESTIF) vs 1/4 (observed) — one analytical test remains.
Beyond that, N-body simulation is required.

---

## Phase 5: Complete The Ground Floor (Dark Energy) 🔄

**Goal:** Bring dark energy from 75% to 100%  
**Prerequisite for:** CMB work (Phase 6)

---

### 5.1 — Regularise Ω_tilt At High Redshift ✅ DONE

**Status:** Hard cutoff applied at z_eff = min(z, 2.0). Model is well-behaved.
All existing SN, BAO, age tests verified unchanged after the cutoff.

---

### 5.2 — Joint SN + BAO + H₀ Fit ✅ DONE

**Results:** α = 0.077–0.089 (bracketed by geometry from both sides).
ALPHA_COSMO = 0.1036 sits within the 2σ range. H₀ and Ωm Gaussian priors
required (SN dominates χ² by 263:1 over BAO). α is geometrically derivable
from x(z) = x₀ × (1+z) × H₀/H(z) — not a free parameter.

---

### 5.3 — DESI DR2 Comparison

**The opportunity:**
DESI DR2 (2024) reported hints of evolving dark energy — w(z) ≠ −1.
ESTIF predicts w ≈ −1.08 at all z < 2. This is a direct, testable
comparison with published data.

**What to do:**
Write `test_desi_comparison.py` that:
1. Downloads or reads the DESI DR2 w(z) measurements
2. Plots ESTIF's predicted w(z) against DESI's observed w(z)
3. Computes χ² for the ESTIF prediction against DESI

**Milestone:** ESTIF w(z) comparison with DESI published, significance reported

---

## Phase 6: The First Floor (CMB / Early Universe) 📋

**Goal:** Extend ESTIF to recombination (z ~ 1100)  
**Prerequisite:** Phase 5 complete (Ω_tilt regularised)  
**Difficulty:** High — requires external tools and possibly collaborators

---

### 6.1 — Understand What CMB Constrains

Before writing any code, understand exactly what the CMB power
spectrum constrains and what ESTIF needs to reproduce:

**The three CMB constraints that matter:**

1. **Sound horizon at recombination** (r_s ~ 147 Mpc)
   Set by baryon-photon plasma oscillations before decoupling.
   Depends on Ωb, Ωm, and the expansion history H(z) at z ~ 1000.
   ESTIF inherits Ωb and Ωm from ΛCDM. The question is whether
   Ω_tilt at z ~ 1000 shifts r_s significantly.

2. **Angular scale of CMB peaks** (l ~ 200 for first peak)
   θ_s = r_s / D_A(z_rec) where D_A is the angular diameter distance
   to recombination. ESTIF modifies D_A via the modified H(z).
   Need to check if the shift in D_A is within Planck's 0.1% precision.

3. **Integrated Sachs-Wolfe effect at late times**
   Dark energy affects CMB at large angles (l < 20) via the ISW effect.
   ESTIF's evolving Ω_tilt predicts a different ISW than ΛCDM.
   This is potentially a unique ESTIF signature in the CMB.

**What to do:**
Write `test_cmb_angle_estimate.py` — a simplified calculation of
the CMB angular scale θ_s under ESTIF. If θ_s differs from ΛCDM
by less than 0.5%, ESTIF is not ruled out by Planck.

**Milestone:** Confirm ESTIF doesn't catastrophically shift CMB peaks

---

### 6.2 — The ISW Effect as a CMB Prediction

**The opportunity:**
The late-time Integrated Sachs-Wolfe (ISW) effect imprints a
specific pattern on CMB temperature at large angular scales.
It is sensitive to exactly the kind of dark energy evolution ESTIF predicts.

ESTIF predicts Ω_tilt increases with z — this is a different ISW
signature than ΛCDM's constant ΩΛ. This could be a genuine,
unique prediction testable against Planck data right now.

**What to do:**
Write `test_isw_prediction.py` that computes the ISW power
spectrum contribution from Ω_tilt(z) and compares to Planck
measurements at low multipoles.

**Milestone:** ESTIF ISW prediction plotted and compared to Planck data

---

### 6.3 — Full CMB Power Spectrum (Requires External Tools)

**The challenge:**
A complete CMB power spectrum calculation requires solving the
Boltzmann equations for photon-baryon perturbations through
recombination. This is not a task for a custom script — it requires
modifying an existing Boltzmann code.

**The tools:**
- **CLASS** (Cosmic Linear Anisotropy Solving System) — open source,
  written in C, well-documented. Can be modified to accept a custom
  dark energy model via a user-defined `w(z)` or `Omega_de(z)`.
- **CAMB** (Code for Anisotropies in the Microwave Background) —
  Python-compatible version available. Similar capability.

**What to do:**
1. Install CLASS: `git clone https://github.com/lesgourg/class_public`
2. Implement ESTIF's Ω_tilt(z) as a custom dark energy fluid
3. Run CLASS with ESTIF parameters and compare power spectrum to Planck

**This is likely where expert collaboration becomes necessary.**
A physicist familiar with CLASS/CAMB could implement this in days
where an independent researcher might spend weeks.

**Milestone:** CMB power spectrum from CLASS with ESTIF dark energy

---

## Phase 7: The Second Floor (Dark Matter) 📋

**Goal:** Address galaxy rotation curves and cluster lensing  
**Prerequisite:** None — can start independently of Phases 5 and 6  
**Difficulty:** Very high — completely new mechanism required

---

### 7.1 — Diagnose The Weak-Field Problem

**The current situation:**
The combined formula gives n ≈ 33 in weak fields (galactic scales).
x^(2n) ≈ 0 for any reasonable galactic curvature. So β ≈ 1 and
Observable ≈ 1 — essentially no correction at all.

ESTIF's strong-field formula is **inert** in galaxy rotation regimes.

**The question:**
Is there a second regime of the tilt geometry — one that activates
at very low curvature rather than high curvature — that could
produce the flat rotation curves?

In MOND (Modified Newtonian Dynamics), gravity is stronger below
a critical acceleration a₀ ≈ 1.2×10⁻¹⁰ m/s². ESTIF would need
something analogous — an additional tilt effect at very low x.

**What to do:**
Write `test_rotation_curve_diagnosis.py` that:
1. Calculates the curvature ratio x at various radii in a typical galaxy
2. Plots what the current ESTIF formula predicts for rotation velocity
3. Shows the gap between prediction and observed flat rotation curves
4. Quantifies what additional effect would be needed to close the gap

**Milestone:** Rotation curve gap quantified — know what needs explaining

---

### 7.2 — The Two Hypotheses

After diagnosing the gap, two distinct hypotheses need testing:

**Hypothesis D1 — A second tilt mode at low curvature:**
Perhaps the tilt geometry has two regimes — a suppression at high x
(what we have) and an enhancement at very low x. Mathematically:

```
n_galactic(x) = N_gal × exp(+B_gal × x)   [grows at small x]
```

This would be a mirror of the strong-field formula, activating at
the opposite extreme. Completely speculative — but testable against
rotation curves if parameters can be constrained.

**Hypothesis D2 — The 4D flow itself creates apparent mass:**
If the inward flow of 3D space through 4D creates apparent density
in regions where the flow is non-uniform, galaxies might experience
an effective dark matter contribution from the flow geometry itself.
This is the deeper conceptual claim — harder to formalise but
more physically motivated.

**What to do:**
Write `test_dark_matter_hypotheses.py` testing D1 against
the Milky Way rotation curve and one external galaxy (NGC 3198
is a standard test case with well-measured flat rotation).

**Milestone:** D1 ruled in or out against two rotation curves

---

### 7.3 — Cluster Lensing and Large Scale Structure

## Phase 7: The Second Floor (Dark Matter) 🔄

**Goal:** Test whether the eddy background can account for dark matter  
**Progress:** Analytical phase complete. Simulation required beyond this point.

---

### 7.1 — Analytical Dark Matter Results ✅ COMPLETE

Five results confirmed through analytical tests:

**Ωm = x₀ identity:**
```
x₀ = R_H / r_universe = 0.310734
Ωm (Planck) =           0.311100
Agreement:               0.12%  — within Planck 1σ

x₀ − Ωb = 0.261734
Ωdm (Planck) = 0.262000
Agreement:     0.10%
```

**Collisionless dynamics:**
- σ(r)/v_escape(r) = 0.5000 exactly — virial condition automatic at every scale
- λ_Jeans(r) = 2.57 × r — self-similar, every scale marginally unstable simultaneously
- Free-fall time at z=10: ~1 Gyr — correct epoch for galaxy formation
- Isothermal collapse (ρ ∝ 1/r²) → flat rotation curves by construction

**Multi-scale observable:**
```
Observable(r) = √β(x_local) × √β(x_galactic) × √β(x_cosmic)
```
At Earth: local and galactic terms = 1.000000 (GR compatible). Cosmic term = 0.830 (dark matter).

**Profile:**
- ESTIF gives v_flat ∝ M^(1/3) — Tully-Fisher exponent off by one step vs observed M^(1/4)

**One analytical test remaining:** Whether the tilt correction to ρ_halo adds
a mass-dependent factor that closes the Tully-Fisher exponent from 1/3 to 1/4.
Script: `test_tully_fisher_correction.py` — fast, on-Mac.

---

### 7.2 — N-Body Simulation (Requires Collaboration) 🔴 BUDGET WALL

**The wall:** v_flat = 220 km/s requires internal halo overdensity
δ ~ 50,000–100,000. This is a simulation output, not derivable analytically.

**What is needed:** A modified N-body code (Gadget-4, AREPO) with the ESTIF
force law ∇(ω²/2) implemented instead of pure Newtonian gravity. Run on
100s of CPU-hours minimum.

**Cannot be done on a Mac Mini.** Requires:
- Access to a university computing cluster, OR
- Collaboration with a computational cosmology group, OR
- Cloud compute allocation (AWS/GCP HPC — ~$500–2000 for a test run)

**What to do now:** Document the analytical predictions precisely so that
when simulation access becomes available, the test is well-defined.
The falsifiable prediction: ESTIF halos should have concentration parameter
c ~ 30–100 (vs NFW c ~ 10–20). If simulation gives c in that range → confirmed.

**Bullet Cluster test:** Cannot be addressed without simulation.
The spatial offset of gas and lensing mass in the Bullet Cluster is very
difficult to explain without particle dark matter. This is the hardest
dark matter test for any alternative model.

---

### 7.3 — Stress-Energy Derivation (Theoretical) 📋

**The open question:** Can ρ_eddy = x₀ × ρ_crit be derived from the 4D
stress-energy tensor projection? This would convert Ωm = x₀ from a
numerical coincidence to a derived result.

**What is needed:** Tensor calculus — the projection of a rotating
hypersurface's kinetic energy onto Tμν. Paper and pencil, no computer.

**This is the most important remaining theoretical problem in the project.**
If it succeeds, ESTIF derives Ωm from first principles. If it fails,
Ωm = x₀ remains a motivated coincidence.

---

## Phase 8: Integration and Complete Model 🔮

**Conditional on:** Phases 5–7 producing positive results

---

### 8.1 — Unified Formula

If Phases 5–7 succeed, the goal is a single geometric framework
that describes gravity, dark energy, and dark matter from one equation.
The current formula handles strong-field gravity and dark energy.
Dark matter requires an extension.

A potential unified form:
```
n(x) = N_strong × exp(−B × x) + N_weak × exp(+B_gal × x)
```

A two-term formula covering both extremes. Whether this is physically
motivated or numerological is the key question for theorists.

---

### 8.2 — CMB Extension

Implement ESTIF dark energy in CLASS and confirm the full CMB
power spectrum is reproduced to within Planck measurement precision.
This would elevate ESTIF from "passes low-redshift tests" to
"passes all known cosmological tests."

---

### 8.3 — Cosmological Simulations

To test large scale structure predictions, N-body simulations with
the ESTIF force law would be required. This is beyond the scope of
current resources — explicitly a collaboration target.

---

## Phase 9: Publication 📄

**When ready:** After at minimum Phase 5 is complete and Phase 6 has
produced at least the CMB angle consistency check.

---

### 9.1 — arXiv Preprint

**Target:** A focused paper on the combined formula and Option A cosmology.
Not waiting for dark matter. Not waiting for full CMB.

**Title candidate:**
"Geometric Dark Energy from 4D Hypersurface Tilt: A Unified Formula
for Strong-Field Gravity and Cosmological Constant"

**Core claims (defensible now):**
1. Combined formula satisfies EHT, Λ, LISA simultaneously
2. GR time dilation is the special case at n = ½
3. Dark energy replaced by Ω_tilt(z) — passes six low-z tests
4. Λ drift prediction: 0.023%/Gyr, approaching EUCLID threshold
5. N_MAX ≈ 5/7 × ln(r_e/l_P) — connection to electron scale

**Honest scope limitations in paper:**
- CMB not addressed (formula diverges at z ~ 1100)
- Dark matter not addressed
- Fractional multipliers 5/7 and 1/3 not derived

**Target journals:**
- Physical Review D
- Classical and Quantum Gravity
- JCAP (Journal of Cosmology and Astroparticle Physics)

---

### 9.2 — Expert Review Before Submission

Three types of expert needed:

1. **Cosmologist** — to evaluate Option A, the Λ drift prediction,
   and the BAO/SN/H₀ results
2. **GR / modified gravity theorist** — to evaluate the tilt formula,
   the n = ½ GR equivalence, and the physical interpretation
3. **Particle physicist** — to evaluate the electron radius connection
   and the N_MAX = 5/7 × ln(r_e/l_P) claim

---

## Open Questions That Drive The Work

These are the questions that will determine what the next theoretical
steps look like. Not all require code.

| Question | Phase | Type |
|---|---|---|
| Why 5/7 and 1/3? | 5 | Paper + pencil |
| Why x = 0.272 geometrically? | 5 | Paper + pencil |
| Can Ω_tilt be derived at high-z rather than cut off? | 5–6 | Theory |
| Does the ISW prediction match Planck? | 6 | Calculation |
| Can a second tilt mode explain rotation curves? | 7 | Numerical |
| Does ESTIF predict a Bullet Cluster analogue? | 7–8 | Conceptual |
| Does the 4D flow create apparent dark matter? | 7–8 | Theory |

---

## Risk Assessment

### If Ω_tilt cannot be regularised naturally (Phase 5.1)

A hard cutoff at z < 2 is an honest fallback. The paper's scope is
limited to low-redshift cosmology, which is still publishable and
scientifically valuable.

### If CMB angle check fails (Phase 6.1)

If ESTIF shifts the CMB acoustic scale by more than 0.5%, the current
Option A formula is ruled out by Planck. The formula would need
fundamental revision. This is a genuine risk — acknowledged upfront.

### If rotation curve tests fail (Phase 7)

Dark matter remains unexplained by ESTIF. The project publishes as
a partial replacement (dark energy only) and explicitly notes dark
matter as an open problem. This is still a significant contribution.

### If LISA detects no signal (~2034)

The strong-field formula is falsified. The cosmological work stands
independently. Science worked correctly.

---

## What Can Be Done Without A Supercomputer

**Doable on a Mac Mini — analytical and scripted:**

| Task | Type | Phase | Est. Time |
|---|---|---|---|
| test_tully_fisher_correction.py | Script | 7.1 | 1 hour |
| test_desi_comparison.py | Script | 5.3 | 1 hour |
| test_cmb_angle_estimate.py | Script | 6.1 | 2 hours |
| test_isw_prediction.py | Script | 6.2 | 4 hours |
| Why 5/7 and 1/3? | Theory | — | Unknown |
| Why x=0.272 geometrically? | Theory | — | Unknown |
| Stress-energy tensor derivation | Theory | 7.3 | Unknown |
| CLASS/CAMB integration (linear only) | Software | 6.3 | 1–2 weeks |

**Requires supercomputer or collaboration:**

| Task | Why Blocked | Path |
|---|---|---|
| v_flat = 220 km/s | δ~50,000 only from N-body | University cluster |
| Tully-Fisher 1/3→1/4 full | Concentration parameter c | N-body simulation |
| Bullet Cluster test | Spatial offset needs DM dynamics | N-body simulation |
| Full CMB power spectrum | Boltzmann equations | CLASS on laptop is OK, but precision fitting needs cluster |
| Large scale structure | 3D matter power spectrum | N-body simulation |

---

## Immediate Next Steps (Priority Order)

1. **Phase 7.1 remaining** — `test_tully_fisher_correction.py`
   Does obs(x_local) at r_virial add M-dependent factor closing 1/3→1/4?

2. **Phase 5.3** — `test_desi_comparison.py`
   ESTIF w = −1.08 vs DESI DR2. Fast script, direct comparison.

3. **Phase 6.1** — `test_cmb_angle_estimate.py`
   Does ESTIF shift the CMB acoustic scale by less than 0.5%?
   If yes → not ruled out by Planck. If no → major problem.

4. **Phase 6.2** — `test_isw_prediction.py`
   Late-ISW from Ω_tilt evolution — unique ESTIF signature in CMB.

5. **Phase 7.3** — Stress-energy tensor derivation (theory)
   If ρ_eddy = x₀ × ρ_crit can be derived, Ωm is no longer a free parameter.

6. **Phase 9.1** — arXiv preprint
   The strong-field + dark energy results are publishable now.
   Dark matter is documented as Phase 7 — acknowledged open problem.

---

## Lessons Learned

### From ESTIF-FD (v1.0) — 2024
The exponential S(t) cosmology was overambitious. Trying to derive
everything simultaneously failed. Lesson: validate components
independently before combining.

### From ESTIF-Gravity (v3.0) — 2025
Accepting ΛCDM as a foundation was correct. It produced clean
strong-field predictions. But the boundary was drawn too conservatively.
The tilt geometry naturally extends to cosmology.

### From ESTIF Option A (v4.0) — 2026
Dark energy replacement works partially. The formula needs a high-z
regularisation. Going straight to CMB without it would repeat the
ESTIF-FD mistake of building on an unstable foundation.

### From Dark Matter Investigation (v6.0) — 2026
The fluid dynamics analogy (Jeans instability) is valid but incomplete.
The correct framework is collisionless dynamics — σ(r) ∝ r, not a single
sound speed. The sand dune analogy works for the hierarchy but space is not
the Sahara: objects orbit, don't collide. The virial condition σ/v_esc = 0.5
is exact and automatic. The simulation wall is real — acknowledge it, document
it, and find collaborators rather than pretending the Mac Mini can do it.

**The pattern:** Validate each floor before building the next.
Know the boundary of what your tools can do.

---

## Version History

| Version | Date | Key Change |
|---|---|---|
| v1.0 | Sep 2024 | Original ESTIF-FD (S(t) cosmology) |
| v1.5 | Jan 2025 | Post-CMB-age-fix updates |
| v2.0 | Oct 2025 | ESTIF-Gravity fork (gravity only) |
| v3.0 | Oct 2025 | β derived from tilt geometry |
| **v4.0** | **Mar 2026** | **Combined formula + Option A cosmology** |
| **v5.0** | **Mar 2026** | **Joint SN+BAO fit, ALPHA_COSMO geometric derivation** |
| **v6.0** | **Mar 2026** | **Gravity=Time=Eddies, Ωm=x₀, collisionless dark matter** |
| **v6.1** | **Mar 2026** | **MOND derived, SPARC validated, DESI DR2 constraint** |
| **v6.2** | **Mar 2026** | **a₀ redshift constancy proved, parameter independence, letter drafted** |

---

## Version History

| Version | Date | Key Change |
|---|---|---|
| v1.0 | Sep 2024 | Original ESTIF-FD (S(t) cosmology) |
| v1.5 | Jan 2025 | Post-CMB-age-fix updates |
| v2.0 | Oct 2025 | ESTIF-Gravity fork (gravity only) |
| v3.0 | Oct 2025 | β derived from tilt geometry |
| **v4.0** | **Mar 2026** | **Combined formula + Option A cosmology** |

---

**End of Roadmap**

---

## ROADMAP UPDATE — v6.1 (18 March 2026)

---

### Phase 5.3 — DESI DR2 w(z) Comparison ✅ DONE — FAILS

**Status:** Complete. Result: FAILS.

ESTIF Ω_tilt(z) tested against DESI DR2 (arXiv:2503.14738, released 19 March 2026).
chi²/N = 10.8 vs DESI DR2 (ΛCDM: 1.9). Pre-existing prediction w_eff ≈ −1.08
falsified at 3.5σ (DESI DR2 w₀ = −0.73 ± 0.10).

Root cause: x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is circular.

Script: `tests/test_desi_wz_consistency.py`
Milestone: The cosmological sector failure is documented. Rework required.

---

### Phase 5.4 — Self-Consistent Ω_tilt(z) (NEW — PRIORITY 1)

**Status:** Not started. Priority 1 — blocking all cosmology claims.

Replace x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) with a self-consistent version:
x(z) = x₀ × (1+z) × H₀/H_ESTIF(z) via iterative solve.

H_ESTIF depends on Ω_tilt, which depends on x(z), which depends on H_ESTIF.
This is a fixed-point iteration — computationally tractable on Mac Mini.

What to do: Write `tests/test_omega_tilt_selfconsistent.py` that:
1. Starts with H_ΛCDM as initial guess for H_ESTIF
2. Iterates until convergence
3. Computes resulting chi²/N vs DESI DR2
4. Compares to current circular version

Milestone: Self-consistent Ω_tilt(z) tested against DESI DR2. If chi²/N < 2.0 → proceed. If still fails → rethink functional form entirely.

---

### Phase 7.0 — MOND Derivation ✅ COMPLETE (v6.1)

**Status:** Complete. Result: SOLID.

Four-step derivation of a₀ = H₀cx₀/√3 from geometry. Zero free parameters.
Confirmed against 87 quality-1 SPARC galaxies: RMS = 15.6%.

Scripts: `tests/derive_mond_from_geometry.py`, `tests/test_sparc_tully_fisher.py`

---

### Phase 7.0b — SPARC Bias Analysis ✅ COMPLETE (v6.1)

**Status:** Complete. Result: Calibration issue, not structural.

Script: `tests/test_sparc_bias_analysis.py`
Conclusion: All bias structure disappears after Υ* correction. Force law is sound.

---

### Phase 7.0c — Multiplier Derivation ✅ PARTIAL (v6.1)

**Status:** Partially complete.

- B = L/3: DERIVED from 3D isotropy (0.69% off)
- N_MAX/5/7: CONDITIONAL on x_c geometric derivation
- x_c = 0.272: OPEN — the remaining theoretical gap

Script: `tests/test_multiplier_derivation.py`

---

### Phase 7.0d — Derive x_c Geometrically (NEW — PRIORITY 2)

**Status:** Not started. Pure theory — no code needed.

x_c = 0.272 is where n(x_c) = 1/2 and ESTIF equals GR time dilation exactly.
It corresponds to r = 3.68 Rs from a black hole — between ISCO (3 Rs) and the photon sphere (1.5 Rs).

Approach: look at Schwarzschild thermodynamics at r = Rs/x_c:
- Hawking temperature: T ∝ 1/Rs — does it connect to x_c?
- Bekenstein-Hawking entropy: S ∝ Rs² — area element at r = Rs/x_c?
- Is x_c the solution to a geometric fixed-point equation?

Milestone: x_c derived from Schwarzschild geometry without calibration input.
If successful: 5/7 follows from 1/3 + x_c → full multiplier derivation complete.

---

### Phase 9 — Publication (NEW)

**Gravity Letter — ACHIEVABLE NOW**

Content: Geometric derivation of a₀ = H₀cx₀/√3 + SPARC validation.
Claim: First derivation of MOND critical acceleration from geometry, zero free parameters.
Confirmed: 87 SPARC quality-1 galaxies, RMS = 15.6%.
Scope: 4–6 pages. No cosmological sector claims.
Venue: MNRAS Letters, ApJL, JCAP.

Checklist:
- ✅ Derivation complete (zero free params)
- ✅ Uniqueness confirmed (1/12 factors qualifies, only with independent justification)
- ✅ SPARC validation (87 galaxies)
- ✅ 1/3 derived from isotropy (same principle as 1/√3 in MOND)
- ⚠️ x_c still observational — stated as open question in paper
- ⚠️ Υ* calibration needs honest caveat (0.85 vs standard 0.50)

---

### Immediate Next Steps (Priority Order)

| Priority | Task | Type | Estimate |
|---|---|---|---|
| 1 | Submit gravity-only letter | Submission | This week |
| 2 | Self-consistent Ω_tilt(z) | Script | 1–2 days |
| 3 | DESI DR2 retest after fix | Script | 0.5 days |
| 4 | Derive x_c geometrically | Theory | Unknown |
| 5 | SPARC with Υ* = 0.65 | Script update | 0.5 days |

