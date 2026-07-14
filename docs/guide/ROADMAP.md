# ESTIF Development Roadmap

**Version:** 6.4.0
**Last Updated:** 12 July 2026
**Status:** v6.4.0 — A1→A1′; Phase 2 null; Fronts 1–3; C-15 (GW sector) DERIVED; Principle P proven non-derivable as a law (RHAC-009). Gravity letter ready.

> ⚠️ **v6.4.0 (12 July 2026):** this roadmap body is the v6.2/v6.3 plan preserved as history and is largely superseded. Authoritative current state: `docs/plan/PATH_ONE_CHECKLIST.md` (item tracker) + `docs/plan/RHAC.md` (RHAC-001…009). Key correction since: "derive Principle P" — listed below and elsewhere as a top target — is **CLOSED** (P is a predictive postulate, not derivable from the axioms; RHAC-009). Redirect foundational effort to x_c = 0.272 or the √3/cH₀ scale.

> ⚠️ **Read the v6.3 update at the bottom of this file first.** Everything between
> here and the *ROADMAP UPDATE — v6.3* heading is the **v6.2 plan, preserved as
> history**. Large parts of it are superseded and should not be worked from:
>
> - **Phases 5.3, 5.4, 6.1, 6.2, 6.3, 8.2** build on the evolving Ω_tilt(z)
>   dark-energy law, which is **retired** (it fits DESI DR2 worse than the plain
>   cosmological constant underneath it — χ²/N 3.35 vs 1.92).
> - The **"Dark Energy — Partial"** section's six passing low-z tests and the
>   **w = −1.08** prediction are retired to an appendix; see
>   `docs/report/VALIDATION_REPORT.md` Part 2.
> - The **Λ drift = 0.023%/Gyr** prediction is no longer load-bearing.
> - The **EHT / Λ / LISA** ✅ marks are now ⚠️ conditional: the derived ESTIF vacuum
>   is exactly Schwarzschild, so any *deviation* from GR must be sourced by the
>   un-derived eddy-stress sector (`CORRECTIONS_v6.3.1.md`, C1).
> - **Ωm = x₀** is a consistency relation, not a prediction (C2).
>
> The v6.2 mission statement, lessons learned, and dark-matter analytics stand.

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
| Age of universe | 13.379 Gyr ✅ (oldest stars ≥ 13.5 Gyr) |
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

**Publication status:** The gravity-only letter is ready for submission now (v6.2). A full paper covering dark energy will require Phase 5 to be complete and Phase 6 to have
produced at least the CMB angle consistency check. See the v6.1 update below for the current publication checklist.

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

---

## ROADMAP UPDATE — v6.3 "THE SPLIT" (8 July 2026)

The mission is unchanged ("replace ΛCDM, don't adjust it"), but the route forks.
Two derivation results since v6.2 restructure the roadmap. See
`MILESTONE_v6.3_THE_SPLIT.md`.

---

### Foundation upgrade — gravity is now DERIVED

The strong-field gravity floor was marked 100% complete in v6.2, but its
foundation was a Schwarzschild match (set n = ½). v6.3 replaces that: the field
equation `rho_eff = m′(r)/(4πr²)` (Poisson in integrated form) is now **forced**
by the flow axioms via Gauss–Codazzi, giving exact Schwarzschild in vacuum and
Newtonian ρ₀ for a source, with no Poisson postulate (Task 4). The foundation is
now derived, not borrowed.

```
Foundation (strong-field gravity):    ████████████  100%  ✅ Complete + DERIVED
Field equation from flow axioms:      ████████████  100%  ✅ NEW (Task 4)
Strong-field pressure/stress sector:  ███░░░░░░░░░   25%  🔄 full T_μν remaining
Cosmology (Path One = frozen eddy):   ████████████  100%  ✅ ties ΛCDM, derived
Cosmology (Path Two = thawing eddy):  █░░░░░░░░░░░   10%  🔬 hard, vorticity T_μν
Dark matter:                          ██░░░░░░░░░░   20%  🔄 analytical done, sim wall
```

---

### The fork replaces the old Phase 5–6 cosmology plan

The prior roadmap's Phase 5.4 (self-consistent Ω_tilt) and Phase 6 (CMB on top
of Ω_tilt) are **superseded**. Phase 5.4 was done (Task 5: 10.8 → 3.35) and
revealed the tilt shape is the problem; Task 6 then showed the frozen-eddy limit
(w = −1) ties ΛCDM and beats the tilt. So the cosmology roadmap becomes:

#### PATH ONE — ESTIF-Core (clean) ✅ recommended default

- **P1.1 Write the axioms into the theory.** ✅ **DONE (v6.3 concept rewrite.)**
  A1–A3 are now stated explicitly in `ESTIF_CONCEPT.md`.
- **P1.2 Resolve the A1 conflict.** ✅ **DONE.** The shrinking-ruler narrative is
  retired in favour of the flow (Painlevé–Gullstrand) picture, and v_flow = cx₀ is
  relabelled as the *sideways component* of a total-c flow.
- **P1.2b Apply the v6.3.1 errata (C1–C6)** across the repo. ⬜ This is now the
  gating item before submission — checklist P-2.
- **P1.3 Rewrite the cosmology sector** as "constant cosmic eddy → cosmological
  constant, χ²/N = 1.92 (ties ΛCDM)." Move Ω_tilt(z), N_MAX, B, the sign-flip and
  the z<2 cutoff to an "explored and set aside" appendix.
- **P1.4 Submit the gravity letter** (strengthened by the derived field equation).

#### PATH TWO — ESTIF-Extended (hard) 🔬 high-risk research

- **P2.1 Vorticity stress tensor.** Derive the leading correction to w = −1 from
  the full rotating-shear / vorticity T_μν (the cosmological half of the T_μν
  work). Target the mild DESI thawing (nominal χ²/N ≈ 0.66 — ⚠️ under-marginalized, see C4; re-derive with rd, H₀, Ωm marginalized before committing). The naive reductions E1
  (stiff, χ²=3232) and E2 (tracker, χ²=754) are already falsified; the full
  off-diagonal tensor is required.
- **P2.2 CMB / ISW** only after P2.1 produces a well-behaved, DESI-consistent
  H(z). On Path One (pure Λ) the CMB check is the standard ΛCDM one.

---

### Immediate Next Steps (v6.3 priority order)

| Priority | Task | Path | Type | Estimate |
|---|---|---|---|---|
| 1 | Apply the v6.3.1 errata (C1–C6) across the repo | One | Writing | days |
| 2 | Preferred-frame / simultaneity section (barrier 5, checklist C-11) | One | Writing | weeks |
| 3 | Submit gravity letter (derived field equation + C1/C6 wording) | One | Submission | after 1–2 |
| 4 | Strong-field pressure/stress sector (full T_μν) | Both | Theory | unknown |
| 5 | Vorticity stress-tensor derivation of leading w(z) | Two | Theory | unknown (hard) |

> This table is a summary. `docs/plan/PATH_ONE_CHECKLIST.md` is the authoritative
> item-level tracker; where the two disagree, the checklist wins.

---

**Roadmap Version:** 6.4.0
**Last Updated:** 12 July 2026
