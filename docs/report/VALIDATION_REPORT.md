# ESTIF Validation Report

**Model Version:** ESTIF v6.3 — "The Split"
**Date:** 8 July 2026
**Status:** Gravity field equation DERIVED (not matched to Schwarzschild). Strong-field complete. MOND derived, SPARC validated. Cosmology reframed: frozen eddy = cosmological constant ties ΛCDM; Ω_tilt(z) retired. Project split into Path One (Core) and Path Two (Extended). Gravity letter ready and strengthened.

> **See also:** `MILESTONE_v6.3_THE_SPLIT.md` for the full v6.3 narrative.

---

## Executive Summary

The headline change since v6.2: **the gravitational field equation is now derived
from the flow axioms** rather than obtained by matching the flow profile to the
Schwarzschild solution. The strong-field tilt formula remains the tool for
deviations from GR (EHT, LISA) and for the MOND acceleration scale, and is
unchanged. The cosmology sector is reframed — its honest best result is a
cosmological constant (the frozen cosmic eddy), which ties ΛCDM on DESI DR2; the
evolving Ω_tilt(z) apparatus is retired.

**Calibrated tilt parameters (strong-field use, unchanged):**
```
N_MAX = 33.265    B = 15.429
n(x)  = N_MAX × exp(−B × x)
β(x)  = √(1 − x^(2n(x)))
Observable = √β(x)
```

---

## Part 0: The Derived Field Equation (v6.3 — new headline result)

Every prior version reproduced Newton's law by an unstated shortcut: the flow
profile was matched to the Schwarzschild solution (tilt exponent set to n = ½ so
that β = √(1 − Rs/r), the GR time-dilation factor), and the force was read off the
matched profile. The fidelity audit (`estif_fidelity_audit.py`) made this explicit.
v6.3 removes it.

### 0.1 The field equation is forced by the flow axioms

Starting from the three flow axioms — A1 (flat 3-slices carried through the bulk),
A2 (everything moves through the bulk at speed c), A3 (empty space is not a source)
— the flow metric is fed through a symbolic Gauss–Codazzi / ADM engine. The engine
**forces**:

```
ρ_eff = m′(r) / (4π r²)      i.e.   dm/dr = 4π r² ρ_eff
```

This is mass continuity — Poisson's equation in integrated form — the standard
(0,0) Einstein equation on flat slices. It is computed from the geometry, not
assumed.

| Check | Result | Status |
|---|---|---|
| ρ_eff = m′(r)/(4π r²) forced from flow metric | matches engine normal projection | ✅ |
| Vacuum (ρ_eff = 0) → v² = 2GM/r uniquely | Birkhoff in flow variables | ✅ |
| v² = 2GM/r → full Einstein tensor = 0 | **exact Schwarzschild**, no approximation | ✅ |
| Uniform-density ball → ρ_eff = ρ₀ | correct Newtonian source | ✅ |
| Volume-conserving ("plughole") flow | gives 1/r⁵ force — **fails** (control) | ✅ ruled out |

**Script:** `tests/estif_task4_field_equation.py` (5/5 checks pass).

### 0.2 Signature and special relativity emerge from the same construction

| Check | Result | Status |
|---|---|---|
| Lorentzian signature (−,+,+,+) from Euclidean bulk + A2 | minus sign produced, not inserted | ✅ |
| Time dilation dτ/dt = √(1 − v²/c²) | exact, zero free parameters | ✅ |
| Photons (dw = 0) null; rest clocks dτ = dt | exact | ✅ |
| Gullstrand–Painlevé flow metric = exact Schwarzschild vacuum | engine-verified | ✅ |

**Scripts:** `estif_flow_signature_dynamics.py` (18/18), `estif_converse_flow_law.py`
(Birkhoff in flow variables), `estif_tmunu_gauss_codazzi.py` (engine validated
against flat FRW and de Sitter).

### 0.3 What remains open in gravity

The strong-field **pressure/stress** sector (a fully relativistic interior with
pressure) requires the complete off-diagonal stress tensor. It is needed for none
of {vacuum, Newton, Schwarzschild}, and is the remaining rigorous gravity step.

**Verdict:** ✅ The former "Poisson postulate" is now a theorem for vacuum, Newton,
and Schwarzschild. The gravity sector rests on a derived foundation.

> **Precision (C6).** What is established is that the three flow axioms uniquely
> *select* the constraint (energy) sector of General Relativity in
> Painlevé–Gullstrand gauge — forcing mass continuity, hence exact Schwarzschild in
> vacuum and the Newtonian source in the weak field — **without matching to the
> Schwarzschild solution**. What is *adopted*, not derived, is the gravitational
> coupling: the identification of the geometric constraint scalar with 8πG × energy
> density. ESTIF does not derive G or the factor 8π from below. The correct claim is
> that the force law is no longer matched to GR's vacuum solution; the coupling to
> matter is still the standard Einstein–Hilbert one.

---

## Part 1: Strong-Field Gravity (deviations from GR — unchanged)

> **Scope note (v6.3):** the tilt formula below is the tool for *strong-field
> deviations from GR* and for the a₀ scale. Its use as a *cosmological dark-energy
> law* is retired (see Part 2). The two uses are independent.

### 1.1 Joint Calibration — Three Simultaneous Tests

#### EHT M87\* Shadow

| Quantity | Value |
|---|---|
| x at photon sphere | 0.6667 |
| n(x) | 0.0011 |
| Observable √β | 0.1741 |
| Shadow predicted | 42.00 μas |
| Shadow observed | 42.0 ± 3.0 μas |
| Tension | **0.00σ** ✅ |

#### Cosmological Constant (as a local-tilt scale)

| Quantity | Value |
|---|---|
| x at cosmic scale | 0.3107 |
| Observable √β | 0.8300 |
| Λ predicted | 1.1056 × 10⁻⁵² m⁻² |
| Λ measured | 1.1056 × 10⁻⁵² m⁻² |
| Ratio | **1.000000** ✅ |

#### LISA GW Delay (65 M☉ merger)

| Quantity | Value |
|---|---|
| x at ISCO | 0.3333 |
| Observable √β | 0.7677 |
| GW delay predicted | 491.7 μs |
| LISA S/N | **49.2σ** ⚠️ conditional — see C1 note below |

**Verdict:** All three are *consistent* with observation simultaneously, with zero
free parameters after calibration.

> ⚠️ **Conditional (C1).** The ESTIF vacuum is exactly Schwarzschild (Part 0), so
> any *deviation* from GR in photon-sphere shadows or in gravitational-wave
> propagation through vacuum must be sourced by the non-vacuum eddy background — a
> sector not yet derived. The EHT and LISA figures above are therefore predictions
> **conditional on** the open eddy-stress derivation, not established results. The
> observations remain consistent with ESTIF; the *deviation from GR* is what awaits
> derivation. The Planck Λ entry is a calibration match, not a vacuum deviation, and
> is unaffected by this caveat.

### 1.2 Gravity = Generalized Time Dilation

`β(x) = √(1 − x^(2n(x)))` reduces to the Schwarzschild factor `τ(x) = √(1 − x)`
when n = ½, occurring naturally at x = 0.2721:
```
β(0.2721) = 0.8531    τ(0.2721) = 0.8532    (4-decimal match) ✅
```
GR time dilation is the special case of ESTIF tilt at n = ½. In v6.3 this is
consistent with the derived vacuum result (exact Schwarzschild); the tilt formula
supplies the deviations from it. See `test_gravity_time_connection.py`.

### 1.3 Natural Scale — Electron Radius

```
ln(r_e/l_P) = 46.608
5/7 × ln(r_e/l_P) = 33.291   vs N_MAX = 33.265   (0.079%) ✅
1/3 × ln(r_e/l_P) = 15.536   vs B     = 15.429   (0.693%) ✅
```
B = L/3 is derived from 3D isotropy; N_MAX = 5/7 × L is conditional on the
still-open geometric derivation of x_c = 0.272. See `test_electron_connection.py`.

### 1.4 GW Delay Mass Dependence (n = 0.05) — conditional (C1)

> ⚠️ Every entry in this table is a *deviation from GR* in vacuum GW propagation.
> Under the exact-Schwarzschild vacuum result (Part 0) these deviations require the
> non-vacuum eddy background as their source. They are conditional on the open
> eddy-stress sector, not established predictions.

| Binary Mass | GW Delay | LISA S/N | Status |
|---|---|---|---|
| 10 M☉ | 32 μs | 3.2σ | ⚠️ conditional |
| 30 M☉ | 95 μs | 9.5σ | ⚠️ conditional |
| 65 M☉ | 207 μs | 20.7σ | ⚠️ conditional |
| 100 M☉ | 318 μs | 31.8σ | ⚠️ conditional |
| 500 M☉ | 1.6 ms | 158σ | ⚠️ conditional |

> **Retired in v6.3:** the Λ-drift prediction (0.023%/Gyr) depended on reading N_MAX
> as a cosmological ratio ln(r_universe/Rs_m87). Under the electron-scale anchoring
> and the cosmology reframe it is no longer load-bearing; it is documented as an
> explored direction.

---

## Part 2: Cosmology — The Reframe (v6.3)

This sector changed the most. The evolving Ω_tilt(z) dark-energy law is retired; the
honest best result is a cosmological constant that comes for free from the derived
gravity.

### 2.1 The frozen cosmic eddy is a cosmological constant

A corollary of the derived field equation (Part 0): a flow whose effective energy
density is *constant* is exact de Sitter — a cosmological constant, w = −1. Tested
against real DESI DR2 BAO data:

```
Frozen eddy (w = −1):   χ²/N = 1.92   —   tied with ΛCDM.
```

Derived, no tilt formula, no fitted parameters beyond matter and the constant.

### 2.2 Why the evolving Ω_tilt(z) apparatus is retired

The old model:
```
H²(z) = H₀² [Ωm(1+z)³ + Ω_tilt(z)],    Ω_tilt(z) = Ω_Λ (obs_now/obs_z)²
```
failed DESI DR2 at χ²/N = 10.8, for two reasons now understood:

| Stage | Model | χ²/N vs DESI DR2 | Note |
|---|---|---|---|
| Original | Ω_tilt(z) circular ruler | 10.80 | x(z) used H_ΛCDM as its own ruler |
| Fix 1 | Ω_tilt(z) self-consistent ruler | 3.35 | circularity removed (Task 5) |
| Reframe | **frozen eddy = Λ (derived)** | **1.92** | ties ΛCDM, beats the tilt |
| Derived eddy E1 (conserved spin) | w = +1 stiff | 3232 | falsified (Task 6) |
| Derived eddy E2 (tracker) | thaws to ~0 | 754 | falsified (Task 6) |
| ΛCDM reference | w = −1 | 1.92 | — |
| Best-fit evolving (CPL, fitted) | w0=−0.85, wa=−0.45 | 0.66 | Path Two target ⚠️ under-marginalized, see C4 |

**Conclusion:** on DESI, the entire Ω_tilt apparatus (dynamic n, N_MAX, B, sign
choice, z<2 cutoff) is a **net negative** relative to the plain cosmological
constant underneath it. The cosmological dark-energy claim is retired to an
"explored and set aside" appendix; the tilt formula's strong-field use (Part 1)
stands. Scripts: `estif_task5_desi_selfconsistent.py`, `estif_task5b_cosmo_eos.py`,
`estif_task6_eddy_eos.py`.

### 2.3 What DESI prefers, and the effective w(z)

DESI DR2 BAO alone prefers a mild *thawing* dark energy (best-fit CPL χ²/N ≈ 0.66,
w rising toward −0.85 today). ESTIF's self-consistent tilt w(z) already leans this
way, tracking DESI's published w(z) to within ~0.05, and its 3.35 is essentially
where DESI's *own* published w0,wa land on the BAO-only subset (3.09).

| model | w_eff(z≈0) |
|---|---|
| ESTIF old (circular Ω_tilt) | −1.358 (superseded) |
| ESTIF self-consistent tilt | −0.80 |
| Frozen eddy (Path One) | −1.00 |
| DESI DR2 published | −0.73 ± 0.10 |

So the honest picture is not "ESTIF fails DESI" — it is "a cosmological constant
ties ΛCDM, and a small *derived* thawing correction could do better if it can be
derived." That is Path Two.

> ⚠️ **Caveat (C4).** The χ²/N ≈ 0.66 CPL target comes from a BAO-alone fit with rd,
> H₀, and Ωm held fixed at Planck values. Fixing nuisance parameters inflates the
> apparent evolving-dark-energy advantage; DESI's own BAO-alone preference for
> evolving DE is considerably milder, and the strong combined-data significance
> comes from adding CMB and supernovae. This target must be re-derived with rd, H₀,
> and Ωm marginalized (or at least profiled) before it is used to justify Path Two.
> The AIC comparison (frozen eddy vs fitted CPL) should be repeated under the same
> marginalization.

> **Retired to appendix:** the Hubble-radius Λ bridge, the Ω_tilt(z) evolution law,
> the six pre-2026 low-z fits (Pantheon+ 2.08–2.33σ, BAO 5/5, age, H₀ tension), and
> the w_eff = −1.358 prediction. Documented as explored; no longer load-bearing.

---

## Part 3: Dark Matter — Analytical Phase (unchanged; now connected to Part 0)

### 3.0 The Ωm Bootstrap (v6.3.2 — conditional on principle P)

Correction C2 established that `Ωm = x₀` is circular: r_universe is the ΛCDM particle
horizon, an integral containing Ωm. The bootstrap asks whether that circularity can
be *solved* rather than removed.

**Principle P (adopted, not derived):** `Ωm = R_H / r_p`, with r_p the particle
horizon. Substituting the horizon integral closes the loop into Ωm · I(Ωm) = 1
which has a **unique** root.

| Case | Measured inputs | Root Ωm | vs Planck (0.3111 ± 0.0056) |
|---|---|---|---|
| Matter + Λ only | none | 0.3043 | −2.2% |
| **With radiation** | **T_CMB, N_eff, h** | **0.31408** | **+0.96%, 0.53σ** |

**Back-prediction:** r_universe = 4.353×10²⁶ m, −1.07% against the 4.4×10²⁶ m import
that C2 flagged as Ωm-dependent. The import is no longer needed.

**Identity:** P is equivalent to the statement that the mean matter pull at the
horizon equals `cH₀/2`. This is the same cH₀ that sets a₀; the ratio computes to
1.00000. Whether that is the mechanism or a coincidence is undetermined.

**Closure — propagating the bootstrap Ωm through the framework:**

| Quantity | With Planck Ωm | With bootstrap Ωm | Note |
|---|---|---|---|
| a₀ | 1.179×10⁻¹⁰ (1.72%) | **1.1920×10⁻¹⁰ (0.66%)** | vs MOND empirical |
| SPARC v_flat | — | ×1.00269 | insensitive |
| DESI DR2 χ²/N | 1.919 (ΛCDM) | **1.618** | fixed-(H₀, rd) |

**Input ledger after adopting P:**
measured  = { H₀, T_CMB, N_eff }
computed  = { Ωm, Ω_Λ, x₀, r_universe, a₀ }
#### Honest flags (all four, none hidden)

1. **Everything above is conditional on P**, which is *not* derived from A1–A3. Part
   B — deriving P — is open. Three candidate routes exist (flow-budget amplitude;
   horizon-acceleration balance; homogeneous field equation); none has been attempted
   in earnest. Without Part B this is a reparametrization of the C2 circularity, not
   an escape from it.
2. **Gaztañaga adjacency.** The causal-universe scale (≈ 0.3176 H₀, reached via
   inflation) sits next to this root. A comparison memo is a **prerequisite for any
   novelty claim** and does not yet exist. No priority is asserted here.
3. **DESI 1.618 is a fixed-ruler result.** Under C4-style marginalization over rd,
   H₀, and Ωm the ordering against ΛCDM could change.
4. **a₀'s empirical target carries ~10% scatter.** Improving 1.72% → 0.66% inside a
   10% band is pleasing, not decisive. It is not a detection.

**Scripts:** `estif_omega_bootstrap.py`, `estif_bootstrap_closure.py`.

**Status:** 🔶 Part A resolved (conditional); Part B open. Until Part B lands and the
Gaztañaga memo clears, §3.1's C2 downgrade stands unchanged — see RHAC Scenario Q.

---

### 3.1 The Ωm = x₀ Consistency Relation

```
x₀ = R_H / r_universe = 0.310734     Ωm (Planck) = 0.311100     (0.12%)
x₀ − Ωb = 0.261734                   Ωdm (Planck) = 0.262000     (0.10%)
```
> **Epistemic status (C2).** This is a **consistency relation**, not an
> Ωm-independent prediction. Here r_universe = 4.4×10²⁶ m is the ΛCDM particle
> horizon — an integral that itself contains Ωm. The 0.12% agreement is therefore a
> self-consistency of the geometric picture with Planck-calibrated values, not a
> derivation of Ωm. An independent prediction requires deriving r_universe from the
> flow framework without the Ωm-dependent horizon integral — an open task (RHAC
> Scenario H). The numerical agreement is real and stands; only its status as a
> "prediction" is withdrawn.

**v6.3 connection:** the derived field equation gives the *local* ρ_eff = m′/4πr².
Whether ρ_eddy = x₀ρ_crit emerges is the *homogeneous* version of the same
calculation — now well-posed. See `test_eddy_dark_matter.py`.

### 3.2 Gravity = Time = Eddies (three-way identity)

At x = 0.272 (n = ½), τ(x) = √β(x) and (ω/H₀)² = x, so `a = −c²∇(ω/H₀)²/2 = GM/r²`.
In v6.3 this agrees with, and is underwritten by, the derived field equation. Multi-
scale observable at Earth: cosmic term 0.830 dominant; local/galactic terms
1.000000 (GR recovered). See `test_eddy_time_gravity.py`, `test_solar_system_eddy.py`.

### 3.3 Collisionless Dynamics

```
σ(r) = r √(2πG ρ_eddy/3)          (grows linearly with r)
σ(r)/v_escape(r) = 0.5000          (exact, every scale) ✅
λ_Jeans(r) = √(2π²/3) r = 2.565 r  (self-similar) ✅
t_ff(z=10) = 1.1 Gyr               (correct epoch) ✅
```
See `test_collisionless_eddy.py`.

### 3.4 The N-Body Wall (unchanged)

`v_flat = 220 km/s` requires internal halo overdensity δ ~ 50,000–100,000 × ρ_eddy —
a simulation output from virialization. Falsifiable prediction: ESTIF halos reach
that δ under the flow force law. Requires a cluster or cloud HPC; a collaboration
target, not a framework limitation.

### 3.5 Tully-Fisher Exponent

ESTIF pure geometry gives v_flat ∝ M^(1/3); the MOND limit `v_flat⁴ = G M_bar a₀`
gives the observed M^(1/4). Resolved via the MOND limit (v6.1). See
`test_tully_fisher_correction.py`.

---

## Part 4: MOND Derivation (v6.1, now on the derived foundation)

```
Step 1 — Force law:  a = −c²/2 ∇(ω/H₀)² = GM/r²   [now underwritten by Part 0]
Step 2 — Cosmic flow (sideways component):  v_flow = c x₀ = c Ωm
Step 3 — 3D isotropic projection:  v_3D = v_flow/√3
Step 4 — Threshold:  a₀ = v_3D H₀ = H₀ c x₀/√3 = 1.179×10⁻¹⁰ m/s²
```
MOND empirical 1.200×10⁻¹⁰; agreement **1.72%**; zero free parameters.

**v6.3 note:** Step 1's force law is the weak-field face of the derived vacuum
result v² = 2GM/r (Part 0), not a Schwarzschild match — closing the "borrowed
recipe" gap that previously sat under Step 1.

**SPARC validation:** 87 quality-1 galaxies, RMS = 15.6% (within BTFR scatter); 82%
within 20%, 97% within 30%; −7.6% bias traces to Υ* calibration. **a₀ redshift
constancy:** H(z) cancels in the comoving frame (deviation 2.22×10⁻¹⁶).
**Parameter independence:** 3,600 H₀/Ωm combinations within ±20% SPARC scatter; 8
published datasets pass. Scripts: `derive_mond_from_geometry.py`,
`test_sparc_tully_fisher.py`, `test_a0_redshift.py`, `test_a0_parameter_independence.py`.

---

## Part 5: The Two Paths (v6.3)

| | Path One — ESTIF-Core (clean) ✅ default | Path Two — ESTIF-Extended (hard) 🔬 |
|---|---|---|
| Gravity | derived field equation (Part 0) | same |
| Cosmology | frozen eddy → Λ, ties ΛCDM (χ²/N = 1.92) | derive small thawing from vorticity T_μν (target 0.66) |
| Ω_tilt(z) | retired to appendix | — |
| Risk | low, publishable now | high, timeline unknown |

---

## Part 6: Historical (unchanged)

- **ESTIF-FD v1.0:** S(t) = exp(−∫H dt), χ² 3.8× worse than ΛCDM → ruled out. Concept kept, equation abandoned.
- **ESTIF-Gravity v3.0:** fixed n could not satisfy EHT and Λ together → resolved by dynamic n.

---

## Part 7: Known Limitations (v6.3)

- **Strong-field pressure/stress sector:** full off-diagonal T_μν remaining (Part 0.3).
- **Path Two cosmology:** leading w(z) correction from the vorticity stress tensor; naive reductions (E1, E2) falsified.
- **x_c = 0.272:** not yet geometrically derived (closes N_MAX = 5/7 × L).
- **ρ_eddy = x₀ρ_crit:** homogeneous version of the derived field equation; open.
- **Dark-matter halos:** N-body simulation (collaboration target).
- **CMB:** on Path One (pure Λ) this is the standard ΛCDM check; on Path Two it follows P2.1.
- Not peer-reviewed.

---

## Summary (v6.3)

| Component | Status | Key Result |
|---|---|---|
| **Gravity field equation** | ✅ **Derived** | mass continuity = Poisson; vacuum → exact Schwarzschild |
| Signature + SR from flow axioms | ✅ Derived | (−,+,+,+) and dτ/dt = √(1−v²/c²) |
| Strong-field formula (EHT+Λ+LISA) | ⚠️ Consistent; *deviations* conditional (C1) | 3 tests simultaneous, 0 free params |
| Gravity = time = eddies | ✅ Confirmed | β = τ at n = ½; Newton from gradient |
| Natural scale (electron radius) | ✅ Identified | N_MAX ≈ 5/7 × L, B = L/3 (derived) |
| MOND a₀ (on derived foundation) | ✅ Derived | 1.72%, zero params, SPARC RMS 15.6% |
| a₀ redshift constancy | ✅ Proved | H(z) cancels (2.22×10⁻¹⁶) |
| Parameter independence | ✅ Confirmed | 3,600 combos within SPARC scatter |
| Ωm = x₀ consistency relation | 🔶 Downgraded (C2) | 0.12%, but r_universe is the Ωm-dependent ΛCDM horizon |
| **Ωm bootstrap** (conditional on P) | 🔶 **NEW v6.3.2** | Ωm·I(Ωm)=1, unique root 0.31408 (0.53σ); a₀ → 0.66%; DESI 1.618 fixed-ruler |
| Collisionless dark matter | ✅ Confirmed | σ/v_esc = 0.5, λ = 2.565r |
| **Cosmology (Path One)** | ✅ **Frozen eddy = Λ** | ties ΛCDM (χ²/N = 1.92) |
| Cosmology Ω_tilt(z) | 🔴 Retired | net negative on DESI vs the Λ limit |
| Cosmology (Path Two) | 🔬 Open | derive thawing from vorticity T_μν |
| Strong-field pressure sector | 🔄 Open | full off-diagonal T_μν |
| v_flat from simulation | 🔴 Budget wall | N-body, cluster/HPC |

---

## Part 8: v6.3 Test Scripts

| Script | Purpose | Status |
|---|---|---|
| `estif_task4_field_equation.py` | Field equation derived (mass continuity = Poisson) | ✅ 5/5 |
| `estif_flow_signature_dynamics.py` | Signature + SR + Newton from flow axioms | ✅ 18/18 |
| `estif_converse_flow_law.py` | Vacuum forces v²=2A/r (Birkhoff in flow variables) | ✅ |
| `estif_tmunu_gauss_codazzi.py` | ADM engine, validated vs FRW + de Sitter | ✅ |
| `estif_task5_desi_selfconsistent.py` | De-circularized DESI DR2 (10.8 → 3.35) | ✅ |
| `estif_task5b_cosmo_eos.py` | DESI-preferred w(z); tilt tracks it within ~0.05 | ✅ |
| `estif_task6_eddy_eos.py` | Frozen-eddy reframe (1.92 ties ΛCDM, beats tilt) | ✅ |
| `estif_fidelity_audit.py` | Axiom-presence audit of the corpus | ✅ |

Prior v6.1/v6.2 test records (MOND, SPARC, multipliers, a₀ constancy, parameter
independence) stand unchanged. The v6.2 DESI entry (`test_desi_wz_consistency.py`,
χ²/N = 10.8) is superseded by Part 2 above.

---

**Validation Report Version:** 6.3.2 / 9 July 2026.
