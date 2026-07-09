# ESTIF Development Status

**Last Updated:** 8 July 2026
**Version:** 6.3 — "The Split"
**Status:** Gravity letter ready (and strengthened). Project split into Path One (ESTIF-Core, clean) and Path Two (ESTIF-Extended, hard). See `MILESTONE_v6.3_THE_SPLIT.md`.

---

## Executive Summary

ESTIF is a geometric model deriving gravity, dark energy, and dark matter from the claim that 3D space is a flat hypersurface carried through a 4D bulk. As of v6.3 the project has been restructured into two tracks after a sequence of derivations (Tasks 4–6, July 2026):

- **Path One — ESTIF-Core (clean, recommended default):** gravity now rests on a *derived* field equation; cosmology is a plain cosmological constant that ties ΛCDM. The tilt apparatus is retired.
- **Path Two — ESTIF-Extended (hard, high-risk):** attempt to derive a small thawing correction to w = −1 from the full vorticity stress tensor.

The single most important change since v6.2: **the gravity field equation is now derived from the flow axioms rather than matched to the Schwarzschild solution.**

> **Precision (C6).** The axioms uniquely *select* the constraint (energy) sector of General Relativity in Painlevé–Gullstrand gauge, forcing mass continuity and hence exact Schwarzschild in vacuum — without matching to the Schwarzschild solution. The gravitational *coupling* (geometric constraint scalar ↔ 8πG × energy density) is adopted, not derived from below. The force law is no longer matched to GR's vacuum solution; the coupling to matter is still Einstein–Hilbert.

---

## Three-Sector Status

### Gravity Sector ✅ SOLID — and now on a DERIVED foundation

| Result | Value | Status |
|---|---|---|
| **Field equation derived (mass continuity = Poisson)** | rho_eff = m′(r)/(4πr²) forced by flow metric | ✅ NEW v6.3 (Task 4) |
| **Vacuum → exact Schwarzschild** | v² = 2GM/r unique, full Einstein tensor = 0 | ✅ NEW v6.3 |
| **Weak-field source → Newtonian ρ₀** | uniform ball returns rho0 exactly | ✅ NEW v6.3 |
| Signature + SR derived from Euclidean bulk + universal-c | (−,+,+,+) and dτ/dt = √(1−v²/c²) | ✅ NEW v6.3 (signature suite) |
| MOND a₀ derived (now on derived field eq) | H₀cx₀/√3 = 1.179×10⁻¹⁰ m/s² (Planck Ωm); **1.1920×10⁻¹⁰ with bootstrap Ωm** | ✅ 1.72% → **0.66%** (conditional on P) |
| SPARC BTFR (87 galaxies, Qual-1) | RMS = 15.6% | ✅ Within observed scatter |
| EHT M87* shadow | 42.0 μas, 0.00σ | ⚠️ consistent; *deviation* conditional (C1) |
| Planck Λ (local tilt, unchanged) | ratio = 1.0000 | ✅ calibration match, not a vacuum deviation |
| LISA GW delay | 491 μs, S/N = 49.2σ | ⚠️ conditional (C1) |
| β = τ at n = ½ (x = 0.272) | GR as special case | ✅ |
| a₀ redshift constancy | H(z) cancels exactly | ✅ v6.2 |
| Parameter independence | 3,600 H₀/Ωm combos within SPARC scatter | ✅ v6.2 |

> ⚠️ **Conditional (C1).** The ESTIF vacuum is exactly Schwarzschild (Task 4), so any *deviation* from GR in photon-sphere shadows or in vacuum GW propagation must be sourced by the non-vacuum eddy background — the sector below. The observations remain *consistent* with ESTIF; the *deviation from GR* is what awaits derivation. The single-speed flow ansatz forces p_r = −ρ, so it can represent only vacuum and Λ-like sources — which is *why* no vacuum deviation is available (`estif_tmunu_task4.py`, limitation L1).

**Open in gravity:** the strong-field pressure/stress sector (relativistic interior with pressure) still needs the full off-diagonal T_μν. Not needed for vacuum, Newton, or Schwarzschild.

### Cosmology Sector 🔄 REFRAMED — frozen eddy ties ΛCDM; tilt retired

| Test | Result | Status |
|---|---|---|
| **Frozen eddy = cosmological constant (DERIVED)** | χ²/N = 1.92 (ties ΛCDM) | ✅ NEW v6.3 — the honest best |
| DESI DR2 — circular Ω_tilt(z) | χ²/N = 10.80 | ❌ superseded |
| DESI DR2 — self-consistent Ω_tilt(z) (Task 5) | χ²/N = 3.35 | 🔄 fixed circularity, still short of ΛCDM |
| Best evolving-w flow (fitted, Task 5b) | χ²/N = 0.66 at w0=−0.85, wa=−0.45 | ⚠️ under-marginalized (C4); rd, H₀, Ωm fixed at Planck — re-derive before use as the Path Two bar |
| Derived eddy E1 (conserved-L spin) | χ²/N = 3232 | ❌ falsified (Task 6) |
| Derived eddy E2 (tracker) | χ²/N = 754 | ❌ falsified (Task 6) |
| Age of universe | 13.379 Gyr | ✅ Stands |
| **Ωm bootstrap (conditional on P)** | Ωm·I(Ωm)=1, unique root **0.31408** (0.96% from Planck, 0.53σ); DESI DR2 χ²/N = **1.618** (fixed-ruler) | 🔶 **NEW v6.3.2** |

**Root cause of the tilt failure (now understood):** the tilt *shape* fits worse than the frozen-eddy limit underneath it. The Ω_tilt apparatus (N_MAX, B, sign-flip, z<2 cutoff) is a net negative on DESI and is **retired** from the cosmology claim under Path One.

### Dark Matter Sector 🟡 Analytical phase complete, simulation wall

Ωm = x₀ (0.12%), Ωdm = x₀ − Ωb (0.10%), σ/v_esc = 0.5 exact, λ_Jeans = 2.565r. v_flat = 220 km/s requires δ ~ 50,000–100,000 — N-body, off-Mac. The stress-energy derivation of ρ_eddy = x₀ρ_crit is now connected to the Task 4 machinery (the local field equation is derived; the homogeneous-source version is the same calculation).

> **Epistemic status of Ωm = x₀ (C2).** This is a **consistency relation**, not an Ωm-independent prediction. r_universe = 4.4×10²⁶ m is the ΛCDM particle horizon, an integral that itself contains Ωm. The 0.12% agreement is a self-consistency of the geometric picture with Planck values. An independent prediction requires deriving r_universe from the flow framework without the Ωm-dependent horizon integral (RHAC Scenario H). The numerical agreement stands; the claim of derivation is withdrawn.

> **The bootstrap (v6.3.2) — a conditional route past C2.** Rather than remove the circularity, solve it. Adopting principle **P** (Ωm = R_H/r_p) closes it into **Ωm·I(Ωm) = 1**, which has a *unique* root: 0.3043 with zero measured inputs; **0.31408** with radiation (inputs = T_CMB, N_eff, h) — **0.96% from Planck, 0.53σ inside its error bar**. It back-predicts r_universe = 4.353×10²⁶ m (−1.07%), so the Ωm-dependent import C2 objected to is no longer needed. Closure: a₀ improves **1.72% → 0.66%**; DESI DR2 gives χ²/N = **1.618** (fixed-ruler). Input ledger: measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ, x₀, r_universe, a₀}. Identity: P ⇔ mean matter pull at the horizon = cH₀/2 — the same cH₀ that sets a₀ (ratio 1.00000).
>
> **This is conditional.** P is *not* derived from A1–A3 (Part B, open — three candidate routes, none attempted). The Gaztañaga comparison (causal-universe scale ≈ 0.3176 H₀) is a prerequisite for any novelty claim and does not yet exist. DESI 1.618 is fixed-ruler and could reorder under C4 marginalization. a₀'s empirical target carries ~10% scatter, so 0.66% is pleasing, not decisive. **Until Part B lands, the C2 downgrade above stands.** See RHAC Scenario Q and `VALIDATION_REPORT.md` §3.0. Scripts: `estif_omega_bootstrap.py`, `estif_bootstrap_closure.py`.

---

## The Two Paths

### Path One — ESTIF-Core (clean) ✅ recommended default
1. Write A2 (universal speed c) and A3 (empty space is not a source) into the theory documents — the fidelity audit found they live only in the derivation scripts.
2. Resolve the A1 conflict: retire the shrinking-ruler narrative in `ESTIF_CONCEPT.md` in favour of the flow (Painlevé–Gullstrand) picture.
3. State cosmology as: constant cosmic eddy → cosmological constant, χ²/N = 1.92 (ties ΛCDM). Retire Ω_tilt(z) to an "explored and set aside" appendix.
4. Submit the gravity letter (unaffected, strengthened by the derived field equation).

### Path Two — ESTIF-Extended (hard) 🔬 high-risk
1. Derive the leading w(z) correction to w = −1 from the full rotating-shear / vorticity stress tensor (the cosmological half of the T_μν work).
2. Target the mild thawing DESI prefers. The two naive reductions (E1, E2) are already proven wrong; the full off-diagonal tensor is required. ⚠️ The nominal χ²/N ≈ 0.66 bar comes from a BAO-alone CPL fit with rd, H₀, Ωm fixed at Planck; fixing nuisance parameters inflates the evolving-DE advantage. Re-derive under marginalization before committing (C4). Note also that a single-speed flow forces p_r = −ρ (L1), so any thawing requires leaving that class — a second metric function or explicit time dependence.

---

## Test Suite Status (new scripts, v6.3)

| Test Script | Purpose | Status |
|---|---|---|
| `estif_tmunu_gauss_codazzi.py` | ADM/Gauss–Codazzi engine, validated vs FRW + de Sitter | ✅ NEW |
| `estif_flow_signature_dynamics.py` | Signature + SR + Newton from Euclidean bulk + universal-c | ✅ NEW (18/18) |
| `estif_converse_flow_law.py` | Vacuum forces v²=2A/r (Birkhoff in flow variables) | ✅ NEW |
| `estif_task4_field_equation.py` | Field equation derived (mass continuity = Poisson) | ✅ NEW (5/5) |
| `estif_task5_desi_selfconsistent.py` | De-circularized DESI DR2 test (10.8 → 3.35) | ✅ NEW |
| `estif_task5b_cosmo_eos.py` | DESI-preferred w(z); tilt tracks it within ~0.05 | ✅ NEW |
| `estif_task6_eddy_eos.py` | Frozen-eddy reframe (1.92 beats tilt 3.35) | ✅ NEW |
| `estif_fidelity_audit.py` | Axiom-presence audit of the corpus | ✅ NEW |

Prior analytical suite (`src/estif_ec_gr_run_simulation.py` → 21/21) unchanged.

---

## Publication Readiness

### Gravity Letter — DRAFTED ✅ and STRENGTHENED
The a₀ derivation now sits on a *derived* field equation (Task 4) rather than a Schwarzschild match. This directly answers the "borrowed recipe" concern. Venue: MNRAS Letters, ApJL, JCAP. Both prior blockers remain resolved.

### Cosmology — RESTRUCTURED
- Path One: publishable now as "constant cosmic eddy → cosmological constant, competitive with ΛCDM." Honest, defensible, no failing tests.
- Path Two: blocked on the vorticity stress-tensor derivation.

---

## Priority Actions

1. ✅ **v6.3.1 errata (C1–C6) applied** (9 July 2026). Script-side C2/C4/C5 remain pending.
2. **Write the Gaztañaga comparison memo** — prerequisite for any Ωm novelty claim (checklist B-4b). Nothing about the bootstrap ships before this.
3. **Attempt Part B**: derive principle P from A1–A3. The single largest available upgrade to Path One's claims — it would convert Ωm, Ω_Λ, r_universe, and a₀ from consistency relations into zero-parameter predictions from three measured inputs.
4. **Preferred-frame / simultaneity section** (barrier 5, checklist C-11) — pre-submission.
5. **Submit the gravity letter** (strengthened by Task 4, wording corrected by C1 + C6).
3. **Rewrite the cosmology sector** around the frozen-eddy = Λ result; move Ω_tilt(z) to an appendix.
4. **(Path Two, optional)** attempt the vorticity stress-tensor derivation of the leading w(z) correction.

---
**Status Document Version:** 6.3.2 / **Last Updated:** 9 July 2026.
