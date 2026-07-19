# ESTIF Development Status

**Last Updated:** 8 July 2026
**Version:** 6.3 — "The Split"
**Status:** Gravity letter ready (and strengthened). Project split into Path One (ESTIF-Core, clean) and Path Two (ESTIF-Extended, hard). See `MILESTONE_v6.3_THE_SPLIT.md`.

> **Target location in repo:** `docs/report/STATUS.md`

---

## Executive Summary

ESTIF is a geometric model deriving gravity, dark energy, and dark matter from the claim that 3D space is a flat hypersurface carried through a 4D bulk. As of v6.3 the project has been restructured into two tracks after a sequence of derivations (Tasks 4–6, July 2026):

- **Path One — ESTIF-Core (clean, recommended default):** gravity now rests on a *derived* field equation; cosmology is a plain cosmological constant that ties ΛCDM. The tilt apparatus is retired.
- **Path Two — ESTIF-Extended (hard, high-risk):** attempt to derive a small thawing correction to w = −1 from the full vorticity stress tensor.

The single most important change since v6.2: **the gravity field equation is now derived from the flow axioms rather than matched to the Schwarzschild solution.**

---

## Three-Sector Status

### Gravity Sector ✅ SOLID — and now on a DERIVED foundation

| Result | Value | Status |
|---|---|---|
| **Field equation derived (mass continuity = Poisson)** | rho_eff = m′(r)/(4πr²) forced by flow metric | ✅ NEW v6.3 (Task 4) |
| **Vacuum → exact Schwarzschild** | v² = 2GM/r unique, full Einstein tensor = 0 | ✅ NEW v6.3 |
| **Weak-field source → Newtonian ρ₀** | uniform ball returns rho0 exactly | ✅ NEW v6.3 |
| Signature + SR derived from Euclidean bulk + universal-c | (−,+,+,+) and dτ/dt = √(1−v²/c²) | ✅ NEW v6.3 (signature suite) |
| MOND a₀ derived (now on derived field eq) | H₀cx₀/√3 = 1.179×10⁻¹⁰ m/s² | ✅ 1.72% match |
| SPARC BTFR (87 galaxies, Qual-1) | RMS = 15.6% | ✅ Within observed scatter |
| EHT M87* shadow | 42.0 μas, 0.00σ | ✅ |
| Planck Λ (local tilt, unchanged) | ratio = 1.0000 | ✅ |
| LISA GW delay | 491 μs, S/N = 49.2σ | ✅ |
| β = τ at n = ½ (x = 0.272) | GR as special case | ✅ |
| a₀ redshift constancy | H(z) cancels exactly | ✅ v6.2 |
| Parameter independence | 3,600 H₀/Ωm combos within SPARC scatter | ✅ v6.2 |

**Open in gravity:** the strong-field pressure/stress sector (relativistic interior with pressure) still needs the full off-diagonal T_μν. Not needed for vacuum, Newton, or Schwarzschild.

### Cosmology Sector 🔄 REFRAMED — frozen eddy ties ΛCDM; tilt retired

| Test | Result | Status |
|---|---|---|
| **Frozen eddy = cosmological constant (DERIVED)** | χ²/N = 1.92 (ties ΛCDM) | ✅ NEW v6.3 — the honest best |
| DESI DR2 — circular Ω_tilt(z) | χ²/N = 10.80 | ❌ superseded |
| DESI DR2 — self-consistent Ω_tilt(z) (Task 5) | χ²/N = 3.35 | 🔄 fixed circularity, still short of ΛCDM |
| Best evolving-w flow (fitted, Task 5b) | χ²/N = 0.66 at w0=−0.85, wa=−0.45 | ℹ️ fitted ceiling / Path Two target |
| Derived eddy E1 (conserved-L spin) | χ²/N = 3232 | ❌ falsified (Task 6) |
| Derived eddy E2 (tracker) | χ²/N = 754 | ❌ falsified (Task 6) |
| Age of universe | 13.379 Gyr | ✅ Stands |

**Root cause of the tilt failure (now understood):** the tilt *shape* fits worse than the frozen-eddy limit underneath it. The Ω_tilt apparatus (N_MAX, B, sign-flip, z<2 cutoff) is a net negative on DESI and is **retired** from the cosmology claim under Path One.

### Dark Matter Sector 🟡 Analytical phase complete, simulation wall

Unchanged from v6.2. Ωm = x₀ (0.12%), Ωdm = x₀ − Ωb (0.10%), σ/v_esc = 0.5 exact, λ_Jeans = 2.565r. v_flat = 220 km/s requires δ ~ 50,000–100,000 — N-body, off-Mac. The stress-energy derivation of ρ_eddy = x₀ρ_crit is now connected to the Task 4 machinery (the local field equation is derived; the homogeneous-source version is the same calculation).

---

## The Two Paths

### Path One — ESTIF-Core (clean) ✅ recommended default
1. Write A2 (universal speed c) and A3 (empty space is not a source) into the theory documents — the fidelity audit found they live only in the derivation scripts.
2. Resolve the A1 conflict: retire the shrinking-ruler narrative in `ESTIF_CONCEPT.md` in favour of the flow (Painlevé–Gullstrand) picture.
3. State cosmology as: constant cosmic eddy → cosmological constant, χ²/N = 1.92 (ties ΛCDM). Retire Ω_tilt(z) to an "explored and set aside" appendix.
4. Submit the gravity letter (unaffected, strengthened by the derived field equation).

### Path Two — ESTIF-Extended (hard) 🔬 high-risk
1. Derive the leading w(z) correction to w = −1 from the full rotating-shear / vorticity stress tensor (the cosmological half of the T_μν work).
2. Target the mild thawing DESI prefers (χ²/N ≈ 0.66). The two naive reductions (E1, E2) are already proven wrong; the full off-diagonal tensor is required.

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

1. **Adopt Path One as the default** and write A2 + A3 into the theory; retire the shrinking-ruler narrative.
2. **Submit the gravity letter** (strengthened by Task 4).
3. **Rewrite the cosmology sector** around the frozen-eddy = Λ result; move Ω_tilt(z) to an appendix.
4. **(Path Two, optional)** attempt the vorticity stress-tensor derivation of the leading w(z) correction.

---

**Status Document Version:** 6.3
**Last Updated:** 8 July 2026
