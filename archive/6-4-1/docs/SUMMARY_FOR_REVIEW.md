# ESTIF v6.3 — Summary for Expert Review

**Author:** Peter Angelov (Independent Researcher, tervion@gmail.com)
**Version:** 6.4.0 (July 2026) — "The Split" + Phase 2 null + A1′
**Repository:** https://github.com/tervion/estif-publication
**Zenodo:** https://zenodo.org/records/17261724
**Validation:** `python3 tests/estif_task4_field_equation.py` | `python3 tests/derive_mond_from_geometry.py` | `python3 src/estif_ec_gr_run_simulation.py` (21/21)

> **See also:** `MILESTONE_v6.3_THE_SPLIT.md`

---

## The Single Claim

3D space is a flat sheet carried through a 4D bulk at the speed of light. The local
tilt and slowing of that flow near mass is gravity; motion through the 4th dimension
is the passage of time. In v6.3 this is made precise by three axioms — flat slices
(A1), universal speed c (A2), empty space sources nothing (A3) — from which gravity
is now *derived*.

---

## What Is New in v6.3

### HEADLINE: The gravitational field equation is DERIVED, not borrowed

Every prior version reproduced Newton's law by matching the flow profile to the
Schwarzschild solution (setting the tilt exponent n = ½ so β = √(1 − Rs/r)) and
differentiating — which presupposes the answer. v6.3 removes this. From the three
flow axioms alone, a symbolic Gauss–Codazzi / ADM engine **forces**

```
ρ_eff = m′(r) / (4π r²)      (mass continuity = Poisson, integrated form)
```

with these engine-verified consequences:
- **vacuum (A3: ρ_eff = 0) → v² = 2GM/r uniquely → exact Schwarzschild** (full
  Einstein tensor vanishes, no weak-field approximation);
- **a uniform matter ball → ρ_eff = ρ₀ exactly** (correct Newtonian source);
- **control:** a volume-conserving ("plughole") flow gives a 1/r⁵ force and fails —
  the vacuum condition, not tuning, selects the Newtonian law.

Lorentzian signature (−,+,+,+) and special-relativistic kinematics emerge from the
same Euclidean-bulk-plus-speed-c construction. The former "Poisson postulate" is now
a theorem for vacuum, Newton, and Schwarzschild.

**Precision on "derived" (C6).** What the axioms establish is that they uniquely
*select* the constraint (energy) sector of General Relativity in Painlevé–Gullstrand
gauge, forcing mass continuity and hence exact Schwarzschild in vacuum — **without
matching to the Schwarzschild solution**, which was the previous gap. What is
*adopted*, not derived from below, is the gravitational coupling: the identification
of the geometric constraint scalar with 8πG × energy density. ESTIF does not derive
G or the factor 8π. This qualifier is the honest form of the headline, and it is the
subject of reviewer question #1.

**Open:** the strong-field pressure/stress sector (full off-diagonal T_μν), needed
for none of {vacuum, Newton, Schwarzschild}.

Scripts: `estif_task4_field_equation.py` (5/5), `estif_flow_signature_dynamics.py`
(18/18), `estif_converse_flow_law.py`, `estif_tmunu_gauss_codazzi.py`.

### Cosmology reframed: frozen eddy = cosmological constant

The evolving Ω_tilt(z) dark-energy law failed DESI DR2 (χ²/N = 10.8). Two results
changed the picture:
- **De-circularization** (Task 5): replacing the ΛCDM ruler in x(z) with ESTIF's own
  H(z) drops the score to **3.35** — a correctness fix, still short of ΛCDM.
- **Frozen-eddy reframe** (Task 6): a constant eddy energy density (derived from the
  field equation → de Sitter, w = −1) scores **χ²/N = 1.92, tying ΛCDM and beating
  the tilt formula's 3.35**. Two attempts to *derive* an evolving eddy EoS failed
  badly (conserved spin w = +1 → χ²/N = 3232; expansion tracker → 754).

**Conclusion:** on DESI, the entire Ω_tilt apparatus is a net negative relative to
the plain cosmological constant underneath it, so the cosmological dark-energy claim
is retired to an appendix. The tilt formula's strong-field use is unaffected.
Scripts: `estif_task5_desi_selfconsistent.py`, `estif_task5b_cosmo_eos.py`,
`estif_task6_eddy_eos.py`.

### The project splits into two paths

- **Path One — ESTIF-Core (clean, recommended):** derived gravity + frozen-eddy
  cosmological constant (ties ΛCDM). Publishable now.
- **Path Two — ESTIF-Extended (hard):** derive the small thawing correction to
  w = −1 that DESI hints at, from the full vorticity stress tensor. The naive
  reductions are already falsified. ⚠️ The nominal target χ²/N ≈ 0.66 comes from a
  BAO-alone CPL fit with rd, H₀, and Ωm held fixed at Planck values; fixing nuisance
  parameters inflates the evolving-DE advantage, so this bar must be re-derived under
  marginalization before it is used to justify Path Two (C4).

### Fidelity audit and consolidation

`estif_fidelity_audit.py` found that axioms A2 and A3 previously lived only in the
derivation scripts, and that `ESTIF_CONCEPT.md` ran a competing "shrinking-ruler"
narrative. v6.3 writes A2 and A3 into the theory, retires the shrinking-ruler
picture in favour of the flow (Painlevé–Gullstrand) picture, and relabels
v_flow = cx₀ ≈ 0.31c as the *sideways component* of a total-c flow (not the full
speed).

---

### The Ωm bootstrap (v6.3.2 — conditional)

Correction C2 showed that `Ωm = x₀` is circular: r_universe is the ΛCDM particle
horizon, an integral containing Ωm. The bootstrap solves the circularity instead of
removing it. Adopting principle **P** (`Ωm = R_H/r_p`) closes it into Ωm · I(Ωm) = 1
which has a **unique** root: 0.3043 with zero measured inputs; **0.31408** with
radiation included (inputs = T_CMB, N_eff, h) — **0.96% from Planck, 0.53σ inside
its error bar**. It back-predicts r_universe = 4.353×10²⁶ m (−1.07%), eliminating
the Ωm-dependent import C2 objected to.

**Closure:** propagating the bootstrap Ωm gives a₀ = 1.1920×10⁻¹⁰ m/s², improving
MOND agreement **1.72% → 0.66%** (SPARC insensitive, v_flat × 1.00269), and DESI DR2
χ²/N = **1.618** vs ΛCDM's 1.919 — *within the fixed-(H₀, rd) test*.

**Input ledger after adopting P:** measured = {H₀, T_CMB, N_eff}; computed =
{Ωm, Ω_Λ, x₀, r_universe, a₀}.

**Identity:** P ⇔ mean matter pull at the horizon = cH₀/2 — the same cH₀ that sets
a₀; the ratio computes to 1.00000.

**This is conditional and is presented as such.** P is *not* derived from A1–A3;
that is Part B, and it is untouched (three candidate routes, none attempted). The
adjacency to Gaztañaga's causal-universe scale (≈ 0.3176 H₀, via inflation) is
unresolved and blocks any novelty claim. The DESI figure is fixed-ruler. a₀'s
empirical target carries ~10% scatter, so 0.66% is pleasing, not decisive. Scripts:
`estif_omega_bootstrap.py`, `estif_bootstrap_closure.py`.

---

## The Tilt Formula (strong-field deviations — unchanged)

```
x        = curvature ratio (Rs/r locally)
n(x)     = 33.265 × exp(−15.429 × x)
β(x)     = √(1 − x^(2n(x)))
Observable = √β(x)
```
Parameters connect to the classical electron radius: N_MAX ≈ 5/7 × ln(r_e/l_P)
(0.08%, conditional on x_c), B ≈ 1/3 × ln(r_e/l_P) (0.69%, derived from 3D isotropy).

---

## The MOND Derivation (v6.1, now on the derived foundation)

```
Step 1 — Force law:  a = −c²/2 ∇(ω/H₀)² = GM/r²    [underwritten by the derived field eq]
Step 2 — Cosmic flow (sideways component):  v_flow = c x₀ = c Ωm
Step 3 — 3D isotropic projection:  v_3D = v_flow/√3
Step 4 — Threshold:  a₀ = H₀ c x₀/√3 = 1.179×10⁻¹⁰ m/s²   (1.72% from empirical, zero params)
                     with bootstrap Ωm: 1.1920×10⁻¹⁰       (0.66%, conditional on P)
```
**SPARC:** 87 quality-1 galaxies, RMS = 15.6% (within BTFR scatter). **a₀ redshift
constancy:** H(z) cancels in the comoving frame (deviation 2.22×10⁻¹⁶).
**Parameter independence:** 3,600 H₀/Ωm combinations within ±20% SPARC scatter.

---

## Structure (v6.3)

| Sector | Status |
|---|---|
| **Gravity — field equation** | ✅ **DERIVED** (vacuum → exact Schwarzschild) |
| Gravity — a₀/MOND/SPARC | ✅ Solid, now on a derived foundation |
| Gravity — EHT/Λ/LISA strong-field | ⚠️ Consistent; *deviation* claims conditional on the un-derived eddy-stress sector (C1) |
| **Cosmology (Path One)** | ✅ Frozen eddy = Λ, ties ΛCDM (χ²/N = 1.92) |
| Cosmology — Ω_tilt(z) | 🔴 Retired (net negative on DESI) |
| Cosmology (Path Two) | 🔬 Open (derive thawing from vorticity T_μν) |
| Dark matter | 🟡 Analytical done; N-body wall |
| **Ωm bootstrap** | 🔶 Conditional on P: unique root 0.31408 (0.53σ); Part B open |

---

## Strong-Field Calibration (numbers unchanged; status qualified — C1)

| Observation | ESTIF Prediction | Result |
|---|---|---|
| EHT M87\* shadow | 42.0 μas | ✅ consistent, 0.00σ |
| Planck Λ (local tilt) | 1.1056 × 10⁻⁵² m⁻² | ✅ ratio = 1.0000 (calibration match) |
| LISA GW delay (65 M☉) | 491 μs | ⚠️ conditional, S/N = 49.2σ |

> ⚠️ **Conditional (C1).** The ESTIF vacuum is exactly Schwarzschild, so any
> *deviation* from GR in photon-sphere shadows or in vacuum GW propagation must be
> sourced by the non-vacuum eddy background — a sector not yet derived. The
> observations remain *consistent* with ESTIF; the predicted *deviation from GR* is
> conditional on the open eddy-stress derivation. The Λ entry is a calibration match,
> not a vacuum deviation, and is unaffected.

---

## What This Is NOT

- Not curve-fitting: gravity's field equation and a₀ are derived.
- Not numerology: the vacuum condition forces v²=2GM/r; the 1/√3 and 1/3 factors
  come from 3D isotropy.
- Not a complete cosmological replacement: the honest cosmology is a cosmological
  constant (ties ΛCDM); the evolving version is retired pending Path Two.

---

## Questions for Expert Review

**Theoretical:**
1. The field equation is derived from the flow axioms via Gauss–Codazzi (mass
   continuity = Poisson; vacuum → exact Schwarzschild). Is the derivation, as
   presented, a genuine first-principles result, or does the identification of the
   flow speed with proper time smuggle in GR structure?
2. The frozen-eddy limit gives w = −1 for free (constant eddy density → de Sitter).
   Is deriving a *small thawing correction* from the rotating-shear/vorticity stress
   tensor (Path Two) a well-posed program, given that conserved-spin and tracker
   reductions both fail against DESI?
3. What is the geometric property of Schwarzschild spacetime at r = 3.68 Rs that
   gives x_c = 0.272? (Closes the N_MAX = 5/7 × L chain.)
4. Can ρ_eddy = x₀ ρ_crit be derived as the homogeneous version of the (now derived)
   local field equation?

**Observational:**
5. On DESI DR2 BAO alone, a cosmological constant ties ΛCDM (1.92) and a best-fit
   evolving w reaches 0.66 *with rd, H₀, Ωm fixed at Planck*; ESTIF's self-consistent
   tilt (3.35) sits near DESI's own published-parameter fit on this subset (3.09). Is
   the frozen-eddy = Λ claim the right one to publish now, with the thawing correction
   as future work — and how much of the 0.66 advantage survives marginalization over
   rd, H₀, and Ωm?
6. The SPARC zero-bias Υ* = 0.85 is above the McGaugh+2014 standard of 0.50. Is 0.85
   within the plausible 3.6 μm mass-to-light range?

**Publication:**
7. Is "derived gravitational field equation from a flowing-hypersurface model +
   geometric derivation of a₀ + SPARC validation" a viable standalone letter for
   MNRAS Letters, ApJL, or JCAP?
8. What framing best distinguishes this from MOND variants, given that a₀ is derived
   from cosmological constants and the field equation from the flow axioms rather
   than by modifying the force law directly?

**On the bootstrap:**
9. ~~Is principle P derivable from A1–A3?~~ **Resolved internally 12 July 2026
   (RHAC-009): NO.** P holds only at our epoch (Ωm(a) = R_H/d_p crosses once, at
   a≈1; →1 vs ½ at high z), so it is not a theorem of the time-symmetric axioms —
   it is a predictive postulate. The cH₀/2 equivalence is that today-condition, not
   a mechanism. Remaining live question for reviewers: how does the postulated root
   (0.31408) relate to Gaztañaga's causal-boundary scale (≈ 0.3176 H₀, via
   inflation)? The comparison memo is still owed (B-4b).

   See docs/plan/GAZTANAGA_COMPARISON.md: the memo concludes there is no defensible novelty claim
   on Ωm (Gaztañaga has priority); ESTIF's distinctive thread is the a₀ = cH₀/2 link.

---

## Honest Assessment

**What is solid:**
- Gravitational field equation derived from the flow axioms: mass continuity =
  Poisson; vacuum → exact Schwarzschild; uniform source → Newtonian ρ₀.
- Lorentzian signature and SR kinematics emerge from the same construction.
- a₀ derived from geometry (1.72%, zero params), 87-galaxy SPARC RMS 15.6%.
- Combined strong-field formula: EHT + Λ + LISA, zero free params.
- GR time dilation as the n = ½ special case; B = L/3 from isotropy.
- Ωm = x₀ (0.12%), Ωdm = x₀ − Ωb (0.10%) — but as a **consistency relation**, not an
  Ωm-independent prediction: r_universe is the ΛCDM particle horizon, which itself
  contains Ωm (C2). The numerical agreement stands; the claim of derivation does not.

**What is conditional (new in v6.3.2):**
- The Ωm bootstrap: unique root 0.31408 (0.53σ from Planck), a₀ → 0.66%, r_universe
  back-predicted to −1.07%, DESI 1.618. **All of it rests on principle P, which is
  not derived — and (RHAC-009, 12 Jul 2026) is not derivable as a law.** §C2's
  downgrade is therefore final: Ωm is fixed by a predictive postulate, not derived.
  P's falsifiable content (Ωm = 0.3141) stands as a postulate.
- The Gaztañaga comparison is unwritten. No priority is claimed on the Ωm result.
- a₀ redshift constancy (algebraic); parameter independence (3,600 combinations).
- **Cosmology:** frozen eddy = cosmological constant ties ΛCDM on DESI DR2 (1.92).

**What was retired (honest):**
- Evolving Ω_tilt(z): net negative on DESI vs the cosmological-constant limit
  (10.8 circular → 3.35 self-consistent → 1.92 frozen). Retired to an appendix.
- The two naive derived eddy EoS models (stiff spin, tracker): falsified against DESI.
- The Λ-drift and Hubble-bridge predictions: no longer load-bearing.

**What remains open:**
- Strong-field pressure/stress sector (full off-diagonal T_μν).
- ~~Path Two: derive the leading w(z) thawing correction~~ — CLOSED (Phase 2, honorable null,
  11 July 2026): the residual dark-energy sector the axioms permit is provably empty; w = −1
  exactly; Λ imported. Any evolving-w must come from outside the flow sector. docs/plan/PHASE2_DECLARATION.md.
- x_c = 0.272 geometric origin; ρ_eddy = x₀ρ_crit homogeneous derivation.
- N-body simulation for halo structure; CMB (standard ΛCDM check on Path One).
- Not peer-reviewed.

---

## Files to Examine

1. `tests/estif_task4_field_equation.py` — field equation derived (key new result)
2. `tests/estif_flow_signature_dynamics.py` — signature + SR + Newton from the axioms
3. `tests/estif_converse_flow_law.py` — vacuum forces v²=2A/r (Birkhoff in flow variables)
4. `tests/estif_tmunu_gauss_codazzi.py` — ADM engine, validated vs FRW + de Sitter
5. `tests/estif_task5_desi_selfconsistent.py` — de-circularized DESI test (10.8 → 3.35)
6. `tests/estif_task6_eddy_eos.py` — frozen-eddy reframe (1.92 ties ΛCDM)
7. `tests/estif_task5b_cosmo_eos.py` — DESI-preferred w(z)
8. `tests/estif_fidelity_audit.py` — axiom-presence audit
9. `tests/derive_mond_from_geometry.py` — the 4-step MOND derivation
10. `tests/test_sparc_tully_fisher.py` — 87-galaxy validation
11. `tests/test_a0_redshift.py`, `tests/test_a0_parameter_independence.py` — a₀ robustness (v6.2)
12. `src/estif_ec_gr_model.py` — core implementation
13. `docs/report/ESTIF_CONCEPT.md` — conceptual foundation (v6.3, axioms A1–A3)
14. `MILESTONE_v6.3_THE_SPLIT.md` — the split and honest status
15. estif_omega_bootstrap.py
16. estif_bootstrap_closure.py

### What is new in v6.4 (11 July 2026)

- **Phase 2 closed — honorable null (strong form).** The residual dark-energy stress permitted
  by A1–A3 was computed in full (pre-registered) and is provably empty; w = −1 exactly; Λ is a
  bare imported constant. ESTIF does not derive dark energy — and it is now proven it cannot,
  from these axioms. docs/plan/PHASE2_DECLARATION.md (RHAC-004).
- **Axiom A1 → A1′.** Slices are even *on average*; matter sources local dents; the mean spatial
  curvature is identically zero at all epochs (registered falsifier, current 0.0007 ± 0.0019).
  This restores the growing density mode (D₊ = GR growth, f(z=0.5) = 0.76) and the radiative
  sector — gravitational waves at speed = c, forced by A2, GW170817-consistent (RHAC-005/006).
  All strict-A1 results survive as the exact-evenness limit.
- **Gaztañaga comparison delivered.** docs/plan/GAZTANAGA_COMPARISON.md — the causal-horizon →
  Ωm ≈ 0.3 / no-DE result is Gaztañaga's (peer-reviewed, 2019+). No ESTIF novelty on Ωm; the
  a₀ link is the only distinctive thread.

**Document Version:** 6.4.0 | **Updated:** 12 July 2026
