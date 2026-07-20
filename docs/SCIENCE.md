# ESTIF — The Framework and the Evidence

**Version:** 6.4.2 · **Last updated:** 21 July 2026
**Author:** Peter Angelov (Independent Researcher) · tervion@gmail.com
**Repository:** https://github.com/tervion/estif-publication · **Zenodo:** https://zenodo.org/records/17261724
**Reproduce the headline results:** `python3 tests/scripts/estif_task4_field_equation.py` (5/5) · `python3 tests/scripts/derive_mond_from_geometry.py` · `python3 src/estif_ec_gr_run_simulation.py` (21/21)

This file consolidates and supersedes ESTIF_CONCEPT.md, VALIDATION_REPORT.md, SUMMARY_FOR_REVIEW.md, PHASE2_DECLARATION.md, GAZTANAGA_COMPARISON.md, JWST_TEST_SPEC.md, and the 13 July doctrine memo. Every claim below carries its epistemic status (derived · adopted · conditional · consistency relation · retired) and its receipt. Project status and the item-level checklist live in `docs/PROJECT.md`; decisions in `docs/RHAC.md`.

---

## 1. The single claim

3D space is an even sheet carried through a 4D bulk at the speed of light. The local tilt and slowing of that flow near mass is gravity; motion through the 4th dimension is the passage of time; the expansion we observe is the geometric projection of the 4D motion onto the sheet. Everything below follows from three axioms about that moving sheet — and, since v6.3, gravity is *derived* from them rather than matched to known solutions.

Two things ESTIF does **not** claim: it does not eliminate dark energy (Λ enters imported, and Phase 2 *proves* it must — §5), and it does not derive Ωm (the bootstrap is conditional on a postulate proven non-derivable — §7).

---

## 2. The three axioms

**A1′ — Space flows, it does not stretch** *(amended 11 July 2026; RHAC-006).*
The 3D sheet is even **on average**: its internal geometry carries no curvature of its own except where matter sources a local deviation (a dent); where nothing sources a dent, exact evenness returns. The spatial average of the curvature is identically zero at all epochs — this is *law*, not an initial condition, and it is the registered kill-shot falsifier (current measurement: 0.0007 ± 0.0019). This is the Painlevé–Gullstrand picture: even spatial slices plus a flow — not space stretching like rubber, and not universal shrinkage (that older picture is retired).

> **Why the amendment.** The strict form (exactly even, everywhere, always) was proven to forbid both the growing density mode — galaxies could not form (no-go theorem, RHAC-005) — and the radiative sector (vs LIGO). A1′ pulls exactly that one over-constraint. All strict-A1 results survive as the exact-evenness limit: vacuum Schwarzschild, the a₀ chain, the Friedmann background, and the Phase 2 null. What A1′ adds: linear growth (f(z=0.5) = 0.76, matching GR and DESI RSD), the gravitational-wave sector, and a container for the ζ = 10⁻⁵ seed (value still imported). Receipts: `tests/scripts/phase2_a1prime/`.

**A2 — Everything moves through the bulk at the speed of light.**
Every object moves through the 4D bulk at exactly *c*. "Sitting still" is moving through the 4th dimension — through time — at full speed; moving through space diverts part of that speed sideways, so a moving clock runs slow. The total motion splits as `dw² + dσ² = c²dt²`, and proper time *is* bulk distance: `c·dτ = dw`.

**A3 — Empty space is not a source.**
Vacuum sources nothing. This is the axiom that forces Birkhoff-type uniqueness in the flow variables and, with the Phase-2 null, forbids any second metric or field.

---

## 3. Gravity is derived, not borrowed *(v6.3 headline; status: derived, with one adopted element)*

Earlier versions reproduced Newton by *matching* the flow profile to the Schwarzschild solution (setting the tilt exponent n = ½) — which presupposes the answer. v6.3 removes that. From the three axioms alone, a symbolic Gauss–Codazzi / ADM engine **forces** the field equation

```
ρ_eff = m′(r) / (4π r²)        (mass continuity = Poisson, integrated form)
```

with engine-verified consequences: **vacuum → exact Schwarzschild** (v² = 2GM/r unique; full Einstein tensor = 0 — Birkhoff in flow variables), and **a uniform matter ball → ρ₀ exactly** (the correct Newtonian source). Lorentzian signature (−,+,+,+) and exact SR kinematics emerge from the same construction (Euclidean bulk + universal-c), with the minus sign *produced*, not inserted.

> **Precision (C6 — the one adopted element).** What is derived is that the axioms uniquely select the constraint (energy) sector of General Relativity in Painlevé–Gullstrand gauge — forcing mass continuity, hence exact Schwarzschild in vacuum and the Newtonian source in the weak field — *without* matching to the Schwarzschild solution. The gravitational **coupling** (geometric constraint scalar ↔ 8πG × energy density) is **adopted** from Einstein–Hilbert, not derived from below.

Gravity, in this picture, is generalized time dilation: static clock rate = local flow speed u(r); gravitational pull = −c²∇u. The strong-field tilt formula (dynamic n, √β projection, electron-scale parameters) remains the tool for *conditional* deviations from GR (§9, C1).

Receipts: `tests/scripts/estif_task4_field_equation.py` (5/5) · `estif_flow_signature_dynamics.py` (18/18) · `estif_converse_flow_law.py` · `estif_tmunu_gauss_codazzi.py` (engine validated vs FRW + de Sitter).

---

## 4. a₀ — the horizon doctrine *(13 July 2026; RHAC-010. Status: SCALE derived; NUMBER open)*

**The doctrine.** a₀ is a horizon-scale acceleration, a₀ ≈ c·H (de Sitter surface gravity). From this alone follow: even (asymptotically constant) rotation curves, the mass-independence of a₀, the Baryonic Tully–Fisher relation v⁴ = G·M·a₀, and the correct order of magnitude (cH/2π ≈ 1.0×10⁻¹⁰ vs observed 1.2×10⁻¹⁰ m/s²). Algebraically, c²√Λ = √3·√Ω_Λ·cH₀ — any "Λ-native" form of a₀ is the Hubble-scale form up to an O(1) factor.

**Retired framings (RHAC-010):** any language stating that the √3 factor, or the full a₀ coefficient, is *derived* from local flow or matter dynamics. The local flow provably cannot make a₀ — the local construction yields (H/2)·v_gal, roughly a thousand times too weak and mass-dependent (`tests/scripts/estif_flow_sim.py`); the horizon background can, and does (`tests/scripts/estif_horizon.py`).

**The working form** a₀ = H₀·c·x₀/√3 = 1.179×10⁻¹⁰ m/s² (1.72% from the empirical 1.2×10⁻¹⁰; 1.1920×10⁻¹⁰ and 0.66% with the bootstrap Ωm — conditional, §7) remains the operational formula; its 1/√3 has the standard equipartition rationale, and Step 1's force law is underwritten by the derived field equation — but under the doctrine the exact prefactor is **constrained to O(1) and known to be horizon-set, not derived**. The scale of a₀ (c·H) is derived; the number is not.

**Evidence for the phenomenology (all with the working form, zero free parameters):**

| Test | Result | Receipt |
|---|---|---|
| SPARC BTFR, 87 quality-1 galaxies (Lelli+2016, Υ* = 0.50) | RMS = 15.6% (within observed scatter); 82% within 20%; −7.6% mean bias traces to the stellar-mass calibration | `tests/scripts/test_sparc_tully_fisher.py` · `derive_mond_from_geometry.py` |
| Weak-lensing BTFR (Mistele+2024, ApJL 969 L3) | predictions ~5–13% low; χ²/N ≈ 1.5 with the 0.1 dex correlated M*/L systematic (5.6–8.0 stat-only) — consistent within systematics, degenerate with the calibration | `tests/scripts/btfr_lensing.py` (11 Jul 2026) |
| Implied-a₀ inversion, systematics fully correlated | +0.7σ (LTG) to +2.0σ (ETG) above the derived value — mild tension, recorded honestly | `tests/scripts/a0_tension_corrected.py` |
| Redshift constancy | H(z) cancels exactly (deviation 2×10⁻¹⁶); consistent with high-z TF at z ≈ 0.75–2.2 | `tests/scripts/test_a0_redshift.py` |
| Parameter independence | 3,600 H₀/Ωm combinations within SPARC scatter; 8 published datasets pass | `tests/scripts/test_a0_parameter_independence.py` |

**The keystone open problem** (the project's central theory target): **(a)** derive the O(1) prefactor from horizon geometry with zero fit freedom — the data-preferred k ≈ 0.128 in a₀ = k·c²√Λ, the cH/2π candidate, and the √3 slot bracket it; a river-model edge calculation lands ~20% high, exactly the spread between candidate acceleration definitions; **(b)** derive the **local mechanism** — what makes the horizon scale c·H act at galaxy radii? Every emergent-gravity approach (Verlinde, Padmanabhan) *assumes* this step; deriving it from the 4D inward flow would be ESTIF's genuine contribution. `tests/scripts/ripple_speed.py` carries the prepared `estif_eom` slot for the equation-of-motion insertion. Receipts for the doctrine: `tests/scripts/a0_horizon_test.py` · `a0_prefactor_derivation.py`.

**Flag — x_c double duty (open).** x_c = 0.272 is used both as the Einstein-mode → galaxy-mode crossover (n(x_c) = ½; r = 3.68 Rs) and, informally, as the a₀-defining edge. Under the simplest edge law a₀ = c·H_Λ·x_c the galaxy data prefer x_c ≈ 0.221, not 0.272 — these may be two different surfaces that have been conflated. Sort by SETTLED vs OPEN before next use (RHAC-010).

**μ is underived.** Differentiating the dynamic-n force law locally sends g_obs/g_N → 0 in the galactic range x ~ 10⁻⁶–10⁻⁸ — the opposite of MOND's rising interpolation; Newton requires n frozen at ½, and no rule fixes where n is evaluated. Standing rule (RHAC-001): no per-gradient formulas — behaviour gradients must emerge as term-dominance of a single equation. Receipt: `tests/scripts/mu_extraction.py`.

---

## 5. Cosmology — the honest reframe, and the Phase 2 proof *(status: Λ imported; residual sector proven empty)*

**The reframe (v6.3).** The evolving Ω_tilt(z) dark-energy law failed DESI DR2 (χ²/N = 10.8 circular); de-circularizing it helped (3.35, Task 5) but revealed the tilt *shape* is the problem — it fits worse than the constant-Λ limit underneath it. That limit — a plain cosmological constant, free from the gravity sector — ties ΛCDM at **χ²/N = 1.92** (ΔAIC = +0.59 with Ωm = x₀; ΔAIC = +18 against the tilt, decisive). The tilt apparatus (N_MAX, B, sign-flip, z < 2 cutoff) is **retired** from the cosmology claim. Receipts: `tests/scripts/estif_task5_desi_selfconsistent.py` · `estif_task5b_cosmo_eos.py` · `estif_task6_eddy_eos.py`.

**Phase 2 — the residual sector is EMPTY (11 July 2026; honorable null, strong form).** Path Two's hope was a small residual stress beyond w = −1 matching DESI's mild thawing hint. Two naive reductions were already falsified (E1 conserved-spin: χ²/N = 3232; E2 tracker: 754). Phase 2 asked the complete question, pre-registered:

> **Binding declaration (fixed before the derivation).** The total residual stress admitted by A1–A3 is to be computed in full and submitted **unmodified**. No sector may be added, dropped, or reinterpreted after seeing the result.

The census, under strict A1 — four sectors:

| Sector | Physical content | Verdict | Why |
|---|---|---|---|
| Shape | anisotropic spatial curvature of the slice | **FORBIDDEN** | excluded by A1 (even slices) — axiomatic, no script needed |
| Rate | time-varying lapse carrying stress | **FORBIDDEN** | A2: the lapse carries no stress | 
| Slosh | bulk spatial flow / momentum density | **ZERO** | pure-divergence theorem — integrates to nothing |
| Swirl | vorticity / rotation of the flow | **FORBIDDEN free; NEGATIVE if sourced** | ⟨δρ⟩ = −⟨ω²⟩/32πG, slaved to matter, ~10¹¹× below ρ_Λ |

**Verdict:** the residual dark-energy sector is **EMPTY**. w = −1 **exactly**; Λ enters as a **bare imported constant**, exactly as in ΛCDM; the 1.92 parity is unaffected; the χ²/N ≈ 0.66 thawing ambition is dead (and stands pre-registered as the bar — `tests/scripts/estif_pathtwo_target_lock.py`, receipt `tests/docs/pathtwo_target_lock_output.txt` — that any *outside-the-axioms* proposal must beat under honest marginalization). This is a strong-form null: not "not found yet" but "provably permitted none" — it converts "we haven't derived dark energy" into "our axioms cannot manufacture dark energy, by construction."

**Retraction.** The "frozen eddy" reading of the constant (a spinning cosmic flow held constant) is retracted: auditor-introduced, author-objected before the derivation, and falsified by the constraint algebra — the swirl sector is exactly the one forbidden free and negative when sourced. The constant survives; the spin reading of it does not.

**The LIGO flag and its resolution.** Strict A1 carries no free radiative modes (no gravitational waves — vs LIGO), the same over-constraint that forbids the growing density mode. Resolution: **A1 → A1′** (§2). The Phase-2 census was re-run under A1′ across **five** sectors (the four above plus the dent/shape-average channel A1′ opens): no Λ-printer appears — **the null stands** (re-census receipt: `tests/scripts/phase2_a1prime/estif_a1prime_recensus.py`, PASS 5/5).

Receipts: door suite `tests/scripts/phase2_a1prime/estif_p2_door1_rate_dial.py`, `_door2_slosh_divergence.py` (the strongest receipt — a theorem), `_door3_swirl_ledger.py`; strict census + A1′ re-audit: `archive/6-4-2 diagnostics/test_UKN.py`, `test_UKN2.py`; growth pair: `estif_growth_nogo_law.py` + `estif_growth_nogo_audit.py` (no-go, five audit attacks bounced) and `estif_a1prime_growth_restored.py` + `estif_a1prime_deepen_exact.py` (restoration, numeric + exact symbolic). Known flags (do not overstate): door 1 asserts the standard lapse-carries-no-stress fact rather than deriving it; door 3's ~10¹¹ suppression is an order-of-magnitude estimate; the LTB correspondence in the audit is literature-backed, not machine-verified.

---

## 6. Growth, gravitational waves, and the four-lock ledger *(fronts 1–3 + C-15; status: derived/computed, honest about degeneracy)*

**Front 1 — growth and JWST (the early-structure test, absorbing the test spec).** The claim under test: massive galaxies at z ≈ 8–13 look too mature for the elapsed time; could ESTIF's growth be faster? The observable is the cumulative comoving number density of massive halos, n(>M, z) for M ~ 10¹¹–10¹² M☉ at z ~ 8–12 (a timing/mass claim — the apparent-size version cannot distinguish, since ESTIF-Core's H(z) is ΛCDM's). Under A1′ the growth factor is D₊ = H·∫da/(aH)³ — GR's linear growth exactly; f(z=0.5) = 0.7603 = the DESI RSD anchor; fσ₈ pulls −0.12σ / +1.03σ. **JWST verdict: outcome 4, honest null** — g(9.1) = 0.998; the too-big-too-early tension is *inherited* from ΛCDM, neither relieved nor worsened. The strict-A1 counterfactual would have been falsified outright (decaying-only growth, f = −0.91) — the A1′ fork rescued this test. Receipts: `tests/scripts/estif_front1_growth_sigma8_jwst.py` · `estif_jwst_growth_spec.py`.

**Front 2 — first black holes.** Fold-back rule ν·σ(M,0)·D(z) = δ_c: the star channel yields ~10² M☉ holes at t ~ 28 Myr (earliest-in-volume) / ~212 Myr (typical 3σ); the direct-collapse channel yields 10⁴–10⁶ M☉ seeds — the only comfortable route to 10⁹ M☉ by z = 7; the primordial channel is closed under the ζ = 10⁻⁵ passport. All numbers are inherited (ΛCDM growth); ESTIF's content is structural — strict A1 forms **no** hole at all. Receipt: `tests/scripts/estif_front2_first_hole_recipe.py`.

**Front 3 — the second lock.** The empty residual sector plus GR-equivalent D₊ force the growth index γ ≈ 0.55 (computed 0.5455–0.5544 over z = 0–5); measured 0.58 ± 0.11 (DESI PV+ShapeFit), pull −0.23σ — PASS. Companion slip lock (Σ, η, μ) = (1, 1, 1) exact. Receipt: `tests/scripts/estif_front3_second_discriminator.py`.

**C-15 — c_gw = c derived (RHAC-008).** Under A1′+A2 the world is one Lorentzian geometry; A3 plus the empty residual sector forbid any second metric or field. TT waves are ripples *of* that geometry and light rides *its* null cones: the vacuum wave operator's principal symbol is g^{μν}k_μk_ν — the light cone — shown symbolically for an arbitrary flow (radial null speed u = v ± c, identical for GW and light). GW170817's |c_gw/c − 1| ≲ 10⁻¹⁵ is passed **structurally**: no dial exists to break it. The same A1′ hinge that opened growth opened radiation — one mechanism, one fork. Receipt: `tests/scripts/estif_C15_gw_sector.py`.

**The four-lock ledger (RHAC-007).** ESTIF-Core forces a *point* where the GR family fits a *region* (0 free dark parameters vs 1–3 fitted): **Ω_k = 0** (kill-shot, exact at all epochs) · **γ ≈ 0.55 and slip = 1** (Front 3) · **w = −1** (Phase 2) · **c_gw = c** (C-15). **Honest finding:** none of the four separates ESTIF-Core from ΛCDM at Ω_k = 0 at linear order — all four separate it from GR's extra freedoms. The genuinely distinguishing content lives off the linear sheet (nonlinear halos — the N-body wall) and in the exactness-as-law structure: ΛCDM *fits* these values; ESTIF-Core *forbids* any others.

---

## 7. Dark matter analytics and the Ωm bootstrap *(status: analytics solid; bootstrap conditional; no Ωm novelty)*

**Analytics (unchanged, connected to the derived machinery):** Ωm = x₀ to 0.12% and Ωdm = x₀ − Ωb to 0.10% — both **consistency relations** (C2, final): r_universe is the ΛCDM particle horizon, which itself contains Ωm. σ/v_esc = 0.5 exact; λ_Jeans = 2.565r; the Tully–Fisher exponent M^(1/4) from the MOND limit. Halo internal structure (v_flat = 220 km/s needs δ ~ 10⁵) is behind the N-body wall.

**The bootstrap (v6.3.2, conditional on P).** Rather than remove the C2 circularity, solve it: adopt **principle P** (Ωm = R_H/r_p, the Hubble-to-particle-horizon ratio) and the circularity closes into **Ωm·I(Ωm) = 1**, which has a **unique root** — 0.3043 with zero measured inputs; **0.31408** with radiation (inputs {H₀, T_CMB, N_eff}): 0.96% from Planck, **0.53σ** inside its error bar. It back-predicts r_universe to −1.07%, so the import C2 objected to is no longer needed. Closure: a₀ → 1.1920×10⁻¹⁰ (0.66% from MOND — a magnitude improvement, *not* a derivation; a₀ is horizon-set, §4); DESI DR2 χ²/N = 1.618 (fixed-ruler). Identity: P ⇔ mean matter pull at the horizon = cH₀/2 — the same cH₀ that sets a₀ (ratio 1.00000). Input ledger after adopting P: measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ, x₀, r_universe, a₀}. Receipts: `tests/scripts/estif_omega_bootstrap.py` · `estif_bootstrap_closure.py`.

**Part B is closed (RHAC-009): P is NOT derivable as a law.** P's equality Ωm(a) = R_H/d_p holds only at a ≈ 1 — the two quantities fall monotonically and cross once, today (high-z: Ωm → 1 while R_H/d_p → ½). A law derivable from time-symmetric axioms must hold at every epoch; P does not. Escape routes closed: an attractor needs coupled dark energy, forbidden by the Phase-2 null; anthropics cannot reproduce 0.3141's precision. **P is a predictive postulate** — honest and falsifiable (Ωm = 0.3141, 0.53σ) — and the "Ωm derived" claim is **retired**. Receipt: `tests/scripts/estif_P_derivation_attempt.py`.

**Honest flags on the bootstrap:** everything is conditional on P; DESI 1.618 is fixed-ruler and could reorder under C4-style marginalization over (rd, H₀, Ωm); a₀'s empirical target carries ~10% scatter, so 0.66% is pleasing, not decisive.

**Prior art — the Gaztañaga verdict (13 July 2026 memo, absorbed; 🔴 no novelty on Ωm).** The causal-horizon-boundary → Ωm ≈ 0.3 / no-dark-energy result is Gaztañaga's, peer-reviewed 2019–2023 (his Ω_Λ ≈ 0.70 ⇒ Ωm ≈ 0.30; causal-universe scale ≈ 0.3176 H₀ via inflation, adjacent to the bootstrap root). **No priority is asserted on Ωm.** Required conduct: cite Gaztañaga prominently; lead with what is original to ESTIF — the a₀ link (the cH₀/2 identity tying the Ωm condition to the acceleration scale) and reaching the number by a different derivation (worthwhile, not novel). Still owed: the numerical cross-check of 0.31408 against his published value (checklist B-4b, item 4).

---

## 8. Expansion as projected inward flow *(status: mechanism as an equation; background inherited)*

The mechanism in one equation: H = c·d(ln b)/dw — expansion is motion through the bulk gradient (engine-verified Friedmann analog). One-flow unification: v² = 2GM/r + H²r², exact Schwarzschild–de Sitter — the gravity eddy and the Hubble flow in a single field. The background CMB is **inherited**: identical H(z) means recombination at z ≈ 1100, ~380,000 yr, and rd = 147.09 Mpc (Planck-calibrated import) pass automatically; the anisotropy *spectrum* as an ESTIF prediction awaits a Boltzmann-code pass (checklist A-5). Age of the universe: 13.379 Gyr.

---

## 9. Known limitations (kept in one place)

- **C1 — conditional deviations.** The ESTIF vacuum is exactly Schwarzschild, so any *deviation* from GR in photon-sphere shadows (EHT: 42.0 μas, 0.00σ) or vacuum GW propagation (LISA delay 491 μs, 49.2σ) must be sourced by the non-vacuum eddy background — a sector not yet derived. The observations are consistent; the deviation claims are conditional.
- **C2 — final.** Ωm = x₀ is a consistency relation; with RHAC-009 the downgrade is permanent (P non-derivable).
- **C4 — under-marginalization.** The χ²/N ≈ 0.66 CPL bar and the bootstrap's DESI 1.618 are fixed-(rd, H₀, Ωm) results; orderings could change under honest marginalization.
- **L1 — single-speed flow.** The flow ansatz forces p_r = −ρ, so it represents only vacuum and Λ-like sources — which is *why* no vacuum deviation is available (C1) and why any evolving-w must leave this class.
- **The coupling (C6)** is adopted, not derived (§3).
- **Preferred frame (barrier 5).** A2 introduces an absolute bulk time; simultaneity is reinterpreted, not derived — an explicit section is owed before submission (checklist C-11).

---

## 10. The honest open questions

1. **The a₀ keystone (RHAC-010):** the O(1) prefactor from horizon geometry, and the local mechanism that makes c·H act at galaxy radii (§4). The scale is derived; the number is not.
2. **x_c = 0.272:** geometric origin (black-hole-edge physics, between ISCO and the photon sphere), plus the double-duty flag (0.272 vs the data-preferred ≈ 0.221).
3. **Strong-field pressure/interior sector** (TOV-level star; full off-diagonal T_μν).
4. **ρ_eddy = x₀ρ_crit from the homogeneous field equation** (well-posed since Task 4; separate from the closed P question).
5. **Halo structure:** N-body simulation (budget wall; collaboration target).
6. **The Bullet Cluster** — the hardest external test for any no-dark-matter framework; unaddressed, and it gates the letter send (checklist B-8).
7. **What is the 4th dimension?** Geometrically well-defined (the direction the flow moves through); no identified counterpart in known field theories.

---

## 11. Where everything lives

Receipt paths cited above are current for the v6.4.2 layout: derivation and test scripts in `tests/scripts/` (Phase-2/A1′ suite in `tests/scripts/phase2_a1prime/`, with its own README of verdicts and expected outputs); cached data and `TEST_INDEX.md` in `tests/docs/`; the strict-A1 UKN pair in `archive/6-4-2 diagnostics/`. Decisions: `docs/RHAC.md` (RHAC-001…011 and Scenarios A–Q). Status and checklist: `docs/PROJECT.md`. History: `CHANGELOG.md`.

---

**Document version:** 6.4.2 · 21 July 2026 · consolidates ESTIF_CONCEPT 6.4.1, VALIDATION_REPORT (v6.3 line), SUMMARY_FOR_REVIEW 6.4.0, PHASE2_DECLARATION 1.0, GAZTANAGA_COMPARISON 1.0, JWST_TEST_SPEC (8 Jul), and the 13 Jul doctrine memo.
