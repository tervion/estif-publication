# ESTIF Changelog

# ESTIF Changelog

All notable changes to this project are documented in this file.

---

## [6.4.1] — 2026-07-13 — Fronts 1–3, four-lock ledger, a₀-horizon doctrine, repo reorg

Housekeeping and one doctrine change. No axiom, derivation, or closed-result
status changes since 6.4.0 — RHAC-004 through RHAC-009 stand as closed. What's
new: three growth/structure fronts executed and filed, a fourth discriminator
lock recorded, the a₀ interpretation is reframed per a dedicated derivation
pass, and the repository layout is cleaned up.

### Added — Fronts 1–3 and the four-lock ledger (RHAC-007)
- Front 1 (growth under A1′): D+ clumping history computed; f(z=0.5)=0.7603
  matches the DESI RSD anchor; JWST verdict is an honest null (tension
  inherited from ΛCDM, neither relieved nor worsened).
- Front 2 (first black holes): fold-back rule ν·σ(M,0)·D(z)=δ_c established;
  star-channel and direct-collapse (no-star) channels characterized;
  primordial channel closed under the ζ=10⁻⁵ passport.
- Front 3 (second discriminator): growth index γ forced to ≈0.55 by the empty
  residual sector plus GR-equivalent D+ (computed 0.5455–0.5544 over z=0–5,
  pull −0.23σ from DESI); companion slip lock (Σ,η,μ)=(1,1,1) confirmed exact.
- **Four-lock ledger** separating ESTIF-Core from GR's extra freedoms: Ω_k=0
  (kill-shot), γ≈0.55 & slip=1 (Front 3), w=−1 (RHAC-004), c_gw=c (RHAC-008).
  None of the four separates ESTIF-Core from flat ΛCDM at linear order —
  distinguishing content lives off the linear sheet (nonlinear halos, N-body
  wall).
- Receipts: `tests/estif_front1_growth_sigma8_jwst.py`,
  `tests/estif_front2_first_hole_recipe.py`,
  `tests/estif_front3_second_discriminator.py`.

### Added — GW sector closed (RHAC-008, C-15)
- c_gw = c derived: under A1′+A2 the world is one Lorentzian geometry; TT
  waves and light share the same null cone for arbitrary flow. GW170817
  passes structurally (predicted deviation = 0 exactly). Receipt:
  `tests/estif_C15_gw_sector.py`.

  ⚠️ Both receipt filenames above are cited in RHAC-007/008 but were not
  located in the repository as of this writing. Locate and commit them, or
  amend the RHAC-007/008 lines to the actual filenames, before treating this
  entry as closed.

### Changed — a₀ reframed as a horizon quantity (13 July 2026 doctrine)
- **a₀ is a horizon-scale acceleration, a₀ ≈ c·H** (de Sitter surface
  gravity) — not a quantity derived from local matter density or from
  principle P. Mass-independence of a₀, flat rotation curves, the Baryonic
  Tully–Fisher relation, and the correct magnitude all follow from this
  alone; the numerical value (1.72% from MOND empirical) is unchanged.
- **Retired:** any language stating the √3 prefactor, or the full numerical
  a₀ coefficient, is *derived*. The prefactor is constrained to O(1) and
  known to be horizon-set, but the precise number is open.
- **Not touched:** the empirical fits, the black-hole exterior solution, or
  the Ωm bootstrap (RHAC Scenario Q) — this changes only the
  *interpretation* of a₀ and any claim that it is derived.
- Receipts: `tests/a0_horizon_test.py`, `tests/a0_prefactor_derivation.py`,
  `tests/estif_flow_sim.py`, `tests/estif_horizon.py`.
- Flagged open: x_c = 0.272 is used in two different roles across the corpus
  (the GR-crossover exponent value, and informally as a horizon-adjacent
  scale); the double duty is unresolved. See
  `docs/plan/ESTIF_document_update_guide.md` §6.
- Record: RHAC-010, `docs/plan/RHAC.md`.

### Changed — src/ suite relabeled to v6.4.1
- `estif_ec_gr_constants.py`, `estif_ec_gr_model.py`,
  `estif_ec_gr_run_simulation.py`: version strings and claim-status comments
  updated to v6.4.1. **No computational change** — all three files remain
  the frozen v6.2 analytical suite; the full run still passes 21/21 with
  identical numbers (x₀=0.3107, a₀=1.1793×10⁻¹⁰). Updated: Ω_tilt(z)
  cosmology flagged RETIRED (historical receipt only, per v6.3 Task 6);
  Ωm = x₀ relabeled a consistency relation (C2); a₀ commentary updated per
  the horizon doctrine above.

### Repository structure
- `PATH_ONE_CHECKLIST.md` moved from `archive/` to `docs/plan/` (the
  ROADMAP v6.4.0 banner already pointed here; the path did not exist until
  this move).
- `docs/LaTeX ` (trailing space in dirname) renamed to `docs/latex`;
  `ESTIF_arXiv_Paper.zip` unzipped in place.
- `tests/Computing growth history against DESI and JWST data/` (a
  doc-staging folder inside `tests/`) emptied and removed:
  `ESTIF_document_update_guide.md` moved to `docs/plan/`;
  `estif_field_dynamics.tex` and its preview PDF moved to `docs/latex/`;
  duplicate copies of `GAZTANAGA_COMPARISON.md` and `PHASE2_DECLARATION.md`
  (byte-identical to the `docs/plan/` originals) deleted.

### Still pending
- Script-side C2 (`estif_pathone_cosmology.py`), C4
  (`estif_pathone_aic_bic.py`), and C5 (`estif_task6_eddy_eos.py`) remain
  unapplied.

---

## [6.4.0] — 2026-07-12 — Phase 2 null, A1 → A1′, C-15 and Part B closed

The residual dark-energy question is settled (against a derived sector), the flow
axiom is amended to admit structure and radiation, and two open items (the GW wave
equation, and the derivability of principle P) are closed. **None of this reverses a
prior decision**; the derived gravity, the tilt retirement, the frozen-limit ΛCDM
parity, and the conditional bootstrap all stand. Terminology: the cosmological term is
now called an **imported cosmological constant**; the "frozen eddy" label is retracted.

### Added / resolved — Phase 2 (RHAC-004)
- **Phase 2 CLOSED — honorable null (strong form).** The *entire* residual stress the
  axioms permit was computed (pre-registered, submitted unmodified) and is **empty**:
  shape and rate FORBIDDEN (A1, A2), slosh ZERO (pure-divergence theorem), swirl
  FORBIDDEN free / NEGATIVE if sourced (⟨δρ⟩ = −⟨ω²⟩/32πG, slaved, ~10¹¹× below ρ_Λ).
  **w = −1 exactly; Λ is a bare imported constant.** ESTIF does not derive dark energy,
  and it is now *proven* it cannot from these axioms. Path Two (derive a thawing
  correction) is closed as a target. (`tests/test_UKN.py`; `docs/plan/PHASE2_DECLARATION.md`)
- **Retracted:** the "frozen eddy" interpretation of the constant Λ (author objected
  pre-derivation; the swirl/spin sector is precisely the forbidden one). The constant
  survives; the "eddy/spin" reading does not. Full rename pending `NAMING.md`.

### Changed — axiom A1 → A1′ (RHAC-005, RHAC-006)
- **Strict-A1 growth no-go theorem (RHAC-005):** under matter, A1–A3 force
  δ̈ = −(3/2)Ωm(a)Hδ with unique solution δ = H(a)/H₀ (decaying); f(0.5) = −0.91 vs
  +0.76 measured. Structure cannot form under strict A1; the strict flow sector also
  carries no free radiative modes (vs LIGO). Audited (5 attacks, all bounced).
  (`tests/test_UKN2.py`)
- **A1 → A1′ (RHAC-006):** slices are even *on average*; matter may source local dents;
  the spatial-curvature average is identically zero at all epochs — **law, and the
  registered falsifier** (current 0.0007 ± 0.0019). Restores (one purchase, three items):
  linear growth D₊ = H·∫da/(aH)³ = GR growth, **f(z=0.5) = 0.76** (DESI RSD-consistent);
  the radiative/GW sector; and a container for the ζ = 10⁻⁵ seed (value still imported).
  **All strict-A1 results survive as the exact-evenness limit** (vacuum Schwarzschild,
  the a₀ chain, the Friedmann background, the Phase 2 null — re-audited under A1′, holds).

### Added — gravitational-wave sector derived (RHAC-008, C-15 closed)
- **c_gw = c derived.** In the single-geometry flow, gravitational waves and light
  share the null cone for arbitrary flow; the GW propagation speed equals c with no
  free parameter, GW170817-consistent. C-15 moves from open to result.

### Changed — Part B closed; the "Ωm derived" claim retired (RHAC-009)
- **Principle P is not derivable as a law.** P (Ωm = R_H/r_p) holds only at our epoch:
  Ωm(a) = R_H/d_p crosses the required value once, near a ≈ 1, and diverges from it at
  high z (→ 1 vs ½), so it cannot follow from the time-symmetric axioms. Escape routes
  closed: a dynamical attractor breaks the Phase-2 null; anthropic selection cannot
  reproduce 0.31408. **P is a predictive postulate, not a theorem; Ωm's status is
  "fixed by P," and the "Ωm derived" language is retired.** The C2 consistency-relation
  downgrade is therefore final. Foundational effort redirected to x_c = 0.272 and the
  √3 / cH₀ scale — both of which underwrite a₀'s *value*, unlike P.

### Added — documents
- `docs/plan/PHASE2_DECLARATION.md` (the honorable-null record; referenced across the repo).
- `docs/plan/GAZTANAGA_COMPARISON.md` — verdict 🔴 **no Ωm novelty**: the causal-horizon →
  Ωm ≈ 0.3 / no-dark-energy result is Gaztañaga's (peer-reviewed, 2019–2023;
  headline Ω_Λ ≈ 0.70 ⇒ Ω_m ≈ 0.30). Cite him prominently; lead with the a₀ link.

---

## [6.3.2] — 2026-07-09 — Ωm bootstrap (Door 2, conditional)

> **Superseded status note (added retroactively, 6.4.0):** the "Part B open" flag
> below was accurate at the time. RHAC-009 (12 July 2026) closed Part B: P is not
> derivable as a law. This entry is preserved as written, dated to when it was
> true.

### Added

- **Ωm bootstrap (Door 2).** Principle P (Ωm = R_H/r_p) inverts the C2 circularity
  into the closed equation **Ωm·I(Ωm) = 1**, with a **unique** root: **0.3043** from
  zero measured inputs; **0.31408** with radiation (inputs = T_CMB, N_eff, h) —
  **0.96% from Planck, 0.53σ** inside its error bar. Back-predicts
  r_universe = 4.353×10²⁶ m (−1.07% vs the import C2 flagged). Identity: P ⇔ mean
  matter pull at the horizon = cH₀/2 — the same cH₀ that sets a₀ (ratio 1.00000).
  (`estif_omega_bootstrap.py`)
- **Closure.** Propagating the bootstrap Ωm: a₀ = 1.1920×10⁻¹⁰ m/s², MOND agreement
  **1.72% → 0.66%** (SPARC insensitive, v_flat ×1.00269); DESI DR2 χ²/N = **1.618**
  vs ΛCDM 1.919, *within the fixed-(H₀, rd) test*. **Input ledger after adopting P:**
  measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ, x₀, r_universe, a₀}.
  (`estif_bootstrap_closure.py`)
- RHAC **Scenario Q**; checklist items **B-4b** (Gaztañaga memo) and **B-4c** (derive P).

### Honest flags

- **Conditional on P**, which is *not* derived from A1–A3. Part B open; three
  candidate routes, none attempted. Without it this reparametrizes C2's circularity
  rather than escaping it.
- **Gaztañaga adjacency** (causal-universe scale ≈ 0.3176 H₀, via inflation).
  Comparison memo is a prerequisite for any novelty claim; it does not exist yet.
- **DESI 1.618 is fixed-ruler.** Could reorder under C4 marginalization.
- **a₀'s empirical target carries ~10% scatter.** 0.66% is pleasing, not decisive.

### Still pending

- Script-side C2 (`estif_pathone_cosmology.py`), C4 (`estif_pathone_aic_bic.py`),
  and C5 (`estif_task6_eddy_eos.py`) remain unapplied.

---

## [6.3.1] — 2026-07-09 — Adversarial errata

Six corrections from an adversarial review pass of the v6.3 documents. **None
reverses a v6.3 decision** — the split, the tilt retirement, the de-circularization,
the frozen-eddy reframe, and the derived constraint sector all stand. C1, C2, and C6
are pre-submission critical. See `CORRECTIONS_v6.3.1.md`.

### Changed — epistemic status of three claims

- **C1 — strong-field deviations are conditional.** The derived ESTIF vacuum is
  *exactly* Schwarzschild (full Einstein tensor = 0), so it produces zero deviation
  from GR in photon-sphere shadows or in vacuum GW propagation. The EHT M87\* and
  LISA figures are therefore demoted from ✅ predictions to ⚠️ predictions
  **conditional on the un-derived eddy-stress sector**. The observations remain
  *consistent* with ESTIF; the *deviation from GR* is what awaits derivation. The
  weak-field a₀/MOND chain is unaffected. Planck Λ is a calibration match, not a
  vacuum deviation, and is unaffected.
- **C2 — Ωm = x₀ is a consistency relation, not a prediction.** Here
  r_universe = 4.4×10²⁶ m is the ΛCDM particle horizon, an integral that itself
  contains Ωm. The 0.12% agreement is a self-consistency of the geometric picture
  with Planck-calibrated values, not an Ωm-independent derivation. An independent
  prediction requires deriving r_universe from the flow framework (RHAC Scenario H).
  The numerical agreement stands; the word "predicted" is withdrawn.
- **C6 — "derived, not borrowed" made precise.** The three flow axioms uniquely
  *select* the constraint (energy) sector of General Relativity in
  Painlevé–Gullstrand gauge, forcing mass continuity and hence exact Schwarzschild in
  vacuum — without matching to the Schwarzschild solution, which was the actual gap.
  The gravitational *coupling* (geometric constraint scalar ↔ 8πG × energy density)
  is **adopted** from Einstein–Hilbert, not derived from below. ESTIF does not derive
  G or the factor 8π.

### Added — a hard constraint and a caveat

- **C3 — σ8/S8 hard filter on the JWST test.** A growth enhancement that persists to
  z = 0 at the 8 Mpc/h scale is grossly excluded: a scale-independent +13% boost
  gives σ8 ≈ 0.917, ~18σ above Planck and in the *wrong direction* for weak-lensing
  S8. D_ESTIF(z) must therefore be **two-sided** — ≳13% enhancement at z ≈ 9,
  decaying to ≈ 1 by z ≲ 2 — or scale-dependent. This makes the derivation target
  more specific, not merely larger. (`docs/report/JWST_TEST_SPEC.md` §5a.)
- **C4 — Path Two target is under-marginalized.** The χ²/N ≈ 0.66 CPL bar comes from
  a BAO-alone fit with rd, H₀, and Ωm held fixed at Planck values. Fixing nuisance
  parameters inflates the evolving-dark-energy advantage; DESI's own BAO-alone
  preference is milder. Re-derive with those parameters marginalized before using
  0.66 to justify Path Two, and repeat the frozen-eddy vs fitted-CPL AIC comparison
  under the same marginalization.

### Fixed

- **C5** — ⬜ *not yet applied.* The Task 6 E1 (conserved-angular-momentum spin)
  explanatory block in `tests/estif_task6_eddy_eos.py` still contains a mid-thought
  fragment and an exponent placeholder. Replacement text is specified in
  `CORRECTIONS_v6.3.1.md`. The numerical result (χ²/N = 3232, falsified) is
  unaffected either way.
- **Documentation writing tasks closed.** The v6.3 concept rewrite is confirmed
  complete: A1–A3 are stated explicitly in `ESTIF_CONCEPT.md`, the shrinking-ruler
  narrative is retired in favour of the flat-slice flow picture, and v_flow = c·x₀
  is relabelled as the *sideways component* of a total-c motion. Four documents that
  still listed this as outstanding work were corrected; `CITATION.cff`'s abstract no
  longer states that A2 and A3 live only in the derivation scripts.

### Repository structure

- v6.3 document suite placed at canonical paths (`MILESTONE_v6.3_THE_SPLIT.md` and
  `CORRECTIONS_v6.3.1.md` to root; `STATUS.md`, `VALIDATION_REPORT.md`, and
  `JWST_TEST_SPEC.md` to `docs/report/`; `PATH_ONE_CHECKLIST.md` to `docs/plan/`).
- `CHANGELOG.md`, `docs/plan/RHAC.md`, and `docs/guide/ROADMAP.md` reconstituted from
  their v6.2 originals plus the v6.3 append/prepend blocks.
- `ROADMAP.md`'s superseded v6.2 phases (5.3, 5.4, 6.1–6.3, 8.2 — all built on the
  retired Ω_tilt law) marked as preserved history, not a working plan.
- `MILESTONE_v6.3_THE_SPLIT.md` kept frozen as the split-day snapshot; its three
  corrected claims carry inline errata notes rather than rewrites.

---

## [6.3.0] — 2026-07-08 — "The Split"

Major milestone. The project forks into **Path One (ESTIF-Core, clean)** and
**Path Two (ESTIF-Extended, hard)**. See `MILESTONE_v6.3_THE_SPLIT.md`.

### Added — gravity field equation is now DERIVED (Task 4)
- The field equation `rho_eff = m′(r)/(4πr²)` (mass continuity = Poisson in
  integrated form) is now **forced** by the flow axioms via the Gauss–Codazzi
  engine — no longer matched to the Schwarzschild solution.
  - Vacuum (rho_eff = 0) → v² = 2GM/r uniquely → **exact Schwarzschild**.
  - Uniform-density ball → rho_eff = rho0 **exactly** (correct Newtonian source).
  - The former "D2 Poisson postulate" is now a theorem for vacuum, Newton, and
    Schwarzschild.
- New scripts: `estif_task4_field_equation.py` (5/5),
  `estif_flow_signature_dynamics.py` (18/18: signature + SR + Newton from a
  Euclidean bulk with a universal speed-c constraint),
  `estif_converse_flow_law.py` (Birkhoff's theorem restated in flow variables),
  `estif_tmunu_gauss_codazzi.py` (ADM engine validated vs flat FRW and de Sitter).

### Changed — cosmology reframed; tilt apparatus retired (Tasks 5, 5b, 6)
- **Circularity fixed:** `x(z)` no longer uses ΛCDM as its own ruler. The
  self-consistent fixed-point solve (`estif_task5_desi_selfconsistent.py`) drops
  the DESI DR2 fit from χ²/N = 10.80 (circular) to **3.35** (self-consistent).
- **Frozen-eddy reframe:** the constant-eddy limit (w = −1, derived from Task 4)
  scores **χ²/N = 1.92, tying ΛCDM and beating the tilt formula's 3.35**. On DESI
  the entire Ω_tilt apparatus (N_MAX, B, sign-flip, z<2 cutoff) is a net negative
  and is **retired** from the cosmology claim under Path One.
  (`estif_task6_eddy_eos.py`.)
- **DESI-preferred w(z):** the self-consistent tilt already thaws toward the
  DESI-preferred curve (within ~0.05); best-fit evolving-w flow reaches χ²/N =
  0.66 (Path Two target). (`estif_task5b_cosmo_eos.py`.)

### Deprecated / retired
- `Ω_tilt(z)` as a cosmology claim (moved to an "explored and set aside"
  appendix under Path One). The tilt formula's local strong-field use (EHT, the
  n(x) dynamic exponent) is unaffected; only the *cosmological* dark-energy claim
  is retired.
- The z < 2 hard cutoff (scaffolding for the retired tilt term) is moot under
  both paths.

### Fixed / flagged (fidelity audit)
- `estif_fidelity_audit.py` found that axioms A2 (universal speed c) and A3
  (empty space is not a source) are present only in the derivation scripts, and
  that `ESTIF_CONCEPT.md` runs a competing "shrinking-ruler" narrative (A1
  conflict). Path One writing tasks: add A2 + A3 to the theory, retire the
  shrinking-ruler narrative, relabel v_flow = cx₀ ≈ 0.31c as the *sideways
  component* of a total-c flow (not the full flow speed).

### Portability
- The DESI test scripts now cache data beside the script
  (`os.path.dirname(__file__)`) instead of a hard-coded absolute path, so they
  run unmodified on any machine.

### Notes
- The gravity letter is unaffected and **strengthened**: the a₀/MOND derivation
now rests on a derived field equation rather than a Schwarzschild match.

---

## [6.2.0] - 2026-03-20 — Pre-publication blockers resolved, letter drafted

**Summary:** Both pre-publication blockers identified by peer review are resolved.
The gravity-only letter is ready for submission. No changes to core physics or test results.

### Added

#### a₀ Redshift Constancy — Algebraic Proof (resolves blocker 1)

The question: if a₀ = H₀cx₀/√3 uses today's H₀, does ESTIF predict a₀ ∝ H(z)?

Answer: No. In the comoving frame (the physically correct frame for galaxy dynamics):
- x(z) = c / [H(z) × r_universe_comoving]
- a₀(z) = H(z) × c × x(z) / √3 = c² / (r_universe_comoving × √3)
- H(z) cancels exactly. Maximum deviation: 2.22×10⁻¹⁶ (floating-point epsilon)
- This is an algebraic identity, not a numerical result

Observational confirmation: Di Teodoro+2021, Übler+2017, Tiley+2019 all consistent
with constant a₀ at z = 0.75–2.2, deviations ≤ 2σ.

Added: `tests/test_a0_redshift.py`  
Added: `results/validated/a0_redshift.png`

#### Parameter Independence — Test 3 (resolves the "just Planck values" criticism)

Tested a₀ = H₀cx₀/√3 across 3,600 combinations of H₀ ∈ [65,75] km/s/Mpc
and Ωm ∈ [0.27, 0.33]:
- 100% of combinations within ±20% SPARC scatter
- 8 published datasets (Planck, WMAP, SH0ES, DES, KiDS, SPT, ACT, H0LiCOW): all pass
- Planck–SH0ES Hubble tension shifts a₀ by only 4.1%

Added: `tests/test_a0_parameter_independence.py`  
Added: `results/validated/a0_parameter_independence.png`

#### 1/√3 Language Reframe (resolves blocker 2)

Changed from "derived from the equipartition theorem" to:
"motivated by 3D spatial isotropy and consistent with the equipartition theorem;
a complete kinetic theory of the eddy background is identified as future theoretical work"

Updated in: `tests/derive_mond_from_geometry.py`, `docs/SUMMARY_FOR_REVIEW.md`,
all letter drafts.

#### Gravity-Only Letter — Draft Complete

`ESTIF_letter_final.docx` — 4–6 page letter ready for submission.
Sections: Abstract, Introduction, Framework, Derivation, Redshift Constancy,
SPARC Validation, Discussion, Conclusions.
Target venues: MNRAS Letters, ApJL, JCAP.

Frame argument added to Section 4: explains why the comoving frame is the
physically correct frame for galaxy dynamics, and why the H(z) cancellation
is not a coordinate trick but a consequence of bound systems decoupling
from expansion.

### Changed

- `src/estif_ec_gr_constants.py`: v6.1 → v6.2
- `src/estif_ec_gr_model.py`: v6.0 → v6.2
- `src/estif_ec_gr_run_simulation.py`: v6.0 → v6.2 (3 strings); main() function added
- `setup.py`: version 2.0.0 → 6.2.0
- `requirements.txt`: updated header
- All documentation files: version and date updated

### Publication Status After v6.2

| Blocker | Status |
|---|---|
| a₀ redshift evolution | ✅ Proved constant — algebraic identity |
| 1/√3 language | ✅ Reframed — motivated by isotropy, formal derivation future work |

---

## [6.1.0] - 2026-03-18 — MOND Derived, SPARC Validated, DESI Constraint, Multipliers

**MILESTONE VERSION**: The MOND critical acceleration a₀ is derived from geometry for the
first time, confirmed against 87 SPARC galaxies. The cosmological sector is tested against
DESI DR2 and found to fail. The 1/3 multiplier for B is derived from first principles.
A clear two-sector picture emerges: gravity sector solid, cosmology sector requires rework.

### Added

#### MOND Derivation — Zero Free Parameters

- Four-step geometric derivation of a₀ = H₀cx₀/√3 = 1.179×10⁻¹⁰ m/s²
- Zero free parameters: H₀ and x₀ = Ωm from Planck 2018; √3 from isotropy
- Agreement: 1.72% from MOND empirical (within MOND measurement uncertainty ~5–10%)
- Uniqueness confirmed: 1/√3 is the only one of 12 candidate factors landing below 5%
  AND having an independent physical derivation
- Added: `tests/derive_mond_from_geometry.py`

#### SPARC Baryonic Tully-Fisher Validation

- Tested against full SPARC catalog (Lelli et al. 2016, AJ 152, 157)
- Quality-1 sample (87 galaxies): RMS = 15.6% — within observed BTFR scatter
- Quality-1+2 sample (129 galaxies): RMS = 18.4%
- Added: `tests/test_sparc_tully_fisher.py`

#### SPARC Bias Analysis

- Mean bias −7.6% traced to stellar mass calibration, not force law structure
- All correlations disappear when Υ* corrected (all p > 0.16)
- Added: `tests/test_sparc_bias_analysis.py`

#### Multiplier Derivation

- B = L/3: DERIVED from 3D isotropic projection (0.69% off)
- N_MAX/5/7: CONDITIONAL on x_c geometric derivation
- Added: `tests/test_multiplier_derivation.py`

#### DESI DR2 Test

- ESTIF Ω_tilt(z) fails DESI DR2 at chi²/N = 10.8 (ΛCDM: 1.9)
- Pre-existing prediction w_eff ≈ −1.358 is 3.5σ from DESI DR2 w₀ = −0.73 ± 0.10
- Root cause: x(z) formula is circular (uses H_ΛCDM as its own ruler)
- Added: `tests/test_desi_wz_consistency.py`

#### Cross-Examination Synthesis

- Added: `tests/cross_examination.py`

### Changed

All documentation files updated to v6.1.

---

*(Earlier versions: v6.0 March 2026, v4.0 March 2026, v2.0 October 2025, v1.0 September 2024 — see git log)*
