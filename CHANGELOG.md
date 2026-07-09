# ESTIF Changelog

All notable changes to this project are documented in this file.

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
  under the same

Major milestone. The project forks into **Path One (ESTIF-Core, clean)** and
**Path Two (ESTIF-Extended, hard)**. See `MILESTONE_v6.3_THE_SPLIT.md`.

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
