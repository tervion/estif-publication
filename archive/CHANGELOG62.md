# ESTIF Changelog

All notable changes to this project are documented in this file.

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
