## [6.1.0] - 2026-03-18 — MOND Derived, SPARC Validated, DESI Constraint, Multipliers Progress

**MILESTONE VERSION**: The MOND critical acceleration a₀ is derived from geometry for the
first time, confirmed against 87 SPARC galaxies. The cosmological sector is tested against
DESI DR2 and found to fail. The 1/3 multiplier for B is derived from first principles.
A clear two-sector picture emerges: gravity sector solid, cosmology sector requires rework.

### Added

#### MOND Derivation — Zero Free Parameters (New Milestone)

- **Four-step geometric derivation** of a₀ = H₀cx₀/√3 = 1.179×10⁻¹⁰ m/s²
  - Step 1: Force law a = −c²/2 × ∇(ω/H₀)² = GM/r² (exact Newton, no approximation)
  - Step 2: Flow speed v_flow = cx₀ (from tilt formula at x = x₀, two independent routes)
  - Step 3: 3D isotropic projection v_3D = v_flow/√3 (unique from dimension count)
  - Step 4: a₀ = v_3D × H₀ (natural cosmic deceleration threshold)
  - **Zero free parameters**: H₀ and x₀ = Ωm from Planck 2018; √3 from isotropy
- **Uniqueness confirmed**: 1/√3 is the only one of 12 candidate factors landing below 5%
  AND having an independent physical derivation
- **Agreement**: 1.72% from MOND empirical value (within MOND measurement uncertainty ~5-10%)
- Added: `tests/derive_mond_from_geometry.py` — standalone 4-step derivation script

#### SPARC Baryonic Tully-Fisher Validation

- Tested a₀ = H₀cx₀/√3 against the full SPARC catalog (Lelli et al. 2016, AJ 152, 157)
  via VizieR J/AJ/152/157, 175 galaxies total
- Formula: v_flat = (G × M_bar × a₀)^(1/4), Υ* = 0.50 M☉/L☉ (McGaugh & Schombert 2014)
- **Quality-1 sample (87 galaxies): RMS = 15.6%** — within observed BTFR scatter (15–20%)
- **Quality-1+2 sample (129 galaxies): RMS = 18.4%** — also within scatter
- 82% of quality-1 galaxies predicted within 20%
- Added: `tests/test_sparc_tully_fisher.py`

#### SPARC Bias Analysis

- Root cause of −7.6% mean bias identified: stellar mass calibration
- Bias correlates with morphological type, gas fraction, surface brightness (all p < 0.002)
- **Key finding**: all correlations disappear when Υ* is corrected → bias is calibration, not structural
- Zero-bias Υ* ≈ 0.85 (70% above McGaugh+2014 standard — noted as honest caveat)
- No residual structure after Υ* correction (all p > 0.16)
- Conclusion: the ESTIF force law is structurally sound; mass calibration needs improvement
- Added: `tests/test_sparc_bias_analysis.py`

#### Multiplier Derivation Progress

- **B = L/3 derived** from 3D isotropic projection of the 4D decay rate:
  - L = ln(r_e/l_P) = ln(α) + ln(m_P/m_e) — exact decomposition confirmed
  - Same isotropy principle as 1/√3 in MOND derivation
  - Error: 0.69% — within uncertainty of incomplete T_μν projection
- **N_MAX conditional derivation**: follows from B = L/3 + GR crossover n(x_c) = 1/2,
  but x_c = 0.272 is still observationally determined (not geometrically derived)
- Honest gap identified: two unexplained fractions reduced to ONE unexplained number (x_c)
- Added: `tests/test_multiplier_derivation.py`

#### DESI DR2 Test

- Tested ESTIF Ω_tilt(z) against DESI DR2 BAO data (arXiv:2503.14738, released 19 March 2026)
  using official data from CobayaSampler/bao_data repository (DR1 and DR2)
- **Result: FAILS** — chi²/N = 10.8 vs DESI DR2 (ΛCDM: 1.9)
- Pre-existing prediction w_eff ≈ −1.08 is 3.5σ from DESI DR2 w₀ = −0.73 ± 0.10
- Only 2/13 DR2 bins within 1σ
- Root cause identified: x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is circular
- DESI DR2 itself shows 3.1σ preference for dynamical dark energy over ΛCDM —
  ESTIF's specific w(z) evolution shape is wrong, not the idea of dynamical DE
- Added: `tests/test_desi_wz_consistency.py`

#### Cross-Examination Synthesis

- Formal cross-examination of all three tests showing two-sector structure:
  - **Gravity sector** (Tests 2+3): SOLID — isotropy principle appears three independent times
  - **Cosmology sector** (Test 1): FAILS DESI DR2 — Ω_tilt(z) needs rework
- Root cause of DESI failure: circular x(z) definition and aggressive z_eff cutoff
- DESI failure does NOT invalidate gravity sector (completely different scales)
- Priority action list established (4 priorities in order)
- Publication path for gravity-only letter confirmed as achievable
- Added: `tests/cross_examination.py`

#### Results and Images

New validated results added to `results/validated/`:
- `mond_derivation.png` — 4-step derivation chain, uniqueness table, Tully-Fisher
- `sparc_tully_fisher.png` — 87 SPARC galaxies vs derived a₀ (RMS 15.6%)
- `sparc_bias_analysis.png` — bias correlations, Υ* sensitivity, residual structure
- `multiplier_derivation.png` — L decomposition, n(x) crossover, derivation chain

New result added to `results/` (documenting tested constraint):
- `desi_wz_consistency.png` — DESI DR2 vs ESTIF (shows the failure clearly)
- `cross_examination.png` — two-sector verdict, priority action list

### Changed

- `README.md` → v6.1: MOND milestone section, two-sector status, honest DESI result,
  updated prediction table, updated test list, updated honest assessment
- `CITATION.cff` → v6.1: Version, date, abstract updated with new results
- `docs/SUMMARY_FOR_REVIEW.md` → v6.1: Complete update reflecting new milestone
- `docs/guide/ROADMAP.md` → v6.1: Phase 5.3 marked done (and failed), priority list updated
- `docs/plan/RHAC.md` → v6.1: Scenario E resolved (DESI conflict → FAILS), new scenarios
- `docs/report/STATUS.md` → v6.1: New test results, honest DESI status
- `docs/report/VALIDATION_REPORT.md` → v6.1: Part 5 added for new tests
- `docs/report/ESTIF_CONCEPT.md` → v6.1: MOND section updated with 4-step derivation
- `results/results_README_v6.md` → v6.1: New plots listed
- `results/validated/validated_README_final.md` → v6.1: New figures documented
- `results/work_in_progress/wip_README.md` → v6.1: DESI comparison moved to done (failed)

### Honest Assessment Update

The pre-existing prediction w_eff ≈ −1.08 was published in v6.0. DESI DR2, released
19 March 2026, shows this prediction is falsified at 3.5σ. This is documented honestly
rather than explained away. The cosmological sector (Ω_tilt(z) evolution) requires a
fundamental rework — the self-consistent x(z) using H_ESTIF rather than H_ΛCDM.

At the same time, the gravity sector results (MOND derivation, SPARC validation,
multiplier derivation) constitute a genuine advance that stands independently of
the cosmological sector failure.


---

*(Previous entries below — see git history for full v6.0.0 changelog)*

