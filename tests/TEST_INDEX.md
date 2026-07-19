# TEST INDEX — provenance and status of every file in `tests/`

**Version:** 6.4.1 · **Created:** 15 July 2026
**Maintenance rule (RHAC discipline):** when a new test is added, add a row here
in the same commit. When a claim is retired, the script is NOT deleted — its
Status changes to HISTORICAL RECEIPT. Never delete receipts.

**Status legend**
- ✅ CURRENT RECEIPT — backs a live claim; cited in current docs
- 🔴 HISTORICAL RECEIPT — the claim it tested is retired/superseded; kept as the
  runnable record (per RHAC "mark, never delete")
- ⚠️ SUPERSEDED — replaced by a later, more correct test of the same question
- 📦 CACHE / DATA — auto-generated or downloaded data; must stay beside its script
- 🔧 RUNNER / TOOL — infrastructure, not a claim receipt
- ❓ UNVERIFIED — version/role inferred; confirm against `git log` date and the
  script's own header docstring, then remove the ❓

Dates below are from the changelog/RHAC record. To cross-check against actual
git history, run the command in the "Verify against git" section at the bottom.

---

## v6.4.1 — 13 July 2026 · a₀-horizon doctrine (RHAC-010)

| File | Role | Status |
|---|---|---|
| `a0_horizon_test.py` | Algebraic equivalence: √Λ-form vs Hubble-form of a₀ (guide §2) | ✅ CURRENT RECEIPT |
| `a0_prefactor_derivation.py` | Prefactor target band; numerology floor; river-model edge 0.26–0.28 (guide §3) | ✅ CURRENT RECEIPT |
| `estif_flow_sim.py` | Local flow cannot produce a₀ (mass-dependence M^(1/3), non-locality; guide §4) | ✅ CURRENT RECEIPT |
| `estif_horizon.py` | Horizon background reproduces flat curves, BTFR, mass-independence, magnitude (guide §5) | ✅ CURRENT RECEIPT |
| `a0_tension_corrected.py` | a₀ tension recomputation from the same session | ✅ ❓ confirm date/role from header |
| `btfr_lensing.py` | Weak-lensing BTFR check (Mistele et al. 2024) | ✅ ❓ confirm date/role from header |

## v6.4.0 — 11–12 July 2026 · Phase 2 null, A1→A1′, Fronts 1–3, C-15, Part B closed

| File | Role | Status |
|---|---|---|
| `test_UKN.py` | Phase 2 residual-stress census under strict A1 — sector EMPTY (RHAC-004) | ✅ CURRENT RECEIPT |
| `test_UKN2.py` | Growth no-go theorem + A1′ re-audit; growth & GW restored (RHAC-005/006) | ✅ CURRENT RECEIPT |
| `phase2_a1prime/README.md` | Index of the A1′ fork receipts | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_a1prime_deepen_exact.py` | A1′: exact dent-deepening solution | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_a1prime_growth_restored.py` | A1′: linear growth D₊ restored, f(0.5)=0.76 | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_growth_nogo_law.py` | Strict-A1 no-go law derivation (RHAC-005) | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_growth_nogo_audit.py` | Five audit attacks on the no-go (all bounced) | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_p2_door1_rate_dial.py` | Phase 2 door 1: rate sector (FORBIDDEN) | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_p2_door2_slosh_divergence.py` | Phase 2 door 2: slosh sector (ZERO) | ✅ CURRENT RECEIPT |
| `phase2_a1prime/estif_p2_door3_swirl_ledger.py` | Phase 2 door 3: swirl sector (forbidden/negative) | ✅ CURRENT RECEIPT |
| `estif_front1_growth_sigma8_jwst.py` | Front 1: D(z), f(z), σ8, fσ8, JWST honest null (RHAC-007) | ✅ CURRENT RECEIPT |
| `estif_front2_first_hole_recipe.py` | Front 2: fold-back rule, first-black-hole channels (RHAC-007) | ✅ CURRENT RECEIPT |
| `estif_front3_second_discriminator.py` | Front 3: growth index γ≈0.55, slip lock (RHAC-007) | ❌ **CITED, NOT IN REPO** — locate & commit, or amend RHAC-007 |
| `estif_C15_gw_sector.py` | C-15: c_gw = c derived, GW170817 structural pass (RHAC-008) | ❌ **CITED, NOT IN REPO** — locate & commit, or amend RHAC-008 |
| `estif_P_derivation_attempt.py` | Principle P non-derivability obstruction (RHAC-009) | ✅ CURRENT RECEIPT |
| `estif_pathtwo_target_lock.py` | Pre-registered Path Two target, committed before derivation | ✅ CURRENT RECEIPT (pre-registration) |
| `pathtwo_target_lock_output.txt` | Timestamped output of the target lock | 📦 CACHE / DATA |

## ~10 July 2026 · Two-machine crossroads (RHAC-001)

| File | Role | Status |
|---|---|---|
| `mu_extraction.py` | Tilt sector read locally → g/g_N → 0; μ(a/a₀) confirmed underived | ✅ CURRENT RECEIPT (documents the open μ gap) |
| `ripple_speed.py` | GW/ripple propagation-speed exploration (RHAC-008 prep) | ❓ confirm date/role from header |

## v6.3.0–v6.3.2 — 8–9 July 2026 · The Split, errata, bootstrap

| File | Role | Status |
|---|---|---|
| `estif_task4_field_equation.py` | Field equation derived: mass continuity = Poisson (5/5) | ✅ CURRENT RECEIPT |
| `estif_flow_signature_dynamics.py` | Signature + SR + Newton from Euclidean bulk + universal c (18/18) | ✅ CURRENT RECEIPT |
| `estif_converse_flow_law.py` | Vacuum forces v²=2A/r — Birkhoff in flow variables | ✅ CURRENT RECEIPT |
| `estif_converse_flow_law2.py` | Variant/extension of the converse flow law | ✅ ❓ confirm relation to v1 from header |
| `estif_tmunu_gauss_codazzi.py` | ADM/Gauss–Codazzi engine, validated vs FRW + de Sitter | ✅ CURRENT RECEIPT |
| `estif_tmunu_task4.py` | Task-4 stress tensor; p_r = −ρ limitation L1 | ✅ CURRENT RECEIPT |
| `estif_task5_desi_selfconsistent.py` | De-circularized DESI DR2 fit: 10.8 → 3.35 | ✅ CURRENT RECEIPT |
| `estif_task5b_cosmo_eos.py` | DESI-preferred w(z); tilt tracks it within ~0.05 | ✅ CURRENT RECEIPT |
| `estif_task6_eddy_eos.py` | E1/E2 eddy EoS falsified; constant-Λ limit ties ΛCDM (1.92) | ✅ CURRENT RECEIPT — ⚠️ C5 text fix still pending in this file |
| `estif_fidelity_audit.py` | Axiom-presence audit of the corpus (71 files) | ✅ CURRENT RECEIPT |
| `estif_omega_bootstrap.py` | Ωm·I(Ωm)=1 unique root 0.31408 (conditional on P; v6.3.2) | ✅ CURRENT RECEIPT |
| `estif_bootstrap_closure.py` | Bootstrap Ωm propagated: a₀ 0.66%, r_u, DESI 1.618 (v6.3.2) | ✅ CURRENT RECEIPT |
| `estif_pathone_cosmology.py` | Path One DESI DR2 cosmology fit | ✅ CURRENT RECEIPT — ⚠️ C2 script-side fix pending |
| `estif_pathone_aic_bic.py` | AIC/BIC comparison, frozen limit vs fitted CPL | ✅ CURRENT RECEIPT — ⚠️ C4 marginalization pending |
| `estif_jwst_growth_spec.py` | JWST growth-test specification (C3 two-sided filter) | ✅ ❓ confirm date from header |
| `dr1_cov.txt` `dr1_mean.txt` `dr2_cov.txt` `dr2_mean.txt` `dr2b_c.txt` `dr2b_m.txt` | DESI DR1/DR2 data caches (fetched beside script by design) | 📦 CACHE — do not move |
| `pathone_dr1_c.txt` `pathone_dr1_m.txt` `pathone_dr2_c.txt` `pathone_dr2_m.txt` | Path One fit caches | 📦 CACHE — do not move |

## v6.1–v6.2 — 18–20 March 2026 · MOND derived, SPARC, blockers resolved

| File | Role | Status |
|---|---|---|
| `derive_mond_from_geometry.py` | Four-step a₀ derivation — **note RHAC-010:** √3 prefactor now a working form, not derived; numbers stand | ✅ CURRENT RECEIPT (reinterpreted) |
| `test_sparc_tully_fisher.py` | 87 quality-1 SPARC galaxies, RMS 15.6% | ✅ CURRENT RECEIPT |
| `test_sparc_bias_analysis.py` | −7.6% bias traced to Υ* calibration | ✅ CURRENT RECEIPT |
| `test_multiplier_derivation.py` | B = L/3 derived; N_MAX = 5/7·L conditional on x_c | ✅ CURRENT RECEIPT |
| `test_mond_sqrt3.py` | 1/√3 uniqueness among 12 candidate factors | ✅ CURRENT RECEIPT (prefactor status per RHAC-010) |
| `test_a0_redshift.py` | a₀ redshift constancy — H(z) cancels (blocker 1) | ✅ CURRENT RECEIPT |
| `test_a0_parameter_independence.py` | 3,600 H₀/Ωm combinations within SPARC scatter (blocker resolution) | ✅ CURRENT RECEIPT |
| `test_tully_fisher_correction.py` | Tilt adds no M-factor; MOND limit gives M^(1/4) (Scenario I) | ✅ CURRENT RECEIPT |
| `test_desi_wz_consistency.py` | Ω_tilt(z) vs DESI DR2: χ²/N = 10.8 FAIL (Scenario E) | ⚠️ SUPERSEDED by `estif_task5_desi_selfconsistent.py`; kept as the failure record |
| `cross_examination.py` | v6.1 cross-examination synthesis | 🔴 HISTORICAL RECEIPT |

## v4.0–v6.0 — March 2026 and earlier · calibration, tilt era, eddy dark matter

**Gravity-sector receipts (claims still live):**

| File | Role | Status |
|---|---|---|
| `test_joint_calibration.py` | EHT + Λ + LISA simultaneous calibration (N_MAX, B) | ✅ CURRENT RECEIPT (C1: deviations conditional) |
| `test_gravity_time_connection.py` | β = τ at n = ½, x = 0.272 | ✅ CURRENT RECEIPT |
| `test_electron_connection.py` | N_MAX, B ↔ ln(r_e/l_P) | ✅ CURRENT RECEIPT |
| `test_eddy_time_gravity.py` | Gravity = time = eddies three-way identity | ✅ CURRENT RECEIPT |
| `test_solar_system_eddy.py` | Formula dormant at solar-system scales (GR recovered) | ✅ CURRENT RECEIPT |
| `test_eddy_dark_matter.py` | ρ_eddy = x₀ρ_crit framework (Ωm = x₀ now C2 consistency relation) | ✅ CURRENT RECEIPT (reinterpreted) |
| `test_collisionless_eddy.py` | σ/v_esc = 0.5 exact; λ_Jeans = 2.565r | ✅ CURRENT RECEIPT |
| `test_virialized_eddy.py` | Virialization of the eddy background | ✅ CURRENT RECEIPT |
| `test_jeans_length_eddy.py` | Self-similar Jeans hierarchy | ✅ CURRENT RECEIPT |
| `test_hierarchical_collapse.py` | Hierarchical collapse / free-fall epoch | ✅ CURRENT RECEIPT |
| `test_x0_matter_term.py` | x₀ as matter-like term | ✅ CURRENT RECEIPT (C2 applies) |

**Tilt-cosmology receipts (claim RETIRED v6.3 — kept per RHAC discipline):**

| File | Role | Status |
|---|---|---|
| `test_estif_cosmology.py` | Option A Ω_tilt(z) cosmology suite | 🔴 HISTORICAL RECEIPT |
| `test_hubble_tilt_cosmology.py` | Hubble-radius Λ bridge | 🔴 HISTORICAL RECEIPT |
| `test_nmax_drift.py` | Λ drift 0.023%/Gyr (no longer load-bearing) | 🔴 HISTORICAL RECEIPT |
| `test_nmax_investigation.py` | N_MAX cosmological-ratio reading | 🔴 HISTORICAL RECEIPT |
| `test_n_gap_hypotheses.py` | H1/H2 observable-projection hypotheses (Scenario 4) | 🔴 HISTORICAL RECEIPT |
| `test_half_power_and_dynamic_n.py` | √β projection + dynamic n | 🔴 HISTORICAL RECEIPT (strong-field use lives in joint calibration) |
| `test_combined_formula.py` | Combined-formula validation | 🔴 HISTORICAL RECEIPT |
| `test_tilt_models.py` | Tilt model comparison | 🔴 HISTORICAL RECEIPT |
| `test_tilt_scan.py` | Tilt parameter scan | 🔴 HISTORICAL RECEIPT |
| `test_lisa_tilt_scan.py` | LISA tilt scan (C1: conditional) | 🔴 HISTORICAL RECEIPT |
| `test_alpha_from_geometry.py` | ALPHA_COSMO geometric bracketing (Phase 5.2) | 🔴 HISTORICAL RECEIPT |
| `test_joint_cosmology_fit.py` | Joint SN+BAO+H₀ fit (Phase 5.2) | 🔴 HISTORICAL RECEIPT |
| `test_fixed_h0_fit.py` | Fixed-H₀ cosmology fit | 🔴 HISTORICAL RECEIPT |
| `test_cosmo_approach_a_fit.py` | Approach-A full fit | 🔴 HISTORICAL RECEIPT |
| `test_cosmo_correction_scaling.py` | Correction-scaling study | 🔴 HISTORICAL RECEIPT |
| `test_cosmological_consistency.py` | Low-z consistency battery | 🔴 HISTORICAL RECEIPT |
| `test_cosmological_replacement.py` | ΩΛ-replacement test | 🔴 HISTORICAL RECEIPT |
| `test_baryons_only_cosmology.py` | Baryons-only control (BAO χ²=409, ruled out) | 🔴 HISTORICAL RECEIPT (control) |
| `test_bao_diagnostic.py` | BAO diagnostic | 🔴 HISTORICAL RECEIPT |
| `test_planck_ruler.py` | Planck acoustic-ruler estimate (Scenario J prep) | 🔴 HISTORICAL RECEIPT |
| `test_two_scale_search.py` | Two-scale formula search | 🔴 HISTORICAL RECEIPT |

**Supernova pipeline (Scenario 7):**

| File | Role | Status |
|---|---|---|
| `test_pantheon_plus_fit.py` | Pantheon+ MU_SH0ES fit (signal suppressed by preprocessing) | 🔴 HISTORICAL RECEIPT |
| `test_pantheon_raw_fit.py` | Raw-magnitude Tripp recovery, 2.09σ | 🔴 HISTORICAL RECEIPT |
| `debug_sn_discrepancy.py` | Pipeline-suppression diagnosis | 🔴 HISTORICAL RECEIPT |

## Infrastructure and outputs

| File | Role | Status |
|---|---|---|
| `run_all_comparisons.py` | ESTIF-Gravity fork runner. Called `observational/compare_{eht_m87,ligo_gw,jwst_galaxies}.py`; targets deleted 21 Mar 2026 (`c971ab8`, archive wipe — collateral, not a retirement decision). Never ran the receipt suite; the old "batch runner" description was wrong. Retired 15 Jul 2026. | 🔴 HISTORICAL RECEIPT |
| `plots/` (57 png) | Figure outputs of the scripts above | 📦 OUTPUT — curated copies live in `results/validated/` |

---

## Verify against git

Run this to get the true add-date of every file under `tests/`, then correct any
❓ rows above and remove the flags:

```bash
cd /Users/peterangelov/estif_publication
git log --diff-filter=A --date=short --pretty=format:'%ad' --name-only -- tests/ | awk '/^[0-9][0-9][0-9][0-9]-/ {d=$0; next} NF {print d, $0}' | sort -u > tests/file_add_dates.txt
open tests/file_add_dates.txt
```

Caveat: if a batch of files entered in one bulk commit, their dates collapse to
the commit date — the header docstrings are then the tiebreaker.

---

**Index Version:** 1.0 · 6.4.1 · 15 July 2026
