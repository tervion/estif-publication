# AUDIT NOTES — Jul 9 Session vs repo state
**Audit date:** 16 July 2026 · **Scope:** `audit/Jul 9 Session/` (31 files) vs `tests/`, `src/`, and doc targets · **Method:** byte-level diff on Claude's machine (all files copied via Filesystem MCP)

**Overall verdict: GOOD, with 4 concrete repo-behind gaps.** The Jul 9 session was mostly integrated correctly. Docs were filed and later improved. But three script-side corrections (C3, C4, part of C2) were written on Jul 9 and never committed — they sit only in the audit folder while the repo carries pre-correction copies, and two repo documents (CORRECTIONS C3 note, TEST_INDEX) currently assert those fixes are in scripts where they are not.

---

## A. Scripts — audit folder vs `tests/` (15 scripts)

### Byte-identical, committed correctly (10)
`estif_bootstrap_closure.py`, `estif_fidelity_audit.py`, `estif_flow_signature_dynamics.py`, `estif_omega_bootstrap.py`, `estif_task4_field_equation.py`, `estif_task5_desi_selfconsistent.py`, `estif_task5b_cosmo_eos.py`, `estif_task6_eddy_eos.py`, `estif_tmunu_gauss_codazzi.py`, `estif_tmunu_task4.py`

Note on `estif_task6_eddy_eos.py`: identical in both places, and **C5 is genuinely unapplied in both** — `Wait --` (line 133) and `a^?` (line 136) are still present. The C5 fix exists only as prescribed replacement text inside `CORRECTIONS_v6.3.1.md`; it was never applied to any copy of the script. TEST_INDEX's "C5 pending" flag is accurate.

### Committed under a different name (1) — resolves a TEST_INDEX ❓
- audit `estif_converse_flow_law.py` **== `tests/estif_converse_flow_law2.py` byte-for-byte.**
- So: **v2 is the Jul 9 rewrite** ((i)(ii)(iii) formalization, engine-identity P0, momentum-constraint check). **v1** (`tests/estif_converse_flow_law.py`) is the earlier variant (A1–A4 assumption set, P1–P9b including matter-interior P9a/P9b, sink-flow control P7, offset control P8). They are **complementary, not redundant** — v1 carries the matter-level Birkhoff content and controls that v2 dropped.
- Suggested TEST_INDEX wording: v1 = "original converse test (A1–A4; P1–P9b incl. matter-level Birkhoff + controls)"; v2 = "Jul 9 rewrite (cleaner formalization + momentum constraint); audit-folder original". Remove the ❓.

### Repo BEHIND the session (2) — corrected scripts written Jul 9, never committed
1. **`estif_pathone_aic_bic.py` (C4).** The audit copy contains the full marginalized comparison block (Ωm and rd-scale freed for every model; marginalized Δχ²; "THE PATH TWO TARGET (honest, marginalized)"), explicitly labeled "Erratum v6.3.1". The repo copy is the pre-C4 version (fixed-parameter Q3; old 0.66-flavored target). TEST_INDEX's "C4 marginalization pending" is technically true of the repo — **but the work is already done and sitting in `audit/Jul 9 Session/`.**
2. **`estif_jwst_growth_spec.py` (C3).** The audit copy contains the σ8/S8 hard-filter section and the two-sided target. The repo copy lacks it entirely. This makes three repo statements currently **false**:
   - `CORRECTIONS_v6.3.1.md` C3: "(The re-issued estif_jwst_growth_spec.py now computes and prints the σ8 exclusion…)" — not true of the committed script.
   - `docs/report/JWST_TEST_SPEC.md` v1.1 names it as companion implementing the filter.
   - TEST_INDEX describes the row as "(C3 two-sided filter)".
   **Cleanup needed before commit:** the audit copy contains the section "[3]" **twice** (a "MANDATORY FILTER" draft immediately followed by the fuller "HARD FILTER" block) — a drafting artifact. Keep the HARD FILTER block, drop the first.

### Diverged both ways (1) — merge, don't overwrite
**`estif_pathone_cosmology.py` (C2).** Neither copy is a superset:
- **Repo copy is the later revision** (cites C2 explicitly, references the Ωm bootstrap): C2 applied in docstring + VERDICT print. ✅
- **But** its closing HONEST SCOPE line says "geometric Omega_m … **zero fitted cosmological parameters**" — an overclaim (H₀, Ωm are imported from Planck) and drops the "consistency relation" qualifier in that one line.
- **The audit copy's closing line is better:** "…CONSISTENCY relation + a physical dark-energy interpretation, with **no fitted DARK-ENERGY parameters (w = −1 fixed)**."
- Pattern (each copy has improvements the other lacks) is consistent with a two-machine fork of the RHAC-001 era. **Merge = repo version + audit's closing two lines.** TEST_INDEX's "C2 pending" is imprecise — C2 is ~90% applied; the residual offender is that one closing line.

### Never committed at all (1)
**`estif_weaklensing_btfr.py`** (B-7, published-fit level: Brouwer 2021 g† + Mistele 2024 transfer argument, tests both a₀ = 1.1793e-10 and bootstrap 1.1920e-10). Superseded by the stronger **data-level** `tests/btfr_lensing.py` (Mistele Table 2 recomputation, χ²/N per subsample, inverted a₀ ± σ) — which also **resolves the TEST_INDEX ❓ on `btfr_lensing.py`'s role**: it is the v6.4.1 data-level successor of Jul 9's B-7 script (date still to confirm from git in the Jul 13 pass).
**Decision needed:** commit B-7 script as a superseded/historical receipt, or leave it in `audit/`. Per "mark, never delete" it deserves at least an index mention if any doc cites B-7.
**Related doc bug found:** repo `docs/plan/PATH_ONE_CHECKLIST.md` line 72 still lists **B-7 as ⬜ unchecked** ("lowest-hanging fruit") — stale; should be ✅ citing `btfr_lensing.py`.

---

## B. Session documents vs repo

| Audit doc | Repo target | State |
|---|---|---|
| `CORRECTIONS_v6.3.1.md` | root | ✅ Filed; identical except dropped "target location" line |
| `PROGRESS_REPORT_2026-07-09.md` | `docs/report/` | ✅ Filed; same trivial normalization |
| `MILESTONE_v6.3_THE_SPLIT.md` | root | ✅ Repo **improved**: frozen-snapshot banner + inline C1/C4/C6 callouts |
| `JWST_TEST_SPEC.md` | `docs/report/` | ✅ Repo **improved**: v1.0 → v1.1, §5a σ8/S8 filter + "where the enhancement cannot come from" (correctly anticipated the strict-A1 no-go) |
| `CHANGELOG_v6.3_PREPEND.md` | `CHANGELOG.md` | ✅ Integrated 7/7 headers |
| `ROADMAP_v6.3_APPEND.md` | `docs/guide/ROADMAP.md` | ✅ Integrated 6/6 headers |
| `RHAC_v6.3_APPEND.md` | `docs/plan/RHAC.md` | ⚠️ 7/9 — "Scenario H (re-update)" and "Summary Statistics (v6.3)" absent. Both **superseded** (H closed by RHAC-009; stats replaced by v6.4.1 table). Content not lost; only the v6.3 intermediate step of Scenario H's trail is missing from the archive. Optional backfill. |
| `ERRATA_v6.3.1.md` | — (none) | ⚠️ Audit-only. This is the **earlier draft** of CORRECTIONS (same six items; draft severities C1 HIGH / C2 MEDIUM / C6 LOW vs filed Critical/Critical/Critical). Its own "target location: project root" was never executed — correctly, since CORRECTIONS supersedes it. Do **not** file to root; decide archive vs leave. |
| `JWST_TIMING_TEST_SPEC.md` | — (none) | ⚠️ Audit-only. Distinct doc: the timing/mass JWST test. Its target (`docs/plan/`) and companion script (`tests/estif_jwst_halo_massfunction.py`) **both absent from repo**. Premise resolved by Front 1's honest null (v6.4.0) → superseded planning doc. Decide: file with superseded banner, or leave. |
| `STATUS.md`, `README.md`, `VALIDATION_REPORT.md`, `ESTIF_CONCEPT.md`, `SUMMARY_FOR_REVIEW.md`, `PATH_ONE_CHECKLIST.md`, `CITATION.cff` | various | Repo versions all stamped **12–13 Jul, v6.4.0/6.4.1** → later descendants of the Jul 9 versions, not forks. Detailed content compliance deferred to the main-documents pass (phase 2). Two previews logged: SUMMARY_FOR_REVIEW stamp lags at 6.4.0; PATH_ONE_CHECKLIST B-7 row stale (above). |

---

## C. Core files (`src/`)

- All three files carry **v6.4.1 headers with the Jul 9 corrections integrated**: C2 in `estif_ec_gr_constants.py` (lines 13, 102, 170) and in `run_simulation`'s final summary ("consistency relation (C2)"); tilt retirement per v6.3/RHAC-004 ("Goal 2 … ✅ Receipt intact (claim RETIRED v6.3 — imported Λ)"); C1 in `estif_ec_gr_model.py` docstring; RHAC-010 prefactor status ("working form, NOT derived"). Computations frozen — consistent with the "frozen v6.2 receipts, relabeled claims" policy.
- **Residual (LOW):** `run_simulation.py`'s CALIBRATION CHECK (EHT+Λ+LISA) section and the final-summary calibration line print bare ✅ with no inline C1 conditional; C1 coverage is only via the module-docstring pointer and model.py. One added print line in `test_calibration()` would close it.
- **`src/files.zip`** = exact duplicate of the three live src files (zipped 15 Jul 10:25). Redundant, grep-invisible, duplicates tracked files. Recommend deleting (transfer artifact).

---

## D. TEST_INDEX corrections arising from this pass

1. `estif_front3_second_discriminator.py` and `estif_C15_gw_sector.py`: index says "❌ CITED, NOT IN REPO" — **both ARE in `tests/` now** (12.86 KB and 10.49 KB). Update rows to ✅ after verifying contents in the Jul 11/12 passes.
2. `test_joint_calibration_derived.py` (11.84 KB) exists in `tests/` with **no index row** — violates the same-commit maintenance rule. Add a row (origin likely identified in a later session folder).
3. `estif_converse_flow_law2.py` ❓ → resolved (section A above).
4. `btfr_lensing.py` ❓ → role resolved (data-level successor of B-7); date to confirm in Jul 13 pass.
5. `estif_jwst_growth_spec.py` row description "(C3 two-sided filter)" is currently false for the committed file; becomes true once the audit version is committed.

---

## E. Proposed actions (Jul 9 scope) — for your decision, nothing executed

1. Commit `audit/Jul 9 Session/estif_pathone_aic_bic.py` → `tests/` (C4 done; closes a "pending" that is actually finished).
2. Clean the duplicated "[3]" section in audit `estif_jwst_growth_spec.py`, then commit → `tests/` (closes C3 script-side; makes CORRECTIONS / JWST_TEST_SPEC / TEST_INDEX statements true).
3. Merge `estif_pathone_cosmology.py`: repo base + audit's closing HONEST-SCOPE lines (kills "zero fitted cosmological parameters").
4. Apply the C5 print-block replacement (text ready in CORRECTIONS) to `tests/estif_task6_eddy_eos.py` — never done anywhere.
5. `PATH_ONE_CHECKLIST.md`: check off B-7, cite `btfr_lensing.py`.
6. Decide disposition: `estif_weaklensing_btfr.py`, `ERRATA_v6.3.1.md`, `JWST_TIMING_TEST_SPEC.md` (commit-as-superseded / archive / leave in audit).
7. TEST_INDEX row updates per section D.
8. Optional: one C1 line in `test_calibration()`; delete `src/files.zip`.

---

## Carry-forward questions for the next session folders

- **Jul 11:** origin of `test_UKN.py` / `test_UKN2.py`, `phase2_a1prime/` suite, `estif_pathtwo_target_lock.py` — check integration + whether `test_joint_calibration_derived.py` originates here.
- **Jul 12:** Fronts 1–3 + C-15 scripts — verify committed copies match session copies (index claimed two were missing; they exist now — were the right versions committed?).
- **Jul 13:** v6.4.1 a₀-horizon suite (`a0_horizon_test.py`, `a0_prefactor_derivation.py`, `estif_flow_sim.py`, `estif_horizon.py`, `a0_tension_corrected.py`, `btfr_lensing.py`) — clears the remaining ❓ date/role flags.

**Notes version:** 1.0 · 16 July 2026
