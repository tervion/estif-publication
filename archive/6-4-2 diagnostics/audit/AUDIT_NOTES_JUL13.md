# AUDIT NOTES — Jul 13 Session vs repo state
**Audit date:** 17 July 2026 · **Scope:** `audit/Jul 13 Session/` (8 files + empty `files/` subdir) vs `tests/`, `tests/_index_view/`, `docs/plan/RHAC.md`, DOC_UPDATE_PLAN targets, TEST_INDEX · **Method:** byte-level diff + unified diffs on Claude's machine (files copied via Filesystem MCP) + runtime re-execution of all 7 scripts and both repo reconstructions (numpy 2.4.4 / scipy 1.17.1 / sympy 1.14.0 / colossus installed; canonical verification remains the Mac mini) + grep-level execution matrix of the doc-update plan

**Overall verdict: THE SESSION THAT EXPLAINS THE CRISIS — AND MOSTLY FIXED ITSELF.** This folder is the origin of the "cited, not in repo" wound: two receipts (front3, C-15) were written on 12 Jul in an ephemeral container and never committed. The repo's answer was honest and unusually good — both were **reconstructed on 15 Jul with declared rebuild notes**, upgraded with PASS/FAIL self-checks, and this pass validates the reconstructions **against the surviving originals** (which turn out to live in this very audit folder): every shared physics number agrees. Two other receipts (`estif_sharpest_audit.py`, `estif_bootstrap_verification.py`) were **never committed at all** and remain the last missing receipts in the corpus. The 21-edit DOC_UPDATE_PLAN is fully executed across all seven target docs.

---

## A. Scripts — audit folder vs `tests/` (7 scripts)

### Byte-identical, committed correctly (3)
`estif_front1_growth_sigma8_jwst.py`, `estif_front2_first_hole_recipe.py`, `estif_P_derivation_attempt.py`.

### Repo carries a declared RECONSTRUCTION; session copy is the surviving original (2)
1. **`estif_C15_gw_sector.py`** — session 6.51 KB (original, 12 Jul) vs repo 10.49 KB (rebuild, 15 Jul). The repo rebuild note states the original "never reached the repository" and was "recovered from the session transcript." **The premise is now outdated: the original survives, byte-preserved, in `audit/Jul 13 Session/`.** Diff verdict: the reconstruction is a strict rigor upgrade — it *machine-verifies* the claim the original only asserted in prose (it constructs □_g, extracts the principal symbol from a plane wave, and proves **P_GW − P_light = 0 for arbitrary v, N**), adds metric-inverse and v±c sanity checks, a GW170817 structural check, a 4-item PASS/FAIL ledger with exit code, and an HONEST SCOPE paragraph per C-11. Same PG metric, same conventions, same verdict narrative.
2. **`estif_front3_second_discriminator.py`** — session 10.00 KB (original) vs repo 12.86 KB (rebuild, 15 Jul, same note). Upgrades: three-lock → **four-lock ledger** (adds C-15), explicit observed γ = 0.58 ± 0.11 (DESI DR1 PV+ShapeFit) comparison, slip algebra Σ = μ(1+η)/2, a μ-deviation growth-ODE demo showing the lock "has teeth," ODE cross-check + 6/11 high-z limit check, PASS/FAIL exit code, tightened WHAT-THIS-IS-NOT scope.

**Cross-validation result (the check the rebuild notes asked for, done against the originals):** all shared numbers agree between original and reconstruction — front3: γ(0) = 0.5544, high-z → 6/11 = 0.5455, spread over z = 0–5, pull **−0.23σ**, identical in both; C-15: same cone, u = v ± c, |c_gw/c − 1| = 0 structural pass, identical in both. Reconstructions self-report **PASS 4/4** and **4/4 — C-15 CLOSED**. Remaining gate: the rebuild notes require a **Mac mini run before RHAC-007/008 are treated as closed** — that run is yours; this pass supplies the original-vs-rebuild agreement the notes could not have.

### Never committed at all (2) — the last missing receipts
**`estif_sharpest_audit.py`** (parameter-honesty ledger) and **`estif_bootstrap_verification.py`** (the numbers behind it). Both are explicitly listed in DOC_UPDATE_PLAN's "Receipts to file first (copy the 7 scripts into `tests/`)" — **that instruction was executed 5/7**. Both run clean and carry load-bearing honest-scope content:
- bootstrap_verification: zero-input root **Ωm = 0.30433 → a₀ off 3.75% (−1.21σ)**; +radiation root **Ωm = 0.31408 → a₀ off 0.66% (+0.53σ)**; naive-vs-repo-exact horizon cross-check agrees to 1.79×10⁻¹³; r_universe back-predicted 4.353×10²⁶ m (−1.07% vs the 4.4×10²⁶ hardcode). **"Sub-percent OR zero-input, not both"** — this file IS the receipt for that sentence.
- sharpest_audit: P tagged as a **fourth postulate** in every relevant line; ~4 imports counted; a₀ named the one irreducible native win; the **Gaztañaga-overlap novelty flag** originates here (SUMMARY_FOR_REVIEW's "Gaztañaga comparison delivered" line answers it).
These are cited nowhere in RHAC (checked), so no archive statement is currently false — but the honesty ledger's headline sentences trace to receipts that exist only in the audit folder. Commit is the obvious disposition.

### Runtime verification (all green; verdicts cross-checked against RHAC)
front1: **honest null** — g(9.1) = 0.9978, JWST tension inherited not worsened, f(0.5) = 0.76 (= DESI), fσ8 pulls −0.12σ / +1.03σ, strict-A1 counterfactual = outcome 5 (falsified), §5a two-sided filter PASS trivially. front2: fold-back recipe computed — star channel typical 3σ → ~212 Myr (rarest ~28 Myr), no-star 10⁴–10⁶ M☉ (rare 5σ → 114 Myr), heavy seed OK / light strained to 10⁹ M☉ by z = 7, primordial closed under ζ = 10⁻⁵. P_derivation: crossing-once obstruction, both escape routes closed, "Ωm derived from the axioms" retired. Every number matches RHAC-007/009 verbatim. Note: **front1/front2 require `colossus`** (installed here via pip) — worth confirming it's in `requirements.txt` (phase 2).

---

## B. Session document — DOC_UPDATE_PLAN_v6.4.0.md (session-only, 20.38 KB)

A verbatim FIND/REPLACE edit script: sections A–G across RHAC.md, STATUS.md, PATH_ONE_CHECKLIST.md, ESTIF_CONCEPT.md, SUMMARY_FOR_REVIEW.md, VALIDATION_REPORT.md, ROADMAP.md, plus "Receipts to file first" and a deferred section H (README, CHANGELOG, CITATION.cff, CORRECTIONS).

**Execution matrix (grep-verified per edit, 21 edits): 21/21 DONE.** 17 verbatim; 4 footers/headers DONE-and-superseded by later v6.4.1 stamps (RHAC header+footer → v6.4.1/13 Jul; STATUS footer → 6.4.1; VALIDATION footer → 6.4.1). Beyond-plan improvements found: RHAC-010 appended, SUMMARY carries the delivered Gaztañaga comparison. Section H status: CHANGELOG has a full **[6.4.1] — 2026-07-13** entry (H's changelog item superseded-done); README + CITATION.cff remain the known open v6.3.2/DOI items (phase 2). The plan's only unexecuted instruction is the receipts line: **5/7 filed** (the two honesty receipts above).

Two stamp nits for phase 2: STATUS footer says "6.4.1 / 12 July" (version bumped, date not); CHANGELOG opens with a duplicated "# ESTIF Changelog" H1.

Disposition: the plan is now a historical execution record — leave in `audit/` (same treatment as ERRATA_v6.3.1.md).

---

## C. `tests/_index_view/` — the materialized index inherits every index error

`_index_view/` (by_status/, by_version/, MISSING.txt) is a file-tree rendering of TEST_INDEX v1.0 — and therefore reproduces its errors: UKN/UKN2 filed under `current/` + `v6.4.0/`, `a0_tension_corrected.py` + `btfr_lensing.py` under `v6.4.1/`, front3 + C-15 absent with a now-stale `MISSING.txt` naming exactly those two. It is **not** covered by `.gitignore`, so ~90 duplicate copies are (apparently) tracked — grep pollution and a second staleness surface. Decision needed: **regenerate after the batched index edit, or gitignore the view** (one-liner check on the mini distinguishes copies from symlinks: `ls -la tests/_index_view/by_status/current/ | head`).

---

## D. TEST_INDEX corrections arising from this pass

1. **`estif_front3_second_discriminator.py`** and **`estif_C15_gw_sector.py`** rows: ❌ "CITED, NOT IN REPO" → ✅ "reconstruction (15 Jul; rebuild note in header; original preserved byte-exact in `audit/Jul 13 Session/`; cross-validated against original 17 Jul)" — **contingent on your Mac mini run** per the rebuild notes.
2. **Add rows** for `estif_sharpest_audit.py` + `estif_bootstrap_verification.py` if committed (proposed roles above).
3. `estif_jwst_growth_spec.py` — **re-verified today: the committed copy is still pre-C3** (no σ8/S8 filter text at all). Jul 9 item stands unchanged.
4. Unchanged ❓/rowless items now confirmed **git-only** (session-folder route fully exhausted): `ripple_speed.py`, `test_joint_calibration_derived.py`, and the four a₀-horizon dates. All four RHAC-010 receipts exist in `tests/` with roles matching the index rows.
5. `files/` subdir in the session folder is empty (.DS_Store only) — flag as intentional-or-not; nothing lost that this audit can detect.

---

## E. Proposed actions (Jul 13 scope) — for your decision, nothing executed

1. **Commit the two honesty receipts** → `tests/` (+ index rows): completes DOC_UPDATE_PLAN's own 7/7 instruction and puts the "sub-percent OR zero-input" sentence on a runnable receipt.
2. **Run the two reconstructions on the Mac mini** (`repo` copies of front3 + C-15; both should exit 0 with PASS 4/4). On green, clear **in one batch**: the two index ❌ rows, RHAC-008's 13 Jul "not found" note, the RHAC header caveat "(receipt filenames unverified…)", and `_index_view/MISSING.txt`.
3. **One line in each rebuild note**: "Original recovered in `audit/Jul 13 Session/` (byte copy); reconstruction cross-validated against it 17 Jul 2026 — all shared numbers agree."
4. **`_index_view/` decision**: regenerate after the consolidated index edit, or gitignore.
5. `requirements.txt`: confirm `colossus` is listed (front1/front2 dependency).
6. Leave DOC_UPDATE_PLAN in `audit/` as a historical execution record.

---

## F. Consolidated cross-pass ledger (audit-folder phase complete)

The four session folders are now fully audited. Standing items by pass, for one batched execution round:

**Jul 9 (still open):** commit C4 `estif_pathone_aic_bic.py`; clean + commit C3 `estif_jwst_growth_spec.py` (re-confirmed pre-C3 today); merge C2 closing lines in `estif_pathone_cosmology.py`; apply C5 print-block to `estif_task6_eddy_eos.py`; B-7/ERRATA/JWST_TIMING dispositions; `src/files.zip` delete; optional C1 line in `test_calibration()`.
**Jul 11:** `btfr_lensing.py` supersession header; letter 6.4 revision before outreach (C1 LISA conditional; honorable-null rewrite of (v); a₀-reference clarification; 1.92-vs-1.965 attribution).
**Jul 12:** RHAC-004/005/006 receipt-line fixes (005 cites the wrong file); PHASE2_DECLARATION same; A1′ re-census receipt decision (write ~20 lines vs reword two citations); UKN draft headers.
**Jul 13:** items E-1…E-5 above.
**TEST_INDEX single batched edit** covering all four passes' section-D items, then `_index_view` regen/gitignore.
**PATH_ONE_CHECKLIST:** B-7 → ✅ citing `btfr_lensing.py` (+ `a0_tension_corrected.py` for the inversion), tally 19/2/12/2 → **20/2/11/2**.
**Git-only datings (one mini one-liner clears all):** UKN, UKN2, `test_joint_calibration_derived.py`, `ripple_speed.py`, a₀-horizon four.

**Next phase (per your original goals):** the main-documents pass — STATUS.md, README.md, CHANGELOG.md, VALIDATION_REPORT.md, ESTIF_CONCEPT.md, RHAC.md, ROADMAP.md in full, plus the known README/CITATION.cff/DOI items, the STATUS "6.4.1/12 July" stamp, the CHANGELOG duplicate H1, and the SUMMARY_FOR_REVIEW 6.4.0 lag. Several were already grep-touched by the plan matrix; the deep content-compliance read is what remains.

**Notes version:** 1.0 · 17 July 2026
