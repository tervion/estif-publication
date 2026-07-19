# AUDIT NOTES — Jul 11 Session vs repo state
**Audit date:** 17 July 2026 · **Scope:** `audit/Jul 11 Session/` (5 files) vs `tests/`, `docs/Letter/`, and TEST_INDEX · **Method:** byte-level diff on Claude's machine (all files copied via Filesystem MCP) + runtime re-execution of all four scripts (numpy 2.4.4 / scipy 1.17.1; canonical verification remains the Mac mini)

**Overall verdict: CLEAN COMMIT, DIRTY METADATA.** Unlike Jul 9, every Jul 11 artifact is committed byte-identical — nothing is missing from the repo. The debt is elsewhere: TEST_INDEX files both a₀ scripts under the wrong session with wrong-era roles, `btfr_lensing.py` carries no inline marker that its inverted-a₀ block is superseded, and the outreach letter (v6.3) still violates C1 on the LISA line despite CORRECTIONS saying "apply before submitting the gravity letter."

The session is one coherent package: **corrected a₀ tension + weak-lensing BTFR + μ-gap receipt + Path Two pre-registration + letter v6.3.** The letter cites all three receipts by filename.

---

## A. Scripts — audit folder vs `tests/` (4 scripts)

### Byte-identical, committed correctly (4 of 4)
`a0_tension_corrected.py`, `btfr_lensing.py`, `estif_pathtwo_target_lock.py`, `mu_extraction.py` — all byte-for-byte identical to the repo copies. No missing commits this session.

### Runtime verification (all four re-executed)
1. **`estif_pathtwo_target_lock.py`** — run against the committed `pathone_dr2_m.txt` / `pathone_dr2_c.txt`: output is **byte-identical to the committed `pathtwo_target_lock_output.txt`**. The pre-registration receipt is fully reproducible: pipeline check χ²/N = 1.965 at (w₀,wₐ)=(−1,0) matches "Path One printed 1.965"; best thawing shape (−0.860, −0.410), χ²/N = 0.657; 1σ box w₀∈[−0.94,−0.78], wₐ∈[−0.80,−0.05]. The NULL branch of its pre-registered pass/fail ("derivation returns exactly w=−1 → Path Two closes") is the branch that later fired (RHAC-004 / honorable null) — the receipt did its job.
2. **`a0_tension_corrected.py`** — prints exactly the letter's §6.3(vii) figures: **+1.17σ kinematic, +1.52σ lensing All R<1000** (plus LTG +0.69σ, ETG +2.04σ), with the footer naming the superseded +2.67σ/+2.12σ run.
3. **`btfr_lensing.py`** — INVERTED block prints **+2.67σ (kinematic), +2.12σ (All R<1000)** — precisely the "prior report" the letter declares superseded (correlated M*/L systematic wrongly beaten down by √N via per-bin weighting).
4. **`mu_extraction.py`** — g_obs/g_N: 6.5 at x≈0.27, 0.68 at x=0.2, ~10⁻¹² at x=0.1, **exactly 0 at galactic x≈10⁻⁶**; analytic and numeric derivatives agree at every row. The μ-gap receipt behaves exactly as letter limitation (ii) states.

### Finding A-1 — `btfr_lensing.py` has no supersession marker (violates "mark, never delete" spirit)
The committed script runs cold and prints +2.67σ/+2.12σ with **zero indication** those numbers are superseded. The supersession exists only in `a0_tension_corrected.py`'s footer, the letter's (vii), and nowhere in the file itself or its TEST_INDEX row. Scope of the error, precisely: only the `best_a0()` inversion is wrong (correlated systematic added per-point before inverse-variance summing). The forward χ²/N test (`report()`) is the still-current data-level successor of Jul 9's B-7; its "chi2/N (stat + 0.1dex)" column also treats the systematic per-bin, but no headline claim cites those numbers — no correction filed or needed there, worth one caveat comment at most.
**Fix:** 3-line header note in `tests/btfr_lensing.py` above `best_a0()`: inverted-a₀ block superseded by `a0_tension_corrected.py` (correlated-systematic fix); forward BTFR test remains current.

---

## B. Session documents vs repo

| Audit doc | Repo target | State |
|---|---|---|
| `ESTIF_letter_v6.3.md` | `docs/Letter/6.3/` | ✅ Filed, byte-identical. Content findings B-1…B-4 below. Note: `docs/Letter/6.2/` has Pages/ + PDF/ exports; 6.3 is markdown-only — no compiled outreach copy exists yet. |

### B-1 — C1 violation in the letter's strong-field box (the significant one)
§2.1 box: "…and the **predicted LISA GW delay (491 µs, S/N = 49σ)** with no free parameters after calibration." CORRECTIONS C1 (filed 8 Jul, marked **pre-submission critical**, "Apply before submitting the gravity letter") demotes exactly this deviation-from-GR claim to *conditional on the open eddy-stress sector*, because the ESTIF vacuum is exactly Schwarzschild. The letter was revised **three days after** C1 was filed and the LISA line got no conditionality. The rest of the box is compliant: EHT "0.00σ tension" is a *consistency* statement (C1 explicitly allows), "satisfies Planck Λ (ratio 1.0000)" is a calibration match not a derivation claim (C1 note + C6 satisfied — the letter never says "derived").
**Fix:** one clause + the C1 footnote on the LISA delay, or drop it from the box.

### B-2 — Limitation (v) is stale post-v6.4
"(v6.3 note: … Path Two thawing derivation is **in progress**.)" — Path Two has since closed with the honorable null (w = −1 exactly, RHAC-004/v6.4). The current truth is *stronger and cleaner* for outreach: the residual sector came back empty under pre-registered pass/fail, no dynamical dark energy is claimed, Path One stands. Also: the letter's "ties ΛCDM on DESI DR2, χ²/N = 1.92" is the **task6 constant-Λ-limit** figure (per its TEST_INDEX row); the Path One DR2 13-bin pipeline figure is **1.965** (runtime-confirmed above). Both are real numbers from different pipelines — a v6.4 letter should pick one and attribute it.

### B-3 — (vii) doesn't say which a₀ the tensions are measured against
The +1.17σ/+1.52σ figures are computed relative to the **bootstrap** a₀ = 1.192×10⁻¹⁰ (both scripts hard-code `A0 = 1.192e-10`), while the letter's headline derivation quotes 1.1793×10⁻¹⁰. Effect is small (kinematic +1.17σ → ≈+1.21σ vs 1.1793) but a referee will compute it. One clarifying phrase closes it.

### B-4 — Compliance items the letter already gets right (for the record)
Prefactor status matches RHAC-010 ("motivated, not formally derived"; "confirmatory test, not a mathematical uniqueness proof"). "Zero free parameters" is scoped to "the MOND derivation itself" with H₀, Ωm declared as imported measured inputs — consistent with the C2 honest-scope standard. Limitation (ii) states the μ gap precisely and cites `mu_extraction.py`. Reference "ESTIF v6.2 (software), doi:…17261724" is *accurately labeled* as the v6.2 DOI — the DOI-mismatch problem lives in `CITATION.cff`, not here.

---

## C. Corrections to the Jul 9 carry-forward map (session folders ≠ what we guessed)

The Jul 9 notes' "carry-forward questions" assigned files to the wrong sessions. Actual layout, now established by listing all three remaining folders:

1. **`btfr_lensing.py` and `a0_tension_corrected.py` originate Jul 11, not Jul 13.** The Jul 9 open question "btfr_lensing date to confirm in Jul 13 pass" is resolved *here*: it is the Jul 11 data-level successor of B-7, and its inverted block was corrected the same session.
2. **The `phase2_a1prime/` suite originates Jul 12, not Jul 11.** `audit/Jul 12 Session/` mirrors `tests/phase2_a1prime/` exactly by name and size (8 files incl. README). Byte-diff belongs to the Jul 12 pass.
3. **No session-folder copy exists anywhere** for: `test_UKN.py`, `test_UKN2.py`, `test_joint_calibration_derived.py`, `ripple_speed.py`, and the four a₀-horizon scripts (`a0_horizon_test.py`, `a0_prefactor_derivation.py`, `estif_flow_sim.py`, `estif_horizon.py`). `audit/Jul 13 Session/files/` is **empty** (.DS_Store only). Their add-dates can only come from git — the TEST_INDEX "Verify against git" one-liner is now the *only* route for these; the session-folder route is closed.
4. **Jul 13 preview:** session copies of `estif_C15_gw_sector.py` (6.51 KB) and `estif_front3_second_discriminator.py` (10.00 KB) are **smaller** than the repo copies (10.49 / 12.86 KB) — repo likely carries later revisions. Which version is the committed one, and whether the session copy is an earlier draft, goes to the Jul 13 pass (this also reframes Jul 9's question "were the right versions committed?").
5. **UKN content characterization** (read from repo, tiny files): `test_UKN.py` = SymPy swirl-sector census (divergence-free swirl; vacuum residual +w/2 ≠ 0 → free tangle forbidden, needs stirrer; isotropic quadratic average; energy ledger **negative**) — the symbolic companion of Jul 12's door-3 swirl ledger, RHAC-004 era. `test_UKN2.py` = SymPy proof that D = H·∫da/(aH)³ solves the growth ODE **exactly** — the symbolic companion of Jul 12's `estif_a1prime_deepen_exact.py`, A1′ era. Content supports ~Jul 11–12 dating; git for the exact day.

---

## D. TEST_INDEX corrections arising from this pass

TEST_INDEX v1.0 is dated **15 Jul** — it predates AUDIT_NOTES_JUL9 (16 Jul), so **none** of the Jul 9 section-D fixes are in it yet (front3/C15 still ❌, `test_joint_calibration_derived.py` still rowless, converse_flow_law2 still ❓, jwst row still asserts the uncommitted C3 filter). This pass adds:

1. **`a0_tension_corrected.py`** — wrong section (listed under v6.4.1 / 13 Jul a₀-horizon doctrine). Move to the 11–12 Jul block; remove ❓. Suggested role: "a₀ inversion, M*/L systematic fully correlated across bins (+1.17σ kin, +1.52σ lens primary); **supersedes the inverted-a₀ block of `btfr_lensing.py`**; cited in letter v6.3 §6.3(vii)."
2. **`btfr_lensing.py`** — same wrong section; remove ❓. Suggested role: "Data-level weak-lensing BTFR vs derived a₀ (Mistele+24 Table 2 recomputation; successor of Jul 9's B-7). ✅ forward test — ⚠️ inverted-a₀ block superseded by `a0_tension_corrected.py`."
3. **`test_UKN2.py` row description is wrong.** Index claims "Growth no-go theorem + A1′ re-audit; growth & GW restored (RHAC-005/006)" — the file contains **only** the exact deepen-mode solution check. The no-go law and its audit live in `phase2_a1prime/estif_growth_nogo_{law,audit}.py`; GW lives in C-15. Reword to: "Exact growing-mode (deepen) solution D = H·I verified symbolically (A1′ era; companion of `estif_a1prime_deepen_exact.py`)."
4. **`test_UKN.py` row** — description ("Phase 2 residual-stress census under strict A1 — sector EMPTY") is close but conflates: the script shows the swirl sector *forbidden in vacuum* and its ledger *negative*, i.e. the door-3 verdict in SymPy. Suggested: "SymPy swirl census: free tangle forbidden in vacuum (residual +w/2), no axis, ledger negative (RHAC-004; companion of `estif_p2_door3_swirl_ledger.py`)." Confirm RHAC attribution wording against RHAC.md in phase 2.
5. **`mu_extraction.py`** — dated "~10 July (RHAC-001 two-machine crossroads)" but its session copy sits in the Jul 11 folder and letter v6.3 (Jul 11) cites it as "(v6.3 update)". Suggest "~10–11 Jul"; git as tiebreaker. Role wording is accurate.
6. **`estif_pathtwo_target_lock.py` / `pathtwo_target_lock_output.txt`** rows — correct as written; add one clause worth having in the index: "output receipt reproduced byte-identical from committed script+data (audited 17 Jul)."
7. Items **1, 2, 3, 5 from Jul 9 section D remain open** and should ride in the same index edit (single-pass update, per working preference).

---

## E. Proposed actions (Jul 11 scope) — for your decision, nothing executed

1. Add the supersession header note in `tests/btfr_lensing.py` (finding A-1). 3 lines, no code change.
2. TEST_INDEX single-pass batch edit: this pass's items D-1…D-6 **plus** the still-unapplied Jul 9 section-D items, in one commit.
3. Letter: create `docs/Letter/6.4/` revision before any outreach send — (a) C1 conditional on the LISA 491 µs line (B-1, the critical one); (b) rewrite limitation (v) with the honorable null — Path Two closed, w = −1 exactly under pre-registered pass/fail (B-2); (c) state the a₀ reference for the (vii) tensions (B-3); (d) pick and attribute one DESI figure — 1.92 (task6 constant-Λ limit) or 1.965 (Path One DR2 13-bin) (B-2).
4. Optional, outreach-ready step: export letter 6.4 to Pages/PDF mirroring the `6.2/` folder structure.
5. No commits needed for the five session artifacts — all already in the repo, byte-identical.

---

## Carry-forward questions for the next session folders

- **Jul 12:** byte-diff the 8-file `phase2_a1prime` mirror (expected identical by name+size — confirm); does its `README.md` correctly index all 7 scripts; do RHAC-005/006 entries cite these filenames correctly; check whether `test_UKN.py`/`test_UKN2.py` should be cross-referenced from `phase2_a1prime/README.md` as symbolic companions.
- **Jul 13:** the two size mismatches — session `estif_C15_gw_sector.py` (6.51 KB) vs repo (10.49 KB), session `estif_front3_second_discriminator.py` (10.00 KB) vs repo (12.86 KB): establish which is the later revision and whether the repo copy supersedes or diverges (two-machine-fork pattern possible, as with C2 on Jul 9). Also: `DOC_UPDATE_PLAN_v6.4.0.md` (20.38 KB, session-only) — check execution status against the repo docs; `estif_bootstrap_verification.py` and `estif_sharpest_audit.py` have no repo counterparts under those names — locate or flag; `files/` subdir is empty — confirm intentional.
- **Git-only datings** (session-folder route closed): `test_UKN.py`, `test_UKN2.py`, `test_joint_calibration_derived.py`, `ripple_speed.py`, a₀-horizon suite. One run of the TEST_INDEX "Verify against git" one-liner on the Mac mini clears all five ❓/rowless items at once.

**Notes version:** 1.0 · 17 July 2026
