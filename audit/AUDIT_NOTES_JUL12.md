# AUDIT NOTES — Jul 12 Session vs repo state
**Audit date:** 17 July 2026 · **Scope:** `audit/Jul 12 Session/` (8 files) vs `tests/phase2_a1prime/`, `tests/test_UKN*.py`, `docs/plan/RHAC.md`, `docs/plan/PHASE2_DECLARATION.md`, TEST_INDEX · **Method:** byte-level diff + AST-normalized comparison on Claude's machine (files copied via Filesystem MCP) + runtime re-execution of all 7 scripts (sympy 1.14.0; canonical verification remains the Mac mini)

**Overall verdict: CLEANEST COMMIT YET, BROKEN RECEIPT WIRING.** All 8 session files are committed byte-identical as `tests/phase2_a1prime/`, and the whole suite runs green against its own README regression block — including the strict-A1 sign flip (f = −0.911) and the A1′ restoration (f = +0.76). But the archive layer points the wrong way: RHAC-004/005/006 and PHASE2_DECLARATION.md cite **only** `test_UKN.py`/`test_UKN2.py` — which this pass proves are docstring-free **draft duplicates** of just two of the seven canonical scripts. RHAC-005 cites the wrong file outright, and one claim (the A1′ five-sector re-audit) has no machine receipt anywhere in the repo.

---

## A. Scripts — audit folder vs `tests/phase2_a1prime/` (7 scripts + README)

### Byte-identical, committed correctly (8 of 8)
`README.md`, `estif_p2_door1_rate_dial.py`, `estif_p2_door2_slosh_divergence.py`, `estif_p2_door3_swirl_ledger.py`, `estif_growth_nogo_law.py`, `estif_growth_nogo_audit.py`, `estif_a1prime_growth_restored.py`, `estif_a1prime_deepen_exact.py` — all byte-for-byte. Nothing missing.

### Runtime verification (all 7 re-executed, README order)
Every script exits 0 and matches the README's "Expected outputs" regression block:
1. **door1** — `sqrt(-2*G*M/r + 1)` / `3*H**2/(8*pi*G)` / `3*H` ✅ (consistency check; ADM-lapse fact asserted, per its own docstring)
2. **door2** — `0` / `[0, 0, 0]` ✅ (pure-divergence theorem for the whole scalar class — the suite's strongest receipt, as the README says)
3. **door3** — `0` / zero momentum residual / isotropic ⟨wᵢwⱼ⟩ / ledger `−(A²+B²+C²)/(32πG)` ✅ (free tangle forbidden; sourced tangle negative-definite)
4. **nogo_law** — `0` / `−0.911` ✅ (δ = H/H₀ exact; strict-A1 sign flip vs DESI +0.76)
5. **nogo_audit** — `[0, 0, 0]` / `0` / `0` / `0` ✅ (gate, slaving, solution, Raychaudhuri closure)
6. **growth_restored** — `0` / `0.76` ✅ (fade mode embedded; deepen-mode f(0.5) numeric, matches DESI RSD and Ωm(a)^0.55)
7. **deepen_exact** — `0` ✅ (exact symbolic proof of D₊ = H·∫da/(aH)³ via the defining-derivative substitution)

Cosmetic only: the README's expected-output block normalizes sympy's spacing/sign formatting (e.g. writes `-(A**2+B**2+C**2)/…` where sympy prints `(-A**2 - B**2 - C**2)/…`). Human-checkable as intended; would trip a naive automated string diff. No action unless regression-checking gets scripted.

The suite's "Known flags" section is a model of the honesty standard — it firewalls the asserted ADM fact, the ~10¹¹ order-of-magnitude estimate, the literature-backed LTB step, and explicitly states script 6 does NOT check GW speed (C-15 separate). The problems in section C below are all cases where documents *outside* the suite fail to match this internal precision.

---

## B. The UKN files resolved — draft duplicates, not companions

AST-normalized comparison (docstrings/comments/whitespace stripped, code structure compared):

- **`tests/test_UKN.py` ≡ `phase2_a1prime/estif_p2_door3_swirl_ledger.py`** — computationally identical.
- **`tests/test_UKN2.py` ≡ `phase2_a1prime/estif_a1prime_deepen_exact.py`** — computationally identical.

So the UKN files are the **docstring-free drafts** of scripts 3 and 7, sitting at `tests/` root under placeholder names. This corrects two prior statements:
1. Jul 9 carry-forward "origin of test_UKN.py / test_UKN2.py" — answered: they are the Phase-2/A1′ era originals. RHAC-004/005/006 (all dated **11 Jul**) cite them because they were the only receipts existing that day; the documented 7-script suite is the **12 Jul** formalization (hence the session-folder naming), and the archive was never re-pointed at it.
2. My own Jul 11 note calling UKN/UKN2 "symbolic companions" of the suite — wrong; they are **duplicates**, carrying zero content the canonical scripts lack.

Consequence: keep both files (mark, never delete — they are the RHAC-cited names), but each needs a 2-line header docstring: "Draft duplicate of `phase2_a1prime/<canonical>.py` (AST-identical); kept as the filename cited in RHAC-004/005/006. Canonical documented copy supersedes."

---

## C. Receipt wiring — RHAC and the Declaration vs the actual receipts

### C-1 — RHAC-004 cites 1 of 3 scripted doors
Receipt line: `tests/test_UKN.py · docs/plan/PHASE2_DECLARATION.md`. The census verdict spans four doors, three of which have scripts — but only the swirl draft is cited. **Door 2's pure-divergence theorem (the strongest receipt in the suite, carrying the "slosh ZERO" verdict RHAC-004 itself states) is uncited**, as is door 1. Also: RHAC-004 states "~1e11 below rho_L" as census fact; the suite README explicitly flags that figure as an order-of-magnitude estimate, not a receipt — one qualifier word would align them.
**Fix:** receipt line → `phase2_a1prime/estif_p2_door{1,2,3}_*.py` (+ keep UKN as historical name) · declaration; add "(order-of-magnitude)" to the 1e11 clause.

### C-2 — RHAC-005 cites the wrong file (the sharp one)
Receipt line: "`tests/test_UKN2.py` (prints [0,0,0], 0, 0, 0)." Two errors in one line:
- The output signature `[0,0,0], 0, 0, 0` belongs to **`estif_growth_nogo_audit.py`** (runtime-confirmed). `test_UKN2.py` prints a single `0`.
- `test_UKN2.py` ≡ the **deepen-mode proof** — an **A1′** receipt. It contains no strict-A1 no-go content at all. RHAC-005's own quoted figure f(0.5) = −0.91 is printed by **`estif_growth_nogo_law.py`**, which goes uncited.
**Fix:** receipt line → `phase2_a1prime/estif_growth_nogo_law.py` (prints 0, −0.911) + `phase2_a1prime/estif_growth_nogo_audit.py` (prints [0,0,0], 0, 0, 0). Optional precision: "(5 attacks bounced; 4/5 machine-verified, LTB literature-backed)" — matching the audit script's own docstring.

### C-3 — RHAC-006: right content, wrong name — and one claim with no receipt anywhere
Cites `tests/test_UKN2.py`. By content that IS the deepen proof, so the citation is right-by-content, wrong-by-canonical-name; the numeric f = 0.76 that RHAC-006 quotes lives in `estif_a1prime_growth_restored.py`, uncited.
Sharper issue: RHAC-006 states "**PHASE-2 RE-AUDIT under A1′: five sectors, no Lambda-printer — null stands. Receipts: tests/test_UKN2.py.**" The deepen proof contains no sector re-audit. **No script in the repo re-runs the census under A1′.** The re-audit claim is currently pencil/discussion-level with a receipt citation pointing at a file that doesn't contain it.
**Decision needed:** (a) write the small A1′ re-census receipt (likely ~20 lines: doors 1–3 verdicts are A2/A3-driven and unchanged; the shape door under A1′ allows sourced dents only, no vacuum Λ-printer — the door-2/door-3 algebra re-run with the A1′ shape term), or (b) reword RHAC-006 + declaration to "pencil-level re-audit; growth-sector receipts: deepen_exact + growth_restored." Option (a) matches the receipt standard better.
Minor, optional: RHAC-006's "gravitational waves (speed = c forced by A2)" anticipates C-15, which the suite README explicitly firewalls as open on 11–12 Jul; RHAC-008 (12 Jul) closes it properly. A cross-ref "(derived: RHAC-008)" would seal the sequence.

### C-4 — PHASE2_DECLARATION.md has the same wiring
Its Receipts line: "`tests/test_UKN.py` (strict A1), `tests/test_UKN2.py` (A1′ re-audit)" — same UKN-only citations, same mislabeled "A1′ re-audit," and its provenance note sends readers to "the test scripts" meaning two undocumented drafts covering 2 of 7 receipts. Same fix batch as C-1…C-3. The declaration is otherwise excellent — the pre-registered binding declaration (§2) is exactly the discipline standard, and the E1/E2 falsification table is a clean record.

### C-5 — RHAC-007…010 are wired correctly (preview)
All later entries cite canonical filenames (`estif_front1/2/3…`, `estif_C15_gw_sector.py`, `estif_P_derivation_attempt.py`, `a0_horizon_test.py`…). The wiring problem is confined to the 11 Jul entries. Note: RHAC-008 carries a 13 Jul note "receipt filename … not found in the repository — locate and commit, or amend" for C-15 and front3. Both files **are** in `tests/` now (Jul 9 notes D-1); the note and the TEST_INDEX ❌ rows clear together once the Jul 13 pass confirms the committed copies are the right versions (session copies are smaller: 6.51 vs 10.49 KB, 10.00 vs 12.86 KB — repo likely carries later revisions).

---

## D. TEST_INDEX corrections arising from this pass

1. **`test_UKN.py` row** — description "Phase 2 residual-stress census under strict A1 — sector EMPTY (RHAC-004)" is wrong twice: the file is the **door-3 swirl draft only**, and "sector EMPTY" is door 2's slosh verdict. Reword: "Draft duplicate (AST-identical) of `phase2_a1prime/estif_p2_door3_swirl_ledger.py`; the filename cited in RHAC-004." Status ⚠️ SUPERSEDED (kept as RHAC-cited name).
2. **`test_UKN2.py` row** — description "Growth no-go theorem + A1′ re-audit; growth & GW restored (RHAC-005/006)" is wrong: the file contains **only** the exact deepen-mode proof — no no-go, no re-audit, no GW. Reword: "Draft duplicate (AST-identical) of `phase2_a1prime/estif_a1prime_deepen_exact.py`; the filename cited in RHAC-005/006." Status ⚠️ SUPERSEDED (kept as RHAC-cited name).
3. **`phase2_a1prime/estif_growth_nogo_audit.py` row** — "Five audit attacks on the no-go (all bounced)": outcome claim fine; add receipt precision "(4/5 machine-verified; LTB literature-backed)" per the script's own docstring.
4. The other six suite rows are **verified accurate** against contents and runtime — no changes.
5. Rides in the standing single-pass index batch together with Jul 9 section-D and Jul 11 section-D items.

---

## E. Proposed actions (Jul 12 scope) — for your decision, nothing executed

1. **RHAC receipt-line edits** (RHAC-004/005/006): append canonical `phase2_a1prime/` filenames; fix 005's wrong-file citation; add "(order-of-magnitude)" to the 1e11 clause; optional C-15 cross-ref in 006. Per RHAC discipline this is an amendment to citations, not to verdicts — nothing about the physics record changes.
2. **PHASE2_DECLARATION.md**: same receipt-line fix + provenance note pointed at `tests/phase2_a1prime/` (and its README) instead of the bare drafts.
3. **A1′ re-census receipt decision** (C-3): write the ~20-line receipt, or reword the two "A1′ re-audit" citations to pencil-level. Your call; (a) is the cleaner precedent.
4. **UKN headers**: 2-line docstring in each draft pointing to its canonical copy (mark, never delete satisfied).
5. **TEST_INDEX**: items D-1…D-3 into the batched edit.
6. No commits needed for the eight session artifacts — all already in the repo, byte-identical.

---

## Carry-forward questions for the next session folder

- **Jul 13:** the two size mismatches (`estif_C15_gw_sector.py` 6.51 KB session vs 10.49 KB repo; `estif_front3_second_discriminator.py` 10.00 vs 12.86 KB) — establish later-revision vs fork; on confirmation, clear in one batch: the two TEST_INDEX ❌ rows AND the stale 13 Jul "not found" note inside RHAC-008. Byte-check the three same-size pairs (`estif_front1…`, `estif_front2…`, `estif_P_derivation_attempt.py`). `estif_bootstrap_verification.py` and `estif_sharpest_audit.py` have no repo counterparts under those names — locate (renamed? e.g. relation to `estif_bootstrap_closure.py`?) or flag for commit/disposition. `DOC_UPDATE_PLAN_v6.4.0.md` (session-only, 20.38 KB) — check execution status against repo docs; likely feeds directly into the phase-2 main-documents pass. `files/` subdir empty — confirm intentional.
- **Git-only datings** (unchanged): `test_UKN.py`, `test_UKN2.py`, `test_joint_calibration_derived.py`, `ripple_speed.py`, a₀-horizon suite — one run of the TEST_INDEX "Verify against git" one-liner on the Mac mini clears them all.
- **Phase-2 preview noted in passing:** `docs/plan/GAZTANAGA_COMPARISON.md` now exists (previously listed as pending Tier C) — include in the main-documents pass.

**Notes version:** 1.0 · 17 July 2026
