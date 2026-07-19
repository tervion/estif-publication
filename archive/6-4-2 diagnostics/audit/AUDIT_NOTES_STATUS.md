# AUDIT NOTES — STATUS.md (main-documents pass, doc 1 of 7)
**Audit date:** 18 July 2026 · **Scope:** `docs/report/STATUS.md` (repo copy, stamped 6.4.1 / 12 Jul) vs the receipt corpus, RHAC record, DOC_UPDATE_PLAN, and findings of the four session-folder passes · **Method:** full read + claim-by-claim cross-check on Claude's machine; pointer targets verified (`PHASE2_DECLARATION.md`, `GAZTANAGA_COMPARISON.md` exist; `VALIDATION_REPORT.md` §3.0 confirmed at line 271); `src/` listing confirmed

**Overall verdict: PLANNED EDITS EXECUTED FAITHFULLY, BUT THE DOCUMENT IS A v6.4.0 BODY WEARING A v6.4.1 STAMP — AND IT INHERITS THREE ERRORS THE AUDITS HAVE SINCE DISPROVEN.** Everything DOC_UPDATE_PLAN scheduled for STATUS.md is in and correct (header block, Phase-2 closure, RHAC-009, C1/C2/C6 wording). The debt is what the plan never scheduled: zero RHAC-010 content under a 6.4.1 stamp, a test table that reproduces the UKN mis-wiring disproven on 17 Jul, a C3 omission in the pending-list, and a "letter ready" headline that the Jul 11 pass gated.

One self-correction first: AUDIT_NOTES_JUL9 §C referred to the runner as `run_simulation.py`; the actual file is **`src/estif_ec_gr_run_simulation.py`** — STATUS.md cites it correctly; my shorthand was the imprecision, not the document.

---

## A. Findings (ordered by severity)

### S-1 — HIGH · v6.4.1 stamp, v6.4.0 content; date also wrong
Header and footer both read **"Version 6.4.1 / Last Updated: 12 July 2026."** v6.4.1 is the 13 Jul a₀-horizon doctrine (RHAC-010) — and STATUS contains **none of it**: no horizon-background statement, no prefactor demotion, no mention of the four receipts (`a0_horizon_test.py` etc.). Two concrete consequences inside the document:
- The gravity table's a₀ row still reads "MOND a₀ **derived** … H₀cx₀/**√3**" with no RHAC-010 qualifier — the exact wording TEST_INDEX already corrects on `derive_mond_from_geometry.py` ("√3 prefactor now a working form, not derived; numbers stand").
- Priority Action 4 ("Derive x_c = 0.272 and/or the √3 factor … closes the N_MAX = 5/7·L chain") is pre-doctrine framing; RHAC-010's target-band / numerology-floor / river-edge structure is the current statement of that problem.
This also explains the Jul 13 stamp nit precisely: the footer version was bumped to 6.4.1 without the content or the date following. **Fix:** either add a short RHAC-010 block (a₀ row qualifier + one paragraph + Action 4 reword) and stamp 6.4.1 / 13 Jul honestly, or restamp the document 6.4.0 / 12 Jul as-is. The first matches the doc's role.

### S-2 — HIGH · Test table reproduces the disproven UKN descriptions (4th document in the fix batch)
The "Test Suite Status" table carries:
- `test_UKN.py` — "Phase 2 residual-stress census (strict A1) — sector empty"
- `test_UKN2.py` — "Phase 2 re-audit under A1′; growth + GW restored"
The Jul 12 pass proved (AST-level) that UKN is the **door-3 swirl draft only** and UKN2 is the **deepen-exact proof only** — no census, no re-audit, no GW, and "sector empty" is door 2's verdict, which UKN doesn't contain. STATUS.md therefore joins RHAC-004/005/006, PHASE2_DECLARATION.md, and TEST_INDEX as the **fourth surface** carrying this wiring; the batched fix (canonical `phase2_a1prime/` citations + draft-duplicate labels) must include these two rows.

### S-3 — HIGH (one-word fix) · Priority Action 1 omits C3 from the pending list
"Script-side **C2/C4/C5** remain pending." Re-verified 17 Jul: the committed `estif_jwst_growth_spec.py` is still pre-C3 (no σ8/S8 filter). The list should read **C2/C3/C4/C5**. As written, STATUS perpetuates the CORRECTIONS_v6.3.1 falsehood family (CORRECTIONS claims the C3 script was re-issued; JWST_TEST_SPEC v1.1 and TEST_INDEX assert the filter is in the file). Fourth member of that family, by omission.

### S-4 — MEDIUM-HIGH · "Gravity letter ready" contradicts the Jul 11 gate
Header status line: "**Gravity letter ready** (and strengthened)"; Priority Action 6: "Submit the gravity letter … **C1 + C6 wording applied**." True of STATUS itself — but the letter that would be sent is v6.3, which **violates C1 on the LISA 491 µs line** (Jul 11 finding B-1, filed as pre-submission critical), plus the stale limitation (v), the unattributed a₀ reference (B-3), and the 1.92/1.965 ambiguity (B-2). Until the `docs/Letter/6.4/` revision exists, "ready" is an overclaim. **Fix:** "ready pending the 6.4 revision (C1 LISA conditional; honorable-null limitation rewrite)" — one clause in each place.

### S-5 — MEDIUM · Executive-summary first sentence: taboo word + strict-A1 phrasing + overclaimed verb
"…deriving gravity, **dark energy**, and dark matter from the claim that 3D space is a **flat** hypersurface…"
- "flat" is the taboo term (project rule: "even" / "unbent") — and worse, it encodes **strict A1**, which RHAC-006 amended; under A1′ the honest phrase is "even on average, with matter sourcing local deviations."
- "deriving dark energy" is no longer the honest verb after the honorable null: Λ is imported, w = −1 exactly; the framework derives gravity and constrains dark energy to a constant. The shop-window sentence should not overclaim what the two bullets under it immediately retract.

### S-6 — LOW · Footer's reconciliation claim is ~90% true; three residues remain
Footer: "'frozen eddy' → 'constant cosmic term (imported Λ)' … reconciled." Surviving residues: "**constant cosmic eddy** → cosmological constant" in Path One item 3 AND in Publication Readiness; test-table label "**Frozen-eddy** reframe" on `estif_task6_eddy_eos.py` (TEST_INDEX's wording — "E1/E2 eddy EoS falsified; constant-Λ limit ties ΛCDM" — is the compliant version). The gravity-sector "non-vacuum eddy background" and DM-sector ρ_eddy usages are the *dark-matter* eddy object and are fine.

### S-7 — LOW / VERIFY IN ESTIF_CONCEPT PASS · Path One items 1–2 likely stale
Item 1 ("write A2/A3 into the theory documents") and item 2 ("resolve the A1 conflict: retire the shrinking-ruler narrative") are v6.3-era to-dos. If DOC_UPDATE_PLAN's ESTIF_CONCEPT edits executed as the matrix says, item 1 may be done (→ ✅) and item 2's "A1 conflict" now collides with the A1→A1′ meaning of that phrase. Deferred to the ESTIF_CONCEPT.md pass; flagged here so the STATUS edit round can close both in one touch.

### S-8 — LOW · χ²/N = 1.92 used three times with no pipeline attribution
Same as letter finding B-2's second half: 1.92 is the **task6 constant-Λ-limit** figure; the Path One DR2 13-bin pipeline prints **1.965** (runtime-reconfirmed 17 Jul). One parenthetical at first use closes it.

### S-9 — STRUCTURAL / DECISION · The test table is a second, stale index
The "Test Suite Status (new scripts, v6.3)" table stops at v6.3.2 + the two UKN drafts: no `phase2_a1prime/` suite, no Fronts 1–3, no C-15, no `estif_P_derivation_attempt.py`, no target-lock, no a₀-horizon four. Since 15 Jul, TEST_INDEX.md is the canonical provenance index — maintaining a parallel table in STATUS is the same disease as `_index_view/` (a second staleness surface). **Proposal:** replace the table with 2–3 headline receipts + "full provenance: `tests/TEST_INDEX.md`," rather than extending it.

### S-10 — LINKAGE · Task 5b row updates when the Jul 9 C4 commit lands
"⚠️ under-marginalized (C4) … re-derive before use as the Path Two bar" is honest today — but the marginalized re-run **exists**, uncommitted, in `audit/Jul 9 Session/estif_pathone_aic_bic.py`. When Jul 9 action 1 executes, this row (and Priority Action 1) should be touched in the same commit.

---

## B. What STATUS.md gets right (for the record)

C1 discipline is **better than the letter's**: EHT row "⚠️ consistent; *deviation* conditional (C1)", LISA row "⚠️ conditional (C1)", Planck-Λ "calibration match, not a vacuum deviation" — exactly the standard the letter's §2.1 box missed. The C6 precision block (constraint sector selected; coupling adopted) is verbatim-faithful. The C2 epistemic-status block and the bootstrap paragraph match the receipts to the digit (0.3043 zero-input / 0.31408 +radiation / 0.96% / 0.53σ / a₀ 1.72% → 0.66% / DESI 1.618 fixed-ruler / r_u −1.07%) with the input ledger declared — and the "sub-percent OR zero-input" tension is implicitly present via the two-root statement (its runnable receipt remains the uncommitted `estif_bootstrap_verification.py`; Jul 13 action E-1). RHAC-009 closure is wired correctly, with a **resolving** pointer to VALIDATION_REPORT §3.0 (verified, line 271). The Ω_k ≡ 0 registered kill-shot is stated. Both doc pointers resolve. Every script filename cited exists in the repo, including `src/estif_ec_gr_run_simulation.py`. Path Two's closure narrative matches RHAC-004 and the target-lock receipt's NULL branch exactly.

---

## C. Proposed actions (STATUS.md scope) — for your decision, nothing executed

1. **S-1:** add the RHAC-010 block (a₀-row qualifier "prefactor: working form, not derived — RHAC-010"; one doctrine paragraph; Action 4 reword) and restamp **6.4.1 / 13 July** — or, if you prefer zero new content, restamp 6.4.0 / 12 July. My recommendation: the former.
2. **S-2:** fold the two test-table rows into the standing UKN wiring batch (RHAC-004/005/006 + DECLARATION + TEST_INDEX + STATUS, one commit).
3. **S-3:** "C2/C4/C5" → "C2/C3/C4/C5" (or apply the C3 script commit from the Jul 9 batch first, then the list shrinks honestly).
4. **S-4:** two clauses gating "ready"/"submit" on the letter 6.4 revision.
5. **S-5:** rewrite the opening sentence — suggested: "…deriving gravity from the claim that 3D space is an even hypersurface (A1′: even on average; matter sources local dents) carried through a 4D bulk at c; cosmology carries an imported cosmological constant (w = −1 exactly)."
6. **S-6:** three wording touches (2× "constant cosmic eddy", 1× "Frozen-eddy reframe" → TEST_INDEX wording).
7. **S-8:** one attribution parenthetical at the first 1.92.
8. **S-9:** your call — replace the table with a TEST_INDEX pointer (recommended) or extend it to v6.4.1.
9. **S-10:** ride the Task 5b row on the C4 commit.

All of these are one edit session; per the single-pass preference I'd fold them into the consolidated round from AUDIT_NOTES_JUL13 §F rather than commit STATUS alone.

---

## D. Carry-forward for the remaining six documents

- **README.md (next):** known open v6.3.2 edits + the CITATION.cff/DOI mismatch orbit it; check whether it repeats the "flat hypersurface" opener, the UKN table, "letter ready," and which version it claims.
- **CHANGELOG.md:** duplicated H1 (known); check the [6.4.1] entry covers RHAC-010 (if yes, it's the model for STATUS's S-1 block); verify [6.4.0] wording vs honorable-null standard.
- **VALIDATION_REPORT.md:** banner (line 7) already reads correctly A1′/RHAC-007/008/009-aware — deep pass to verify Parts 5/7/Summary "read as closed" instruction is tolerable vs needing inline edits; §3.0 Part-B closure confirmed present.
- **ESTIF_CONCEPT.md:** resolves S-7 (A2/A3 present? shrinking-ruler retired? A1′ stated?); taboo-term sweep.
- **RHAC.md:** the receipt-line fixes are already specified (Jul 12 pass C-1…C-3); this pass adds a full-content read of entries 001–010 vs verdicts.
- **ROADMAP.md:** ROADMAP_v6.3_APPEND integrated 6/6 (Jul 9); check post-v6.4 items (fronts, RHAC-010 redirect) exist.

**Notes version:** 1.0 · 18 July 2026
