# AUDIT NOTES — CHANGELOG.md (main-documents pass, doc 3 of 7)
**Audit date:** 18 July 2026 · **Scope:** `CHANGELOG.md` (repo copy, latest entry [6.4.1] — 2026-07-13) vs `README.md` (the pass reference), the receipt corpus, RHAC record, and the four session-folder passes · **Method:** full read + entry-by-entry cross-check on Claude's machine; the five [6.4.1] repository-structure claims verified by live listings (`docs/`, `docs/plan/`, `docs/latex/`, `tests/`); all receipt filenames cited by any entry checked against the `tests/` listing. Per instruction, `CITATION.cff`, `SUMMARY_FOR_REVIEW.md`, and `docs/Letter/` were not opened.

**Overall verdict: THE MOST ACCURATE LEDGER IN THE REPO — BUT IT CARRIES THE UKN WIRING INTO A FIFTH DOCUMENT AND SKIPS TWO DAYS OF HISTORY.** The [6.4.1] entry is confirmed as the second model (with README) for the STATUS S-1 fix: full RHAC-010 coverage, the complete four-lock ledger including lock 2, and five repository-structure claims that all verify true — one of which retires a standing ledger item (`docs/LaTeX ` trailing space: fixed 13 Jul, the ledger is stale). The debt: [6.4.0]'s receipt lines cite only the UKN drafts (including RHAC-005's wrong-file error, reproduced verbatim), the [6.4.1] "Still pending" list omits C3 (fifth member of that falsehood family), the 10–11 Jul session — including a published-σ-figure supersession, exactly what changelogs exist to record — has no entry at all, and the ⚠️ receipt warning in [6.4.1] is now stale in the good direction (both files exist as validated reconstructions).

Resolution of the README-pass carry-forward questions for this document: duplicated H1 — **confirmed** (CL-10.1). [6.4.1] covers RHAC-010 fully — **yes, verified; second model for S-1** (§B). [6.4.0] vs the honorable-null standard — **wording meets the standard; receipt wiring does not** (CL-1, CL-5). 15 Jul events recorded — **no, the log stops at 13 Jul** (CL-3, CL-6). 1.92 attribution — **correct at origin: [6.3.0] ties it to `estif_task6_eddy_eos.py`; this is the attribution model for S-8/R-8** (§B).

One correction to the README pass first: AUDIT_NOTES_README §B stated "the wiring disease stops at four documents." That was checked against README only and is now falsified — CHANGELOG [6.4.0] is the fifth surface (CL-1).

---

## A. Findings (ordered by severity)

### CL-1 — HIGH · [6.4.0] is the fifth UKN wiring surface; RHAC-005's wrong-file citation reproduced verbatim
Three receipt problems inside one entry:
1. **Phase-2 receipt line** reads "(`tests/test_UKN.py`; `docs/plan/PHASE2_DECLARATION.md`)" — citing the door-3 swirl draft for the *full four-door census*. Door 2's pure-divergence theorem — the very "slosh ZERO" verdict quoted in the same sentence, and the suite's strongest receipt — is uncited, as is door 1. Identical to the RHAC-004 finding (Jul 12 pass C-1).
2. **The A1→A1′ section** attaches "(`tests/test_UKN2.py`)" to the strict-A1 no-go paragraph (RHAC-005 content: δ = H/H₀, f(0.5) = −0.91, five attacks). Same wrong-file error as RHAC-005 itself: UKN2 is AST-identical to the **deepen-exact proof** — an A1′ receipt containing zero no-go content. The −0.91 is printed by `phase2_a1prime/estif_growth_nogo_law.py`; the attack audit is `estif_growth_nogo_audit.py`. Both uncited.
3. **The A1′ restoration content** (D₊ = H·∫da/(aH)³ = GR growth, f(z=0.5) = 0.76) cites **no receipt at all** — canonical: `estif_a1prime_growth_restored.py` + `estif_a1prime_deepen_exact.py`.

The standing wiring batch (RHAC-004/005/006 + PHASE2_DECLARATION + TEST_INDEX + STATUS test table) must therefore include CHANGELOG [6.4.0]. Mechanism is the changelog's own [6.3.2] precedent: amend citations inline or add a retroactive note — this is an amendment to citations, not to verdicts; the physics record is unchanged. Historical nuance worth preserving in the fix wording: on 11–12 Jul the UKN drafts *were* the only receipts, so the citations were honest when written; the 7-script suite is the 12 Jul formalization the archive was never re-pointed at.

### CL-2 — HIGH (one-word fix) · [6.4.1] "Still pending" omits C3
"Script-side C2 (…), C4 (…), and C5 (…) remain unapplied." Re-verified 17 Jul: the committed `estif_jwst_growth_spec.py` is still pre-C3 (no σ8/S8 filter). The list should read **C2/C3/C4/C5** — fifth member of the C3 falsehood family (CORRECTIONS_v6.3.1, JWST_TEST_SPEC v1.1, TEST_INDEX, STATUS S-3, now CHANGELOG), this one by omission. Note: [6.3.2] carries the same three-item list, but that entry is frozen and inherited CORRECTIONS' (false) re-issue claim in good faith on 9 Jul — the fix surface is [6.4.1]'s live list. Cleanest resolution: land the Jul 9 C3 commit first, then the list shrinks honestly in the same touch (rides with STATUS S-3).

### CL-3 — MEDIUM-HIGH · The ⚠️ receipt warning in [6.4.1] is stale (in the good direction)
"Both receipt filenames above are cited in RHAC-007/008 but were not located in the repository as of this writing. Locate and commit them, or amend…" — honest discipline on 13 Jul; outdated now. Both files exist as the **15 Jul reconstructions** with declared rebuild notes, cross-validated 17 Jul against the surviving originals in `audit/Jul 13 Session/` (all shared numbers agree; reconstructions self-report PASS 4/4). The warning clears in the same batch as its twins — RHAC-008's "not found" note, the two TEST_INDEX ❌ rows, `_index_view/MISSING.txt` — **gated on the Mac mini run (Jul 13 action E-2)**. The replacement text should record the resolution, not delete the history: reconstruction committed 15 Jul; original preserved byte-exact in `audit/Jul 13 Session/`; cross-validated 17 Jul.

### CL-4 — MEDIUM-HIGH · The 10–11 Jul window has no entry; a published-figure supersession is unrecorded
Nothing sits between [6.3.2] (9 Jul) and [6.4.0] (12 Jul). Unrecorded notable changes, all committed 10–11 Jul:
- **The a₀-tension supersession**: +2.67σ/+2.12σ (inverted-a₀ block of `btfr_lensing.py`) superseded by **+1.17σ kinematic / +1.52σ lensing** (`a0_tension_corrected.py`, correlated-M*/L-systematic fix). A correction of previously reported σ figures is precisely what a changelog exists to record; today it lives only in a script footer and letter v6.3 §6.3(vii).
- **`estif_pathtwo_target_lock.py` + `pathtwo_target_lock_output.txt`** — the pre-registration receipt (χ²/N = 0.657 thawing bar; the NULL branch that later fired). The anti-curve-fitting discipline artifact has no changelog line.
- `btfr_lensing.py` (data-level BTFR, successor of Jul 9's B-7), `mu_extraction.py` (μ-gap receipt), letter v6.3 revision.

"All notable changes to this project are documented in this file" is currently false for that window. Options: **(a)** a retro-dated **[6.3.3] — 2026-07-11** entry explicitly marked "(added retroactively, [date])" per the [6.3.2] precedent — recommended, because the σ-supersession belongs where a reader would look for it chronologically; **(b)** fold into the forthcoming consolidated entry (CL-6). Your call.

### CL-5 — MEDIUM · The receiptless "re-audited under A1′" claim appears here too
[6.4.0]: "the Phase 2 null — **re-audited under A1′, holds**." Third surface of the Jul 12 pass's C-3 finding (RHAC-006, PHASE2_DECLARATION, now CHANGELOG): no script in the repo re-runs the census under A1′. The standing decision — write the ~20-line A1′ re-census receipt, or reword all citations to pencil-level — now covers three documents; whichever branch you pick, this line rides it.

### CL-6 — MEDIUM / DECISION · The log stops at 13 Jul; 15–18 Jul structural events unrecorded
Not yet logged: the two receipt reconstructions (15 Jul), **TEST_INDEX.md creation** (15 Jul — the canonical provenance index, a structural change to how the repo is navigated), `_index_view/` materialization, and the audit itself (16–18 Jul). Proposal consistent with single-pass discipline: **one new entry when the consolidated fix round lands** ([6.4.2] or [Unreleased] → dated at commit), covering the reconstructions + TEST_INDEX + the batched corrections, and absorbing the CL-3 resolution text (and CL-4, if option (b) is chosen).

### CL-7 — LOW · "~10¹¹× below ρ_Λ" needs the order-of-magnitude qualifier
[6.4.0] states the figure as census fact; the suite README explicitly flags it as an order-of-magnitude estimate, not a receipt. Same one-word "(order-of-magnitude)" fix as RHAC-004 (Jul 12 pass C-1) — ride the same touch.

### CL-8 — LOW · Notation drift inside a single paragraph
[6.4.0] RHAC-009 section writes "P (Ωm = R_H/**r_p**)" and "Ωm(a) = R_H/**d_p** crosses…" three lines apart; [6.3.2] uses r_p throughout. Canonical form per the receipts and RHAC-009: Ωm(a) = R_H(a)/d_p(a). Rides README R-10.1 (unify or gloss once, one decision for all surfaces).

### CL-9 — LOW · The unsourced 0.12% originates here
[6.3.1] C2: "The 0.12% agreement is a self-consistency…" — the origin of the figure README quotes (R-10.3 flagged: no receipt named this pass prints 0.12%). [6.3.1] is a frozen entry, so the fix surface is README's Goal 3, not the changelog — but if the sourcing hunt (name the script or replace with the receipted 0.96%/−1.07%) finds the figure unreproducible, a retroactive note lands here too. Flag only.

### CL-10 — LOW · Cosmetics / optional-retroactive cluster
1. **Duplicated H1** — "# ESTIF Changelog" twice at the top (known, Jul 13 pass §B). Delete one line.
2. **Date-separator inconsistency**: [6.4.x]/[6.3.x] use "—", [6.2.0]/[6.1.0] use "-". Trivial; touch only if editing anyway.
3. **Optional [6.3.2]-style retroactive notes** on frozen entries whose headline status later changed: [6.1.0] "MOND Derived — Zero Free Parameters" and its "√3 from isotropy" derivation language (retired by RHAC-010); [6.2.0]'s 1/√3 reframe (superseded by the horizon doctrine); [6.3.0]'s "χ²/N = 0.66 (Path Two target)" line (target closed by C4 + RHAC-004). Reverse-chronological reading already covers all three via the entries above them — decision item, defensible either way.
4. **"Audited (5 attacks, all bounced)"** — optional precision "(4/5 machine-verified; LTB literature-backed)" per the audit script's own docstring; same wording item as RHAC-005/TEST_INDEX (Jul 12 D-3).
5. **"Full rename pending `NAMING.md`"** — same dangling Tier-C forward reference as README R-10.5; rides the NAMING.md decision.

### CL-11 — LINKAGE · Three changelog statements ride on pending actions
1. [6.3.1] "Documentation writing tasks closed… A1–A3 are stated explicitly in `ESTIF_CONCEPT.md`, the shrinking-ruler narrative is retired…" — now the **third surface** staking this on the ESTIF_CONCEPT pass (with STATUS S-7 and README R-11.3). Verification lands at doc 5.
2. [6.4.1] four-lock and front rows word the front3/C-15 results as closed; the repo receipts are the 15 Jul reconstructions — unconditional wording is gated on the E-2 mini run (same gate as README R-11.1 and CL-3).
3. The [6.4.1] "Still pending" list shrinks when the Jul 9 script batch (C2 merge, C3+C4 commits, C5 application) lands — update in the same commit (with CL-2).

---

## B. What CHANGELOG.md gets right (for the record)

**All five [6.4.1] repository-structure claims verify true — and one retires a standing ledger item.** Live listings this pass confirm: `docs/latex` exists lowercase with no trailing space; `ESTIF_arXiv_Paper.zip` is unzipped in place (`arxiv_paper_unzipped/` present, grep-visibility restored); `estif_field_dynamics.tex` + preview PDF are at `docs/latex/`; the `tests/Computing growth history…` staging folder is gone; `PATH_ONE_CHECKLIST.md` and `ESTIF_document_update_guide.md` are at `docs/plan/`. **The standing ledger item "`docs/LaTeX ` trailing space + zip invisible to grep" is obsolete — fixed 13 Jul.** Second stale ledger item found by the main-documents pass (first: README's v6.3.2 edits).

**[6.4.1] is confirmed as the second model for the STATUS S-1 fix.** Full RHAC-010 coverage: the doctrine statement (a₀ ≈ c·H, de Sitter surface gravity), the retirement verbs on the √3/full-coefficient "derived" language, the explicit not-touched scope (empirical fits, BH exterior, Ωm bootstrap), all four a₀-horizon receipts, and the RHAC-010 record pointer.

**The x_c double-duty flag names its document** — `docs/plan/ESTIF_document_update_guide.md` §6, verified to exist — which resolves README R-10.2's dangling "(§6)": README should copy the name from here.

**The four-lock ledger is complete, including lock 2** (γ ≈ 0.55 & slip = 1) — the lock README is missing (R-3) — with the honest C-11 clause attached ("none of the four separates ESTIF-Core from flat ΛCDM at linear order; distinguishing content lives off the linear sheet"). CHANGELOG is the model for the R-3 addition.

**1.92 is receipt-attributed at its origin**: [6.3.0] ties the figure to `estif_task6_eddy_eos.py` in the same sentence — the attribution model for STATUS S-8 / README R-8. No bare 1.92 appears in the later entries.

**The [6.3.2] retroactive-supersession note is the mechanism precedent** for every fix this pass proposes on frozen entries (CL-1, CL-3, CL-4, CL-10.3): annotate, never rewrite; date the annotation.

**The numbers are exact.** [6.3.2] bootstrap block: 0.3043 zero-input / 0.31408 with {T_CMB, N_eff, h} / 0.96% / 0.53σ / r_universe 4.353×10²⁶ (−1.07%) / a₀ 1.1920×10⁻¹⁰ (1.72% → 0.66%) / DESI 1.618 fixed-ruler vs 1.919 — all match the receipts to the digit, and the entry demonstrates the "sub-percent OR zero-input" discipline *by construction* (sub-percent claimed only with the input ledger named). [6.4.1] fronts: f(0.5) = 0.7603, γ 0.5455–0.5544, pull −0.23σ, ζ = 10⁻⁵ — match the originals and the reconstructions (cross-validated 17 Jul). [6.4.0]: f(0.5) = −0.91 vs +0.76, Ω_k 0.0007 ± 0.0019, ⟨δρ⟩ = −⟨ω²⟩/32πG — all match. src relabel: x₀ = 0.3107, a₀ = 1.1793×10⁻¹⁰, 21/21, "no computational change" — consistent with the Jul 9 §C src findings.

**Honesty wiring elsewhere is faithful**: C1/C2/C6 blocks verbatim-faithful to CORRECTIONS; C5 honestly "⬜ not yet applied" (still true — never applied anywhere); [6.4.0]'s honorable-null wording meets the standard (pre-registered, submitted unmodified, strong form, "does not derive dark energy, and it is now proven it cannot", frozen-eddy retraction); [6.4.1]'s header framing is honest ("No axiom, derivation, or closed-result status changes since 6.4.0"); the ⚠️ warning itself was the log doing its job on 13 Jul. Every receipt filename cited by any entry exists in `tests/` (verified against the listing); the changelog cites neither of the two still-uncommitted honesty receipts, so no citation is false.

**Taboo-sweep rulings (recorded so later sweeps don't false-positive):** "flat ΛCDM" ([6.4.1]), "flat FRW" ([6.3.0]), and "flat rotation curves" ([6.4.1]) are domain-standard cosmology/astrophysics usage, compliant — extends the README R-10.6 ruling. "Flat-slice flow picture" ([6.3.1]) *is* grid-ontology usage, but in a frozen entry that was accurate on 9 Jul (pre-A1′), with the superseding [6.4.0] A1→A1′ entry sitting directly above it in reverse-chronological reading — no action.

---

## C. Proposed actions (CHANGELOG.md scope) — for your decision, nothing executed

1. **CL-1:** fold [6.4.0]'s receipt lines into the standing UKN wiring batch — now **six surfaces in one commit**: RHAC-004/005/006, PHASE2_DECLARATION, TEST_INDEX, STATUS test table, CHANGELOG [6.4.0] (+ `_index_view` regen/gitignore after). Amend citations to the canonical `phase2_a1prime/` filenames, keep the UKN names as the historical citation, add the two missing restoration receipts.
2. **CL-2:** add C3 to the [6.4.1] pending list — or land the Jul 9 C3 commit first and update the list in that commit (preferred; same touch as STATUS S-3).
3. **CL-3:** after the E-2 Mac mini run is green, replace the ⚠️ box with the resolution note (reconstruction 15 Jul; original preserved in `audit/Jul 13 Session/`; cross-validated 17 Jul) — same batch as RHAC-008's note, the TEST_INDEX ❌ rows, and MISSING.txt.
4. **CL-4:** decision — retro-dated **[6.3.3] — 2026-07-11** marked "(added retroactively)" covering the σ-supersession, the target-lock pre-registration, btfr_lensing, mu_extraction, and letter v6.3 (recommended), or fold into the CL-6 entry.
5. **CL-5:** rides the standing A1′ re-census decision (write the ~20-line receipt vs reword to pencil-level) — three surfaces now.
6. **CL-6:** decision — one consolidated new entry when the fix round lands, covering the 15 Jul reconstructions + TEST_INDEX creation + the batch itself.
7. **CL-7:** "(order-of-magnitude)" qualifier — ride the RHAC-004 touch.
8. **CL-8 / CL-9 / CL-10:** ride their existing decisions (R-10.1 notation, R-10.3 sourcing, H1 deletion, optional retroactive notes, D-3 precision wording, NAMING.md).
9. **Ledger maintenance:** strike "`docs/LaTeX ` trailing space / zipped paper" from the outstanding-items list — verified fixed 13 Jul.

Per the single-pass preference, all of the above folds into the consolidated execution round from AUDIT_NOTES_JUL13 §F rather than a CHANGELOG-only commit.

---

## D. Carry-forward for the remaining four documents

- **VALIDATION_REPORT.md (next, doc 4):** banner already A1′/007/008/009-aware per the Jul 13 matrix; deep pass on whether Parts 5/7/Summary "read as closed" is tolerable or needs inline edits; §3.0 Part-B closure confirmed at line 271; check for UKN mentions (would make a sixth wiring surface), 1.92/1.965 attribution, C1 conditionality on the EHT/LISA rows, and the footer stamp (6.4.1 per the matrix — verify content matches the stamp, the S-1 disease).
- **ESTIF_CONCEPT.md (doc 5):** now carries **four dependents** — STATUS S-7, README R-11.3, CL-11.1, and R-1's need for a canonical A1′ statement to point at. Verify A2/A3 present, shrinking-ruler retired, A1′/RHAC-006 stated, v_flow = cx₀ relabelled as the sideways component; full taboo-term sweep.
- **RHAC.md (doc 6):** receipt-line fixes already specified (Jul 12 C-1…C-3); this pass adds the full-content read of entries 001–010 vs verdicts; the stale 008 "not found" note clears with E-2 (same batch as CL-3); the ~1e11 qualifier lands with CL-7.
- **ROADMAP.md (doc 7):** 6/6 append integrated (Jul 9); check post-v6.4 items exist (Fronts 1–3, the RHAC-010 redirect of the old √3 item); letter-gating language consistent with the R-2/S-4 formulation.

**Notes version:** 1.0 · 18 July 2026
