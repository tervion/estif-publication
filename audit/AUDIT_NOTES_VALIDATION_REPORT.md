# AUDIT NOTES — VALIDATION_REPORT.md (main-documents pass, doc 4 of 7)
**Audit date:** 18 July 2026 · **Scope:** `docs/report/VALIDATION_REPORT.md` (repo copy; header "v6.4.1 / 8 July 2026", footer "6.4.1 / 12 July 2026") vs the receipt corpus, RHAC record, the four session-folder passes, and the STATUS/README/CHANGELOG passes · **Method:** full read + claim-by-claim cross-check on Claude's machine; every cited script filename checked against the live `tests/` listing (all exist); the §3.1 "0.12%" figure recomputed from its own printed values; the banner's coverage instruction mapped against the actual state of Parts 2, 5, 7, and the Summary. Per instruction, `CITATION.cff`, `SUMMARY_FOR_REVIEW.md`, and `docs/Letter/` were not opened.

**Overall verdict: THE STRONGEST HONESTY CORE IN THE REPO, UNDER THE WORST HEADER IN THE REPO.** The body is a model of the discipline standard — C1 conditionality is wired better here than anywhere else, §3.0's bootstrap block is the best statement of the "sub-percent OR zero-input" honesty in any document, all ~22 cited receipts exist, and no UKN wiring reaches this file. But the header block is a **three-era chimera**: a v6.4.1 version number on v6.3's nickname, v6.3's date (8 July), and v6.3's status text — which uses the retracted "frozen eddy" term and declares the "gravity letter ready" in direct contradiction of the retraction banner sitting ten lines below it. And the 6.4.1 stamp is false in the STATUS S-1 sense: the document contains zero RHAC-010 content while Part 4 still presents the √3 prefactor as derived.

Resolution of the CHANGELOG-pass carry-forward questions for this document: banner A1′/007/008/009-aware — **confirmed; RHAC-010-unaware** (V-2). Parts 5/7/Summary "read as closed" tolerable? — **moot in the stated form: all three were already edited inline to CLOSED; the banner's caveat now points at a state that no longer exists, while missing the two surfaces that still read open** (V-3). §3.0 Part-B closure — **present and correctly wired**. UKN mentions — **none; the wiring disease stays at five surfaces**. 1.92/1.965 — **1.92 quoted for the Path One rows; the Path One pipeline itself (1.965) is uncited anywhere in the report** (V-5). C1 on EHT/LISA — **exemplary; this is the model surface** (§B). Footer stamp vs content — **S-1 disease confirmed, sharpest instance yet** (V-1, V-2).

One README-pass update first: **R-10.3 is resolved by this pass.** The unsourced "0.12%" traces to §3.1 here, and the arithmetic verifies from the section's own printed values: |0.310734 − 0.3111| / 0.3111 = 0.118% ≈ 0.12% (and the companion 0.10% = |0.261734 − 0.262| / 0.262). The figure is sourced, not phantom — README's fix is to cite the §3.1 pair (or state the ratio), not to hunt a script or replace the number. The CL-9 "retroactive note if unreproducible" contingency is void.

---

## A. Findings (ordered by severity)

### V-1 — HIGH · The header block is a three-era chimera that contradicts its own banner
Four fields, three eras:
- **"Model Version: ESTIF v6.4.1 — 'The Split'"** — a 6.4.1 number wearing v6.3's milestone nickname (`MILESTONE_v6.3_THE_SPLIT.md`).
- **"Date: 8 July 2026"** — the v6.3.0 date; the footer says 12 July; the 6.4.1 claim implies 13 July. Three dates in one document.
- **"Status: … frozen eddy = cosmological constant ties ΛCDM … Project split into Path One … and Path Two … Gravity letter ready and strengthened."** — v6.3's status text, verbatim in spirit: it uses **"frozen eddy," the exact label the v6.4.0 banner below retracts** (RHAC-004), announces the split as news, presents Path Two as newly opened (pre-null framing), and carries the "letter ready" overclaim (V-4).
The live status field of the validation report — the first thing a referee reads — asserts terminology and readiness that the document's own banner withdraws ten lines later. This is the sharpest instance of the S-1 family: STATUS had a wrong stamp; VALIDATION has a wrong stamp *and* a self-contradicting header. **Fix:** rewrite the header block wholesale — honest version + date (per the V-2 decision), post-null status line (imported Λ, w = −1 exactly; Path Two closed 11 Jul; letter gated per the agreed formulation), no retracted terms.

### V-2 — HIGH · 6.4.1 stamp, zero RHAC-010 content (STATUS S-1 twin)
The 12 Jul banner predates the a₀-horizon doctrine and the footer bump to 6.4.1 brought no content with it (the Jul 13 execution matrix recorded exactly this: footer stamp superseded-done, nothing else). Concrete consequences inside the document:
- **Part 4** presents the four-step MOND chain with "Step 3 — 3D isotropic projection: v_flow/√3" as a derivation, no qualifier. RHAC-010 demoted the √3 prefactor to a working form; TEST_INDEX's `derive_mond_from_geometry.py` row and README's Goal-3 paragraph carry the retirement verbs; the validation report — the document a referee checks the derivation against — does not.
- **Summary row "MOND a₀ … ✅ Derived"** — unqualified; needs "(prefactor: working form, RHAC-010)".
- **Part 7's x_c item** ("not yet geometrically derived (closes N_MAX = 5/7 × L)") is pre-doctrine framing — RHAC-010's target-band / numerology-floor / river-edge structure is the current statement (same as STATUS Action 4).
**Fix (same decision as S-1, and the two documents should move together):** either restamp 6.4.0 / 12 July as-is, or add RHAC-010 (one banner clause + Part 4 qualifier + Summary a₀-row qualifier + Part 7 reword) and stamp 6.4.1 / 13 July honestly. If the content branch is chosen, one added Summary row for the four locks (Ω_k ≡ 0, γ ≈ 0.55 & slip = 1, w = −1, c_gw = c) would close the Summary's silence on the 6.4.0 additions in the same touch.

### V-3 — MEDIUM-HIGH · The banner's "read as closed" instruction is doubly stale
Banner: "where Parts 5, 7, and the Summary still show it as an open derivation, read it as closed."
- (a) **Parts 5, 7, and the Summary no longer show it open** — all three were edited inline (Part 5's Path Two column 🟢 CLOSED; Part 7's struck-through item; Summary's Closed row). The DOC_UPDATE_PLAN edits executed; the caveat describes a pre-edit state.
- (b) The surfaces that **do** still read open are not in the banner's list: **Part 2.2's** stage-table row "0.66 | Path Two target ⚠️ under-marginalized", **Part 2.3's** closing "…a small *derived* thawing correction could do better if it can be derived. **That is Path Two.**", and the C4 caveat's forward-looking "before it is used to justify Path Two."
**Fix:** reword the banner clause to name 2.2/2.3 — or make the two small inline touches ("was Path Two — closed 11 Jul 2026, RHAC-004"; "the target is retained as the historical bar the null was tested against") and delete the clause. Either way the banner stops pointing at ghosts.

### V-4 — MEDIUM · "Gravity letter ready and strengthened" — third surface of the letter contradiction
STATUS says ready (S-4, overclaim — v6.3 letter violates C1 on the LISA line); README says blocked (R-2, stale in the other direction); VALIDATION's header says ready. Three documents, three positions. The agreed single formulation ("v6.3 exists; submission gated on the 6.4 revision — C1 LISA conditional, honorable-null limitation rewrite; relativistic-consistency section unwritten, receipts in hand") lands here as part of the V-1 header rewrite.

### V-5 — MEDIUM · Path One rows quote the task6 figure; the Path One pipeline is uncited in the entire report
Part 5's Path One cosmology cell and the Summary's Path One row both read "ties ΛCDM (χ²/N = 1.92)" — the `estif_task6_eddy_eos.py` constant-Λ-limit figure (correctly receipt-attributed at its Part 2 origin, per the CHANGELOG [6.3.0] model). But the actual Path One DR2 13-bin pipeline — `estif_pathone_cosmology.py`, printing **1.965** (runtime-reconfirmed 17 Jul) — appears **nowhere in the validation report**: not cited in Part 2, Part 5, Part 8, or the Summary. The report validates Path One's cosmology with a figure from a different pipeline while omitting Path One's own receipt. **Fix:** attribution parenthetical at the Part 5 first use (S-8/R-8 twin) + add `estif_pathone_cosmology.py` (and its 1.965) to the Part 2 script line or Part 8. Linkage: `estif_pathone_aic_bic.py` is likewise uncited — defensible while C4 is pending, but the marginalized rerun exists uncommitted in `audit/Jul 9 Session/` (S-10); when that commit lands, 2.2's target row and 2.3's C4 caveat update in the same touch, and the CPL row's v6.3-era rounding (w₀=−0.85, wₐ=−0.45, 0.66) reconciles against the target-lock print (−0.860, −0.410, 0.657).

### V-6 — MEDIUM · Part 7 (Known Limitations) omits the μ(a/a₀) gap
The referee-facing limitations list has the pressure sector, x_c, ρ_eddy, halos, CMB, peer review — and no interpolation-function gap. The letter carries it as limitation (ii) with `mu_extraction.py` as receipt (g_obs/g_N → exactly 0 at galactic x ≈ 10⁻⁶, runtime-reconfirmed 17 Jul); README lacks it too (R-6). The validation report is the second-worst place for the omission after README. **Fix:** one Part 7 line citing the receipt.

### V-7 — LOW · §3.1's C2 block still points at "an open task (RHAC Scenario H)"
Scenario H (derive r_universe independently) was closed by RHAC-009: the P route is proven impossible, and §3.0's status line plus the banner both say so — but §3.1's epistemic-status block, three screens up from §3.0's closure, still reads "an open task (RHAC Scenario H)." One clause: "(route closed: RHAC-009; see §3.0)."

### V-8 — LOW · Retracted "frozen eddy" — six instances (S-6/R-9 family, largest count)
Header status line (rides V-1), executive summary ("the frozen cosmic eddy"), 2.1 code-block label, 2.2 table row "frozen eddy = Λ (derived)", 2.3 table "Frozen eddy (Path One)", Part 8 task6 row "Frozen-eddy reframe". The banner's retraction covers the reading; the header instance is the only live-field offender. Body instances are frozen-part tolerable; if touched, TEST_INDEX's task6 wording ("E1/E2 eddy EoS falsified; constant-Λ limit ties ΛCDM") is the compliant model. Note also 2.1/2.2's "(derived)" flag on the constant: the corollary (constant ρ → de Sitter) is legitimately derived; the *value* is imported — the later-edited surfaces (Part 5, Summary) already say "imported" correctly.

### V-9 — LOW · r_p / d_p drift inside §3.0 itself
The P definition line uses r_p ("Ωm = R_H / r_p, with r_p the particle horizon"); honest flag 1 uses d_p ("Ωm(a) = R_H/d_p crosses once"). Same section, both notations. Rides the standing R-10.1/CL-8 unification decision (canonical per receipts and RHAC-009: d_p).

### V-10 — LOW · Part 8 nits
Task6 row shows bare ✅ without the C5-pending flag TEST_INDEX carries (C5 remains unapplied anywhere — still true). Part 8 is honestly scoped "v6.3 Test Scripts" + "prior records stand," so it is not a shadow index in the S-9 sense — but with TEST_INDEX now canonical, one closing pointer line ("full provenance: `tests/TEST_INDEX.md`") is the cheap S-9-consistent upgrade.

### V-11 — LINKAGE · Banner claims riding on pending actions
"The GW sector c_gw = c is DERIVED (RHAC-008)" and "the four exactness locks (RHAC-007)" are worded unconditionally; the repo receipts are the 15 Jul reconstructions — cross-validated 17 Jul against the surviving originals (all shared numbers agree), **gated on the Mac mini run** (Jul 13 action E-2), same gate as README R-11.1 and CL-3/CL-11.2. No wording change needed if E-2 lands first; otherwise "(reconstruction; mini-run pending)". The V-5 pipeline additions ride the C4 commit as noted.

---

## B. What VALIDATION_REPORT.md gets right (for the record)

**C1 discipline is the best in the corpus — the model surface.** The 1.1 LISA row is flagged conditional at the value itself; the full C1 block is verbatim-faithful to CORRECTIONS (consistent-vs-deviation split; Planck-Λ explicitly exempted as a calibration match); **every row** of the 1.4 mass table carries ⚠️ conditional plus a preamble stating why; the Summary row splits "Consistent; *deviations* conditional (C1)". The letter's §2.1 box should be fixed by copying this document.

**§3.0 is the strongest bootstrap statement in the repo.** Every number matches the receipts to the digit (0.3043 zero-input / 0.31408 with {T_CMB, N_eff, h} / +0.96%, 0.53σ / r_universe 4.353×10²⁶ m, −1.07% / a₀ 1.1920×10⁻¹⁰, 1.72% → 0.66% / SPARC ×1.00269 / DESI 1.618 fixed-ruler). The four honest flags are all present, none hidden: conditionality on P **with the RHAC-009 closure wired in** ("Part B is CLOSED"), Gaztañaga no-novelty with the memo pointer (`docs/plan/GAZTANAGA_COMPARISON.md` — exists, verified), the fixed-ruler caveat on 1.618, the 10%-scatter deflation of the 0.66%. The input ledger is explicit, so "sub-percent OR zero-input" is avoided by construction. Both cited scripts (`estif_omega_bootstrap.py`, `estif_bootstrap_closure.py`) are committed — §3.0 cites nothing that lives only in `audit/`.

**§3.1's C2 block is the honest-downgrade model** ("the numerical agreement is real and stands; only its status as a 'prediction' is withdrawn") — and it settles R-10.3: the 0.12% is the recomputable ratio of its own two printed values.

**The Path Two closure edits executed.** Part 5's column, Part 7's struck-through item, and the Summary row all read CLOSED with the honorable-null wording standard (residual sector empty; w = −1 exactly; Λ imported; RHAC-004). The banner's staleness (V-3) is a caveat problem, not a content problem.

**No UKN anywhere.** The wiring disease stays at five surfaces; the validation report is clean. Every cited receipt exists in `tests/` (checked against the live listing, ~22 filenames including the Part 8 table); no citation is false; the two uncommitted honesty receipts are cited nowhere here.

**The failure record is exemplary.** Part 2.2's stage table preserves the full descent (10.80 circular → 3.35 self-consistent → 1.92 constant limit; E1 = 3232, E2 = 754 falsified) with the C4 caveat honestly attached to the 0.66 target; Part 1.4's Λ-drift retirement note names what was retired and why; Part 6 keeps the two dead ancestors on the record. The C6 precision block is verbatim-faithful (constraint sector selected; 8πG coupling adopted; "ESTIF does not derive G").

**Taboo-sweep rulings (recorded so later sweeps don't false-positive):** "flat FRW" (0.2, Part 8) and "flat ΛCDM"-style comparator usage are domain-standard, compliant (extends the R-10.6/CL rulings). "A1 (flat 3-slices)" in Part 0.1 and "flat slices" in the 0.1 table are grid-ontology usage — but **banner-covered by design**: the banner's global instruction ("Throughout, 'A1 (flat slices)' now reads A1′…") is this document's chosen mechanism, and it is stated before any body use. Tolerable as-is; they convert automatically if the V-2 content branch triggers a body-edit round.

---

## C. Proposed actions (VALIDATION_REPORT.md scope) — for your decision, nothing executed

1. **V-1:** rewrite the header block — honest version + single date (per the V-2 branch), post-null status line with no retracted terms, letter line in the agreed gated formulation. This is the one non-negotiable touch regardless of the V-2 decision.
2. **V-2:** same decision as STATUS S-1, taken jointly for both documents: restamp 6.4.0 / 12 Jul, **or** add RHAC-010 (banner clause + Part 4 qualifier + Summary a₀-row qualifier + Part 7 x_c reword) and stamp 6.4.1 / 13 Jul. Recommendation: the latter, matching the S-1 recommendation — plus the optional four-lock Summary row.
3. **V-3:** reword the banner's read-as-closed clause to name 2.2/2.3, or make the two inline touches and drop the clause.
4. **V-4:** header letter line — the same single formulation as STATUS S-4 / README R-2, one wording across all three documents.
5. **V-5:** attribution parenthetical at Part 5's first 1.92 (task6 constant-Λ limit); add `estif_pathone_cosmology.py` (1.965) to Part 2's script line or Part 8; ride the C4 commit for the 2.2/2.3 target-row updates.
6. **V-6:** one μ-gap line in Part 7 citing `mu_extraction.py`.
7. **V-7:** "(route closed: RHAC-009; see §3.0)" clause in §3.1.
8. **V-8/V-9/V-10:** ride their standing decisions (S-6 wording touches if the body is opened; R-10.1 notation; optional TEST_INDEX pointer in Part 8).
9. **V-11:** no text action; clears with the E-2 mini run (same batch as CL-3, RHAC-008's note, the TEST_INDEX ❌ rows, MISSING.txt).

Per the single-pass preference, all of the above folds into the consolidated execution round from AUDIT_NOTES_JUL13 §F rather than a VALIDATION-only commit. Note the pairing discipline: V-1/V-2 and STATUS S-1 are one decision executed on two documents.

---

## D. Carry-forward for the remaining three documents

- **ESTIF_CONCEPT.md (next, doc 5):** carries **four dependents** — STATUS S-7, README R-11.3, CL-11.1 (all staking "A2/A3 stated, shrinking-ruler retired" on this pass), and R-1's need for a canonical A1′ statement to point at. Verify: A1′/RHAC-006 stated; A2/A3 present; shrinking-ruler retired; v_flow = cx₀ labelled as the sideways component (VALIDATION Part 4 has "(sideways component)" — CONCEPT should match); full taboo-term sweep; check whether it uses "frozen eddy" or post-null Λ language; stamp-vs-content check (the S-1/V-2 disease is now 2-for-2 on `docs/report/`).
- **RHAC.md (doc 6):** receipt-line fixes already specified (Jul 12 C-1…C-3); full-content read of 001–010 vs verdicts; the 008 "not found" note clears with E-2; the ~1e11 "(order-of-magnitude)" qualifier lands with CL-7; **add from this pass:** verify Scenario H's closure is recorded on H itself (V-7's origin — a reader landing on H should meet the RHAC-009 pointer), and confirm Scenario Q wording matches §3.0's citation of it.
- **ROADMAP.md (doc 7):** 6/6 append integrated (Jul 9); check post-v6.4 items exist (Fronts 1–3 follow-ups, the RHAC-010 redirect of the old √3 item); letter-gating language consistent with the V-4/S-4/R-2 formulation; check for a Path Two item still listed as live (V-3's disease in roadmap form).

**Notes version:** 1.0 · 18 July 2026
