# AUDIT NOTES — README.md (main-documents pass, doc 2 of 7)
**Audit date:** 18 July 2026 · **Scope:** `README.md` (repo copy, stamped 6.4.1 / 13 Jul) vs the receipt corpus, RHAC record, the four session-folder passes, and AUDIT_NOTES_STATUS carry-forward items · **Method:** full read + claim-by-claim cross-check on Claude's machine; all 13 reproduce-paths verified against repo listings; `requirements.txt` read and one receipt header opened to confirm a dependency claim (`tests/estif_task4_field_equation.py`). Per instruction, `CITATION.cff`, `SUMMARY_FOR_REVIEW.md`, and `docs/Letter/` were **not opened**; DOI statements below rely on the standing ledger and the Jul 11 letter-pass record only.

**Overall verdict: THE HONEST STAMP THE OTHERS LACK — BUT THE SHOP WINDOW STILL SELLS STRICT A1.** README is the best-integrated main document audited so far: the 6.4.1 label sits on genuine 6.4.1 content (RHAC-010 at three sites, honorable null twice, RHAC-009 retirement verbs, Gaztañaga no-novelty disclosed, bootstrap numbers to the digit). The standing "README v6.3.2 edits outstanding" ledger item is **obsolete** — those edits are in. The debt is concentrated at the top and at the edges: the opening sentence and the derivation paragraph state the pre-amendment axiom ("flat hypersurface", "flat 3-slices") while "(A1′)" appears twice in the status table without ever being defined; the OPEN PROBLEMS block contradicts STATUS on letter readiness; lock 2 of the four-lock ledger (γ ≈ 0.55) is absent; and the reproduce instructions fail on a clean install because `requirements.txt` never lists sympy.

Resolution of the STATUS-pass carry-forward questions for this document: "flat hypersurface" opener — **yes, repeated** (R-1). UKN table — **absent** (clean; no fifth wiring surface). "Letter ready" — **inverted**: README says *blocked* (R-2). Version claimed — **6.4.1, and honestly so** (§B).

---

## A. Findings (ordered by severity)

### R-1 — HIGH · Opener and derivation paragraph state strict A1; "A1′" used twice, defined nowhere
Line 1 of the body: "3D space is a **flat** hypersurface carried through a 4D bulk." The milestone section repeats it as an axiom statement: "Starting only from the flow axioms — **flat 3-slices**, everything moving through the bulk at speed c…". Two problems in one wording, same as STATUS S-5:
- "flat" is the taboo term (project rule: "even" / "unbent");
- it encodes **strict A1**, amended by RHAC-006. Under A1′ the honest phrase is "even on average; matter sources local dents; unsourced evenness returns."

Compounding it: the status table labels two rows "(A1′)" — growth and GW — and README **never defines A1′ or mentions the amendment**. A first-time reader meets an undefined symbol in the table while the opener asserts the superseded axiom. The document that most needs the RHAC-006 one-liner is the one visitors read first.

Nuance worth stating in the fix (this is precision, not damage control): the gravity derivation is **unaffected** by A1→A1′ — matter-sourced local dents are exactly what A1′ permits, so Task 4's chain survives verbatim. One clause saying so prevents a referee from reading the amendment as undermining the headline result.

Receipt-side note: `estif_task4_field_equation.py`'s frozen docstring carries the same pre-amendment wording ("A1 flat 3-slices…"). Per the frozen-receipts/relabeled-claims policy the **document is the fix surface**, not the receipt; an optional 2-line header note in the receipt (btfr_lensing-style) is a decision item, not a requirement.

### R-2 — HIGH · Letter readiness: README contradicts STATUS, in the opposite direction
OPEN PROBLEMS item 6: "Relativistic-consistency section (does 'flowing space' break GR?) — **blocks the gravity letter**." STATUS says "Gravity letter **ready** (and strengthened)." Both are wrong, differently:
- STATUS overclaims (S-4): the sendable letter is v6.3, which violates C1 on the LISA line — submission is gated on the 6.4 revision (Jul 11 findings B-1…B-3).
- README's blocker is stale in the other direction: the *physics* for the relativistic-consistency section exists as receipts — `estif_flow_signature_dynamics.py` (18/18: signature + SR + Newton from Euclidean bulk + universal c) and C-15 (GW/light share the null cone). What's missing is the **written section**, not the derivation.

**Fix:** one consistent formulation across both documents: "letter v6.3 exists; submission gated on the 6.4 revision (C1 LISA conditional; honorable-null limitation rewrite); relativistic-consistency section unwritten — receipts in hand (flow_signature 18/18, C-15)."

### R-3 — MEDIUM-HIGH · Lock 2 of the four-lock ledger is missing; Fronts 2 and the JWST null unrepresented
README carries lock 1 (Ω_k ≡ 0 kill-shot, with current 0.0007 ± 0.0019), lock 3 (w = −1, twice), and lock 4 (c_gw = c row) — but **not lock 2**: the growth index γ ≈ 0.55 (0.5455–0.5544 over z = 0–5, pull −0.23σ vs DESI) and the slip lock (Σ, η, μ) = (1, 1, 1). The growth row stops at f(0.5) = 0.76. Also absent: Front 2 (fold-back rule, first-black-hole channels) and Front 1's **JWST honest null** — the latter being exactly the kind of disclosure README's honesty register is built on. One extended row (or one added row) covers it. Stamp linkage: front3's repo copy is the 15 Jul reconstruction — see R-11 before wording the addition as unconditional.

### R-4 — MEDIUM-HIGH · Reproduce block fails on a clean install: `requirements.txt` lists no sympy (verified), no colossus
README: "Every script below runs on a laptop" after the standard `pip install -r requirements.txt`. Verified this pass: the **first command** of the reproduce block, `tests/estif_task4_field_equation.py`, opens with `import sympy as sp` — and `requirements.txt` (numpy, scipy, matplotlib, astropy only) never lists sympy. Same for the other symbolic receipts in the block (flow_signature, converse, tmunu engine — "Deps: sympy" in their own headers). `colossus` (front1/front2 dependency, Jul 13 action E-5) is likewise absent. Additional `requirements.txt` observations while it was open:
- header stamped "**ESTIF v6.1**" — two doctrine generations stale;
- the `numpy<2.0.0` pin is motivated by the src model code, but every audit re-execution ran green on numpy 2.4.4 — the pin comment deserves a re-check or a scope note ("pin applies to src/, receipts verified on 2.x");
- trailing line `#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2` — unexplained marker, Oct 2025 vintage. Flagged for disposition, **not corrected** (per the intentional-naming rule).

### R-5 — MEDIUM · DOI badge + BibTeX label v6.4.1 under the v6.2-era DOI
Badge and citation block both point `10.5281/zenodo.17261724` at "ESTIF … v6.4.1". Per the standing ledger item and the Jul 11 letter-pass record (which found the letter *accurately* labels that DOI "v6.2 (software)"), 17261724 is the v6.2 deposition. README is therefore the second surface of the known CITATION.cff mismatch: a reader following the DOI lands on a v6.2 archive while citing v6.4.1. Two clean resolutions: mint a new Zenodo version for v6.4.x and update badge + BibTeX, or (if 17261724 is being used as the concept DOI) say so explicitly in the citation block. `CITATION.cff` itself untouched this pass per instruction; it rides the same decision.

### R-6 — MEDIUM · μ(a/a₀) gap not disclosed
The letter carries it as limitation (ii) with `mu_extraction.py` as receipt (runtime-reconfirmed Jul 11: g_obs/g_N → exactly 0 at galactic x ≈ 10⁻⁶). README's open lists never state it. OPEN PROBLEMS item 4 (the local mechanism for the horizon scale) is *related but not the same claim* — a referee cross-reading letter and README will ask where the repo's front page discloses the interpolation-function gap. One line, citing the receipt, either in "What is still open" or as OPEN PROBLEMS 4a.

### R-7 — MEDIUM-LOW · Goal 2 presents Path Two as live; the same document closes it twice
"…the evolving-dark-energy version **is** Path Two." The Path table and "What is still open" both say CLOSED (11 Jul, honorable null, RHAC-004). One clause: "…**was** Path Two — closed 11 Jul 2026 (honorable null; w = −1 exactly)."

### R-8 — LOW · χ²/N = 1.92 three times, no pipeline attribution (S-8 twin)
Path table, status row, and the reproduce comment all quote 1.92 bare. It is the task6 **constant-Λ-limit** figure; the Path One DR2 13-bin pipeline prints **1.965** (runtime-reconfirmed 17 Jul). One parenthetical at first use.

### R-9 — LOW · "frozen eddy" wording survives in the reproduce block (S-6 family)
Section comment "the frozen-eddy reframe" and "frozen eddy (1.92) beats tilt (3.35)". RHAC-004's sanctioned term: "constant cosmic term (imported Λ)". Two touches; TEST_INDEX's task6 row wording is the compliant model.

### R-10 — LOW · Cosmetics / verify cluster
1. **Notation drift in Goal 3:** `r_universe`, `r_p`, and `R_H` all appear; the receipts and RHAC-009 state P as Ωm(a) = R_H(a)/**d_p**(a). Unify (or gloss once).
2. **"(§6)" in OPEN PROBLEMS 3** has no named document — the section reference dangles (the a₀-horizon guide sections cited in TEST_INDEX are the likely referent). Name the doc.
3. **"holds to 0.12%"** (Goal 3 consistency relation): no receipt named this pass prints 0.12%. Either name the script that does, or replace with the receipted bootstrap figures (0.96% / −1.07%).
4. **`estif_converse_flow_law.py` cited** (v1). Fine — v1 carries the matter-level Birkhoff content — but per the Jul 9 pass v2 is the cleaner Jul 9 formalization; optional dual cite.
5. **"(nickname pending NAMING.md)"** — forward reference to a Tier-C doc that doesn't exist yet; honest as phrased ("pending"). No action unless NAMING.md is being dropped.
6. **"flat rotation curves" / "flat curves"** (Goal 3, `estif_horizon.py` comment): domain-standard astrophysics usage, **ruled compliant** — the taboo targets grid/hypersurface ontology, not rotation-curve phenomenology. Recorded so the terminology sweep doesn't false-positive.

### R-11 — LINKAGE · Three README claims ride on pending actions from the session passes
1. **"c_gw = c derived (C-15 closed, RHAC-008)"** — the repo receipt is the 15 Jul reconstruction; original-vs-rebuild cross-validation done (Jul 13 pass, all shared numbers agree), but the rebuild notes gate closure on a **Mac mini run** (action E-2). The row's wording becomes unconditional the moment that run is green.
2. **Gaztañaga no-novelty + "fixed by the predictive postulate P"** — this content originates in `estif_sharpest_audit.py` / `estif_bootstrap_verification.py`, both **still uncommitted** (action E-1). No README statement is false, but its honesty headline traces to receipts that exist only in `audit/`.
3. **"Writing tasks — done in v6.3. A2 and A3 are now stated in ESTIF_CONCEPT.md; the shrinking-ruler narrative is retired"** — pre-answers STATUS S-7. Two documents now stake this on the ESTIF_CONCEPT pass; verification lands there.

---

## B. What README.md gets right (for the record)

**The stamp is honest — the mirror image of STATUS S-1.** 6.4.1 label, 6.4.1 content: RHAC-010 integrated at three sites (status-table a₀ row with "exact prefactor open"; the Goal-3 doctrine paragraph with the correct retirement verbs — "√3 prefactor… **retired**", "exact O(1) prefactor remains open"; OPEN PROBLEMS 1–4 restating the doctrine's problem structure). README and the CHANGELOG [6.4.1] entry are the two models for the STATUS S-1 fix.

**The v6.3.2 ledger item for README is obsolete.** The bootstrap block is present and matches the receipts to the digit: root 0.31408, 0.96% from Planck, 0.53σ, a₀ → 1.1920×10⁻¹⁰ (0.66%), r_universe −1.07%, DESI DR2 1.618 fixed-ruler — with the input ledger {H₀, T_CMB, N_eff} declared, so the "sub-percent OR zero-input" trap is avoided by construction (sub-percent is claimed *with* the three inputs named, never as zero-input).

**Honesty wiring is correct throughout:** C1 note verbatim-faithful (consistent-vs-deviation split; Planck-Λ as calibration match); C6 precision block faithful (constraint sector derived, 8πG coupling adopted); C2 stated as a consistency relation with the circularity named; RHAC-009 Part B closed with "Ωm derived" **retired**; Gaztañaga no-novelty disclosed prominently, with "cite him prominently" and the peer-review dates; honorable null wired twice with the 11 Jul date, RHAC-004, and a resolving pointer to `PHASE2_DECLARATION.md` (exists, verified); the opening bullet list claims gravity/time/expansion/dark-matter and — correctly — **no dark-energy bullet**; Ω_k kill-shot stated with the current measurement.

**Structural hygiene:** no UKN surface anywhere (the wiring disease stops at four documents); no parallel test-status table (the S-9 disease absent — the reproduce list is a curated set of headline receipts, not a shadow index); all 13 reproduce paths exist in the repo, including `src/estif_ec_gr_run_simulation.py` (21/21); DESI cache-beside-script description matches the TEST_INDEX cache policy; SPARC (RMS 15.6%, 87), f(0.5) = 0.76, 0.00σ / 1.0000 / 49.2σ all match receipts.

**No TEST_INDEX corrections arise from this pass** — a first.

---

## C. Proposed actions (README.md scope) — for your decision, nothing executed

1. **R-1:** rewrite the opener — suggested: "3D space is an **even** hypersurface (A1′: even on average; matter sources local dents) carried through a 4D bulk at c." Add one A1′/RHAC-006 definition line in the milestone section, reword "flat 3-slices" in the derivation paragraph, and add the one-clause precision that the gravity chain is unaffected by the amendment. Optional: 2-line header note in the task4 receipt.
2. **R-2:** reconcile the letter line with STATUS in one formulation (gated on the 6.4 revision); reword OPEN PROBLEMS 6 from "blocks" to "section unwritten; receipts in hand."
3. **R-3:** add lock 2 (γ ≈ 0.55, −0.23σ; slip (1,1,1)) to the status table; optional clauses for Front 2 and the JWST honest null. Word as unconditional only after the E-2 mini run (or stamp "(reconstruction; mini-run pending)").
4. **R-4:** `requirements.txt` commit: add `sympy`, add `colossus`, restamp the header, scope-note or revisit the numpy pin, decide the trailing APPROVED-FORK marker. Merges with Jul 13 action E-5.
5. **R-5:** Zenodo decision — mint a v6.4.x version DOI and update badge + BibTeX, or annotate 17261724 as the concept DOI. Same decision governs `CITATION.cff` (excluded this pass).
6. **R-6:** one μ-gap disclosure line citing `mu_extraction.py`.
7. **R-7:** Goal-2 clause: "was Path Two — closed 11 Jul 2026."
8. **R-8:** attribution parenthetical at the first 1.92 (task6 constant-Λ limit; Path One DR2 13-bin = 1.965).
9. **R-9:** two "frozen eddy" → "constant cosmic term (imported Λ)" touches.
10. **R-10:** notation unification (d_p), name the "(§6)" document, source or replace "0.12%", optional converse v2 dual-cite.

Per the single-pass preference, all of the above folds into the consolidated execution round from AUDIT_NOTES_JUL13 §F rather than a README-only commit.

---

## D. Carry-forward for the remaining five documents

- **CHANGELOG.md (next):** duplicated H1 (known); does [6.4.1] cover RHAC-010 fully (it is the presumed second model for STATUS S-1 — verify); [6.4.0] wording vs the honorable-null standard; whether anything records the **15 Jul events** (front3/C-15 reconstructions, TEST_INDEX creation) or the log stops at 13 Jul; 1.92 attribution check.
- **VALIDATION_REPORT.md:** banner already A1′/007/008/009-aware (Jul 13 matrix); deep pass on Parts 5/7/Summary "read as closed" tolerability; §3.0 confirmed at line 271; check for UKN mentions and 1.92/1.965 attribution.
- **ESTIF_CONCEPT.md:** now carries **three dependents** — STATUS S-7, README's "done in v6.3" claim (R-11.3), and R-1's need for a canonical A1′ statement to point at. Verify A2/A3 present, shrinking-ruler retired, A1′/RHAC-006 stated; full taboo-term sweep.
- **RHAC.md:** receipt-line fixes already specified (Jul 12 pass C-1…C-3); this pass adds the full-content read of entries 001–010 vs verdicts; RHAC-008's stale "not found" note clears with E-2.
- **ROADMAP.md:** 6/6 append integrated (Jul 9); check post-v6.4 items exist (fronts, RHAC-010 redirect of the old √3 item); letter-gating language consistent with the R-2 formulation.

**Notes version:** 1.0 · 18 July 2026
