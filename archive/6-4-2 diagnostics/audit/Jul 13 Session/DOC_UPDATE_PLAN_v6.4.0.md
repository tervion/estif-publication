# DOC UPDATE PLAN — session of 12 July 2026 (→ v6.4.0)

**What this pass records:** (1) C-15 closed — GW sector `c_gw = c` derived;
(2) Principle P proven **not derivable as a law** — priority redirect;
(3) Fronts 1–3 + the four exactness-lock ledger; (4) version/date hygiene on the
two docs that predate A1′ (`VALIDATION_REPORT.md`, `ROADMAP.md`).

**How to run it:** apply every edit below in one sitting, then a single commit
(`docs: consolidate to v6.4.0 — C-15 closed, P non-derivable, Fronts 1–3`).
Each `FIND:` must match the file **verbatim** (unicode included). If any FIND
does not match, stop and flag it — don't guess.

**Receipts to file first (copy the 7 scripts into `tests/`):**
`estif_front1_growth_sigma8_jwst.py`, `estif_front2_first_hole_recipe.py`,
`estif_front3_second_discriminator.py`, `estif_sharpest_audit.py`,
`estif_bootstrap_verification.py`, `estif_P_derivation_attempt.py`,
`estif_C15_gw_sector.py` → `/Users/peterangelov/estif_publication/tests/`.

---

## A. RHAC.md — append three new records

Insert the following **after** the `## AXIOM AMENDMENT · 2026-07-11 · A1 → A1′`
block's closing `---`, immediately **before** `### Updated Summary Statistics (v6.3)`.

```markdown
## RHAC-007 · 2026-07-12 · Fronts 1–3 executed — growth, first holes, second lock
FRONT 1 (growth under A1'): D+ clumping history computed; f(0.5)=0.7603 = DESI
  RSD anchor; fsigma8 pulls -0.12 / +1.03 sigma (PV z=0.07; FS BGS z=0.295).
  JWST verdict = OUTCOME 4 (honest null): g(9.1)=0.998, too-big-too-early tension
  INHERITED from LCDM, not relieved, not worsened. Strict-A1 counterfactual would
  have scored outcome 5 (falsified) — the fork rescued this test.
FRONT 2 (first-black-hole recipe): fold-back rule nu*sigma(M,0)*D(z)=delta_c.
  Star channel -> ~1e2 Msun hole at t~28 Myr (earliest-in-volume) / ~212 Myr
  (typical 3-sigma). Direct-collapse (no-star) channel -> 1e4-1e6 Msun holes;
  the only comfortable route to 1e9 Msun by z=7 (heavy seed OK, light strained).
  Primordial channel CLOSED under the zeta=1e-5 passport. All numbers inherited
  (LCDM growth); ESTIF content is structural (strict A1 forms NO hole).
FRONT 3 (second differing number): growth index gamma forced to ~0.55 by the
  empty residual sector + GR-equivalent D+; computed 0.5455-0.5544 over z=0-5;
  measured 0.58+/-0.11 (DESI PV+ShapeFit), pull -0.23 sigma, PASS. Companion
  slip lock (Sigma,eta,mu)=(1,1,1) exact.
FOUR-LOCK LEDGER (ESTIF-Core forces a POINT; GR family fits a REGION, 0 free
  dark params vs 1-3 fitted): Omega_k=0 (kill-shot) . gamma~0.55 & slip=1
  (Front 3) . w=-1 (RHAC-004) . c_gw=c (RHAC-008). HONEST: none separates
  ESTIF-Core from FLAT LCDM at linear order; all separate it from GR's extra
  freedoms. Distinguishing content lives OFF the linear sheet (nonlinear halos,
  N-body wall) + the exactness-as-law structure.
Receipts: tests/estif_front1_growth_sigma8_jwst.py,
  estif_front2_first_hole_recipe.py, estif_front3_second_discriminator.py.

## RHAC-008 · 2026-07-12 · C-15 closed — GW sector derived (c_gw = c)
CLAIM (asserted RHAC-006, now DERIVED): gravitational waves propagate at c.
DERIVATION: under A1'+A2 the world is ONE Lorentzian geometry (even-on-average
  slices + a flow); A3 + the empty residual sector (RHAC-004) forbid any second
  metric/field. TT waves are ripples OF that geometry; light rides null cones OF
  that geometry; the vacuum wave operator's PRINCIPAL SYMBOL is g^{mu nu} k_mu
  k_nu = the light cone. Shown symbolically for an ARBITRARY flow v: radial null
  speed u = v +/- c, identical for GW and light (flow tilts both cones the same).
  Flow friction (~H) and dent curvature are lower-order and cannot move the
  characteristic speed.
GW170817: predicts |c_gw/c - 1| = 0 exactly (measured bound ~1e-15). PASS
  STRUCTURALLY -- no dial exists to break it (the extra field every c_gw != c
  theory needs is forbidden by A3 + the empty residual sector).
A1' HINGE: strict A1 froze the TT sector (a ripple = a local deviation from exact
  evenness, forbidden) -- the SAME over-constraint that killed the growing mode
  (RHAC-005). A1' opened growth and radiation together, one mechanism, one fork.
Receipt: tests/estif_C15_gw_sector.py. Checklist item C-15 -> DONE.

## RHAC-009 · 2026-07-12 · Principle P is NOT derivable as a law (obstruction proven)
QUESTION (Scenario Q Part B / checklist B-4c / SUMMARY reviewer Q9): can P
  (Omega_m = R_H/r_p) be derived from A1-A3? Prior status: OPEN, billed as "the
  single largest available upgrade to Path One's claims."
RESULT: NO -- structural obstruction, not a cleverness gap. P's core equality
  (instantaneous Omega_m(a) = R_H(a)/d_p(a)) holds ONLY at a~=1; both quantities
  fall monotonically through cosmic history and CROSS ONCE, at today. High-z
  limit: Omega_m -> 1 while R_H/d_p -> 1/2 (EdS Schwarzschild-horizon value). A
  law derivable from time-symmetric axioms must hold at EVERY epoch; P does not;
  therefore P is not a theorem of A1-A3. The repo's own Route (ii) (g_horizon =
  cH0/2) IS this today-condition, re-expressed.
ESCAPE ROUTES (both closed): (a) ATTRACTOR -- fixing a preferred Omega_m ratio
  needs dynamical/coupled dark energy, FORBIDDEN by the Phase-2 null (RHAC-004:
  residual sector empty, w=-1 exactly); pursuing it breaks a locked result.
  (b) ANTHROPIC -- explains the O(1) coincidence in the weak sense but CANNOT
  reproduce P's 0.3141 precision (broad window), and is not ESTIF-specific (same
  move available to LCDM). It restates the coincidence, does not solve it.
STATUS: P is a PREDICTIVE POSTULATE (route c) -- honest, publishable as such,
  with the falsifiable number Omega_m = 0.3141 (0.53 sigma from Planck). The
  claim "Omega_m derived from the axioms" is RETIRED. Priority redirect: derive-P
  is CLOSED; the tractable foundational targets that actually underwrite a0's
  VALUE are x_c = 0.272 (pure Schwarzschild geometry) and the sqrt(3)/cH0
  acceleration scale.
Receipt: tests/estif_P_derivation_attempt.py. Downgrades Scenario Q Part B and
  Scenario H from "open upgrade" to "closed as a derivation."
```

Then bump the RHAC header and footer:

**FIND** (header):
```
**Last Updated:** 8 July 2026
**Status:** Active — updated to v6.3 "The Split". Gravity letter ready. Path One / Path Two fork recorded below.
```
**REPLACE:**
```
**Last Updated:** 12 July 2026
**Status:** Active — v6.4.0. Gravity letter ready. Fronts 1–3 + C-15 filed (RHAC-007/008); Principle P closed as non-derivable (RHAC-009).
```

**FIND** (footer):
```
**Document Version:** 6.4.0 (v6.3 "The Split" + errata + bootstrap + Phase 2 null + A1′ fork)
**Last Updated:** 11 July 2026
```
**REPLACE:**
```
**Document Version:** 6.4.0 (+ Fronts 1–3, C-15 GW sector derived, P non-derivable)
**Last Updated:** 12 July 2026
```

**FIND** (Scenario Q status):
```
**Status:** 🔶 Part A resolved (conditional); Part B open — the single largest
available upgrade to Path One's claims.
```
**REPLACE:**
```
**Status:** 🔶 Part A resolved (conditional); **Part B CLOSED 12 Jul 2026 (RHAC-009) — P is not derivable as a law (epoch-dependence obstruction). P is a predictive postulate; the "Ωm derived" path is retired.**
```

**FIND** (Scenario H re-update status):
```
**Status:** 🟡 Active, sharpened → now reduces to Part B of Scenario Q.
```
**REPLACE:**
```
**Status:** 🔴 CLOSED as a derivation (RHAC-009): it reduces to Scenario Q Part B, now proven non-derivable. Ωm = x₀ remains a consistency relation / predictive postulate, never a theorem.
```

> Manual follow-up (low value, skip if short on time): recount the two "Updated
> Summary Statistics" tables to reflect RHAC-007/008/009 and the two closures.

---

## B. STATUS.md

**FIND** (priority #3):
```
3. **Attempt Part B**: derive principle P from A1–A3. The single largest available upgrade to Path One's claims — it would convert Ωm, Ω_Λ, r_universe, and a₀ from consistency relations into zero-parameter predictions from three measured inputs.
```
**REPLACE:**
```
3. ~~Attempt Part B (derive principle P from A1–A3).~~ **CLOSED 12 July 2026 (RHAC-009):** deriving P as a law is impossible — P holds only at our epoch (Ωm(a) = R_H/d_p crosses once, at a≈1; diverges to 1 vs ½ at high z), so it cannot follow from the time-symmetric axioms. Escape routes closed: an attractor breaks the Phase-2 null; anthropics can't reproduce 0.3141. P is a predictive postulate; the "Ωm derived" claim is retired. Redirect foundational effort to x_c = 0.272 or the √3/cH₀ scale — both underwrite a₀'s value, unlike P.
```

**FIND** (end of the §3.0 bootstrap blockquote):
```
**Until Part B lands, the C2 downgrade above stands.** See RHAC Scenario Q and `VALIDATION_REPORT.md` §3.0.
```
**REPLACE:**
```
**Part B is now closed (RHAC-009): P is not derivable as a law, so the C2 downgrade is final — Ωm's status is "fixed by predictive postulate P", not "derived".** See RHAC-009 and `VALIDATION_REPORT.md` §3.0.
```

**FIND** (add C-15/GW to the Gravity Sector table — insert a new row after the "MOND a₀ derived" row; place it right before the "SPARC BTFR" row):
```
| SPARC BTFR (87 galaxies, Qual-1) | RMS = 15.6% | ✅ Within observed scatter |
```
**REPLACE:**
```
| **GW sector: c_gw = c derived (C-15)** | single-geometry; GW & light share the null cone (arbitrary flow); GW170817 passed structurally | ✅ NEW 12 Jul 2026 (RHAC-008) |
| SPARC BTFR (87 galaxies, Qual-1) | RMS = 15.6% | ✅ Within observed scatter |
```

**FIND** (footer):
```
**Status Document Version:** 6.3.2 / **Last Updated:** 9 July 2026.
```
**REPLACE:**
```
**Status Document Version:** 6.4.0 / **Last Updated:** 12 July 2026.
```

> Note: the "Priority Actions" list has a pre-existing numbering bug (two "3."s).
> Optional to fix while here; not required.

---

## C. PATH_ONE_CHECKLIST.md

**FIND** (status date):
```
**Status date:** 11 July 2026 (v6.4.0 — Phase 2 null + A1′ fork)
```
**REPLACE:**
```
**Status date:** 12 July 2026 (v6.4.0 — Fronts 1–3, C-15 closed, P non-derivable)
```

**FIND** (Two highest-value open derivations, item 2):
```
2. **Principle P** (B-4c). If derived, it converts the Ωm bootstrap from conditional
   to a genuine zero-parameter prediction of Ωm, Ω_Λ, r_universe, and a₀ from three
   measured inputs — the single largest available upgrade to Path One's claims.
```
**REPLACE:**
```
2. ~~Principle P (B-4c).~~ **CLOSED 12 July 2026 (RHAC-009):** P is not derivable
   as a law — it holds only at our epoch. Ωm's status is "fixed by predictive
   postulate", not "derived". The tractable foundational targets are now
   x_c = 0.272 (Schwarzschild geometry) and the √3/cH₀ acceleration scale.
```

**FIND** (C-15 row):
```
| C-15 | GW sector under A1′: derive the dent wave equation from the flow axioms (speed = c forced by A2; GW170817-consistent). Opened by RHAC-005/006 | ⬜ |
```
**REPLACE:**
```
| C-15 | GW sector under A1′: c_gw = c DERIVED — one Lorentzian geometry, GW & light share the null cone (principal symbol = light cone), shown for arbitrary flow; GW170817 passed structurally. RHAC-008 | ✅ 12 Jul 2026 |
```

**FIND** (B-4c row):
```
| B-4c | Part B: derive principle P from the flow axioms. Routes: (i) flow-budget amplitude; (ii) horizon-acceleration balance (P ⇔ g_horizon = cH₀/2 — the same cH₀ that sets a₀, ratio computes to 1.00000); (iii) homogeneous field equation (needs the vorticity attachment shared with Path Two). | ⬜ |
```
**REPLACE:**
```
| B-4c | Part B: derive P from the axioms — **CLOSED 12 Jul 2026 (RHAC-009).** Obstruction: P's equality Ωm(a)=R_H/d_p holds only at a≈1 (crosses once; →1 vs ½ at high z), so no time-symmetric-axiom derivation exists. Route (ii)'s g_horizon=cH₀/2 IS that today-condition. Escapes closed (attractor breaks Phase-2; anthropic can't hit 0.3141). P = predictive postulate. | ✅ closed |
```

**FIND** (A-6 row):
```
| A-6 | **The distinguishing prediction.** Path Two thawing route CLOSED 11 July 2026 (honorable null, RHAC-004). Live routes now: (i) mean spatial curvature ≡ 0 exactly at all epochs — registered kill-shot (current 0.0007 ± 0.0019); (ii) C-11 discriminator hunt (A1′ reproduces GR at linear order) | ⬜ | thawing (after the C4 re-derived target) | ⬜ |
```
**REPLACE:**
```
| A-6 | **The distinguishing prediction.** FOUR exactness locks catalogued (RHAC-007): Ωm/Ω_k=0 (kill-shot), γ≈0.55 & slip=1 (Front 3), w=−1 (RHAC-004), c_gw=c (RHAC-008). Honest finding: none separates ESTIF-Core from FLAT ΛCDM — all separate it from GR's extra freedoms. Genuinely-distinct content lives OFF the linear sheet (nonlinear halos, N-body wall) + the exactness-as-law structure. | 🔶 catalogued |
```

**FIND** (tally row):
```
| **17** | **1** | **15** | **2** |
```
**REPLACE:**
```
| **19** | **2** | **12** | **2** |
```
(Check: 19 + 2 + 12 + 2 = 35, matches "35 tracked items".)

---

## D. CONCEPT — ESTIF_CONCEPT.md (one edit)

**FIND** (Honest Open Questions → "Why Ωm = x₀?"):
```
- **Why Ωm = x₀?** The homogeneous version of the Task 4 calculation — whether
  `ρ_eddy = x₀ρ_crit` emerges — is the central dark-matter target, now well-posed.
```
**REPLACE:**
```
- **Why Ωm = x₀?** Reframed 12 July 2026 (RHAC-009): the bootstrap principle P
  (Ωm = R_H/r_p) that would fix Ωm is a **predictive postulate, not derivable**
  from the axioms — it holds only at our epoch. Ωm = x₀ stays a consistency
  relation. Whether ρ_eddy = x₀ρ_crit emerges from the homogeneous field equation
  is a separate, still-open question.
```

> Optional larger CONCEPT pass (a later decision, not this commit): add a short
> "Gravitational Waves" paragraph under the gravity section (c_gw = c, one
> geometry, GW170817) and a one-line four-lock summary. Flagged, not scripted.

---

## E. SUMMARY_FOR_REVIEW.md

**FIND** (reviewer Q9):
```
9. Is principle P (Ωm = R_H/r_p) derivable from the flow axioms A1–A3? It is
   equivalent to requiring that the mean matter pull at the particle horizon equal
   cH₀/2 — the same cH₀ that sets a₀ (ratio 1.00000). Is that equivalence a mechanism
   or a coincidence? And how does the resulting root (0.31408) relate to Gaztañaga's
   causal-boundary scale (≈ 0.3176 H₀), which is reached by a different route
   (inflation)? Is this the same result in different notation?
```
**REPLACE:**
```
9. ~~Is principle P derivable from A1–A3?~~ **Resolved internally 12 July 2026
   (RHAC-009): NO.** P holds only at our epoch (Ωm(a) = R_H/d_p crosses once, at
   a≈1; →1 vs ½ at high z), so it is not a theorem of the time-symmetric axioms —
   it is a predictive postulate. The cH₀/2 equivalence is that today-condition, not
   a mechanism. Remaining live question for reviewers: how does the postulated root
   (0.31408) relate to Gaztañaga's causal-boundary scale (≈ 0.3176 H₀, via
   inflation)? The comparison memo is still owed (B-4b).
```

**FIND** ("What is conditional" bootstrap bullet):
```
- The Ωm bootstrap: unique root 0.31408 (0.53σ from Planck), a₀ → 0.66%, r_universe
  back-predicted to −1.07%, DESI 1.618. **All of it rests on principle P, which is
  not derived.** If P falls, the bootstrap falls with it and §C2's downgrade is the
  final word. Part B is the test.
```
**REPLACE:**
```
- The Ωm bootstrap: unique root 0.31408 (0.53σ from Planck), a₀ → 0.66%, r_universe
  back-predicted to −1.07%, DESI 1.618. **All of it rests on principle P, which is
  not derived — and (RHAC-009, 12 Jul 2026) is not derivable as a law.** §C2's
  downgrade is therefore final: Ωm is fixed by a predictive postulate, not derived.
  P's falsifiable content (Ωm = 0.3141) stands as a postulate.
```

**FIND** (footer):
```
**Document Version:** 6.3.2 | **Updated:** 9 July 2026
```
**REPLACE:**
```
**Document Version:** 6.4.0 | **Updated:** 12 July 2026
```

> Optional: add the 7 session scripts to the "Files to Examine" list.

---

## F. VALIDATION_REPORT.md (predates A1′ — banner + P fix + footer)

**FIND** (the Status line at top):
```
**Status:** Gravity field equation DERIVED (not matched to Schwarzschild). Strong-field complete. MOND derived, SPARC validated. Cosmology reframed: frozen eddy = cosmological constant ties ΛCDM; Ω_tilt(z) retired. Project split into Path One (Core) and Path Two (Extended). Gravity letter ready and strengthened.
```
**REPLACE:**
```
**Status:** Gravity field equation DERIVED (not matched to Schwarzschild). Strong-field complete. MOND derived, SPARC validated. Cosmology reframed: frozen eddy = cosmological constant ties ΛCDM; Ω_tilt(z) retired. Project split into Path One (Core) and Path Two (Extended). Gravity letter ready and strengthened.

> ⚠️ **v6.4.0 banner (12 July 2026) — this report predates the A1′ amendment.** Throughout, "A1 (flat slices)" now reads **A1′ (even on average; matter dents locally; ⟨curvature⟩ ≡ 0 as law)** — RHAC-006. All strict-A1 results below survive as the exact-evenness limit. Added since: linear growth D₊ and the radiative sector (opened by A1′); **the GW sector c_gw = c is DERIVED** (RHAC-008); the four exactness locks (RHAC-007); and **Principle P proven NOT derivable as a law** (RHAC-009) — the §3.0 bootstrap's Part B is closed and its "Ωm derived" path is retired.
```

**FIND** (§3.0 honest flag #1):
```
1. **Everything above is conditional on P**, which is *not* derived from A1–A3. Part
   B — deriving P — is open. Three candidate routes exist (flow-budget amplitude;
   horizon-acceleration balance; homogeneous field equation); none has been attempted
   in earnest. Without Part B this is a reparametrization of the C2 circularity, not
   an escape from it.
```
**REPLACE:**
```
1. **Everything above is conditional on P**, which is *not* derived from A1–A3 —
   and as of 12 July 2026 (RHAC-009) **cannot be**: P holds only at our epoch
   (Ωm(a) = R_H/d_p crosses once, at a≈1), so no time-symmetric-axiom derivation
   exists. Part B is CLOSED. The bootstrap is therefore a reparametrization of the
   C2 circularity resting on a predictive postulate, not an escape from it.
```

**FIND** (footer):
```
**Validation Report Version:** 6.3.2 / 9 July 2026.
```
**REPLACE:**
```
**Validation Report Version:** 6.4.0 / 12 July 2026 (A1′ banner; P closed; GW sector added).
```

---

## G. ROADMAP.md (predates A1′ — banner + footer)

**FIND** (version block near top):
```
**Version:** 6.3.1
**Last Updated:** 9 July 2026
**Status:** Project split into Path One (ESTIF-Core) and Path Two (ESTIF-Extended). Gravity letter ready. Ω_tilt cosmology retired.
```
**REPLACE:**
```
**Version:** 6.4.0
**Last Updated:** 12 July 2026
**Status:** v6.4.0 — A1→A1′; Phase 2 null; Fronts 1–3; C-15 (GW sector) DERIVED; Principle P proven non-derivable as a law (RHAC-009). Gravity letter ready.

> ⚠️ **v6.4.0 (12 July 2026):** this roadmap body is the v6.2/v6.3 plan preserved as history and is largely superseded. Authoritative current state: `docs/plan/PATH_ONE_CHECKLIST.md` (item tracker) + `docs/plan/RHAC.md` (RHAC-001…009). Key correction since: "derive Principle P" — listed below and elsewhere as a top target — is **CLOSED** (P is a predictive postulate, not derivable from the axioms; RHAC-009). Redirect foundational effort to x_c = 0.272 or the √3/cH₀ scale.
```

**FIND** (footer):
```
**Roadmap Version:** 6.3 — "The Split"
**Last Updated:** 8 July 2026
```
**REPLACE:**
```
**Roadmap Version:** 6.4.0
**Last Updated:** 12 July 2026
```

---

## H. Phase 2 — needs a read pass before exact edits (do NOT guess)

These four weren't in the audit read; give me the go-ahead and I'll read them and
produce exact FIND/REPLACE the same way:

- **README.md** — apply the pending v6.3.2 bootstrap text (checklist P-4), add the
  C-15/GW result and the four-lock line, bump to v6.4.0.
- **CHANGELOG.md** — add a v6.4.0 entry (Fronts 1–3, C-15 derived, P non-derivable).
- **CITATION.cff** — bump version to 6.4.0 and **fix the DOI mismatch** (it points a
  6.3.2 at a 6.2 DOI; reconcile against the canonical Zenodo record 17261724).
- **CORRECTIONS_v6.3.1.md** — likely no change (frozen errata doc); confirm.

**Out of scope for this commit (separate tasks):** the lineage/NAMING.md tag still
owed by RHAC-006; the `docs/LaTeX ` trailing-space folder; the LaTeX/Pages letter
sources. Flag when you want those.

---

## Apply order (one pass, one commit)

1. Copy the 7 scripts into `tests/`.
2. RHAC.md: append RHAC-007/008/009, then the four header/footer/status edits (A).
3. STATUS.md (B) · PATH_ONE_CHECKLIST.md (C) · ESTIF_CONCEPT.md (D) ·
   SUMMARY_FOR_REVIEW.md (E) · VALIDATION_REPORT.md (F) · ROADMAP.md (G).
4. Single commit. Phase 2 (H) after a read pass, as a second commit.
