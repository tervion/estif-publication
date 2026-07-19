# ERRATA — v6.3.1 (adversarial-review corrections)

**Date:** 8 July 2026
**Applies to:** all v6.3 documents and scripts
**Purpose:** six corrections found in a fresh critical pass. None reverses a v6.3
decision (the split, the tilt retirement, the frozen-eddy reframe all stand).
Each hardens the position against peer review. Numbers below are verified by
`/tmp/verify.py`-style recomputation (marginalized BAO fit; σ8 scaling).

> **Target location in repo:** project root `ERRATA_v6.3.1.md`, and apply the
> per-file edits listed under each item.

---

## Correction 1 — Schwarzschild-vs-deviations contradiction (SEVERITY: HIGH)

**The problem.** Task 4 proves the ESTIF *vacuum* is **exactly** Schwarzschild
(full Einstein tensor = 0). But the docs simultaneously list, as ✅ predictions:
EHT M87\* shadow deviation (√β = 0.174) and a 491 μs LISA GW propagation delay.
**These cannot both be true as stated.** Exact Schwarzschild vacuum ⇒ the shadow
is the GR shadow and GW propagation is the GR (null, luminal) result — *zero*
deviation. A nonzero shadow/GW deviation can only arise from the **non-vacuum eddy
background** (a real stress-energy filling space), which is exactly the sector
ESTIF has **not** yet derived.

**The fix (wording).** Demote the EHT and LISA deviation claims from "✅ prediction"
to **"conditional on the (open) eddy-stress sector."** They are *not* vacuum
predictions. The Planck-Λ tie is a *scale* statement, not a propagation deviation,
and may stay — but relabel it as a calibration coincidence, not a test passed.

**Per-file edits:**
- `docs/report/VALIDATION_REPORT.md` — Part 1.1: add a header line to the EHT and
  LISA tables: *"CONDITIONAL — sourced by the non-vacuum eddy background (open
  sector); the derived vacuum is exact Schwarzschild with zero shadow/GW
  deviation."* Change their ✅ to 🔬 in the summary table.
- `docs/report/STATUS.md` — Gravity Sector table: append "(conditional on eddy-
  stress sector)" to the EHT and LISA rows; change ✅ to 🔬.
- `README.md` — "Status at a glance": same annotation on the EHT/LISA line.
- `docs/SUMMARY_FOR_REVIEW.md` — "Strong-Field Calibration" table: same annotation;
  and add to referee question set: *"Do the EHT/LISA deviations survive once the
  vacuum is required to be exactly Schwarzschild (Task 4)? They must be sourced by
  the eddy background, which is not yet derived."*
- `docs/report/ESTIF_CONCEPT.md` — "The Combined Formula (strong-field deviations)":
  add a sentence: *"These deviations are predicted by the tilt formula but are NOT
  yet reconciled with the derived exact-Schwarzschild vacuum; they require the eddy
  background to source them, which is open. Treat as conditional."*

**Unaffected:** the a₀/MOND weak-field chain (it does not use the shadow/GW
deviations). The gravity *letter* stands.

---

## Correction 2 — Ωm = x₀ is a consistency relation, not an independent prediction (SEVERITY: MEDIUM)

**The problem.** x₀ = (c/H₀)/r_universe, and r_universe = 4.4×10²⁶ m is the ΛCDM
particle horizon — an integral that **contains Ωm**. So "Ωm predicted geometrically
from x₀" is circular: Ωm appears on both sides once r_universe is unpacked. It is a
**consistency relation among Planck-calibrated quantities**, not a parameter-free
prediction, until ESTIF derives r_universe internally.

**The fix (wording).** Replace "Ωm predicted geometrically" / "predicted, not
fitted" with **"Ωm = x₀ is a consistency relation (0.12%); it becomes a genuine
prediction only once r_universe is derived within ESTIF rather than taken as the
ΛCDM horizon."** RHAC Scenario H had this right; restore that caution.

**Per-file edits:**
- `tests/estif_pathone_cosmology.py` — patched (see delivered file): the VERDICT no
  longer says "predicted not fitted"; it says "consistency relation, pending an
  internal derivation of r_universe."
- `README.md`, `docs/report/STATUS.md`, `docs/report/VALIDATION_REPORT.md`,
  `docs/SUMMARY_FOR_REVIEW.md`, `CITATION.cff`, `MILESTONE_v6.3_THE_SPLIT.md`,
  `docs/report/ESTIF_CONCEPT.md` — wherever "Ωm = x₀" is called a prediction,
  append "(consistency relation; r_universe still taken as the ΛCDM horizon)".

---

## Correction 3 — The JWST growth boost must be redshift-localized (σ8/S8 filter) (SEVERITY: HIGH for the JWST route)

**The problem.** The JWST spec says a ~13% growth enhancement at z≈9 relieves the
tension. But **if that enhancement persists to z=0, σ8 → 0.81 × 1.13 = 0.92**
(verified), against Planck σ8 = 0.81 and weak-lensing S8 ≈ 0.76–0.78. That is
grossly excluded (>5σ) **and in the wrong direction** — low-z lensing mildly
prefers *suppressed*, not enhanced, growth (the S8 tension).

**The consequence.** D_ESTIF(z) cannot be a constant boost. It must **switch on at
high z (strong at z≈9) and switch off by z ≲ 2** (D_ESTIF/D_LCDM → 1 today). This is
a much sharper, harder requirement than "boost growth 13%": the mechanism must be
transient. This does not kill the JWST route, but it is a hard filter the derived
D_ESTIF(z) must pass, and it strongly constrains the allowed modified-gravity form.

**The fix.** `docs/report/JWST_TEST_SPEC.md` and `tests/estif_jwst_growth_spec.py`
updated to add the σ8/S8 filter as a mandatory constraint (see delivered files).

---

## Correction 4 — The Path Two target (χ²/N ≈ 0.66, "4σ evolving DE") is a fixed-parameter artifact (SEVERITY: HIGH for Path Two strategy)

**The problem.** The AIC/BIC script and several docs cite: best-fit CPL χ²/N ≈ 0.66,
evolving DE preferred at ~4.1σ, frozen eddy "ties ΛCDM at 1.92." **All of these hold
H₀, Ωm, and rd fixed to Planck.** BAO measures DM/rd, DH/rd, DV/rd — degenerate with
rd — and is sensitive to Ωm. Fixing them inflates the apparent significance.

**Verified, under proper marginalization** (Ωm and the rd-scale free for every
model, BAO-alone):

| model | χ²/N (Planck-fixed) | χ²/N (marginalized) |
|---|---|---|
| ΛCDM | 1.92 | **0.79** (Ωm→0.297) |
| best-fit CPL | 0.66 | **0.43** (Ωm→0.387, w0=−0.18, wa=−2.72) |

Marginalized **Δχ²(ΛCDM − CPL) = 4.65 for 2 params → ~2.2σ**, not ~4.1σ. And the
marginalized CPL best-fit runs to a non-physical corner — BAO-alone does **not**
robustly prefer a sensible thawing; the strong, sensible evolving-DE signal comes
from **CMB + SNe combinations**, which ESTIF is not yet tested against.

**The consequences for strategy:**
1. The honest Path Two bar is **not** "reach χ²/N ≈ 0.66." Against marginalized ΛCDM
   (0.79), a derived thawing must clear a **~2.2σ**, not ~4σ, gap on BAO-alone — an
   *easier* bar, but also a *smaller* prize than advertised.
2. Any real Path Two claim must eventually be tested against **DESI + CMB + SNe with
   full marginalization** (CLASS/Cobaya-level), which is beyond the current tooling.
   BAO-alone, fixed-parameter χ² is not sufficient evidence either way.
3. The Path One conclusion is **unaffected**: frozen eddy and ΛCDM are the *same
   model*, so marginalization shifts both identically — the ΔAIC ≈ 0.6 tie is robust.
   Only the *absolute* numbers (1.92) are convention-dependent.

**The fix.** `tests/estif_pathone_aic_bic.py` updated to run the marginalized
comparison and report the ~2.2σ figure alongside the fixed-parameter one, with the
Path Two target restated honestly (see delivered file). Docs citing "4σ" / "0.66 as
the Path Two target" should reference the marginalized numbers and the CMB+SNe caveat:
- `MILESTONE_v6.3_THE_SPLIT.md`, `docs/report/STATUS.md`,
  `docs/guide/ROADMAP.md`, `docs/report/VALIDATION_REPORT.md`,
  `docs/SUMMARY_FOR_REVIEW.md` — replace "target χ²/N ≈ 0.66" with "target: clear the
  marginalized ~2.2σ BAO gap (and, eventually, survive DESI+CMB+SNe with full
  marginalization)."

---

## Correction 5 — Task 6 E1 derivation is sloppy (SEVERITY: LOW — conclusion robust)

**The problem.** `estif_task6_eddy_eos.py`, model E1, wavers between ρ ∝ a⁻⁵ and a⁻⁶
(w = 2/3 vs +1) and the delivered file contains mid-thought "Wait —" prose. The
*conclusion* (a conserved-spin eddy blueshifts and fails as dark energy) is robust
for any w > 0, but the text is not publication-clean.

**The fix.** `tests/estif_task6_eddy_eos.py` E1 comment block rewritten: state the
robust claim ("conserved specific angular momentum ⇒ w > 0 ⇒ blueshifts ⇒ cannot be
dark energy; the exact index between +2/3 and +1 does not matter for the verdict"),
remove the stream-of-consciousness text (see delivered file).

---

## Correction 6 — "Derived, not borrowed" needs precise wording (SEVERITY: LOW — framing)

**The problem.** The headline "gravity derived, not borrowed" is defensible but
loose. What Task 4 rigorously shows: the flow axioms uniquely select the
**Painlevé–Gullstrand-gauge form of GR's constraint (energy) sector**, and the
vacuum condition forces exact Schwarzschild. What is **imported**, not derived: the
coupling itself — that a geometric scalar equals 8πG × (energy density). ESTIF
assumes the Einstein–Hilbert coupling; it does not derive G or the field-equation
form from a deeper principle.

**The fix (wording).** In headline claims, prefer: **"the flow axioms force the
vacuum solution and the mass-continuity (Poisson) constraint without a Schwarzschild
match; the gravitational coupling 8πG is assumed (Einstein–Hilbert), not derived."**
The SUMMARY referee-question #1 already anticipates this; make the headline language
match its caution across README, STATUS, MILESTONE, CONCEPT.

---

## Net effect

| Item | Before | After (verified) |
|---|---|---|
| EHT/LISA deviations | ✅ predictions | 🔬 conditional on open eddy-stress sector |
| Ωm = x₀ | "predicted" | consistency relation (r_universe still ΛCDM) |
| JWST growth boost | "13% relieves tension" | 13% **and must vanish by z≲2** (σ8 filter) |
| Evolving-DE preference | ~4.1σ, target 0.66 | **~2.2σ marginalized**; needs CMB+SNe |
| E1 scaling | a⁻⁵/a⁻⁶ waver | w > 0, verdict robust; text cleaned |
| "Derived gravity" | headline | "vacuum + constraint derived; 8πG assumed" |

**Decisions unchanged:** the split, the tilt retirement, the frozen-eddy = ΛCDM
result (the ΔAIC tie is marginalization-robust), and the JWST-is-gravity-not-
expansion finding. These corrections make the claims survivable, not weaker.

**Errata version:** 1.0 · **Date:** 8 July 2026
