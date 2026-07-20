# CORRECTIONS — v6.3.1 (errata to the v6.3 documentation)

**Date:** 8 July 2026
**Purpose:** Six corrections identified in an adversarial review pass of the v6.3
documents. None reverses a v6.3 decision (the split, the tilt retirement, the
de-circularization, the frozen-eddy reframe all stand). These harden the position
against peer review. Apply before submitting the gravity letter.

> Priority: **C1, C2, C6 are pre-submission critical.** C3, C4, C5 are honesty/
> future-work fixes.

---

## C1 — Exact-Schwarzschild vs strong-field deviations (CRITICAL)

**The problem.** Task 4 proves the ESTIF **vacuum** is *exactly* Schwarzschild
(full Einstein tensor = 0). But the validation documents still list the EHT M87\*
shadow deviation and the 491 μs LISA GW delay as ✅ predictions. **These cannot
coexist as written:** an exact-Schwarzschild vacuum produces *zero* deviation in
photon-sphere shadows and in gravitational-wave propagation through vacuum. A sharp
referee catches this immediately.

**The resolution.** The strong-field deviations can only exist if they are sourced
by the **non-vacuum eddy background** (the tilt/eddy stress that is NOT yet
derived). So every strong-field *deviation* claim must be demoted from a firm
prediction to **conditional on the open eddy-stress sector**. The weak-field a₀/MOND
chain is unaffected (it never invokes the vacuum-deviation claim). The EHT shadow
*consistency* (ESTIF is consistent with the observed 42 μas) is fine; it is the
claim of a *predicted deviation from GR* that must be qualified.

**Apply everywhere the tilt formula's EHT/Λ/LISA results appear as ✅ predictions.**
Global rule: change the status from ✅ to ⚠️ and attach this footnote:

> ⚠️ Conditional. The ESTIF vacuum is exactly Schwarzschild (Task 4), so any
> deviation from GR in shadows or GW propagation must be sourced by the non-vacuum
> eddy background — a sector not yet derived. These strong-field deviation figures
> are therefore predictions *conditional on* the open eddy-stress derivation, not
> established results. The observations remain *consistent* with ESTIF; the
> *deviation from GR* is what awaits derivation.

Specific instances:
- `docs/report/VALIDATION_REPORT.md` §1.1 (EHT/Λ/LISA table), §1.4 removed already,
  §1.5 (GW mass table), Summary table rows "Strong-field formula", "Λ drift".
- `docs/SUMMARY_FOR_REVIEW.md` "Strong-Field Calibration" table.
- `docs/report/ESTIF_CONCEPT.md` "The Combined Formula (strong-field deviations)"
  table — add the footnote to the EHT/Λ/LISA row.
- `docs/report/STATUS.md` gravity-sector rows "EHT M87* shadow", "LISA GW delay".
- `README.md` "EHT / Λ / LISA (local tilt)" status row.
- `MILESTONE_v6.3_THE_SPLIT.md` — the line "the a₀/MOND derivation sits on a derived
  field equation" is fine; add a sentence that the strong-field *deviation* claims
  are conditional on the eddy-stress sector.

**Note on Planck Λ:** the "Λ = 1.0000 ratio" is a calibration match, not a vacuum
deviation, so it is not affected by C1 — but see C6 on calling it "derived."

---

## C2 — Ωm = x₀ is a consistency relation, not an independent prediction

**The problem.** x₀ = (c/H₀)/r_universe, and r_universe = 4.4×10²⁶ m is the ΛCDM
particle horizon — an integral that **itself contains Ωm**. So "Ωm predicted
geometrically as x₀" is circular at the level that matters: it is a *consistency
relation among Planck-calibrated quantities*, not a prediction independent of Ωm.
RHAC Scenario H had the right suspicion; the v6.3 docs drifted past it.

**The resolution.** Downgrade the language from "predicted/derived" to "consistency
relation," everywhere it appears, until ESTIF derives r_universe internally (without
an Ωm-dependent horizon integral). Replacement text:

> Ωm = x₀ = (c/H₀)/r_universe holds to 0.12%. This is a **consistency relation**:
> r_universe is here the ΛCDM particle horizon, which itself depends on Ωm, so the
> agreement is a self-consistency of the geometric picture with Planck values, not
> an Ωm-independent prediction. An independent prediction requires deriving
> r_universe from the flow framework without the Ωm-dependent horizon integral —
> an open task (RHAC Scenario H).

Specific instances:
- `tests/estif_pathone_cosmology.py`: in the docstring and the VERDICT print, change
  "Omega_m predicted geometrically as x0 (0.12%), not fitted" →
  "Omega_m = x0 is a CONSISTENCY RELATION (r_universe is the LCDM horizon, which
  depends on Omega_m); not an Omega_m-independent prediction." Keep the χ² result;
  only the interpretive claim changes.
- `README.md`, `docs/report/STATUS.md`, `docs/report/ESTIF_CONCEPT.md`,
  `docs/report/VALIDATION_REPORT.md §3.1`, `docs/SUMMARY_FOR_REVIEW.md`: replace
  "predicted/geometric Omega_m" phrasing with the consistency-relation text above.
- `CITATION.cff` abstract: "Omega_m = x0 predicted" → "Omega_m = x0 as a consistency
  relation."

**This weakens one marketing claim but removes a circularity a referee would
exploit. The 0.12% numerical agreement is real and stays; only its epistemic status
is corrected.**

---

## C3 — The JWST spec is missing the σ8/S8 filter (HARD constraint)

**The problem.** `JWST_TEST_SPEC.md` states that a ~13% growth enhancement at z≈9
relieves the JWST tension, but omits the low-redshift structure constraint. A
scale-independent 13% growth boost that **persists to z=0** raises σ8 from 0.811 to
**0.917 — ~18σ above Planck**, and in the *wrong* direction relative to weak-lensing
S8 (which mildly prefers *suppressed* low-z growth). So the naive "boost growth 13%"
is grossly excluded by present-day data.

**The resolution.** The D_ESTIF(z) enhancement must be **transient (concentrated at
high z, decaying to ≈1 by z ≲ 2)** and/or **scale-dependent (enhancing small-scale/
high-k halo formation while leaving the 8 Mpc/h scale that sets σ8 essentially
untouched)**. This sharpens the derivation target substantially. Add this section to
the JWST spec and script (both re-issued in v6.3.1 with it included):

> **HARD FILTER — σ8/S8.** The growth enhancement must NOT persist to z=0 at the
> 8 Mpc/h scale. A persistent scale-independent 13% boost gives σ8 ≈ 0.92 (~18σ
> above Planck; wrong sign for lensing). Therefore D_ESTIF(z) must be either
> transient (high-z only, → 1 by z ≲ 2) or scale-dependent (small-scale only). This
> is a two-sided constraint: enough early growth for JWST, ≈ standard late growth
> for σ8. It makes the D_ESTIF derivation target more specific, not merely "larger."

*(The re-issued `estif_jwst_growth_spec.py` now computes and prints the σ8 exclusion
and states the two-sided requirement.)*

---

## C4 — The Path Two target (χ²/N ≈ 0.66) needs re-derivation with marginalization

**The problem.** The session's BAO-alone CPL fit found ~4σ evolving-DE preference
and χ²/N ≈ 0.66 for the best CPL, **holding rd, H₀, and Ωm fixed at Planck values.**
DESI's *own* BAO-alone preference for evolving DE is much milder; the strong
combined-data significance comes from adding CMB and supernovae. Fixing nuisance
parameters **inflates** Δχ². So 0.66 is an optimistic, under-marginalized target.

**The resolution.** Before committing to Path Two, re-derive the target with rd, H₀,
Ωm **marginalized** (or at least profiled) over their DESI/Planck uncertainties. The
true bar a derived w(z) must clear may be **easier** (if fixed-parameter fitting
inflated the CPL advantage) **or** the payoff **smaller** (if marginalization shrinks
the evolving-DE preference). Add this caveat wherever 0.66 is cited as the Path Two
target:

> ⚠️ The χ²/N ≈ 0.66 CPL target is from a BAO-alone fit with rd, H₀, Ωm fixed at
> Planck. Fixing nuisance parameters inflates the evolving-DE advantage; DESI's own
> BAO-alone evolving-DE preference is milder. Re-derive this target with those
> parameters marginalized before using it to justify Path Two. The AIC comparison
> (frozen-eddy vs fitted CPL) should be repeated under marginalization.

Specific instances: `MILESTONE_v6.3_THE_SPLIT.md`, `docs/report/STATUS.md`,
`docs/report/VALIDATION_REPORT.md §2.2/§2.3`, `docs/SUMMARY_FOR_REVIEW.md`,
`tests/estif_pathone_aic_bic.py` (the "Path Two target" block), and
`JWST_TEST_SPEC.md` if it cites 0.66.

---

## C5 — Task 6 E1 derivation is sloppy; clean the file

**The problem.** In `tests/estif_task6_eddy_eos.py`, the E1 (conserved-angular-
momentum spin) block wavers between ρ ∝ a⁻⁵ and a⁻⁶ (w = 2/3 vs w = +1) and the
delivered file contains a mid-thought "Wait --" and an "a^?" placeholder. The
*conclusion* (the spin component blueshifts and fails as dark energy) is robust for
any exponent ≥ 5 (all give w > 0, all excluded), but the text is not publication-
clean.

**The resolution.** Replace the E1 explanatory `print(...)` block with:

```python
print("""  Physical picture: the cosmic eddy is bulk rotation. A rotating patch of
  comoving size ~a has moment of inertia I ~ M a^2; conserving angular momentum
  L = I*omega per comoving patch gives omega ~ a^-2. The rotational energy density
  rho_rot ~ I omega^2 / a^3 then scales as a steep NEGATIVE power of a -- a
  "stiff"/blueshifting component (w > 0). The exact exponent depends on how the
  comoving mass and volume factors are booked, but for ANY such conserved-L
  reduction the exponent is >= 5, i.e. w >= 2/3: the component GROWS toward the
  past and is negligible today. It behaves as extra early matter/stiff fluid,
  NOT as dark energy. We take the representative stiff case w = +1 below.""")
```

The numerical result (χ²/N = 3232, falsified) is unchanged.

---

## C6 — "Derived, not borrowed" must be worded precisely

**The problem.** The headline "gravity is derived, not borrowed" is *mostly* right
but overclaims if unqualified. What Task 4 rigorously establishes: the flow axioms
(flat slices + universal-c + vacuum-sources-nothing) **uniquely select the
Painlevé–Gullstrand-gauge form of the constraint (energy) sector of GR**, forcing
mass continuity → exact Schwarzschild in vacuum. What is **imported**, not derived,
is the **coupling constant and form** — the identification of the geometric
constraint scalar with 8πG × (energy density). ESTIF does not derive Newton's G or
the factor 8π from below; it adopts the Einstein–Hilbert coupling.

**The resolution.** Use this precise wording in the paper and headlines:

> ESTIF *derives* that its three flow axioms uniquely select the constraint (energy)
> sector of General Relativity in Painlevé–Gullstrand gauge — forcing mass
> continuity (Poisson in integrated form), hence exact Schwarzschild in vacuum and
> the Newtonian source in the weak field — **without matching to the Schwarzschild
> solution**. The gravitational *coupling* (geometric scalar ↔ 8πG × energy
> density) is adopted, not derived from below. So the correct claim is: the force
> law is no longer *matched to GR's vacuum solution* (the previous gap), but the
> *coupling to matter* is still the standard Einstein–Hilbert one.

This matches the caution already present in `SUMMARY_FOR_REVIEW.md` reviewer
question #1. Update the stronger "derived, not borrowed" phrasings in
`MILESTONE_v6.3_THE_SPLIT.md`, `README.md`, `docs/report/ESTIF_CONCEPT.md`
("Gravity Is Now Derived, Not Borrowed" — keep the title but add this precision
paragraph), and `docs/report/VALIDATION_REPORT.md §0` to carry this qualifier.

---

## Summary of what changed

| # | Correction | Severity | Nature |
|---|---|---|---|
| C1 | Exact-Schwarzschild vs EHT/LISA deviations | **Critical** | deviation claims → conditional on eddy-stress |
| C2 | Ωm = x₀ is a consistency relation | **Critical** | remove hidden Ωm-circularity |
| C3 | JWST spec missing σ8/S8 filter | High | growth boost must be transient/scale-dependent |
| C4 | Path Two target under-marginalized | High | re-derive 0.66 with rd/H₀/Ωm marginalized |
| C5 | Task 6 E1 text sloppy | Low | file cleanup, result unchanged |
| C6 | "Derived, not borrowed" wording | **Critical** | coupling is adopted, not derived |

**None of these reverses a v6.3 decision.** The split, the tilt retirement, the
de-circularization (10.8 → 3.35), the frozen-eddy reframe (ties ΛCDM), and the
derived constraint sector all stand. C1, C2, and C6 must be applied before the
gravity letter is submitted; C3 and C4 before Path Two or the JWST test is pursued;
C5 whenever convenient.

**Version:** 6.3.1 · **Date:** 8 July 2026
