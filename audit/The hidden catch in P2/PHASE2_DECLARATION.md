# PHASE 2 DECLARATION — The Residual Dark-Energy Sector Is Empty

**Date:** 11 July 2026
**Version:** 6.4.0
**Status:** CLOSED — honorable null (strong form). w = −1 exactly; Λ is a bare imported constant.
**Record of decision:** RHAC-004 (execution), RHAC-005 (strict-A1 growth no-go), RHAC-006 (A1 → A1′ fork).
**Receipts:** `tests/test_UKN.py` (strict A1), `tests/test_UKN2.py` (A1′ re-audit).

> ⚠️ **Provenance note.** This document records the *result and the method* of Phase 2 as
> logged in RHAC-004/005/006. The underlying constraint-algebra derivations (the
> pure-divergence theorem for the slosh sector, the ⟨δρ⟩ = −⟨ω²⟩/32πG slaving of the swirl
> sector) live in the test scripts. Anyone extending this should re-run those scripts and
> confirm the sector verdicts against them before citing the numbers downstream.

---

## 1. What Phase 2 asked

ESTIF-Extended (Path Two) was premised on a hope: that the cosmic flow, beyond the frozen
`w = −1` limit that Path One already ties ΛCDM with, hides a *small residual stress* — a
near-frozen thawing of the dark-energy equation of state — that would match the mild
evolution DESI DR2 BAO hints at (nominal target χ²/N ≈ 0.66, itself under-marginalized; C4).

Two naive reductions had already been tried and falsified against DESI (Task 6):

| Reduction | Derived w | χ²/N vs DESI DR2 | Verdict |
|---|---|---|---|
| E1 — conserved-angular-momentum spin | +1 (stiff) | 3232 | falsified |
| E2 — expansion tracker | thaws to ~0 | 754 | falsified |

Phase 2 was the honest, complete version of that question: **not** another guessed reduction,
but the *entire* residual stress the axioms permit, computed and submitted without cherry-picking.

---

## 2. The binding declaration (pre-registered, before the derivation)

To remove any suspicion of post-hoc sector selection, the following was fixed **before** the
constraint algebra was run:

> **Binding declaration.** The total residual stress admitted by axioms A1–A3 is to be
> computed in full and submitted **unmodified**. No sector may be added, dropped, or
> reinterpreted after seeing the result. Whatever the constraint algebra returns is the verdict,
> null or otherwise.

This is the methodological core of Phase 2. The result below is what that pre-registered
procedure returned.

---

## 3. The census — the four residual sectors and their verdicts

Under strict A1 (exactly even slices) + A2 (universal speed c) + A3 (vacuum sources nothing),
the residual stress decomposes into four candidate channels. Each was evaluated:

| Sector | Physical content | Verdict | Why |
|---|---|---|---|
| **Shape** | Anisotropic spatial curvature of the slice | **FORBIDDEN** | Excluded by A1 (even slices) |
| **Rate** | Time-varying lapse carrying stress | **FORBIDDEN** | A2: the lapse carries no stress; motion is at c |
| **Slosh** | Bulk spatial flow / momentum density | **ZERO** | Pure-divergence theorem — integrates to nothing |
| **Swirl** | Vorticity / rotation of the flow | **FORBIDDEN free; NEGATIVE if sourced** | ⟨δρ⟩ = −⟨ω²⟩/32πG, slaved to matter, and ~10¹¹× below ρ_Λ |

**No channel yields a free, positive, evolving dark-energy stress.** Shape and rate are
forbidden by the axioms; slosh vanishes identically; swirl, the one channel that survives, is
not free — it is slaved to matter, it carries the *wrong sign* (negative energy contribution),
and its magnitude is eleven orders below the observed Λ.

---

## 4. Verdict — honorable null, strong form

> **The residual dark-energy sector is EMPTY.**
> `w = −1` exactly. Λ enters the theory as a **bare imported constant**, exactly as in ΛCDM.
> The frozen-limit χ²/N = 1.92 parity with ΛCDM (Path One) is unaffected.
> **The Path Two ambition of a derived thawing correction (χ²/N ≈ 0.66) is dead.**

This is a *strong-form* null, not a "not-yet": it is not that a thawing term was sought and not
found, but that the axioms **provably permit none**. That is a more valuable result than a weak
null — it converts "we haven't derived dark energy yet" into "our axioms cannot manufacture
dark energy, by construction." ESTIF does **not** eliminate dark energy; it imports Λ, and it
is now *proven* that it must.

---

## 5. Retraction — the "frozen eddy" label

The interpretation that the constant-Λ limit corresponds to a *frozen eddy* (a spinning cosmic
flow held constant) is **retracted**. It was auditor-introduced; the author objected before the
derivation ("whatever motion it is, it is most certainly not spinning"); and the constraint
algebra falsified it directly — the swirl (spin) sector is precisely the one that is forbidden
free and negative when sourced. **The constant survives; the "eddy"/"spin" reading of it does
not.** Documents and scripts should refer to the cosmological term as an imported cosmological
constant, not a frozen eddy. (Full rename batch pending `NAMING.md`.)

---

## 6. The LIGO flag, and its resolution under A1′

Phase 2 (strict A1) surfaced a genuine problem, flagged honestly: **the strict flow sector
carries no free radiative modes** — no propagating gravitational waves — which is in tension
with LIGO/Virgo's direct detections. This is the same over-constraint that RHAC-005 proved also
forbids the growing density mode (galaxies could not form).

**Resolution (RHAC-006): A1 → A1′.** The axiom was amended: the slice is even *on average*,
with matter permitted to source local dents; the spatial average of the curvature is identically
zero at all epochs (this is law and is the registered falsifier — current 0.0007 ± 0.0019).
A1′ restores (one purchase, three items): linear growth (D₊ = H·∫da/(aH)³ = GR growth,
f(z=0.5) = 0.76, DESI-consistent), the radiative/GW sector (speed = c forced by A2,
GW170817-consistent), and a container for the ζ = 10⁻⁵ seed.

**Phase 2 re-audit under A1′ (RHAC-006):** the five-sector census was re-run under the amended
axiom. No Λ-printer appears; **the null stands.** A1′ adds the growing mode and the radiative
sector without opening a residual dark-energy channel. All strict-A1 results survive as the
exact-evenness limit.

---

## 7. Consequences for the project

- **Path One (ESTIF-Core) is unaffected and, if anything, reinforced:** its cosmology is a
  cosmological constant (χ²/N = 1.92, ties ΛCDM), and it is now proven this is the *honest
  best* the axioms permit — not a placeholder awaiting a better derivation.
- **Path Two (ESTIF-Extended) is closed** as a derivation target. Any future evolving-w must
  come from *outside* the A1–A3/A1′ flow sector (a second metric function or explicit time
  dependence, leaving the single-speed-flow class; note L1: a single-speed flow forces
  p_r = −ρ).
- **The far-future / vacuum-product speculation** (RHAC-003) is now conditional on exactly this
  result: w = −1 forever and no new vacuum physics. The mainstream cousin (de Sitter horizon
  temperature, Gibbons–Hawking ~10⁻³⁰ K) is noted; ESTIF Path One is itself already a
  vacuum-flow-product claim (the constant Λ).

---

## 8. One-line summary

*Phase 2 asked whether the flow hides dark energy. Pre-registered and computed in full, the
answer is a proven no: the residual sector is empty, w = −1 exactly, Λ is imported. ESTIF does
not derive dark energy — and now it is proven that it cannot, from these axioms.*

---

**Declaration version:** 1.0 · 11 July 2026 · Records RHAC-004/005/006. Verify sector algebra
against `tests/test_UKN.py` and `tests/test_UKN2.py` before downstream citation.
