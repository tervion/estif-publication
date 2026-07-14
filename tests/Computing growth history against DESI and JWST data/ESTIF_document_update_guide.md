# ESTIF — Document Update Guide (a₀ = horizon quantity)

**Date:** 2026-07-13
**Purpose:** This is the Path-1 consolidation. One session established a single consistent
result about a₀; this file tells you exactly what to change in your existing documents so
they all tell that one story. Format: **FIND** (what to search for) → **WHY** → **REPLACE
WITH** (paste-ready text).

> Do the replacements in order. Nothing here touches your empirical fits or your black-hole
> exterior solution — only the *interpretation* of a₀ and any claim that it is *derived*.

---

## 0. The one-line result (the thing being stamped everywhere)

> **a₀ is a horizon quantity: a₀ ≈ c·H (the de Sitter surface-gravity scale).**
> Mass-independence of a₀, flat rotation curves, the Baryonic Tully–Fisher relation
> (v⁴ ∝ M), and the observed magnitude (~1.2×10⁻¹⁰ m/s²) all *follow* from this.
> **Still open:** the exact O(1) prefactor (cH vs cH/2π vs the value the data prefers,
> k ≈ 0.128 in a₀ = k·c²√Λ). It is constrained to O(1) and known to be horizon-set,
> but the precise number is **not** derived.

---

## 1. Master status block — paste at the TOP of the main paper / notebook page 1

**COPY THIS verbatim into the top of your primary ESTIF document, replacing any existing
"summary of a₀" or abstract-level claim about where a₀ comes from:**

```
STATUS OF a0 (current):
  DERIVED / SETTLED:
    - a0 is a horizon-scale acceleration: a0 ≈ c·H  (de Sitter surface gravity).
    - From this ALONE follow: flat rotation curves, mass-independent a0,
      the Baryonic Tully-Fisher relation v^4 = G·M·a0, and the correct
      magnitude (cH/2pi = 1.0e-10 vs observed 1.2e-10 m/s^2).
    - Equivalence: c^2·sqrt(Lambda) = sqrt(3)·sqrt(Omega_Lambda)·c·H0, so any
      "Lambda-native" form for a0 is algebraically the Hubble-scale form up to
      an O(1) factor. Matter (Omega_m) is absorbed into H0 via Friedmann, not
      eliminated.
  OPEN / NOT DERIVED:
    - The O(1) prefactor. The 87-galaxy fit prefers k ≈ 0.128 in a0 = k·c^2·sqrt(Lambda);
      this is NOT 1/(2pi)=0.159, sqrt(3)/(2pi), or any clean constant yet identified.
    - The mechanism by which the horizon scale c·H acts LOCALLY (see §6).
```

---

## 2. Replace: "a₀ derived from matter density / the 31% coincidence (P)"

**FIND** — any passage that:
- derives a₀ from the matter density Ω_m ≈ 0.31, or the "matter budget," or
- treats **P** (the Ω_m ≈ size-ratio ≈ 0.31 coincidence) as load-bearing / a premise, or
- says a₀ "comes from the amount of matter."

**WHY** — Writing `c²√Λ = √3·√Ω_Λ·c·H₀` shows the "Λ-native" and Hubble forms are the *same
expression*; the matter dependence hides inside H₀, it is not removed. So a₀ is anchored to
c·H, and P is a downstream curiosity, not a foundation. (This is the P-2 rewrite, now confirmed.)

**REPLACE WITH:**

```
a0 does not derive from the matter density. Using the flat-universe Friedmann relation,
c^2·sqrt(Lambda) = sqrt(3)·sqrt(Omega_Lambda)·c·H0, so any expression for a0 written in terms
of Lambda is algebraically identical, up to an O(1) factor, to the cosmic-horizon form c·H0.
The apparent matter dependence is carried by H0 and is not eliminated by switching horizons.
a0 is therefore a horizon-scale acceleration, and the numerical coincidence
Omega_m ≈ (size ratio) ≈ 0.31 (previously labelled "P") is a downstream observation, not a premise.
```

---

## 3. Replace: "the √3 factor / the a₀ prefactor is derived"

**FIND** — any passage presenting the √3 "3-D space" factor, or the full numerical prefactor of
a₀, as derived from first principles.

**WHY** — It is not derived. The bare thermodynamic value 1/(2π)=0.159 overshoots the target
band; the fit prefers k≈0.128, which matches no clean constant yet; the river-model edge lands
~20% high, and that 20% is exactly the ambiguity between acceleration definitions.

**REPLACE WITH:**

```
The O(1) prefactor of a0 (previously justified by a sqrt(3) "3-D space" argument) is NOT yet
derived. Findings:
  - The horizon-thermodynamic value 1/(2pi) = 0.159 overshoots the target band [0.118, 0.133].
  - The 87-galaxy fit prefers k ≈ 0.128 in a0 = k·c^2·sqrt(Lambda); this equals neither
    1/(2pi), sqrt(3)/(2pi) = 0.276, nor any clean geometric constant identified so far.
  - A river-model "edge" calculation (acceleration at the surface where inflow = x_c·c) gives
    0.26–0.28·c·H_Lambda, ~20% high — and that ~20% is precisely the difference between the
    candidate acceleration definitions (proper acceleration vs redshifted surface gravity vs
    flow gradient).
Deriving this prefactor from horizon geometry with zero fit freedom is the central open problem.
The SCALE of a0 (c·H) is derived; the NUMBER is not.
```

---

## 4. Replace: "the inward flow produces a₀ from local galaxy dynamics"

**FIND** — any claim that the inflow/river dynamics yields a₀ from local (galaxy-scale) physics,
or that a₀ is a local property of the flow near a galaxy.

**WHY** — Direct simulation (both combination rules) proves the local flow makes (H/2)·v_gal,
which is mass-dependent (∝ M^{1/3}) and ~1000× too small. a₀ ≈ c·H needs the speed of light,
which enters the flow only at the cosmic horizon.

**REPLACE WITH:**

```
a0 is NON-LOCAL. A direct simulation of the inward-flow river model, testing BOTH combination
rules (velocities-add, then a=(v.grad)v; and accelerations-add), with a local de Sitter
background (v_bg = H·r), produces an extra inward acceleration of order (H/2)·v_gal — H times
the galaxy's own inflow speed (tens of km/s). This is:
  - the correct MOND-like DIRECTION (inward), but only ~1% of Newtonian pull inside galaxy radii,
  - MASS-DEPENDENT: a0 ∝ M^{1/3} (observed a0 is mass-independent), and
  - ~1000x too small in absolute scale.
Root cause: a0 ≈ c·H requires the speed of light c, which appears in the flow only at the cosmic
horizon (r = c/H), not at galaxy radii. a0 therefore cannot be a purely local galaxy calculation.
[Reproduce: estif_flow_sim.py]
```

---

## 5. ADD (new result — currently has no home in any document)

**Insert as a new subsection** (e.g., "a₀ as the de Sitter horizon acceleration") in the main
paper and/or the gravity letter.

**COPY THIS:**

```
HORIZON-REFERENCED RESULT.
Replacing the local background (velocity H·r) with the horizon's own acceleration scale
a_H = k·c·H (de Sitter surface gravity) — the one scale the local flow provably cannot generate
— reproduces the full observed phenomenology from that single input:
  - flat rotation curves;
  - the Baryonic Tully-Fisher relation v^4 = G·M·a_H, with correct normalization
    (~180 km/s at a Milky-Way baryonic mass);
  - mass-independent a0 (∝ M^0), in contrast to the local model's M^{1/3};
  - correct magnitude: cH/2pi = 1.0e-10 vs observed 1.2e-10 m/s^2.
Because the local flow cannot manufacture c·H (it makes H·v_gal instead), and injecting ONLY c·H
turns on all of the above, this is direct evidence that a0 IS the de Sitter horizon acceleration
— a cosmic/boundary quantity. The prefactor k remains undetermined (see §3).
[Reproduce: estif_horizon.py]
```

---

## 6. FLAG as open — the x_c = 0.272 double-duty question

**FIND** — anywhere x_c = 0.272 is used simultaneously as (a) the Einstein-mode → galaxy-mode
crossover AND (b) the scale that sets a₀, without noting they might differ.

**WHY** — Under the simplest edge law a₀ = c·H_Λ·x_c, the galaxy data prefer x_c ≈ 0.221, not
0.272. These may be two different surfaces that have been conflated.

**REPLACE WITH (add this note at the point of first use of x_c):**

```
OPEN: the mode-crossover x_c = 0.272 (Einstein-mode -> galaxy-mode) and the a0-defining edge
may be DIFFERENT surfaces. Under a0 = c·H_Lambda·x_c the data prefer x_c ≈ 0.221. Resolve whether
these are one surface (requiring the flow law to supply the reconciling ~0.8 factor) or two
distinct scales that have been conflated.
```

---

## 7. Do NOT touch (these survived — leave as written)

- **Phase-2 result: Λ is a dead constant** (no dynamical dark energy). Not only intact — it is
  exactly what the horizon result *needs*: a constant Λ ⇒ fixed de Sitter horizon ⇒ constant a₀.
  If anything, add a one-line note that Phase-2 now *supports* the a₀ = cH picture.
- **The 87-galaxy rotation-curve fit itself.** The empirical fit stands; only its *interpretation*
  (a₀ = horizon, not matter) is updated.
- **The black-hole exterior / vacuum solution.** Untouched by any of this.

---

## 8. Open-problems register (paste as an appendix; keep it visible)

```
OPEN PROBLEMS (a0 / gravity):
  1. The a0 prefactor — the O(1) number. CENTRAL. Needs derivation from horizon geometry
     with zero fit freedom (the former "sqrt(3)" slot).
  2. x_c = 0.272 derivation from black-hole-edge physics (clean strong-field problem).
  3. Whether the x_c-edge and the a0-edge are the same surface (§6).
  4. The LOCAL mechanism: what makes the horizon scale c·H act at galaxy radii? Every
     emergent-gravity theory (Verlinde, Padmanabhan) ASSUMES this; deriving it from the
     4D inward flow would be ESTIF's genuine contribution.
  5. Interior / pressure solution (gravity inside a star) — still a hole.
  6. Relativistic-consistency section (does "flowing space" break GR?) — blocks the gravity letter.
```

---

## 9. File-by-file checklist

| File | Actions |
|---|---|
| **Zenodo preprint** | Apply §1 status block; replacements §2, §3, §4; add §5; add §6 flag; add §8 register. |
| **GitHub README** | Rewrite the one-paragraph abstract to the §0 one-liner; link the four scripts (§10) as reproducibility receipts; add §8 register. |
| **Gravity letter (McGaugh draft)** | Reframe a₀ as horizon quantity (§2); add §5 (BTFR / mass-independence fall out) as a supporting result. Do NOT send until: (a) relativistic-consistency section written, (b) Bullet Cluster answer prepared, (c) Gaztañaga overlap checked — all pre-existing blockers from the six-path plan. |
| **Working notebook** | Stamp §1 on page 1; resolve the "Monday flat / Friday flat-on-average" contradictions by sorting each claim into the SETTLED vs OPEN split above. |

---

## 10. Provenance — receipts backing each claim

Keep these in the repo; each figure/claim above is reproducible:

```
a0_horizon_test.py        -> §2  (sqrt(Lambda) vs Hubble form; the algebraic equivalence)
a0_prefactor_derivation.py-> §3  (target band, numerology floor, river-model edge 0.26-0.28)
estif_flow_sim.py         -> §4  (local flow, BOTH rules; mass-dependence M^{1/3}, non-locality)
estif_horizon.py          -> §5  (horizon background: flat curves, BTFR, mass-indep, magnitude)
```

---

**Bottom line for the notebook:** the flow picture cannot make a₀ locally, but the horizon
can — and does. The only thing left underived is the pure number, which is where the
investigation began (the √3 / 0.128-vs-0.159 prefactor). Settle that, and a₀ is closed.
