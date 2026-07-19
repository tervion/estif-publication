# MILESTONE v6.3 — THE SPLIT

**Date:** 8 July 2026
**Supersedes headline status in:** STATUS.md, ROADMAP.md, README.md, CHANGELOG.md

> ⚠️ **Frozen snapshot, with errata.** This document records the state of the
> project on the day of the split and is not updated as work proceeds. Three of its
> claims were corrected the following day by `CORRECTIONS_v6.3.1.md`: the precise
> sense of "derived" (C6), the conditional status of strong-field deviation claims
> (C1), and the under-marginalized Path Two target (C4). The corrections are marked
> inline below. Later developments — including the Ωm bootstrap — deliberately do
> **not** appear here; see `CHANGELOG.md` and `docs/report/STATUS.md`.

---

## One-sentence statement

As of v6.3, ESTIF splits into two tracks — **Path One (ESTIF-Core, the clean
path)** and **Path Two (ESTIF-Extended, the hard path)** — because a sequence
of derivations (Tasks 4–6, July 2026) proved that the gravity sector can now
stand on a *derived* field equation rather than a borrowed one, and that the
cosmology sector's honest best result is a plain cosmological constant, not the
tilt-based Ω_tilt(z) apparatus that fails DESI DR2.

---

## Why this milestone exists

Between March and July 2026 the project's two sectors moved in opposite
directions, and the move was large enough to force a structural decision rather
than another incremental patch. Three results drove it.

### Result 1 — The gravity field equation is now DERIVED, not postulated (Task 4)

Every prior version obtained Newtonian gravity by an unstated move: the flow
profile was matched to the Schwarzschild solution (set n = ½ so β = √(1−x) = the
GR time-dilation factor), and the "force law" then differentiated that matched
profile. Feeding in the Schwarzschild potential and getting Newton back is
differentiating a potential — it presupposes the answer.

Task 4 removed that. Starting only from the flow axioms (flat 3-slices +
universal speed c + "empty space is not a source"), the Gauss–Codazzi engine
**forces**

```
rho_eff = m'(r) / (4 pi r^2)      i.e.   dm/dr = 4 pi r^2 rho_eff
```

— mass continuity, which is Poisson's equation in integrated form, the standard
(0,0) Einstein equation on flat slices. This was computed from the metric, not
assumed. Consequences, engine-verified:
- vacuum (rho_eff = 0) → v² = 2GM/r uniquely → **exact Schwarzschild**;
- a uniform-density ball → rho_eff = rho0 **exactly** (correct Newtonian source).

The old "D2 postulate" that every previous script had to type in by hand is now
a theorem for vacuum, Newton, and Schwarzschild.

> **Corrected 9 July (C6).** "Derived, not borrowed" is precise only in this sense:
> the three flow axioms uniquely *select* the constraint (energy) sector of General
> Relativity in Painlevé–Gullstrand gauge, forcing mass continuity and hence exact
> Schwarzschild in vacuum — **without matching to the Schwarzschild solution**, which
> was the actual gap. What is *adopted*, not derived from below, is the gravitational
> coupling: the identification of the geometric constraint scalar with 8πG × energy
> density. ESTIF does not derive G or the factor 8π.

Script:
`tests/estif_task4_field_equation.py` (5/5 checks pass). Supporting derivations:
`estif_flow_signature_dynamics.py` (18/18), `estif_converse_flow_law.py`
(Birkhoff's theorem in flow variables), `estif_tmunu_gauss_codazzi.py`
(engine validated against flat FRW and de Sitter).

**Remaining gap (honest):** the strong-field *pressure/stress* sector (a
relativistic interior with pressure) still needs the full off-diagonal T_μν.
It is needed for none of {vacuum, Newton, Schwarzschild}.

**Prerequisite for claiming this in the paper:** the fidelity audit
(`tests/estif_fidelity_audit.py`) confirmed that two of the three axioms —
universal speed c (A2) and "empty space is not a source" (A3) — are currently
present only in the derivation scripts, not in the theory documents, and that
the concept document still runs a competing "shrinking-ruler" picture (A1
conflict). Path One's first job is to write A2 + A3 into the theory and retire
the shrinking-ruler language in favour of the flow (Painlevé–Gullstrand)
picture. The physical content was confirmed by the author as intended; it is a
writing-and-consolidation task, not new physics.

### Result 2 — The cosmology circularity is fixed, but the tilt shape is the real problem (Task 5 / 5b)

The Ω_tilt(z) sector failed DESI DR2 at χ²/N = 10.8. The root cause was
circularity: `x(z) = x₀(1+z)H₀/H_ΛCDM(z)` used ΛCDM as its own ruler. Task 5
replaced that ruler with ESTIF's own H(z) via a fixed-point solve (no ΛCDM
anywhere). Result on real DESI DR2:

| model | χ²/N |
|---|---|
| ΛCDM | 1.92 |
| ESTIF old (circular ruler) | 10.80 |
| **ESTIF new (self-consistent ruler)** | **3.35** |

A large correctness improvement, but still short of ΛCDM. Task 5b localized the
residual: it is the *shape* of the tilt term at low-to-mid z, not the ruler.
Notably, ESTIF's self-consistent tilt produces a w(z) that already thaws in the
direction DESI prefers (tracking DESI's published w(z) within ~0.05), and its
3.35 is essentially where DESI's *own* published w0,wa land on the BAO-only
subset (3.09). Script: `tests/estif_task5_desi_selfconsistent.py`,
`tests/estif_task5b_cosmo_eos.py`.

### Result 3 — The frozen-eddy reframe: the simple derived answer beats the fitted one (Task 6)

The attempt to *derive* the eddy's equation of state from first principles was
run and, honestly, the two natural reductions failed badly:

| eddy model | derived w | χ²/N vs DESI DR2 |
|---|---|---|
| E1 conserved-angular-momentum spin | +1 (stiff) | 3232 |
| E2 expansion tracker | thaws to ~0 | 754 |

Both make the eddy dilute or blueshift the wrong way; DESI wants a component
that stays near w = −1 with only mild thawing. But lining the models up exposed
the decisive fact:

> The **frozen-eddy limit** — a constant cosmic eddy energy density, which Task 4
> derives for free (constant density → de Sitter → w = −1) — scores **χ²/N = 1.92**,
> tying ΛCDM, and **beats the project's own tilt formula (3.35)**.

On DESI, the entire tilt apparatus (N_MAX, B, the Ω_tilt sign-flip, the z < 2
cutoff) is a **net negative** relative to the simpler, derived, cosmological-
constant limit underneath it. Script: `tests/estif_task6_eddy_eos.py`.

---

## The decision: split into two tracks

Applying the project's own RHAC decision filter (Step 2, filter 4: *"Does it
improve agreement with observations?"*), the tilt cosmology fails against its own
frozen-eddy limit and is marked 🔴 pivot, not 🟡 active. Rather than delete the
tilt line outright or keep defending it, the project forks:

### Path One — ESTIF-Core (the clean path) ✅ recommended default

**Thesis:** the flow principle gives *derived* gravity, and its simplest,
honest cosmology is a cosmological constant that ties ΛCDM.

Scope:
- Gravity sector as-is, now strengthened: the a₀/MOND derivation sits on a
  *derived* field equation (Task 4), not a Schwarzschild match.
  > **Corrected 9 July (C1).** Because the ESTIF vacuum is *exactly* Schwarzschild,
  > the strong-field *deviation* claims (EHT shadow offset, LISA GW delay) cannot be
  > sourced by the vacuum. They are conditional on the un-derived eddy-stress sector.
  > The observations remain consistent with ESTIF; the deviation from GR awaits
  > derivation. The a₀/MOND chain is unaffected — it never invokes a vacuum deviation.
- Cosmology = frozen cosmic eddy → cosmological constant (χ²/N = 1.92, derived,
  zero tilt parameters).
- The tilt apparatus (Ω_tilt(z), N_MAX, B, sign-flip, z<2 cutoff) is **retired**
  from the cosmology claim and moved to a documented "explored and set aside"
  appendix.
- Writing tasks: put A2 (universal c) and A3 (vacuum sources nothing) into the
  theory documents; resolve the A1 shrinking-ruler conflict in favour of flow.

Deliverable: the gravity letter (already drafted, unaffected — in fact
strengthened) plus a short, honest cosmology statement. Defensible today.

### Path Two — ESTIF-Extended (the hard path) 🔬 high-risk research

**Thesis:** the cosmic eddy produces a *small, near-frozen thawing* on top of
w = −1 that matches the newest DESI preference (nominal target χ²/N ≈ 0.66).

> **Corrected 9 July (C4).** The 0.66 target comes from a BAO-alone CPL fit with rd,
> H₀, and Ωm held fixed at Planck values. Fixing nuisance parameters inflates the
> apparent evolving-dark-energy advantage; DESI's own BAO-alone preference is much
> milder. The target must be re-derived under marginalization before it is used to
> justify Path Two.

Scope:
- Start from the derived w = −1 frozen limit and compute the **leading
  correction** from the full rotating-shear / vorticity stress tensor — the
  cosmological half of the T_μν work, aimed correctly (a small perturbation, not
  a new dominant term).
- The two naive reductions (E1, E2) are already proven wrong; this requires the
  complete off-diagonal stress tensor, genuinely hard, timeline unknown.

Deliverable (if it succeeds): a derived evolving dark energy competitive with —
or better than — ΛCDM on DESI. If it fails: fall back to Path One, no worse off.

---

## Status summary after this milestone

| Sector | v6.2 status | v6.3 status |
|---|---|---|
| Gravity — field equation | matched to Schwarzschild (borrowed) | **DERIVED** (Task 4); pressure sector open |
| Gravity — a₀ / MOND / SPARC | solid | solid, now on a derived foundation |
| Cosmology — circularity | circular (χ²/N 10.8) | **fixed** (self-consistent, 3.35) |
| Cosmology — honest best | tilt Ω_tilt(z), fails DESI | **frozen eddy = Λ, ties ΛCDM (1.92)** |
| Cosmology — tilt apparatus | defended | **retired to appendix (Path One)** |
| Project structure | single track | **split: Core (clean) + Extended (hard)** |

---

## Provenance

All numbers above are runtime output of scripts in `tests/`, each of which
fetches the real DESI DR2 data (CobayaSampler) and/or runs the validated
Gauss–Codazzi engine. Reproduce with:

```bash
python3 tests/estif_task4_field_equation.py         # field equation derived
python3 tests/estif_flow_signature_dynamics.py      # signature + SR + Newton
python3 tests/estif_converse_flow_law.py            # Birkhoff in flow variables
python3 tests/estif_task5_desi_selfconsistent.py    # de-circularized DESI test
python3 tests/estif_task5b_cosmo_eos.py             # DESI-preferred w(z)
python3 tests/estif_task6_eddy_eos.py               # frozen-eddy reframe
python3 tests/estif_fidelity_audit.py               # axiom presence audit
```

**Milestone version:** 6.3 — "The Split"
**Date:** 8 July 2026
