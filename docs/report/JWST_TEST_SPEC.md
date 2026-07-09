# SPECIFICATION — The JWST Early-Structure-Formation Test

**Status:** Ready specification. Baseline computed; ESTIF input identified; target quantified.
**Companion script:** `tests/estif_jwst_growth_spec.py` (runs the LCDM baseline + target on a laptop)
**Date:** 8 July 2026

> **Target location in repo:** `docs/report/JWST_TEST_SPEC.md`

---

## 1. The claim being tested

JWST finds galaxies at z ≈ 8–13 that look too massive and too mature for how
little cosmic time had passed — a live, unresolved tension in ΛCDM. ESTIF's
hypothesis: early structure assembled **faster** than standard gravity allows,
because the flow's effective gravity for density perturbations is stronger than
Newtonian, so massive galaxies exist early **without** violating anything.

This is a *timing/mass* claim, not a *size* claim. (The apparent-size version
fails: angular size is fixed by H(z), and ESTIF-Core's H(z) is identical to
ΛCDM's, so it predicts the same sizes. Only the growth-rate/timing version can
distinguish ESTIF.)

---

## 2. The observable to predict

The **cumulative comoving number density of massive halos** at high redshift,

```
n(>M_halo, z)   for M_halo ~ 1e11–1e12 M_sun,  z ~ 8–12
```

equivalently the **maximum available stellar-mass density** `ρ*(>M, z) = ε · f_b ·
ρ_collapsed(>M, z)`. This is what the JWST massive-galaxy candidates constrain:
too many massive galaxies too early means too many massive halos too early.

---

## 3. The ΛCDM baseline (computed, Planck18, Sheth–Tormen)

Real numbers from the companion script (physical Mpc⁻³):

| z | n(>10¹⁰) | n(>10¹¹) | n(>3×10¹¹) | n(>10¹²) |
|---|---|---|---|---|
| 8.0 | 1.1×10⁻² | 1.6×10⁻⁴ | 7.5×10⁻⁶ | 1.1×10⁻⁷ |
| **9.1** | 3.9×10⁻³ | **2.3×10⁻⁵** | **8.0×10⁻⁷** | 6.2×10⁻⁹ |
| 10.0 | 2.2×10⁻³ | 8.6×10⁻⁶ | 2.1×10⁻⁷ | — |
| 12.0 | 4.6×10⁻⁴ | 7.9×10⁻⁷ | — | — |

**The number to beat (z = 9.1):** n(>10¹¹ M_sun) = 2.3×10⁻⁵ Mpc⁻³;
n(>3×10¹¹ M_sun) = 8.0×10⁻⁷ Mpc⁻³.

**The tension (Boylan-Kolchin 2023 / Labbé 2023):** the two most massive z ≈ 7.5–9.1
JWST candidates (M\* ~ 10^10.5–10^11) require converting a very large fraction of
available baryons to stars under ΛCDM — a star-formation efficiency pushed toward
or past plausible values. ΛCDM does not have enough massive halos early enough.

---

## 4. The one ESTIF input required

Everything above is standard. The single missing ESTIF piece is:

> **D_ESTIF(z) — the linear growth factor under the ESTIF field equation,**
> linearized on an FRW background.

This is the **natural extension of the derived field equation (Task 4) to linear
perturbations**: perturb the flow metric on an expanding background, run the same
Gauss–Codazzi machinery, and read off the growth of the density contrast δ. It is
an **analytic ODE calculation — laptop-tractable, not N-body**.

The physical question it answers: does ESTIF's flow give perturbations an
**effective gravitational source stronger than 4πGρ** (an effective G_eff > G, or an
extra eddy-clustering term)? If yes, growth is enhanced and structure forms earlier.

---

## 5. The target (quantified)

Massive-halo abundance is **exponentially** sensitive to growth. From the companion
script, at z = 9.1:

| growth enhancement D_ESTIF/D_LCDM | boost in n(>10¹¹) | boost in n(>3×10¹¹) |
|---|---|---|
| +5% | 1.95× | 2.4× |
| **+13%** | **4.8×** | **7.8×** |
| +20% | 9.0× | 18× |

Inverse (reservoir boost → required growth enhancement):

```
 2× reservoir  ←  +5%  growth enhancement at z = 9.1
 5× reservoir  ←  +13% growth enhancement at z = 9.1
10× reservoir  ←  +20% growth enhancement at z = 9.1
```

**Target: D_ESTIF(9.1) ≳ 1.13 × D_LCDM(9.1)** relieves the tension (5× reservoir).
This is a **small** modification — but it is subject to the σ8/S8 hard filter in
§5a: the boost must be *transient or scale-dependent*, not persistent. See §5a.

---

## 5a. HARD FILTER — σ8/S8 (the low-redshift structure constraint)

The 13% target above is **necessary but not sufficient**. A growth enhancement
that *persists to z = 0* at the 8 Mpc/h scale overproduces present-day structure:

| growth boost persisting to z=0 | σ8(z=0) | vs Planck (0.811 ± 0.006) |
|---|---|---|
| +5% | 0.852 | +7σ |
| +13% | 0.917 | +18σ |
| +20% | 0.973 | +27σ |

A persistent 13% boost gives σ8 ≈ 0.92 — ~18σ above Planck, and the **wrong sign**
relative to weak-lensing S8 (KiDS/DES mildly prefer *suppressed* low-z growth).

**Therefore D_ESTIF(z) must be TWO-SIDED:**
- **enough early growth** for JWST (≳13% at z ≈ 9), AND
- **≈ standard late growth** for σ8 (enhancement → 1 by z ≲ 2).

This is satisfied only if the enhancement is **transient** (concentrated at high z)
and/or **scale-dependent** (stronger at small scales / high k, leaving the 8 Mpc/h
scale that sets σ8 essentially untouched — a natural signature of a scale-dependent
effective G). This makes the derivation target **more specific, not merely larger**:
D_ESTIF/D_LCDM must be a *localised* early or small-scale boost, not a monotonic one.

---

## 6. What is laptop-tractable vs blocked

| Piece | Method | Status |
|---|---|---|
| ΛCDM baseline abundance | Sheth–Tormen (colossus) | ✅ done (this spec) |
| Growth → abundance mapping | Press–Schechter sensitivity | ✅ done (this spec) |
| **D_ESTIF(z)** | perturbed Gauss–Codazzi (ODE) | 🔬 **the derivation to do — laptop** |
| Collapse threshold δ_c | spherical collapse in ESTIF | 🔬 laptop (secondary) |
| ESTIF abundance prediction | ST mass function with D_ESTIF | ✅ trivial once D_ESTIF known |
| Nonlinear internal halo structure (rotation curves, δ~10⁵) | N-body | 🔴 **blocked — but SEPARATE test** |

**Key point:** the population-level JWST test is **not** behind the N-body wall. The
wall only blocks the *internal* structure of individual halos (the dark-matter
rotation-curve question). The *abundance/timing* question — which is what JWST
actually constrains — needs only the linear growth derivation, which is analytic.

---

## 7. Path assignment (corrected)

At z ≈ 9 the universe is matter-dominated (dark energy ≈ 0.2% of the density), so
**both** ESTIF cosmologies give essentially ΛCDM's H(z) there:
- Path One (frozen eddy): H(z) = ΛCDM exactly.
- Path Two (thawing): the difference is at *low* z; negligible at z ≈ 9.

Therefore the growth boost **cannot come from the expansion history**. It must come
from ESTIF's **modified gravity** — a stronger effective perturbation source — which
is tied to the **derived field equation**, i.e. the **gravity sector (Path One's
core)**, not the dark-energy sector.

*(This corrects an earlier framing that placed JWST relief in Path Two's "different
clock." The clock difference is dynamically irrelevant at z ≈ 9; the lever is
gravity/growth.)*

---

## 8. The go/no-go, stated plainly

1. Derive D_ESTIF(z) by linearizing the field equation on FRW (perturbed
   Gauss–Codazzi). Analytic, laptop.
2. Compute D_ESTIF(9.1)/D_LCDM(9.1).
3. **If ≳ 1.13 AND transient/scale-dependent (passes §5a):** ESTIF relieves the
   JWST tension with a derived, zero-free-parameter growth enhancement — a genuine,
   distinguishing, falsifiable win exactly where ΛCDM strains, without breaking σ8.
   This would be a flagship result. (A boost that relieves JWST but persists to z=0
   fails σ8 and is NOT a win — the two-sided filter is mandatory.)
4. **If ≈ 1.0:** ESTIF's growth matches ΛCDM and does *not* resolve JWST; the tension
   stays a shared problem. Honest null, no worse than ΛCDM.
5. **If < 1.0:** ESTIF makes it worse — a potential falsification of the modified-
   growth claim.

Run `tests/estif_jwst_growth_spec.py` to reproduce the baseline and the target; it
takes the growth enhancement as a knob until D_ESTIF(z) is derived.

---

**Specification version:** 1.0 · **Date:** 8 July 2026 · **Companion:** `estif_jwst_growth_spec.py`
