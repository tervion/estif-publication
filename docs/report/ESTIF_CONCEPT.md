# ESTIF: The Conceptual Foundation

**Document type:** Conceptual framework — not a mathematical derivation
**Purpose:** Explains the physical intuition behind ESTIF before any equations
**Version:** 6.3 (July 2026) — "The Split"

>
> **What changed in v6.3 (read this first):**
> 1. **The gravity field equation is now derived, not borrowed.** Earlier versions
>    reproduced Newton by matching the flow profile to the Schwarzschild solution.
>    v6.3 derives it from the flow axioms directly (see *Gravity Is Now Derived*).
> 2. **The three axioms are now stated explicitly.** A fidelity audit found that two
>    of them lived only in the derivation scripts. They are written into the theory
>    below (see *The Three Axioms*).
> 3. **The shrinking-ruler picture is retired.** It conflicted with the flat-slice
>    flow picture the derivations require. Cosmic expansion is now described as
>    projected inward flow, not universal shrinkage.
> 4. **Cosmology is reframed.** The honest best result is a frozen cosmic eddy =
>    cosmological constant (ties ΛCDM on DESI DR2). The evolving Ω_tilt(z) apparatus
>    is retired to an "explored and set aside" status (see *The Cosmology Reframe*).
> 5. **The project splits into two paths** (see *The Two Paths*).

---

## The Central Idea

3D space is a flat sheet carried steadily through a 4th spatial dimension we
cannot perceive, measure, or point to directly. Everything the sheet contains is
carried with it.

Two things we observe are, in this picture, shadows of that motion:
- **Gravity** is the local tilt and slowing of the flow near mass.
- **Cosmic expansion** is the ordinary growth of our 3D sheet, carried along by
  the flow — the geometric projection of 4D motion onto the surface we inhabit.

The expansion we observe is real, measurable, and consistent with all
cosmological observations. What ESTIF questions is its *cause*: not a mysterious
outward push, but the geometry of a sheet in motion.

---

## The Three Axioms

Everything in ESTIF follows from three statements about the moving sheet. In
v6.3 they are stated explicitly, because the gravity derivation (below) rests on
exactly these and nothing else.

**A1 — Space flows, it does not stretch.**
The 3D sheet is *flat* (its internal geometry has no curvature of its own) and is
carried bodily through the 4D bulk. This is the Painlevé–Gullstrand picture: flat
spatial slices plus a flow. It is *not* a picture of space stretching like rubber,
nor of everything shrinking together (that older idea is retired in v6.3 — see the
note under *Cosmic Expansion*).

**A2 — Everything moves through the bulk at the speed of light.**
Every object — you, a planet, a photon — moves through the 4D bulk at exactly one
speed: *c*. What we call "sitting still" is moving through the 4th dimension (through
*time*) at full speed. What we call "moving through space" is stealing part of that
speed and pointing it sideways within the sheet, which leaves less pointing
timeward — so a moving clock runs slow. In symbols, the total motion splits as

```
dw² + dσ² = c² dt²
```

where `dw` is the step through the 4th dimension (the passage of time) and `dσ` is
the step within our 3D space. Proper time *is* the distance travelled through the
4th dimension: `c dτ = dw`.

**A3 — Empty space is not a source.**
Where there is no matter, the flow passes through undisturbed and originates
nothing of its own. Gravity begins only at mass, exactly as an electromagnetic
field begins only on its charges and merely crosses the vacuum between them. Empty
space carries the flow but sources no gravity.

Everything below — Lorentzian time, special relativity, Newton's law, and exact
Schwarzschild geometry — is a consequence of A1 + A2 + A3.

---

## Why We Can't Perceive It: The Ant on the Soccer Ball

Imagine drawing meridian lines on a soccer ball — exactly like longitude lines on a
globe. An ant starts walking from the equator toward one of the poles, following one
of these lines.

To the ant, the path is perfectly straight. It has no reason to think otherwise. It
cannot perceive the curvature of the ball's surface. It has no instrument that can
detect that its straight line is, from a higher-dimensional perspective, curving
toward every other meridian and eventually converging at the pole.

This is exactly the situation described by General Relativity. The ant is us. The
ball is spacetime. The lines converging at the poles are geodesics — paths that feel
straight but are shaped by the curvature of the surface they exist on.

This is well-established physics. ESTIF does not dispute it. **ESTIF lives one level
deeper:** it asks what the surface is moving through, and how fast.

---

## Cosmic Expansion as Projected Inward Flow

Now the soccer ball (3D space) is moving through a room (4D space) in the one
direction the ant cannot perceive — straight through the room, orthogonal to its
surface. The ant, walking its meridians on the surface, experiences:

1. Its path curving due to the ball's curvature → **gravity** (GR, already known);
2. The ball moving through the room → **cosmic flow through 4D** (the ESTIF claim);
3. The ordinary growth of the surface as it is carried along → **the expansion we
   measure**.

The ant cannot directly observe #2 or #3. It observes only their *shadows* — the
projections of these 4D effects onto the 3D surface it inhabits. Cosmic expansion is
the shadow of the flow.

> **Retired in v6.3 — the shrinking-ruler picture.**
> Earlier versions of this document explained expansion as *everything shrinking
> together* (your ruler shrinks, the thing you measure shrinks, so distant objects
> look farther away). That mechanism is inconsistent with axiom A1 (flat slices
> carried through the bulk), which the gravity derivation requires — flat-slice flow
> and universal shrinkage are two different geometries and cannot both be the
> foundation. v6.3 keeps the flow and drops the shrinkage. Under the flow picture,
> the 3D sheet expands as in standard cosmology, driven by matter plus a
> cosmological constant that ESTIF identifies as the frozen cosmic eddy (see *The
> Cosmology Reframe*). The scale-and-speed intuition that motivated the old analogy
> is preserved — a smaller sheet moving at the same bulk speed *c* covers proper
> distance differently — but the mechanism is projection of motion, not shrinking of
> rulers.

---

## Time Is Motion Through the Fourth Dimension

Axiom A2 is the heart of the framework, so it deserves its own picture.

A clock far from all mass, sitting still, spends its entire light-speed budget
moving through the 4th dimension — moving through *time* — and so ticks as fast as
anything can tick. Speed it up through space, and part of that budget is redirected
sideways within the sheet; less points timeward, so it ticks slower. That is
special relativity, recovered exactly:

```
dτ/dt = √(1 − v²/c²)
```

Park the clock next to a mass, and the *flow itself* carries it partly sideways in
the 4th dimension; again less of the motion points timeward, so it ticks slower.
That is gravitational time dilation. One rule — motion through the bulk at *c* —
produces both kinds of slowing, and they are the same event seen twice: the sideways
part is what we feel as gravity, the leftover timeward part is what the clock reads.

**A note on the cosmic flow speed.** In the cosmology sector you will see the value
`v_flow = c·x₀ ≈ 0.31c`. This is *not* the full speed of the flow — under A2 nothing
moves through the bulk at 0.31c; everything moves at *c*. It is the **sideways
(within-sheet) component** of the cosmic flow at the largest scale; the remaining
budget points timeward. Reading 0.31c as "the flow speed" is a mislabel corrected in
v6.3: it is one component of a motion whose total is always *c*.

---

## What Gravity Is Under This Framework

In GR, mass curves spacetime. ESTIF asks: curves it *into what*?

The answer: mass tilts the 3D sheet into the 4th dimension and slows the flow's
advance in that direction. Near a massive object the sheet is no longer perpendicular
to the 4th dimension — it leans toward it, and clocks carried by the flow there run
slow (the A2 picture above). A photon travelling near the mass traverses a path
extended in 4D; a 3D observer measures only the 3D projection of that path.

This gives two complementary descriptions of the same thing:
- **The flow description** (primary, and now *derived* — see the next section): near
  mass the flow acquires an inward 3-velocity `v_r`, and gravity is the gradient of
  that flow.
- **The tilt description** (used for strong-field *deviations* from GR): how steeply
  the sheet leans, `sin θ = (Rs/r)^n`, and how much of the 4D correction remains
  visible in 3D, `β = cos θ = √(1 − (Rs/r)^{2n})`.

The two meet at a specific curvature (see *Gravity As Generalized Time Dilation*).

---

## Gravity Is Now Derived, Not Borrowed (v6.3 — the central upgrade)

For every version before v6.3, this document obtained Newton's law by an unstated
shortcut: the flow profile was *matched* to the Schwarzschild solution (the tilt
exponent was set to n = ½ so that β = √(1 − Rs/r), the GR time-dilation factor),
and the force was then read off from that matched profile. Feeding in the
Schwarzschild potential and getting Newton back is differentiating a potential — it
presupposes the answer. The fidelity audit made this explicit, and it was the single
weakest point in the framework.

**v6.3 removes the shortcut.** Starting from A1 + A2 + A3 and nothing else, the flow
metric (flat slices, a flow at speed *c*) is fed through a symbolic Gauss–Codazzi /
ADM engine, and the engine *forces* the effective energy density to be

```
ρ_eff = m′(r) / (4π r²)        equivalently        dm/dr = 4π r² ρ_eff
```

This is **mass continuity — Poisson's equation in integrated form** — the standard
(0,0) Einstein equation on flat slices. It was *computed from the geometry*, not
assumed. The consequences, all engine-verified:

- **Vacuum (axiom A3: ρ_eff = 0)** forces the flow speed to `v² = 2GM/r` uniquely —
  the escape-velocity ("river") profile — and this makes the *full* Einstein tensor
  vanish: it is **exact Schwarzschild**, not a weak-field approximation.
- **Inside matter**, a uniform-density ball returns `ρ_eff = ρ₀` **exactly** — the
  correct Newtonian source.

So the old "Poisson postulate" is now a **theorem** for vacuum, Newton, and
Schwarzschild.

> **What "derived" means precisely (C6).** The axioms uniquely *select* the
> constraint (energy) sector of General Relativity in Painlevé–Gullstrand gauge.
> That is what forces mass continuity, hence exact Schwarzschild in vacuum and the
> Newtonian source in the weak field — **without matching to the Schwarzschild
> solution**, which was the actual gap. What is *adopted*, not derived from below, is
> the gravitational *coupling*: the identification of the geometric constraint scalar
> with 8πG × energy density. ESTIF does not derive Newton's G or the factor 8π. The
> honest claim is that the force law is no longer matched to GR's vacuum solution;
> the coupling to matter is still the standard Einstein–Hilbert one.

The naive alternative — space draining like water down a plughole
(volume-conserving flow) — was tested and *fails*, giving a 1/r⁵ force; the vacuum
condition instead conserves the free-fall energetics that give exactly Newton. The
axiom picks the right law automatically; it is not tuned.

Supporting results, all from the same engine:
- Lorentzian signature `(−,+,+,+)` **emerges** from a Euclidean bulk plus the A2
  speed-*c* constraint — the minus sign is produced by solving the constraint, not
  inserted.
- Special-relativistic kinematics (time dilation, null photons, rest-frame clocks)
  follow with zero free parameters.
- The vacuum condition uniquely forcing `v² = 2A/r` is **Birkhoff's theorem restated
  in flow variables**.

**What remains open (honest).** The strong-field *pressure/stress* sector — a fully
relativistic interior with pressure — needs the complete off-diagonal stress tensor.
It is needed for *none* of {vacuum, Newton, Schwarzschild}, but it is the remaining
rigorous step for gravity.

Scripts: `estif_task4_field_equation.py` (5/5 checks), `estif_flow_signature_dynamics.py`
(18/18), `estif_converse_flow_law.py`, `estif_tmunu_gauss_codazzi.py` (engine validated
against flat FRW and de Sitter).

---

## Gravity As Generalized Time Dilation

The tilt description and the flow description meet at one curvature. The Schwarzschild
time-dilation factor — how much slower clocks tick near mass — is

```
τ(x) = √(1 − x)     where x = Rs/r
```

The ESTIF tilt suppression is `β(x) = √(1 − x^{2n(x)})`. These are identical when
`x^{2n} = x`, i.e. when **n = ½**. The dynamic exponent `n(x) = 33.265 × exp(−15.429 x)`
passes through n = ½ naturally at **x = 0.272** — a value that emerged from
calibration, not choice.

| Regime | n(x) | β vs τ | Meaning |
|---|---|---|---|
| Flat space (x → 0) | 33.3 | β = τ = 1 | No gravity, no time dilation |
| Crossover (x = 0.272) | **0.500** | **β = τ exactly** | ESTIF = GR time dilation |
| Cosmological (x = 0.311) | 0.275 | β < τ | ESTIF weaker than GR |
| M87\* photon sphere (x = 0.667) | 0.001 | β ≪ τ | Strongly suppressed |
| Horizon (x → 1) | → 0 | β → 0 | Time stops |

GR time dilation is thus *one special case* of the ESTIF tilt family (n = ½). Below
the crossover ESTIF predicts a stronger effect than GR; above it, weaker. In v6.3
this sits on firmer ground: the vacuum flow derivation already gives exact
Schwarzschild, and the tilt formula is the tool for the *deviations* from it that
EHT and LISA could detect. See `test_gravity_time_connection.py`.

---

## The Combined Formula (strong-field deviations)

```
x          = curvature ratio (Rs/r locally)
n(x)       = 33.265 × exp(−15.429 x)     ← dynamic tilt exponent
β(x)       = √(1 − x^{2n(x)})
Observable = √β(x)                        ← square-root (amplitude) projection
```

Jointly calibrated, this reproduces three independent strong-field observations with
no free parameters after calibration:

| Test | Result | Status |
|---|---|---|
| EHT M87\* shadow | 0.00σ tension, shadow = 42.0 μas | ⚠️ consistent; *deviation* conditional |
| Planck Λ (as a local-tilt scale) | ratio = 1.0000 | ✅ calibration match |
| LISA GW delay (65 M☉) | 491 μs, S/N = 49σ | ⚠️ conditional |

> ⚠️ **Conditional (C1).** The derived vacuum is *exactly* Schwarzschild, so a
> deviation from GR in shadows or in vacuum GW propagation cannot come from the
> vacuum. It must be sourced by the non-vacuum eddy background — a sector not yet
> derived. The observations remain *consistent* with ESTIF; the predicted *deviation
> from GR* awaits derivation. A single-speed flow forces p_r = −ρ, which is precisely
> why the ansatz has no vacuum deviation to offer. The Λ entry is a calibration match,
> not a vacuum deviation, and is unaffected.

**Why √β:** you measure amplitude, but energy scales as amplitude squared. The 4D
correction has amplitude β; the 3D measurement captures √β. **Why n is dynamic:** n
varies with the local environment exactly as *g* varies planet to planet in Newton's
law — the trajectory equation is universal, the input changes. Strong gravity
suppresses the tilt exponent; weak gravity lets it grow. The formula is dormant in
flat space and activates only in extreme curvature. See `test_joint_calibration.py`.

> **Scope note (v6.3).** The combined formula is the framework's tool for
> *strong-field deviations from GR* and for the a₀ scale (below). Its use as a
> *cosmological dark-energy law* (Ω_tilt(z)) is retired — see *The Cosmology Reframe*.
> The two uses are independent; retiring the cosmological one leaves the strong-field
> one untouched.

---

## Gravity = Eddies = Time

All three descriptions of gravity — time dilation, geometric tilt, and 4D eddy spin
— are the same phenomenon. Define the eddy spin rate `ω(x) = H₀ x^{n(x)}`. Then:

| Description | Formula |
|---|---|
| GR time dilation | τ(x) = √(1 − x) |
| ESTIF tilt | √β(x) = √(1 − x^{2n(x)}) |
| Eddy spin energy | (ω/H₀)² = x^{2n(x)} |

At the crossover x = 0.272 (n = ½), `(ω/H₀)² = x`, so gravitational acceleration is
the spatial gradient of the squared eddy spin, `a = −c²∇(ω/H₀)²/2`, which returns
`GM/r²`. In v6.3 this is consistent with — and underwritten by — the derived field
equation: the eddy description and the flow description agree, and the flow
description is the one that is now derived from first principles.

---

## The Electron As Natural Scale

The formula's two calibrated parameters trace to one fixed scale — the classical
electron radius in Planck units:

```
N_MAX ≈ 5/7 × ln(r_e / l_P) = 33.291   (0.079% off)
B     ≈ 1/3 × ln(r_e / l_P) = 15.536   (0.693% off)
```

`r_e` is the scale where electromagnetic self-energy equals rest-mass energy — the
boundary between electromagnetism and gravity — built from fundamental constants
only and identical everywhere at all times. This connects the tilt parameters to a
rigid, universe-independent ruler. The **B = L/3** multiplier is derived from 3D
isotropy (the same argument as the 1/√3 in the MOND section); **N_MAX = 5/7 × L**
remains conditional on the still-open geometric derivation of the crossover x_c =
0.272. See `test_electron_connection.py`, `test_multiplier_derivation.py`.

---

## The MOND Derivation — First Principles

The MOND critical acceleration a₀ ≈ 1.2×10⁻¹⁰ m/s² has predicted galactic rotation
curves for 40 years, but its origin has always been *fitted*. ESTIF derives it from
geometry with zero free parameters:

```
Step 1 — Force law:  a = −c²/2 ∇(ω/H₀)² = GM/r²
Step 2 — Cosmic flow (sideways component):  v_flow = c·x₀ = c·Ωm
Step 3 — 3D isotropic projection:  v_3D = v_flow/√3
Step 4 — Threshold:  a₀ = v_3D × H₀ = H₀·c·x₀/√3
```

```
Result:  a₀ = 1.179×10⁻¹⁰ m/s²    MOND empirical: 1.200×10⁻¹⁰    Agreement: 1.72%
```

**v6.3 note on Step 1.** Step 1's force law is now *underwritten by the derived field
equation* (the *Gravity Is Now Derived* section): `a = −c²∇(ω/H₀)²/2` with
`(ω/H₀)² = Rs/r` is the weak-field face of the vacuum result `v² = 2GM/r`, which is
forced by the flow axioms rather than matched to Schwarzschild. This closes the
"borrowed recipe" gap that previously sat under Step 1.

**Why 1/√3 is forced, not chosen.** For isotropic motion in N dimensions,
equipartition gives each axis 1/N of the kinetic energy, so `v_1D = v_rms/√3` in 3D.
This is standard kinetic theory, applied before any comparison with MOND; a full
kinetic theory of the eddy background is identified as future work. Of 12 candidate
factors, only 1/√3 gives < 5% agreement, and it is the only one with an independent
physical basis. The same factor recurs in the sound speed `c_s = v_rms/√3`, the Jeans
criterion, the virial theorem, and the B = L/3 multiplier.

**Tully–Fisher — validated against SPARC.** `v_flat⁴ = G·M_bar·a₀` gives `v_flat ∝
M^{1/4}` — the observed exponent. Against 87 quality-1 SPARC galaxies (Lelli et al.
2016, Υ* = 0.50): **RMS = 15.6%**, within the observed BTFR scatter; 82% within 20%,
97% within 30%. The −7.6% mean bias traces entirely to the stellar-mass calibration.
See `derive_mond_from_geometry.py`, `test_sparc_tully_fisher.py`.

**a₀ redshift constancy.** In the comoving frame appropriate for bound galaxies, H(z)
cancels exactly from a₀ = H₀cx₀/√3, giving `a₀ = c²/(r_universe,comoving·√3) =
constant` (deviation 2.22×10⁻¹⁶, machine epsilon). Consistent with high-z
Tully–Fisher at z ≈ 0.75–2.2 (Di Teodoro+2021, Übler+2017, Tiley+2019, all ≤ 2σ).
**Parameter independence:** across 3,600 combinations of H₀ ∈ [65,75] and Ωm ∈
[0.27,0.33], every case stays within ±20% SPARC scatter; 8 published datasets pass.

---

## The Cosmology Reframe (v6.3)

This is the sector that changed most. The short version: **ESTIF's honest best
cosmology is a plain cosmological constant, and it comes for free from the derived
gravity.**

### The frozen cosmic eddy is a cosmological constant

The derived field equation has a clean corollary: a flow whose effective energy
density is *constant* is exact de Sitter space — a cosmological constant, `w = −1`.
The cosmic eddy, if it does not dilute, *is* dark energy as a cosmological constant.
Tested against real DESI DR2 BAO data, this frozen-eddy limit scores

```
χ²/N = 1.92   —   tied with ΛCDM.
```

It is derived (from constant eddy density), needs no tilt formula, and no fitted
parameters beyond matter and the constant itself.

### Why the evolving Ω_tilt(z) apparatus is retired

Earlier versions replaced the cosmological constant with an *evolving* geometric
term `Ω_tilt(z)` built from the tilt formula:

```
H²(z) = H₀² [Ωm(1+z)³ + Ω_tilt(z)],    Ω_tilt(z) = Ω_Λ (obs_now/obs_z)²
```

Against DESI DR2 this failed badly (χ²/N = 10.8). The failure had two layers, both
now understood:

1. **Circularity.** The curvature `x(z) = x₀(1+z)H₀/H_ΛCDM(z)` used ΛCDM as its own
   ruler. Fixing this — solving for `x(z)` with ESTIF's *own* H(z) self-consistently
   — dropped the score from **10.8 → 3.35**. A real correctness fix (`estif_task5_
   desi_selfconsistent.py`).
2. **The tilt shape itself.** Even de-circularized, 3.35 is worse than the frozen
   limit's 1.92. Two natural attempts to *derive* how the eddy energy evolves both
   failed against DESI: a conserved-spin eddy gives a stiff `w = +1` (χ²/N = 3232);
   an expansion-tracking eddy dilutes to matter-like (χ²/N = 754). DESI wants a
   component that stays *near* `w = −1` with only mild thawing — which is close to
   what the tilt already does, but the frozen limit does it better.

The conclusion the data forces: **on DESI, the entire Ω_tilt apparatus (the dynamic
n, N_MAX, B, the sign choice, the z < 2 cutoff) is a net negative relative to the
plain cosmological constant underneath it.** So the cosmological dark-energy claim is
retired to an "explored and set aside" appendix. The tilt formula's *strong-field*
use (EHT, LISA, the a₀ scale) is unaffected — that is a different application at a
different scale.

### What DESI actually prefers (and where Path Two aims)

DESI DR2 BAO alone is best fit by a mild *thawing* dark energy (w rising toward −0.85
today from more negative in the past; best-fit CPL χ²/N ≈ 0.66). ESTIF's
self-consistent tilt already leans this way — its w(z) tracks DESI's published curve
to within ~0.05 — and its 3.35 is essentially where DESI's *own* published w0,wa land
on this BAO-only subset (3.09). So the honest picture is not "ESTIF fails DESI"; it
is "a plain cosmological constant ties ΛCDM, and a small derived thawing correction
*could* do better — if it can be derived." That derivation is Path Two.
Scripts: `estif_task5b_cosmo_eos.py`, `estif_task6_eddy_eos.py`.

> **Retired to appendix:** the Hubble-radius Λ bridge, the Λ-drift (0.023%/Gyr)
> prediction, the Ω_tilt(z) evolution law, and the pre-2026 six-low-z-test fits.
> These are documented as explored directions that the frozen-eddy result and the
> DESI comparison superseded. They are not deleted; they are no longer load-bearing.

---

## Dark Matter — Analytical Results (unchanged), Now Connected to the Derived Machinery

The eddy dark-matter results stand, and v6.3 connects them to the derived field
equation.

**The numerical identity:**

```
x₀ = R_H/r_universe = 0.310734     Ωm (Planck) = 0.311100     (0.12%)
x₀ − Ωb = 0.261734                 Ωdm (Planck) = 0.262000     (0.10%)
```

> **This is a consistency relation, not a prediction (C2).** Here r_universe is the
> ΛCDM particle horizon — an integral that itself contains Ωm. The 0.12% agreement is
> a self-consistency of the geometric picture with Planck-calibrated values, not an
> Ωm-independent derivation. Making it a prediction requires deriving r_universe from
> the flow framework without the Ωm-dependent horizon integral. The numbers stand;
> the word "predicted" does not.

**Collisionless dynamics (not a fluid):**velocity dispersion `σ(r) = r√(2πGρ_eddy/3)`
grows linearly with scale; the virial ratio `σ/v_escape = 0.5` holds *exactly* at
every scale (bound orbits are generic); the Jeans length `λ_Jeans = √(2π²/3) r =
2.565 r` is self-similar (every scale marginally unstable at once).

**Solar-system confirmation:** the multi-scale observable `√β(x_local)·√β(x_galactic)·
√β(x_cosmic)` is 1.000000 for the local and galactic terms and 0.830 for the cosmic
term — the formula is correctly dormant at solar-system scales, recovering GR
exactly.

**The connection to Task 4 (v6.3).** The derived field equation gives the *local*
effective density `ρ_eff = m′(r)/4πr²` from the flow. The open dark-matter question —
whether `ρ_eddy = x₀·ρ_crit` emerges from the 4D kinetic energy of the rotating
hypersurface — is the *homogeneous* version of exactly this calculation. It is now
well-posed rather than abstract.

**The N-body wall (unchanged):** `v_flat = 220 km/s` requires internal halo
overdensity `δ ~ 50,000–100,000 × ρ_eddy`, which comes only from virialization in a
simulation. The falsifiable prediction — ESTIF halos reach that δ — is a collaboration
target, not a limitation of the framework. See `test_collisionless_eddy.py`,
`test_solar_system_eddy.py`.

---

## The Two Paths (v6.3)

The results above split the project into two tracks. See
`MILESTONE_v6.3_THE_SPLIT.md`.

### Path One — ESTIF-Core (clean) ✅ recommended default
- Gravity on the *derived* field equation (Task 4).
- Cosmology = frozen cosmic eddy → cosmological constant (ties ΛCDM, χ²/N = 1.92).
- The Ω_tilt(z) apparatus retired to an appendix.
- Remaining writing already done in v6.3: axioms A2 and A3 written into this
  document; the shrinking-ruler narrative retired; v_flow = cx₀ relabelled as the
  sideways component of a total-*c* flow.
- Deliverable: the gravity letter (unaffected, and strengthened by the derived
  field equation), plus a short honest cosmology statement. Defensible today.

### Path Two — ESTIF-Extended (hard) 🔬 high-risk research
- Derive the *leading correction* to `w = −1` from the full rotating-shear /
  vorticity stress tensor — the cosmological half of the T_μν work, aimed correctly
  (a small perturbation on top of the frozen limit, not a new dominant term).
- Target the mild DESI thawing (χ²/N ≈ 0.66) ⚠️ under-marginalized, see C4. The two naive reductions are already 
  falsified; the full off-diagonal stress tensor is required.

---

## The Honest Open Questions

**Answered in v6.3:**
- *Is the gravity force law derived or borrowed?* → **Derived, with one qualifier.**
  The vacuum condition forces `v² = 2GM/r` (exact Schwarzschild); mass continuity =
  Poisson follows from the flow metric. No Schwarzschild match. The *coupling* to
  matter (8πG) is adopted from Einstein–Hilbert, not derived (C6).
- *Do the flow axioms hold Lorentzian time and SR?* → Yes; signature and SR kinematics
  emerge from a Euclidean bulk plus the speed-*c* constraint.
- *What is ESTIF's honest cosmology?* → A frozen cosmic eddy = cosmological constant,
  tying ΛCDM on DESI DR2.

**Still open:**
- **Strong-field pressure/stress sector** of gravity (full off-diagonal T_μν). Needed
  for none of vacuum/Newton/Schwarzschild, but the remaining rigorous gravity step.
- **Path Two cosmology:** derive the leading `w(z)` correction from the vorticity
  stress tensor. The naive reductions (stiff spin; expansion tracker) are falsified.
- **Why Ωm = x₀?** The homogeneous version of the Task 4 calculation — whether
  `ρ_eddy = x₀ρ_crit` emerges — is the central dark-matter target, now well-posed.
- **Why x_c = 0.272?** The geometric origin of the GR crossover; closing it completes
  the N_MAX = 5/7 × L chain.
- **What is the 4th dimension?** Geometrically well-defined (the direction the flow
  moves through), with no identified counterpart in known field theories.
- **Dark-matter halos:** N-body simulation (collaboration target).

---

## Summary

ESTIF begins from three statements about a moving sheet: 3D space is flat and carried
through a 4D bulk (A1); everything moves through the bulk at the speed of light, so
that time is motion through the 4th dimension (A2); and empty space sources nothing
(A3).

From these, v6.3 *derives* gravity rather than borrowing it: the flow metric forces
mass continuity — Poisson's equation — with the vacuum condition giving exact
Schwarzschild and a matter source giving exactly Newton. Lorentzian time and special
relativity fall out of the same construction. This closes the framework's longest-
standing gap, the previously unstated matching of the force law to the Schwarzschild
solution.

The strong-field tilt formula — dynamic n, the √β projection, the electron-scale
parameters — remains the tool for deviations from GR (EHT, LISA) and for the MOND
acceleration a₀ = H₀cx₀/√3 (1.72% from empirical, validated against 87 SPARC
galaxies at RMS 15.6%), now underwritten by the derived force law.

Cosmology is reframed honestly. The evolving Ω_tilt(z) dark-energy law failed DESI
DR2 (χ²/N = 10.8); de-circularizing it helped (3.35) but did not reach ΛCDM, and the
frozen-eddy limit — a plain cosmological constant, derived for free from the gravity
sector — ties ΛCDM at 1.92 and beats the tilt. So the tilt cosmology is set aside,
and the project splits: **Path One** publishes the derived gravity plus the frozen-
eddy cosmological constant; **Path Two** attempts to derive the small thawing
correction DESI hints at, from the full vorticity stress tensor.

The ant knows about the ball. ESTIF is about the room the ball moves through — and in
v6.3, about the one fact that ties the whole picture together: we move through that
room at the speed of light, and that motion is what we feel as time, and as gravity.

---

**Related documents:**
- `MILESTONE_v6.3_THE_SPLIT.md` — the v6.3 split and honest status
- `README.md` — project overview and key results
- `docs/report/STATUS.md` — current validation status
- `docs/SUMMARY_FOR_REVIEW.md` — one-page expert summary

**Key v6.3 derivation scripts:**
- `tests/estif_task4_field_equation.py` — field equation derived (mass continuity = Poisson)
- `tests/estif_flow_signature_dynamics.py` — signature + SR + Newton from the flow axioms
- `tests/estif_converse_flow_law.py` — vacuum forces v²=2A/r (Birkhoff in flow variables)
- `tests/estif_tmunu_gauss_codazzi.py` — ADM engine, validated vs flat FRW + de Sitter
- `tests/estif_task5_desi_selfconsistent.py` — de-circularized DESI test (10.8 → 3.35)
- `tests/estif_task5b_cosmo_eos.py` — DESI-preferred w(z); tilt tracks it
- `tests/estif_task6_eddy_eos.py` — frozen-eddy reframe (1.92 ties ΛCDM, beats tilt)
- `tests/estif_fidelity_audit.py` — axiom-presence audit of the corpus

**Prior test scripts** (gravity sector, unchanged): `derive_mond_from_geometry.py`,
`test_sparc_tully_fisher.py`, `test_joint_calibration.py`, `test_gravity_time_connection.py`,
`test_electron_connection.py`, `test_collisionless_eddy.py`, `test_solar_system_eddy.py`,
`test_a0_redshift.py`, `test_a0_parameter_independence.py`, and the cosmology-sector
scripts now superseded by the reframe (`test_desi_wz_consistency.py`,
`test_hubble_tilt_cosmology.py`, `test_nmax_drift.py`, the Pantheon+ fits).
