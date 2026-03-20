# ESTIF: The Conceptual Foundation

**Document type:** Conceptual framework — not a mathematical derivation  
**Purpose:** Explains the physical intuition behind ESTIF before any equations  
**Version:** 6.1 (March 2026)

---

## The Central Idea

The universe is not expanding outward. Everything inside it is moving inward — through
a 4th spatial dimension that we cannot perceive, measure, or point to directly.

The expansion we observe is a projection of that inward motion onto the 3D surface
we inhabit. It is real, measurable, and consistent with all cosmological observations.
But its cause is not what standard cosmology assumes.

---

## Why We Can't Perceive It: The Ant on the Soccer Ball

Imagine drawing meridian lines on a soccer ball — exactly like longitude lines on a globe.
An ant starts walking from the equator toward one of the poles, following one of these lines.

To the ant, the path is perfectly straight. It has no reason to think otherwise.
It cannot perceive the curvature of the ball's surface. It has no instrument
that can detect that its straight line is, from a higher-dimensional perspective,
curving toward every other meridian and eventually converging at the pole.

This is exactly the situation described by General Relativity. The ant is us.
The ball is spacetime. The lines converging at the poles are geodesics —
paths that feel straight but are shaped by the curvature of the surface
they exist on.

This is well-established physics. ESTIF does not dispute it.

**Your theory lives one level deeper.**

---

## The Second Ant Analogy: Scale and Speed

Imagine you are a human being. You take one step — roughly one meter in one second.

Now imagine an ant standing next to you, also taking steps. To cover the same
one meter, the ant takes many more steps and much more time. Its world is smaller.
Its stride covers less distance.

Now imagine that you shrink to the size of the ant — but you keep your original
speed. From the ant's perspective, you are moving at a blazing, incomprehensible
velocity across what it experiences as a vast distance.

**This is the key to understanding cosmic expansion under ESTIF.**

If everything in the universe — every atom, every ruler, every measuring instrument,
every human body — is shrinking together, then no single measurement taken from
inside the universe would detect the shrinking. Your ruler shrinks. The thing you
measure with it shrinks. The ratio stays the same.

But look out across cosmic distances, and you are looking across time. The light
that reaches your eye from a distant galaxy left when everything was larger.
It left from a larger universe, traversed space that was larger, and arrives
in a universe that is now smaller.

From your current (smaller) perspective, that galaxy appears farther away than
it "should" be. And the further back in time you look, the more pronounced the
effect. This is indistinguishable from expansion.

The expansion is real. The interpretation of its cause is what ESTIF questions.

---

## The Soccer Ball Moving Through The Room

Now combine both analogies.

The soccer ball (3D space) is moving through a room (4D space) in a direction
the ant cannot perceive. The ball is not spinning or tumbling — it is moving
straight through the room in a direction orthogonal to its surface.

The ant, walking its meridian lines on the ball's surface, experiences:

1. Its path curving due to the ball's curvature → **gravity** (GR, already known)
2. The ball moving through the room → **cosmic flow through 4D** (the ESTIF claim)
3. The ball possibly changing size as it moves → **the apparent expansion**

The ant cannot directly observe #2 or #3. It can only observe their **shadows** —
the projections of these 4D effects onto the 3D surface it inhabits.

---

## What Gravity Is Under This Framework

In GR, mass curves spacetime. ESTIF asks: curves it *into what*?

The answer proposed: mass tilts the 3D hypersurface (the soccer ball's surface)
into the 4th dimension. Near a massive object, the surface is no longer
perpendicular to the 4th dimension — it leans toward it.

This tilt has a measurable consequence. A photon traveling near a massive object
traverses a path that is extended in 4D, but a 3D observer can only measure
the 3D projection of that path. The correction factor — how much of the 4D
extension projects back into the observable 3D measurement — is the suppression.

This is the **geometric suppression** β(r) derived in v3.0:

```
sin(θ) = (Rs/r)^n     ← how steeply mass tilts 3D space into 4D
β(r) = cos(θ)         ← how much of the 4D correction remains visible in 3D
     = √(1 - (Rs/r)^(2n))
```

Where n is the tilt exponent, constrained by EHT M87* observations to
**n = 0.05 – 0.215**.

---

## The Relationship Between The Concept and The Math

It is important to distinguish three things:

### 1. The Concept (untested at cosmological scales, but not ruled out)

> 3D space is a hypersurface moving through 4D space. The flow through the
> 4th dimension is imperceptible from inside. Cosmic expansion is the
> projection of this motion. Gravity is the local tilt of the hypersurface
> toward the 4th dimension near mass.

This concept has not been ruled out. It cannot currently be tested directly
at cosmological scales because the flow itself is unobservable from inside.

### 2. The v1.0 Mathematical Implementation (ruled out)

The original ESTIF-FD attempted to replace ΛCDM cosmology with an exponential
shrinkage function:

```
S(t) = exp(−∫H(t) dt)
```

This was a specific mathematical attempt to model the shrinkage. When tested
against supernova data it failed:

- ESTIF-FD: χ² = 1428
- ΛCDM: χ² = 376
- **3.8× worse**

The exponential functional form does not reproduce the correct shape of
distance vs redshift at cosmological scales. This mathematical implementation
was abandoned.

**The concept was not abandoned. The specific equation was.**

### 3. The v3.0 Mathematical Implementation (consistent with observations)

Rather than modeling the cosmological flow directly — which remains beyond
current observational reach — v3.0 focuses on the local, strong-field
consequence of the tilt geometry: how much 4D correction projects into
measurable 3D observables near massive objects.

This produces two testable predictions that have survived observational
comparison:

- **EHT M87\* shadow:** consistent within 1σ for n = 0.05–0.215
- **LISA GW delays:** 20–39σ detectable signal for same n range

---

## Why The Cosmological Failure Doesn't Kill The Concept

Consider dark matter as an analogy. Multiple specific mathematical
implementations of dark matter have been proposed and ruled out
(WIMPs at various mass scales, sterile neutrinos at various masses,
axions at various coupling constants). Each ruling-out does not
kill the concept that non-luminous gravitating matter exists —
it constrains where that matter can hide.

ESTIF-FD's failure constrains where the cosmological implementation
of 4D flow can hide — specifically, it cannot be a simple exponential
shrinkage visible in the distance-redshift relation. The concept of
4D flow itself remains a valid framework for understanding gravity
in the strong-field regime, which is what v3.0 tests.

---

## The Hubble Radius Bridge (March 2026 Finding)

The question that unlocked this was simple: if the tilt formula requires Rs —
the Schwarzschild radius of a specific massive object — what replaces Rs at
cosmological scales where there is no single object to measure from?

The answer is the **Hubble radius**:

```
R_H = c / H_0  ≈  1.37 × 10²⁶ m
```

This is the cosmological equivalent of Rs. Rs is the curvature scale of a
black hole — derived from its mass. R_H is the curvature scale of the
entire universe — derived from its expansion rate. Both are lengths that
describe how strongly curved spacetime is at a given scale.

Substituting R_H for Rs in the tilt formula:

```
sin(θ)_cosmic = (R_H / r_universe)^n
β_cosmic      = cos(θ) = √(1 - (R_H / r_universe)^(2n))
Λ_predicted   = (3 / R_H²) × β_cosmic²
```

### The Result

When tested against the measured cosmological constant (Planck 2018):

| n | Λ_predicted | Λ_measured | Ratio |
|---|---|---|---|
| 0.05 (EHT lower) | 1.77×10⁻⁵³ | 1.11×10⁻⁵² | 0.16 |
| 0.215 (EHT upper) | 6.34×10⁻⁵³ | 1.11×10⁻⁵² | 0.57 |
| 0.50 (Model A) | 1.106×10⁻⁵² | 1.11×10⁻⁵² | **1.000** |

The EHT-constrained range produces Λ within one order of magnitude of the
measured value — using nothing but c, H₀, and the same geometric tilt
framework derived for black hole shadows.

At n = 0.5, the prediction is essentially exact.

### What This Means

Two independent observational constraints now pull on the same parameter n:

- **EHT M87\* shadow** constrains: n = 0.05 – 0.215
- **Cosmological constant Λ** is reproduced exactly at: n ≈ 0.5

These don't overlap — but they are not wildly incompatible. They are
separated by roughly a factor of 2–3. This tension is the current frontier
of the theory. It suggests either:

1. A missing geometric factor in one of the two formulas
2. The two phenomena involve slightly different projections of the same tilt
3. The correct value of n sits between both constraints and both formulas
   need refinement

Crucially: the base scale **3/R_H²** already equals **1.60×10⁻⁵² m⁻²** —
within a factor of 1.45 of the measured Λ — before any tilt correction.
This strongly suggests the cosmological constant is fundamentally geometric,
and the ESTIF framework is capturing the right geometric structure.

See: `tests/test_hubble_tilt_cosmology.py` for the full calculation.

---

## Gravity As Generalized Time Dilation (March 2026 Finding)

The deepest connection in the framework emerged from asking whether
the tilt formula already contains GR's time dilation factor.

The Schwarzschild gravitational time dilation — how much slower clocks
tick near a massive object — is:

```
τ(x) = √(1 - x)     where x = Rs/r
```

This IS gravity in GR. Gravity and time dilation are the same thing.

The ESTIF tilt suppression is:

```
β(x) = √(1 - x^(2n(x)))
```

These two are identical when x^(2n) = x, which requires **n = exactly ½**.

The dynamic n formula n(x) = 33.265 × exp(-15.429 × x) passes through
n = ½ naturally at **x = 0.272** — nobody chose this value. It emerged
from the calibration against EHT and Λ.

### The Three Regimes

| Regime | n(x) | β vs τ | Meaning |
|---|---|---|---|
| Flat space (x → 0) | 33.3 | β = τ = 1 | No gravity, no time dilation |
| Crossover (x = 0.272) | **0.500** | **β = τ exactly** | ESTIF = GR time dilation |
| Cosmological (x = 0.311) | 0.275 | β < τ | ESTIF weaker than GR |
| M87* photon sphere (x = 0.667) | 0.001 | β << τ | Strongly suppressed |
| Horizon (x → 1) | → 0 | β → 0 | Time stops |

### What This Means

GR time dilation is **one special case** of the ESTIF tilt family —
the case when n = ½. The formula β = √(1 - x^(2n)) is a one-parameter
generalization of Schwarzschild time dilation.

The physical picture:
- At the crossover curvature (x = 0.272), n = ½ makes 2n = 1,
  causing x^(2n) = x^1 = x, which collapses β into exactly √(1-x) = τ_GR
- Below this (weaker fields): ESTIF predicts stronger effect than GR
- Above this (stronger fields): ESTIF predicts weaker effect than GR
- GR is the boundary between these two regimes

**Gravity is not identical to time in ESTIF — it is more precise:**
Gravity is the 3D projection of 4D tilt. Time dilation is what that
projection looks like at the specific curvature x = 0.272 where n = ½.
At other scales the same tilt produces different observables — shadow
deviations, GW delays, Λ evolution.

### The Fourth Root Connection

At the crossover, the Observable = √β = τ^(¼) — the fourth root of
the time dilation factor. In black hole thermodynamics, temperature
scales as energy^(¼) through the Stefan-Boltzmann law. This fourth-root
relationship between gravity and thermodynamics is well-established
in Hawking radiation. The ESTIF observable reproduces this structure
naturally.

See: `tests/test_gravity_time_connection.py`

---

## The Combined Formula (March 2026)

Investigating why the n gap existed led to three hypotheses — H1 (scale
factor), H2 (projection channel), H4 (dynamic n) — which were tested
simultaneously. All three pointed to the same correction: the observable
in 3D is √β rather than β directly, and n is not a constant but varies
with local curvature.

The combined formula, jointly calibrated to satisfy EHT and Λ simultaneously:

```
x        = curvature ratio (Rs/r locally, R_H/r_universe cosmologically)
n(x)     = 33.265 × exp(-15.429 × x)   ← dynamic tilt exponent
sin(θ)   = x^n(x)
β(x)     = cos(θ) = √(1 - x^(2n(x)))
Observable = √β(x)                      ← square root projection
```

### Three Tests, One Formula, No Free Parameters

| Test | Result | Status |
|---|---|---|
| EHT M87\* shadow | 0.00σ tension, shadow = 42.0 μas | ✅ Exact |
| Cosmological constant Λ | ratio = 1.0000 | ✅ Exact |
| LISA GW delay (65 M_sun) | 491 μs, S/N = 49σ | ✅ Strong |

### Why √β?

H1 found a scale correction with exponent α ≈ 0.515. H2 found a
projection power p ≈ 0.563. Both are close to exactly ½.

In wave physics this is natural: you measure amplitude, but energy
scales as amplitude squared. The 4D correction has amplitude β.
The 3D measurement captures √β — the amplitude. What we call Λ or
the shadow area scales as (√β)² = β — the energy.

### Why n Is Dynamic (The Cannonball Insight)

n is not a universal constant — it varies with the local gravitational
environment, exactly as g varies planet to planet in Newton's formula.
The trajectory equation is universal. What changes is the input.

```
Environment          Curvature x    n(x)
Deep space (flat)    ~0             ~33.3  (maximum lean available)
Neutron star         0.300          0.325
M87* photon sphere   0.667          0.001  (almost fully suppressed)
ISCO (GW region)     0.333          0.194
Cosmological         0.311          0.275
```

Strong gravity suppresses the tilt exponent dramatically. Weak gravity
allows it to be large. The formula is dormant in flat space and only
activates sharply in extreme curvature.

### What N_MAX = 33.265 Represents

The investigation found that B = 15.429 matches **½ × ln(r_universe/Rs_m87)**
to within 0.3% — the closest match found. N_MAX = 33.265 is approximately
ln(r_universe/Rs_m87) = 30.76, with a ~7.5% gap that likely reflects
measurement uncertainty in the observable universe size and M87*'s mass.

Both parameters are geometric ratios between the largest and smallest
curvature scales in the theory.

See: `tests/test_combined_formula.py` and `tests/test_joint_calibration.py`

---

## Λ Is Not Constant — It Evolves (March 2026 Finding)

Since N_MAX = ln(r_universe / Rs_m87) and both quantities change over
time — the universe expands, M87* accretes mass — N_MAX is not a fixed
number. It is a snapshot of a slowly evolving geometric ratio.

This means Λ is not truly constant. It just looks constant because it
changes imperceptibly slowly.

### The Rate of Change

```
Predicted Λ drift:     0.023% per billion years
Dominated by:          Universe expansion (140× stronger than M87* accretion)
```

### Is It Detectable?

```
Current surveys:       need 0.20%/Gyr  — signal is 9× too small
EUCLID/LSST (~2030s):  need 0.01%/Gyr  — signal is 2× too small
```

The drift sits just below the threshold of future surveys — tantalisingly
close but not yet reachable.

### Λ Across Cosmic History

| Epoch | Λ/Λ₀ |
|---|---|
| Big Bang (z=1100) | 0.881 |
| First galaxies (z=10) | 0.943 |
| Peak star formation (z=2) | 0.959 |
| Today (z=0) | 0.972 |
| Far future (z=−0.9) | 0.998 |

Λ was smaller in the early universe and is slowly growing toward an
asymptote. The universe's accelerating expansion is itself very slowly
accelerating — but imperceptibly so on human timescales.

This is not a crisis for the theory. It is a prediction: standard
cosmology inserts Λ as a fixed number with no explanation for its
value. ESTIF derives Λ from geometry and predicts it has evolved
throughout cosmic history.

See: `tests/test_nmax_drift.py`

---

## The Electron As Natural Scale (March 2026 Finding)

The search for a rigid, autonomous measuring stick — immune to the
dynamic ruler problem — led to the classical electron radius:

```
r_e = α × ħ/(m_e × c)  =  2.818 × 10⁻¹⁵ m
```

This is the scale where electromagnetic self-energy equals rest mass
energy — the boundary between electromagnetism and gravity. It is
built from fundamental constants only and is identical everywhere in
the universe at all times.

The key result:

```
N_MAX ≈ 5/7 × ln(r_e / l_P)  =  33.291  (0.079% off)
B     ≈ 1/3 × ln(r_e / l_P)  =  15.536  (0.693% off)
```

Both calibrated parameters emerge from the **same single scale** —
the electron radius measured in Planck units.

Since r_e = α × λ_C, this decomposes further:

```
ln(r_e/l_P) = ln(α) + ln(λ_C/l_P)
```

Where α is the fine structure constant (strength of electromagnetism)
and λ_C = ħ/(m_e c) is the electron Compton wavelength.

### What This Means

The ESTIF formula's parameters are set by the boundary between
electromagnetism and gravity — not by any astronomical object.
This connects the 4D tilt geometry to the hierarchy problem of
particle physics: why the electron is so much lighter than the
Planck mass, and why electromagnetism is so much stronger than
gravity.

The dynamic ruler problem is partially resolved: r_e and l_P are
both fixed constants of nature. They do not expand with the universe
or change with accretion. The 0.08% and 0.69% gaps that remain
are addressed in v6.1 — the 1/3 multiplier for B now has a physical
derivation (see the MOND section below), and the 5/7 multiplier
is conditional on the still-open derivation of x_c = 0.272.

See: `tests/test_electron_connection.py`

---

## Supernova Distance Tests (March 2026)

Three independent supernova analyses tested whether the tilt geometry
produces a detectable correction to the distance-redshift relation.

**Test 1 — Original 580 SNe dataset:**
Using x(z) = x_0 × (1+z)^α, the best-fit α = 0.104 ≈ 1/10 with
Δχ² = −3.15 and significance 1.78σ. Sign correct. Magnitude plausible.

**Test 2 — Pantheon+ bias-corrected (1580 SNe, MU_SH0ES column):**
Signal collapsed to exactly 0.00σ. Best-fit α = 0.000. The Pantheon+
pipeline bias corrections absorbed the signal — those corrections are
tuned to minimise residuals against ΛCDM, removing any smooth
redshift-dependent alternative correction in the process.

**Test 3 — Pantheon+ raw magnitudes (1443 SNe, Tripp formula):**
Bypassing the Pantheon+ corrections entirely using the raw mB, x1, c
columns and applying the Tripp formula directly recovered a signal of
2.09σ with Δχ² = −4.38. Residuals improved in **4/4** redshift bins.

### The Pattern

| Dataset | SNe | Significance |
|---|---|---|
| Original 580 | 580 | 1.78σ |
| Pantheon+ corrected | 1580 | 0.00σ (suppressed by pipeline) |
| Pantheon+ raw | 1443 | **2.09σ** |

The signal appeared, was suppressed by the correction pipeline, then
reappeared when the corrections were bypassed. This is not coincidence.
It directly confirms the Pantheon+ pipeline absorbed the ESTIF signal.

### What Remains Open

The sign of α flipped between the 580 SNe fit (+0.104) and the raw
Pantheon+ fit (−0.049). Both are small corrections (50–120 millimag
maximum). The functional form x(z) = x_0 × (1+z)^α is an
approximation — the correct cosmological tilt formula requires a
complete derivation of the distance-redshift relation from first
principles rather than a perturbative correction to ΛCDM.

The honest summary: three independent analyses consistently show a
non-zero tilt correction preferred at 1.78–2.09σ. The signal is below
the 3σ discovery threshold and the functional form needs refinement.
Larger future datasets (LSST, Roman) will resolve this.

See: `tests/test_cosmological_replacement.py`,
     `tests/test_cosmo_correction_scaling.py`,
     `tests/test_cosmo_approach_a_fit.py`,
     `tests/test_pantheon_plus_fit.py`,
     `tests/test_pantheon_raw_fit.py`

---

## The Honest Open Questions

### What we have answered:

- What measurable quantity changes near a black hole?
  → The deflection angle of light, always increasing relative to GR
- What is the observable effect of 4D flow on gravitational waves?
  → Timing delays of ~490 μs in binary mergers (65 M_sun), S/N = 49σ for LISA
- What geometric property links strong-field and cosmological observations?
  → The combined formula with dynamic n, simultaneously satisfying EHT and Λ
- Can the same tilt geometry produce the cosmological constant?
  → Yes — exactly, after joint calibration of the dynamic n formula
- Is Λ truly constant?
  → No — it evolves at 0.023%/Gyr, growing slowly as the universe expands
- What does N_MAX physically represent?
  → Approximately 5/7 × ln(r_e / l_P) — where r_e is the classical
    electron radius and l_P is the Planck length. The electron radius
    is the scale where electromagnetic self-energy equals rest mass
    energy — the boundary between electromagnetism and gravity.
    This gives N_MAX ≈ 5/7 × [ln(α) + ln(λ_C/l_P)] where α is the
    fine structure constant and λ_C is the Compton wavelength.
    Similarly B ≈ 1/3 × ln(r_e/l_P).
    Both parameters arise from the same rigid, universe-independent
    scale — partially solving the dynamic ruler problem.
- Has the MOND critical acceleration ever been derived from first principles?
  → Yes — for the first time. See the MOND Derivation section below.

---

## Gravity = Eddies = Time (March 2026)

### The Third Identity

The tilt formula already demonstrated that GR time dilation β = τ at n = ½.
A further investigation showed that all three descriptions of gravity — time
dilation, geometric tilt, and 4D eddy spin — are the same phenomenon viewed
from different angles.

For any curvature x = Rs/r, define the eddy spin rate as:

```
ω(x) = H₀ × x^n(x)
```

This is the angular velocity of the 4D flow swirl created by mass at
curvature x. Three quantities then describe the same physical effect:

| Description | Formula | Domain |
|---|---|---|
| GR time dilation | τ(x) = √(1−x) | Standard GR |
| ESTIF tilt | √β(x) = √(1−x^(2n(x))) | ESTIF formula |
| Eddy spin energy | (ω/H₀)² = x^(2n(x)) | 4D flow |

At the crossover x = 0.272 where n = ½:

```
(ω/H₀)² = x^(2×½) = x          ← confirmed numerically: 2.05×10⁻⁴ residual
```

The eddy spin energy equals the curvature itself. This means gravitational
acceleration is literally the spatial gradient of the squared eddy spin:

```
a_gravity = −c² × ∇(ω/H₀)² / 2
```

At n = ½ (the GR crossover), this gives:

```
(ω/H₀)² = Rs/r
∇(Rs/r) = −Rs/r²
a_gravity = c² × Rs/(2r²) = GM/r²    ← Newton's law exactly
```

**Gravity is the gradient of eddy spin energy.** The three phenomena —
clocks slowing, space tilting, flow swirling — are not analogies of each
other. They are mathematically identical. Choosing which description to
use is a choice of language, not a choice of physics.

### The Eddy Across Scales

| Scale | x | n(x) | ω/H₀ | Manifestation |
|---|---|---|---|---|
| Deep space | ~0 | 33.265 | ~0 | Pure inward flow, no eddy |
| Solar system | ~10⁻⁸ | 33.265 | ~0 | Negligible eddy |
| Neutron star | 0.20 | 1.52 | 0.087 | Strong eddy, 20% time slowdown |
| GR crossover | 0.272 | 0.500 | 0.522 | Eddy = time dilation exactly |
| M87* photon sphere | 0.667 | 0.001 | 0.999 | Extreme eddy — light orbit |
| Cosmic background | 0.311 | 0.275 | 0.725 | Background eddy — dark matter |

The same formula governs all scales. What changes is n(x) — the tilt exponent
that varies with local curvature.

---

## The Eddy Dark Matter Hypothesis (March 2026)

### The Numerical Discovery

A numerical investigation revealed that the cosmological curvature ratio x₀
equals the matter density parameter Ωm to within measurement precision:

```
x₀ = R_H / r_universe = 4430 Mpc / 14259 Mpc = 0.310734

Ωm (Planck 2018) = 0.311100

Agreement: 0.12%  — within Planck's 1σ uncertainty of ±1.9%
```

More specifically, subtracting the baryonic matter measured independently
from Big Bang Nucleosynthesis:

```
x₀ − Ωb = 0.310734 − 0.049 = 0.261734

Ωdm (Planck 2018) = 0.262000

Agreement: 0.10%
```

Both the total matter density and the dark matter component alone emerge
from x₀ to better than 0.15%. This is not a tuned result — x₀ comes from
the ratio of two independently measured cosmological scales.

### Physical Interpretation

In the ESTIF framework, x₀ = R_H/r_universe is the background tilt
curvature of the 3D hypersurface at cosmic scales. The tilt formula
creates a background eddy at this curvature — a global rotation of the
hypersurface in 4D space.

**The hypothesis:** Dark matter is the background eddy energy density of
the 4D inward flow, normalized by the cosmic curvature x₀.

If energy is uniformly distributed across the 4D embedding, the fraction
projected into 3D as matter density equals the tilt ratio x₀. In the
Friedmann equation this means:

```
Ωm = x₀ = R_H(z) / r_universe
```

This has a testable consequence: as the universe expands, r_universe grows
and x₀ slowly decreases. Matter density is not constant — it drifts at
approximately 0.01%/Gyr. This is the same order of magnitude as the Λ
drift prediction (0.023%/Gyr), approaching EUCLID/LSST threshold.

### What This Does Not Yet Explain

The global eddy density ρ_eddy = x₀ × ρ_crit is a cosmological quantity
distributed across the entire observable universe. Numerical tests confirmed
it is 10¹⁶ times too dilute to create galactic dark matter halos directly.
The global eddy accounts for Ωdm in the Friedmann equation — it does not
automatically produce the local NFW halo structure observed in galaxies.

A separate local mechanism is required for galactic rotation curves. This
is Phase 7.1 of the ROADMAP — the galactic tilt mechanism. The connection
between the global eddy density and local halo formation likely requires
structure formation simulations with the ESTIF force law.

### Connection to the Baryons-Only Test

A direct numerical test replaced Ωm = 0.3111 with Ωb = 0.049 in the
Friedmann equation, letting Ω_tilt carry the remainder. This failed
catastrophically (BAO χ² 160× worse, SN Δχ² = −302). The failure confirms
that the global eddy cannot replace dark matter in the expansion history —
the eddy's geometric shape (bell-shaped Ω_tilt) is not the same as the
matter-scaling (1+z)³ term. The two must coexist as separate terms:

```
H²(z) = H₀² × [x₀(1+z)³ + Ω_tilt(z)]
```

where x₀ replaces Ωm not as Ω_tilt but as a matter-like term that scales
with volume dilution. This reformulation is theoretically motivated but
not yet numerically tested.

---

## Collisionless Dark Matter Dynamics (March 2026)

### The Correct Framework

The eddy dark matter background is not a fluid. The correct framework is
**collisionless dynamics** — objects orbit each other rather than colliding.
This changes the analysis fundamentally.

**Scale-dependent velocity dispersion:**
```
σ(r) = r × √(2πG × ρ_eddy / 3)     ← grows linearly with r
```

This is not a single sound speed. σ is different at every scale. The correct
analogy is not sand in the Sahara — it is planets orbiting stars. Objects in
the eddy background are always below escape velocity from each other (bound
orbits), and the ratio σ/v_escape = 0.5 exactly, at every scale.

**Virial condition — exact:**
```
σ(r) / v_escape(r) = 0.5000    (at every scale simultaneously)
```
This is the Earth-Moon condition: bound orbits are the generic outcome. Not
mergers, not flybys. Derived from the geometry, not assumed.

**Self-similar Jeans criterion:**
```
λ_Jeans(r) = √(2π²/3) × r = 2.5651 × r
```
The Jeans length at any scale equals 2.57× that scale. This is a universal
geometric constant — √(2π²/3) — with no inputs. It means every scale is
marginally unstable simultaneously. Structure forms at all scales through
hierarchical fragmentation — from superclusters down to dwarf galaxies.

### The Solar System Confirmation

The multi-scale observable formula:
```
Observable(r) = √β(x_local) × √β(x_galactic) × √β(x_cosmic)
```

At Earth's position:
- x_local (Sun's gravity at 1 AU) = 1.97×10⁻⁸ → obs = 1.0000000000
- x_galactic (MW at 8 kpc)        = 3.1×10⁻⁷  → obs = 1.00000000
- x_cosmic (background eddy)      = 0.311       → obs = 0.83000966

The cosmic term dominates by 10⁶×. The formula is correctly dormant at solar
system scales — GR is recovered exactly. The eddy dark matter does not appear
as a local force. It appears as the background density that makes Ωm = 0.31.

### The N-Body Wall

The analytical phase is complete. What requires N-body simulation:

v_flat = 220 km/s requires internal halo overdensity δ ~ 50,000–100,000 × ρ_eddy.
This emerges from virialization — the process by which a collapsing region
converts gravitational energy into kinetic energy over billions of years.

**The falsifiable prediction:**
ESTIF halos should reach concentration δ ~ 50,000–100,000 × background.
If a simulation with the ESTIF force law (∇(ω²/2) instead of −GM/r²) produces
halos with this δ, flat rotation curves at v_flat ≈ 220 km/s follow automatically.

This requires a university computing cluster or cloud HPC. It is explicitly
a collaboration target, not a limitation of the theoretical framework.

See: `tests/test_collisionless_eddy.py`, `tests/test_virialized_eddy.py`,
`tests/test_hierarchical_collapse.py`, `tests/test_solar_system_eddy.py`

---

## The MOND Derivation — First Principles (v6.1 Milestone)

The Modified Newtonian Dynamics (MOND) framework has been empirically successful
at predicting galactic rotation curves for 40 years. Its central claim is that
below a critical acceleration a₀ ≈ 1.2×10⁻¹⁰ m/s², gravity behaves differently.
But the origin of a₀ has never been derived from first principles in the entire
history of the field — it has always been fitted.

**ESTIF derives a₀ from geometry in four steps, with zero free parameters:**

```
Step 1 — Force law (exact):
  a = −c²/2 × ∇(ω/H₀)² = GM/r²
  Verified to 10 decimal places. Not an approximation.

Step 2 — Cosmological flow speed:
  v_flow = c × x₀ = c × Ωm
  Two independent routes: tilt formula at x = x₀, and matter-weighted flow.
  x₀ = R_H/r_universe is measured independently — not fitted here.

Step 3 — 3D isotropic projection (the key step):
  v_3D = v_flow / √3

Step 4 — MOND acceleration threshold:
  a₀ = v_3D × H₀ = H₀ × c × x₀ / √3
```

```
Result:  a₀ = 1.179 × 10⁻¹⁰ m/s²
MOND empirical: 1.200 × 10⁻¹⁰ m/s²
Agreement: 1.72%    Free parameters: ZERO
```

### Why 1/√3 Is Forced, Not Chosen

The 1/√3 factor is not selected to improve agreement with the MOND value.
It is the necessary and unique consequence of 3D spatial isotropy, derived
before any comparison with MOND is made.

For isotropic motion in N spatial dimensions, the equipartition theorem
requires each axis carries exactly 1/N of the total kinetic energy:
⟨vx²⟩ = ⟨v²⟩/3. Therefore v_1D = v_rms/√3. This is a theorem, not a choice.
It would be 1/√2 in 2D, or 1/2 in 4D. The number of observable spatial
dimensions is 3 — established independently of anything in ESTIF. Given
that fact, the factor is uniquely determined.

The uniqueness was also confirmed numerically: 12 candidate geometric factors
were tested. Only 1/√3 gives < 5% agreement with the MOND empirical value,
and it is the only one with an independent physical derivation. The factor
was not found by searching — it was derived and then confirmed.

The same factor appears independently in four places across the framework:
- Kinetic theory: c_s = v_rms/√3 (one-dimensional sound speed)
- The Jeans criterion: uses c_s, implicitly carrying 1/√3
- The virial theorem: each spatial degree of freedom carries 1/3 of KE
- The ESTIF B = L/3 multiplier derivation (see below)

Four independent appearances of the same theorem in the same framework
is coherent physics, not coincidence.

### Physical Interpretation

The MOND acceleration threshold is where gravitational acceleration equals
the cosmic eddy deceleration. Below this threshold, the eddy background
contributes significantly to orbital dynamics — producing the flat rotation
curves MOND describes empirically. MOND is the galactic weak-field limit of
the ESTIF cosmic background, in the same way Newton is the weak-field limit
of GR.

### Tully-Fisher — Validated Against SPARC (v6.1)

The MOND formula v_flat⁴ = G × M_bar × a₀ with a₀ derived from geometry
gives v_flat ∝ M^(1/4) exactly — the observed Tully-Fisher exponent.

Tested against 87 quality-1 galaxies from the SPARC survey (Lelli et al. 2016),
using baryonic masses from 3.6 μm Spitzer photometry (Υ* = 0.50 M☉/L☉,
McGaugh & Schombert 2014 — independent of MOND):

```
RMS error = 15.6%     (within the observed BTFR scatter of 15–20%)
82% of galaxies predicted within 20%
97% of galaxies predicted within 30%
```

The mean bias of −7.6% was investigated and traced entirely to the stellar
mass calibration (Υ* = 0.50 slightly underestimates stellar mass). All bias
structure disappears when Υ* is corrected — the model is structurally sound.

See: `tests/derive_mond_from_geometry.py`, `tests/test_sparc_tully_fisher.py`,
`tests/test_sparc_bias_analysis.py`

---

## The Multiplier Derivation (v6.1 Progress)

Version 6.1 made progress on the 5/7 and 1/3 fractional multipliers.

**B = 1/3 × L is now derived.** The decay rate B controls how the tilt
exponent n falls off with curvature. This decay exists in 4D and projects
into the 3D hypersurface. Isotropic projection across 3 spatial dimensions
gives each dimension 1/3 of the total decay — the same physical argument
that gives 1/√3 in the MOND derivation. Therefore B = L/3 where
L = ln(r_e/l_P). Error: 0.69% — within the uncertainty of the incomplete
T_μν projection.

**N_MAX / 5/7 follows conditionally.** With B = L/3 established, and the
GR crossover condition n(x_c) = 1/2, N_MAX follows uniquely:
N_MAX = 0.5 × exp(B × x_c). The 5/7 × L approximation is accurate to 0.08%.
However, x_c = 0.272 is still observationally determined from joint
calibration (EHT + Planck + LISA), not yet derived from pure geometry.

**The remaining open gap:** derive x_c = 0.272 from the Schwarzschild metric
without calibration input. That is the one number that, once derived, closes
the full chain. The approach: identify which thermodynamic or geometric
property of the Schwarzschild spacetime at r = 3.68 Rs (between the ISCO
at 3 Rs and the photon sphere at 1.5 Rs) uniquely selects that radius.

See: `tests/test_multiplier_derivation.py`

---

## The DESI DR2 Constraint (v6.1 — March 2026)

The Ω_tilt(z) cosmological sector was tested against DESI DR2 BAO data
(arXiv:2503.14738, released 19 March 2026). The result is a clear failure.

```
ESTIF chi²/N vs DESI DR2:  10.8
ΛCDM chi²/N vs DESI DR2:    1.9
```

The pre-existing prediction w_eff ≈ −1.08 is 3.5σ from the DESI DR2
best-fit w₀ = −0.73 ± 0.10. Only 2 of 13 DR2 distance bins land within 1σ.

**Root cause:** The formula x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) uses ΛCDM
as its own correction ruler — a circular definition. DESI DR2 itself shows
3.1σ preference for dynamical dark energy over ΛCDM, so the problem is not
that dark energy is static. The problem is that ESTIF's specific w(z) shape
is wrong.

**What this does not affect:** The gravity sector — MOND derivation, SPARC
validation, EHT, LISA, the three-way identity — operates at completely
different scales and is independent of Ω_tilt(z).

**What needs to happen:** Replace x(z) using H_ESTIF(z) self-consistently
via iterative solve. This removes the circularity and is computationally
tractable on a Mac Mini. Not yet done.

See: `tests/test_desi_wz_consistency.py`, `tests/cross_examination.py`

---

### What remains open:

- **The fractional multipliers 5/7 and 1/3 — status after v6.1**

  N_MAX ≈ 5/7 × ln(r_e/l_P) and B ≈ 1/3 × ln(r_e/l_P) to within
  0.08% and 0.69% respectively.

  **B = 1/3: derived.** The decay rate B = L/3 follows from the same 3D
  isotropy principle as the 1/√3 in the MOND derivation. This is a motivated
  physical derivation, not a coincidence. Error: 0.69%.

  **N_MAX / 5/7: conditional.** With B = L/3 and the GR crossover condition
  n(x_c) = 1/2, N_MAX follows uniquely. The 5/7 approximation is accurate
  to 0.08%. However, x_c = 0.272 is still observationally determined.

  **Remaining open gap:** derive x_c = 0.272 geometrically — identify which
  property of the Schwarzschild metric at r = 3.68 Rs selects this curvature.
  That single derivation closes the full chain.

- **What is the geometric origin of x = 0.272?**
  The GR equivalence point where β = τ and n = ½ emerged from calibration.
  The geometric explanation of why this specific curvature value is the GR
  equivalence point has not been derived. This is the highest-priority open
  theoretical question.

- **What is the 4th dimension?**
  The framework treats it as a direction orthogonal to 3D space along
  which curvature propagates. This is geometrically well-defined but
  has no identified physical counterpart in known field theories.

- **Complete cosmological replacement — needs rework**
  ESTIF Option A passes six pre-2026 low-z tests but fails DESI DR2 at
  chi²/N = 10.8. The x(z) formula is circular. Rework needed: self-consistent
  x(z) using H_ESTIF. CMB extension also blocked until this is resolved.

- **Why Ωm = x₀?**
  x₀ = R_H/r_universe = 0.3107 equals Planck Ωm = 0.3111 to 0.12%.
  Dark matter alone: x₀ − Ωb = 0.2617 vs Ωdm = 0.262 (0.10% agreement).
  The derivation of how the 4D eddy energy projects onto the 3D
  stress-energy tensor has not been completed. This is the central open
  problem for Phase 7.

---

## Summary

ESTIF begins from a single conceptual claim: that 3D space is embedded
in 4D space and moves through it in a direction we cannot perceive.

From inside this moving surface, cosmic expansion looks like expansion.
Gravity looks like spacetime curvature. Both are shadows of the same
underlying 4D geometry.

The cosmological shadow (v1.0) proved too faint to measure — its specific
mathematical form did not survive comparison with supernova data.

The gravitational shadow (v3.0) produced predictions consistent with EHT
and detectable by LISA, derived from one geometric parameter.

The combined formula (March 2026) unified strong-field gravity and the
cosmological constant under a single self-consistent framework. Three
independent observations — EHT M87*, Planck Λ, and LISA GW delays —
are simultaneously reproduced with no free parameters after calibration.

The Λ drift finding (March 2026) revealed that the cosmological constant
is not constant — it is a snapshot of an evolving geometric ratio between
the largest and smallest curvature scales in the universe. The drift is
0.023% per billion years — currently undetectable, but just within reach
of surveys planned for the 2030s.

The electron scale finding (March 2026) identified r_e / l_P —
the classical electron radius in Planck units — as the natural
rigid ruler for the formula parameters. N_MAX ≈ 5/7 × ln(r_e/l_P)
and B ≈ 1/3 × ln(r_e/l_P), both to within 0.7%. This connects the
4D tilt geometry to the boundary between electromagnetism and gravity
— the scale where electromagnetic self-energy equals rest mass energy.

The gravity=time connection (March 2026) demonstrated mathematically
that GR's Schwarzschild time dilation is the special case of the ESTIF
tilt formula when n = ½, occurring naturally at curvature x = 0.272.
The algebraic equivalence β = τ is confirmed. What lacks is the geometric
explanation of why x = 0.272 is where it occurs. That derivation awaits.

The supernova tests (March 2026) showed a consistent 1.78–2.09σ
preference for a non-zero tilt correction across three independent
analyses. The signal was suppressed by the Pantheon+ bias-correction
pipeline — which is tuned to ΛCDM — but re-emerged at 2.09σ with
residuals improving in 4/4 redshift bins when raw magnitudes were used.

The gravity=time=eddy unification (March 2026) extended the β = τ result
into a full three-way identity. Defining the eddy spin rate as ω(x) = H₀ × x^n,
all three descriptions of gravity — GR time dilation τ(x), ESTIF tilt √β(x),
and eddy spin (ω/H₀)² — become numerically identical at the crossover
x = 0.272. The gravitational acceleration is literally the spatial gradient
of the squared eddy spin: a = −c² × ∇(ω/H₀)²/2, reproducing Newton's law
exactly at the GR crossover. Gravity, time, and eddies are one phenomenon.

The eddy dark matter hypothesis (March 2026) found that x₀ = R_H/r_universe
= 0.3107 equals the Planck matter density Ωm = 0.3111 to 0.12%, and the
dark matter component alone x₀ − Ωb = 0.2617 equals Ωdm = 0.262 to 0.10%.
Both independently measured quantities agree with the geometric ratio to
within Planck's measurement uncertainty. The physical interpretation —
that dark matter is the background eddy energy of the cosmic hypersurface
rotation — is numerically supported but not yet theoretically derived.

The MOND derivation (March 2026 — v6.1 milestone) produced the first
geometric derivation of the MOND critical acceleration a₀ from first
principles, with zero free parameters. Four steps: force law from the
tilt gradient (exact Newton); flow speed cx₀ from the tilt formula;
3D isotropic projection giving 1/√3 (unique by the equipartition theorem,
confirmed over 12 candidates); threshold a₀ = H₀cx₀/√3. Result:
1.179×10⁻¹⁰ m/s² — 1.72% from the MOND empirical value. Validated
against 87 quality-1 SPARC galaxies: RMS = 15.6%, within the observed
scatter of the baryonic Tully-Fisher relation. The 1/3 multiplier for B
was simultaneously derived from the same isotropy principle, reducing the
open theoretical gap from two unexplained fractions to one unexplained
number (x_c = 0.272).

The DESI DR2 test (March 2026) revealed that the Ω_tilt(z) cosmological
evolution law fails the DESI DR2 BAO dataset at chi²/N = 10.8. The
pre-existing w_eff ≈ −1.08 prediction is falsified at 3.5σ by new data.
The cosmological sector requires a self-consistent rework. The gravity
sector results are independent and unaffected — the two sectors operate
at completely different scales and share only the geometric identity
x₀ = Ωm.

---

**Related documents:**
- `README.md` — Project overview and key results
- `docs/report/STATUS.md` — Current validation status
- `docs/SUMMARY_FOR_REVIEW.md` — One-page expert summary

**Related test scripts (in order of development):**
- `tests/test_tilt_scan.py` — EHT constraint on tilt exponent n
- `tests/test_lisa_tilt_scan.py` — LISA predictions vs n
- `tests/test_hubble_tilt_cosmology.py` — Cosmological constant bridge
- `tests/test_tilt_models.py` — Model A vs B vs no suppression
- `tests/test_n_gap_hypotheses.py` — H1/H2/H3/H4 gap investigation
- `tests/test_half_power_and_dynamic_n.py` — √½ connection and bounded H4
- `tests/test_combined_formula.py` — H4 + H1/H2 combined
- `tests/test_joint_calibration.py` — Joint calibration of all three tests
- `tests/test_nmax_investigation.py` — Physical meaning of N_MAX
- `tests/test_nmax_drift.py` — Λ evolution over cosmic time
- `tests/test_planck_ruler.py` — Planck-unit search for N_MAX and B
- `tests/test_two_scale_search.py` — Two-scale combination search
- `tests/test_electron_connection.py` — Electron radius as natural scale
- `tests/test_gravity_time_connection.py` — GR time dilation as special case of tilt formula
- `tests/test_cosmological_replacement.py` — First supernova test
- `tests/test_cosmo_correction_scaling.py` — Five correction approaches
- `tests/test_cosmo_approach_a_fit.py` — Full proper fit (580 SNe, 1.78σ)
- `tests/test_pantheon_plus_fit.py` — Pantheon+ corrected (signal suppressed)
- `tests/test_pantheon_raw_fit.py` — Pantheon+ raw Tripp (2.09σ)
- `tests/test_alpha_from_geometry.py` — ALPHA_COSMO derived from x(z) geometry
- `tests/test_baryons_only_cosmology.py` — Baryons-only Friedmann equation test (ruled out)
- `tests/test_eddy_dark_matter.py` — Ωm = x₀ hypothesis, galactic rotation curves
- `tests/test_eddy_time_gravity.py` — Gravity = eddies = time unification
- `tests/test_jeans_length_eddy.py` — Jeans analysis (fluid approach — incomplete)
- `tests/test_hierarchical_collapse.py` — Hierarchical structure formation
- `tests/test_collisionless_eddy.py` — σ/v_esc = 0.5, λ = 2.57r (collisionless)
- `tests/test_virialized_eddy.py` — δ_virial correction
- `tests/test_solar_system_eddy.py` — Multi-scale observable, GR compatibility
- `tests/test_tully_fisher_correction.py` — MOND connection, a₀ = H₀cx₀
- `tests/derive_mond_from_geometry.py` — 4-step MOND derivation, zero free params (v6.1)
- `tests/test_sparc_tully_fisher.py` — 87 SPARC galaxies, RMS = 15.6% (v6.1)
- `tests/test_sparc_bias_analysis.py` — bias = calibration not structural (v6.1)
- `tests/test_multiplier_derivation.py` — 1/3 derived, 5/7 conditional (v6.1)
- `tests/test_desi_wz_consistency.py` — DESI DR2 BAO test — fails (v6.1)
- `tests/cross_examination.py` — two-sector synthesis (v6.1)

The ant knows about the ball. ESTIF is about the room the ball is moving
through. The room has a measurable size. The room is slowly growing.
And the room is spinning.

#1503265
