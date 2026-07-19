# SPEC — The JWST Early-Structure-Formation Test

**Purpose:** the precise specification of the test that could distinguish ESTIF
from ΛCDM using JWST's early massive galaxies — the "timing/mass" argument, not
the "apparent size" argument (which dies to the identical-H(z) tie).
**Status:** specification + laptop-runnable ΛCDM baseline. The distinguishing
ESTIF calculation is gated on one open physics question (below).
**Target location in repo:** `docs/plan/JWST_TIMING_TEST_SPEC.md`
**Companion script:** `tests/estif_jwst_halo_massfunction.py`

---

## 1. Why the apparent-size version fails, and the timing version doesn't

**Size version (dead):** "early galaxies look too big, but they were normal size
— there's just more space now." Apparent angular size vs distance is fixed by the
angular-diameter distance, which depends only on H(z). ESTIF-Core has the
*identical* H(z) to ΛCDM (frozen eddy = cosmological constant), so it predicts the
*same* apparent sizes. Reinterpreting expansion changes the words, not the sizes.
This version cannot break the tie.

**Timing/mass version (live):** the JWST tension is not about sizes — it is that
galaxies are too *massive* too *early*. Under ΛCDM the dark-matter halo mass
function sets an absolute ceiling on the stellar mass that can exist by z ≈ 9,
because M⋆ ≤ ε · f_b · M_halo with ε ≤ 1. The observed masses push ε to
implausible values. This is a claim about how fast structure assembles — a regime
where ESTIF can, in principle, differ from ΛCDM. This is the version worth pursuing.

---

## 2. The exact observable

**Primary:** the cumulative comoving stellar-mass density in massive galaxies at
high redshift,

```
rho_star(> M_star, z)   at z ≈ 7.5 and z ≈ 9.1
```

**Equivalent (cleaner, baryon-physics-independent):** the required star-formation
efficiency ε that a given cosmology needs to produce the observed rho_star:

```
epsilon_required(z) = rho_star_observed(> M_star, z) / [ f_b · rho_halo(> M_halo, z) ]
```

where `rho_halo(> M_halo, z)` is the mass density in halos above the threshold
`M_halo = M_star / (f_b · ε)`, computed from the halo mass function, and
`f_b = Omega_b / Omega_m ≈ 0.157` is the cosmic baryon fraction.

A cosmology is in tension if it requires ε implausibly close to 1. A cosmology
*eases* the tension if it raises `rho_halo(> M_halo, z)` (more massive halos
earlier), lowering ε_required.

---

## 3. The ΛCDM number to beat (the target)

From Boylan-Kolchin (2023, Nature Astronomy 7, 731) using Labbé et al. (2023) data
and Planck 2020:

| Redshift | ΛCDM required ε | Plausible? |
|---|---|---|
| z ≈ 7.5 | **> 0.57** | strained (normal ε ~ 0.1–0.2) |
| z ≈ 9.1 | **≈ 0.99** | implausible (ε = 1 is the hard ceiling) |

**The benchmark alternative (Early Dark Energy):** boosts the z ≈ 9.1 baryonic
reservoir by a factor **3.3×**, dropping required ε from 0.99 to **0.72**.

**Therefore the ESTIF target is explicit:**
> To materially resolve the tension, ESTIF must boost the high-z massive-halo
> reservoir by MORE than EDE's 3.3× — enough to bring ε_required(z ≈ 9.1) down to
> a plausible value (say ε < 0.3, i.e. a reservoir boost of roughly ≥ 5× over
> ΛCDM). Matching EDE (3.3×, ε → 0.72) would be interesting but not decisive;
> beating it would be a genuine, distinguishing win.

**Honesty caveat (must be stated in any writeup):** the tension is contested. Some
2025 simulations (e.g. COLIBRE) find ΛCDM consistent with JWST once baryon physics
is modeled fully, and the whole argument depends on stellar-mass estimates and the
IMF. So a successful ESTIF result would be "ESTIF also resolves it, more naturally,"
not "ESTIF rescues a definitively broken ΛCDM."

---

## 4. The ESTIF inputs — what the calculation needs

The abundance of rare massive halos at high z is exponentially sensitive to the
amplitude of density fluctuations `sigma(M, z)`. The halo mass function needs three
inputs; the first is standard, the second and third are where ESTIF could differ:

**(I) The linear power spectrum P(k)** — standard, from the transfer function
(Eisenstein–Hu). Identical for ESTIF-Core and ΛCDM (same early universe, same
Omega_m = x0 ≈ Omega_m).

**(II) The linear growth factor D(z)** — how a small overdensity grows with time.
Governed by
```
D'' + 2 H(z) D' − 4 π G rho_eff D = 0
```
- Under **ESTIF-Core (strict)**: H(z) = ΛCDM and gravity = standard (Task 4 gives
  the ordinary Poisson equation). Then D(z) is IDENTICAL to ΛCDM. **No boost. No
  resolution.** — see the gate in §5.
- A boost requires an EXTRA gravitational source or a modified growth term (§5).

**(III) The collapse threshold delta_c** — the linear overdensity at which a region
collapses (≈ 1.686 in ΛCDM, from spherical collapse). If ESTIF's collapse dynamics
differ (e.g. the eddy background accelerates collapse), delta_c is lower, which
raises massive-halo abundance. The ESTIF spherical-collapse calculation with the
eddy background is the concrete sub-task here — and it is laptop-analytic (a 1D ODE),
not N-body.

---

## 5. THE GATE — the one open physics question this test rests on

**Strict Path One does not resolve the JWST tension.** If ESTIF-Core is exactly
ΛCDM expansion + standard gravity, then D(z), delta_c, and the halo mass function
are all identical to ΛCDM, and ε_required is unchanged. Reinterpreting dark energy
as a frozen eddy changes nothing about structure growth.

**A resolution requires ESTIF to deviate from ΛCDM in the growth sector**, via one
of three candidate mechanisms — each of which is a real, statable ESTIF claim that
must be derived:

1. **Enhanced eddy clustering.** If the cosmic eddy is not exactly pressureless
   cold dark matter but clusters more efficiently at high z (e.g. it has a
   head-start from being a coherent flow rather than a particle gas), it boosts the
   effective source in the growth equation. This is the *homogeneous-source* face
   of the same physics as the rotation-curve/halo problem — the same N-body wall.
2. **Lower collapse threshold.** If the eddy background lowers delta_c (faster
   collapse — consistent with the ~1 Gyr free-fall time at z = 10 already noted in
   the ESTIF docs), massive-halo abundance rises exponentially. This is the most
   accessible lever: a spherical-collapse ODE with the eddy term. Laptop-doable.
3. **A different early clock (Path Two).** If the vorticity derivation gives a
   genuinely different H(z) at high z (the thawing), the growth history and the
   available cosmic time both change. This is the only mechanism that is *uniquely*
   the falling-through-4D picture rather than generic modified gravity.

**Where each mechanism lives:**

| Mechanism | Which path | Computability |
|---|---|---|
| (2) lower delta_c via eddy collapse | Path One-adjacent (analytic) | **Laptop** — spherical-collapse ODE |
| (1) enhanced eddy clustering | dark-matter sector (same wall) | N-body for full result; analytic estimate possible |
| (3) different high-z H(z) | Path Two | needs the vorticity w(z) derivation first |

**Honest summary of the gate:** the JWST timing test is *not* free under Path One.
It requires ESTIF structure growth to differ from ΛCDM, which is a claim beyond
"matches ΛCDM." The most tractable version — mechanism (2), a lower collapse
threshold from the eddy background — is a laptop calculation and is the right first
step. If it produces a ≥ 5× reservoir boost, ESTIF resolves JWST more naturally
than EDE. If it produces < 3.3×, EDE already does better and the argument is weak.

---

## 6. The decision criterion

Run the ΛCDM baseline (companion script) to fix ε_required(z) under standard
assumptions and confirm it lands near the Boylan-Kolchin benchmark (ε ~ 0.57 at
z ≈ 7.5, ~1 at z ≈ 9.1). Then compute the ESTIF halo mass function with the
mechanism-(2) collapse threshold (and, later, mechanisms 1 and 3), and read off:

```
boost = rho_halo_ESTIF(> M_halo, z≈9.1) / rho_halo_LCDM(> M_halo, z≈9.1)
epsilon_ESTIF(z≈9.1) = epsilon_LCDM(z≈9.1) / boost   (approx, for the reservoir)
```

- **boost ≥ 5× (ε → < 0.3):** decisive — ESTIF resolves JWST, beats EDE. Flagship result.
- **3.3× ≤ boost < 5×:** ESTIF matches/edges EDE. Interesting, publishable as
  "ESTIF naturally eases the tension."
- **boost < 3.3×:** weak — EDE and others do better; do not lead with this.
- **boost ≈ 1× (mechanism (2) gives nothing):** the tension needs mechanism (1) or
  (3); the argument moves fully to Path Two + N-body.

---

## 7. What the companion script does (and does not) do

**Does (laptop, now):**
- Computes the ΛCDM halo mass function (Press–Schechter / Sheth–Tormen) with an
  Eisenstein–Hu transfer function, Planck parameters.
- Computes the available baryon reservoir and ε_required(z) at z ≈ 7.5, 9.1;
  checks it lands near the Boylan-Kolchin benchmark.
- Runs the SENSITIVITY analysis: how large a growth/abundance boost is needed to
  bring ε down to plausible values — i.e. quantifies the ESTIF target precisely.
- Provides a clearly-marked slot for ESTIF's `delta_c` and `D(z)` so the
  distinguishing calculation drops in once mechanism (2) or (3) is derived.

**Does not (gated):**
- It does not compute the ESTIF boost itself — that needs the eddy spherical-
  collapse delta_c (mechanism 2, next laptop step) or the Path Two H(z), neither of
  which is yet derived. The script isolates exactly what those must deliver.

---

**Spec version:** 1.0 (8 July 2026) — for ESTIF v6.3
**References:** Boylan-Kolchin 2023 (Nat. Astron. 7, 731); Labbé et al. 2023
(Nature 616, 266); Menci et al. 2022 (ApJL 938, L5, EDE constraint); COLIBRE /
Chaikin et al. 2025 (ΛCDM-consistent counterview).
