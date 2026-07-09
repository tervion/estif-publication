# THE PATH ONE CHECKLIST

**Scope:** ESTIF-Core (Path One) — goals, sub-goals, milestones, and current status.
**Status date:** 8 July 2026 (v6.3 + v6.3.1 errata)
**Companions:** `MILESTONE_v6.3_THE_SPLIT.md` · `CORRECTIONS_v6.3.1.md` · `docs/plan/RHAC.md`

> Maintenance rule (RHAC discipline): update status here when a scenario closes;
> never delete items — mark them ✅ / ❌ and date them.

**Legend:** ✅ complete · 🔶 partial / downgraded · ⬜ open · ❌ externally blocked

---

## The six barriers (from the adversarial review), with effort class

| # | Barrier | Fix | Effort |
|---|---|---|---|
| 1 | Internal contradiction: exact-Schwarzschild vacuum vs EHT/LISA deviation claims | Apply C1 (deviations → conditional on eddy-stress sector) | days (written; needs applying) |
| 2 | "Derived, not borrowed" overclaims the coupling | Apply C6 (axioms select GR's constraint sector in PG gauge; coupling adopted) | days (written; needs applying) |
| 3 | Ωm = x₀ hides an Ωm-dependent horizon integral | Apply C2 (downgrade to consistency relation) | days (written; needs applying) |
| 4 | Pressure/interior sector (relativistic star, TOV-level) | New derivation | weeks |
| 5 | Preferred-frame nature of A2 (absolute bulk time; simultaneity reinterpreted, not derived) | New explicit section in the paper | weeks (writing) |
| 6 | Peer review — the only external test yet faced | Submit the gravity letter | process |

**Two highest-value open derivations:**

1. **D_ESTIF(z)**, the linear growth factor (perturbed field equation on FRW).
   Laptop-analytic. Feeds three checklist items at once: the dark-matter mass
   function (B6), the CMB anisotropy prediction (A5), and the JWST distinguishing
   test (A6). Must pass the two-sided σ8 filter (C3). ⚠️ Note that Task 4 returns
   ρ_eff = ρ₀ *exactly* for a uniform source — the derived field equation is standard
   Poisson, so G_eff > G cannot come from it as it stands. Any enhancement must arise
   from the generalized (non-single-speed) flow.
2. **Principle P** (B-4c). If derived, it converts the Ωm bootstrap from conditional
   to a genuine zero-parameter prediction of Ωm, Ω_Λ, r_universe, and a₀ from three
   measured inputs — the single largest available upgrade to Path One's claims.

---

## GOAL C — Gravity = Time  *(most advanced)*

| # | Item | Status |
|---|---|---|
| C-1 | A2 stated as law: dw² + dσ² = c²dt²; proper time c·dτ = dw | ✅ |
| C-2 | Lorentzian signature (−,+,+,+) derived from Euclidean bulk + A2 (minus sign produced, not inserted) | ✅ |
| C-3 | Exact SR recovered: time dilation, null photons, rest clocks — zero free parameters | ✅ |
| C-4 | Static clock rate = local flow speed u(r); gravitational pull = −c²∇u (engine, Christoffels) | ✅ |
| C-5 | Field equation derived: ρ_eff = m′(r)/(4πr²) — mass continuity = Poisson (Task 4, 5/5) | ✅ |
| C-6 | Vacuum (A3) ⇒ v² = 2GM/r unique ⇒ **exact Schwarzschild** (Birkhoff in flow variables) | ✅ |
| C-7 | Uniform-density ball ⇒ ρ_eff = ρ₀ exactly (correct Newtonian source) | ✅ |
| C-8 | Axioms A1–A3 written into the theory; shrinking-ruler narrative retired (v6.3 concept rewrite) | ✅ |
| C-9 | Apply C1 in repo: strong-field deviation claims (EHT/LISA) → conditional on eddy-stress sector | ⬜ |
| C-10 | Apply C6 in repo: precise "derived" wording (coupling to matter is adopted, not derived) | ⬜ |
| C-11 | Preferred-frame / relativity-of-simultaneity section (barrier 5) | ⬜ |
| C-12 | Pressure/interior sector — TOV-level relativistic star (barrier 4; not needed for vacuum claims) | ⬜ |
| C-13 | Moving-observer bulk identity inside eddies (flagged caveat from the signature suite) | ⬜ |
| C-14 | Peer review of the gravity letter (barrier 6 — the external test) | ⬜ |

## GOAL B — No dark matter  *(split down the middle)*

| # | Item | Status |
|---|---|---|
| B-1 | a₀ = H₀cx₀/√3 derived (1.72% from empirical, zero free parameters); Step 1 now underwritten by the derived field equation | ✅ |
| B-2 | SPARC validation: 87 quality-1 galaxies, RMS = 15.6% (within BTFR scatter) | ✅ |
| B-3 | a₀ redshift constancy (H(z) cancels; 2.22×10⁻¹⁶) + parameter independence (3,600 H₀/Ωm combos; 8 datasets) | ✅ |
| B-4 | Ωm = x₀ — **downgraded (C2)** to a consistency relation, then **partially recovered (Door 2)**: principle P ("Ωm = R_H/r_p") inverts the circularity into the closed equation Ωm·I(Ωm) = 1 with a unique root — 0.3043 with zero inputs, **0.31408** with radiation (measured T_CMB, N_eff, h): **0.96% from Planck, 0.53σ inside its error bar**. Closure: a₀ improves 1.72% → **0.66%**; r_universe back-predicted to −1.07%; DESI DR2 χ²/N = **1.618** (fixed-ruler). Ledger after adopting P: measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ, x₀, r_universe, a₀}. **Conditional on P — Part B (deriving P from A1–A3) is open; three candidate routes, none attempted.** | 🔶 |
| B-4b | Gaztañaga comparison memo — the causal-universe scale (≈ 0.3176 H₀, via inflation) is adjacent to the bootstrap root. **Prerequisite for any Ωm novelty claim.** | ⬜ |
| B-4c | Part B: derive principle P from the flow axioms. Routes: (i) flow-budget amplitude; (ii) horizon-acceleration balance (P ⇔ g_horizon = cH₀/2 — the same cH₀ that sets a₀, ratio computes to 1.00000); (iii) homogeneous field equation (needs the vorticity attachment shared with Path Two). | ⬜ |
| B-5 | ρ_eddy = x₀ρ_crit from the homogeneous field equation (well-posed since Task 4) | ⬜ |
| B-6 | D_ESTIF(z) linear growth — laptop-analytic; must pass the two-sided σ8/S8 filter (C3: transient or scale-dependent) | ⬜ |
| B-7 | Weak-lensing BTFR check (KiDS flat circular velocities to ~1 Mpc) against the derived a₀ — **lowest-hanging fruit on the tree** | ⬜ |
| B-8 | Bullet Cluster response — hardest external test for any no-dark-matter framework; unaddressed | ⬜ |
| B-9 | Halo internal structure (v_flat = 220 km/s, δ ~ 10⁵) — N-body simulation | ❌ blocked (compute/collaborator) |

## GOAL A — "Expansion is really 4D infall"  *(mechanism done; proof structurally hard)*

| # | Item | Status |
|---|---|---|
| A-1 | Mechanism as an equation: H = c·d(ln b)/dw — expansion is motion through the bulk gradient (engine-verified Friedmann analog) | ✅ |
| A-2 | Frozen eddy ⇒ w = −1 ⇒ ties ΛCDM on DESI DR2 (χ²/N = 1.92; with Ωm = x₀: 1.97; ΔAIC = +0.59, indistinguishable) | ✅ |
| A-3 | One-flow unification: v² = 2GM/r + H²r² — exact Schwarzschild–de Sitter, gravity eddy + Hubble flow in a single field | ✅ |
| A-4 | Circularity removed (10.8 → 3.35, self-consistent x(z)); tilt cosmology retired (ΔAIC = +18 vs frozen eddy, decisive) | ✅ |
| A-5 | CMB anisotropy spectrum as an ESTIF prediction (needs D_ESTIF(z) + Boltzmann code). Background CMB is **inherited** — identical H(z) ⇒ recombination z ≈ 1100, ~380,000 yr, rd = 147.09 Mpc (Planck-calibrated import) pass automatically | ⬜ |
| A-6 | **The distinguishing prediction.** Identical H(z) ⇒ infall vs expansion cannot be separated on background data. Live routes: two-sided growth signature (JWST spec, after C3) or a Path Two thawing (after the C4 re-derived target) | ⬜ |
| A-R | ~97,300-year recombination age — **retired with exponential shrinkage; EC-fork only; not a Path One prediction** | ❌ retired by design |

## Cross-cutting — publication track

| # | Item | Status |
|---|---|---|
| P-1 | v6.3 documentation suite + v6.3.1 errata + canonical Path One scripts (DESI parity; AIC/BIC) delivered | ✅ |
| P-2 | Apply the errata (C1–C6) across the repo; commit as v6.3.1 | ✅ 9 July 2026 |
| P-4 | Apply the bootstrap amendments across the repo; commit as v6.3.2 | ⬜ |
| P-3 | Gravity-letter final pass (C1/C6 wording + preferred-frame section) → submit (barrier 6) | ⬜ |

---

## Tally (34 tracked items)

| ✅ complete | 🔶 partial | ⬜ open | ❌ blocked/retired |
|---|---|---|---|
| **16** | **1** | **15** | **2** |

**Pattern:** the derivation core of Path One is finished (all of Goal C's physics,
the a₀ chain, the frozen-eddy parity). What remains is documentation hygiene
(P-2), two referee-facing sections (C-11, C-12), a handful of cheap undone tests
(B-7 is the cheapest), the growth derivation that unlocks three items at once
(B-6/A-5/A-6), and the external gauntlet (C-14).

---

## Lineage note

v6.3 is **not** part of the FD or EC forks. Family tree: ESTIF-FD v1.0 (2024,
exponential shrinkage; ruled out by SNe) → sibling branches: **EC fork** (June
2025; shrinkage kept, variable H(t), 97,300-yr CMB embraced) and the **main line**
v3.0 → v6.x → **v6.3** (flat-slice flow; shrinkage formally retired; field
equation derived). EC's vortex late-acceleration mechanism corresponds to Task 6's
model E1, falsified against DESI DR2 (χ²/N = 3232). EC's "Sink Law" intuition is
the ancestor of the Task 4 derivation, with the velocity/acceleration conflation
corrected (velocity ∝ r^(−1/2), not 1/r²).

**Checklist version:** 1.0 · **Date:** 8 July 2026
