# THE PATH ONE CHECKLIST

**Scope:** ESTIF-Core (Path One) — goals, sub-goals, milestones, and current status.
**Status date:** 12 July 2026 (v6.4.0 — Fronts 1–3, C-15 closed, P non-derivable)
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

1. **D_ESTIF(z) — RESOLVED (11 July 2026), two-stage.** Under strict A1 the growth
   factor is forced to δ = H(a)/H₀ — decaying only (no-go theorem, RHAC-005).
   Under A1′ (RHAC-006) the deepen mode is restored: D₊ = H·∫da/(aH)³, which
   reproduces GR linear growth exactly — f(z=0.5) = 0.76, matching DESI RSD.
   Consequence for C3/JWST: no ESTIF-specific enhancement exists at linear order;
   the JWST tension is inherited from ΛCDM, not relieved. B-6 closed; A-5 unblocked.
2. ~~Principle P (B-4c).~~ **CLOSED 12 July 2026 (RHAC-009):** P is not derivable
   as a law — it holds only at our epoch. Ωm's status is "fixed by predictive
   postulate", not "derived". The tractable foundational targets are now
   x_c = 0.272 (Schwarzschild geometry) and the √3/cH₀ acceleration scale.

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
| C-15 | GW sector under A1′: c_gw = c DERIVED — one Lorentzian geometry, GW & light share the null cone (principal symbol = light cone), shown for arbitrary flow; GW170817 passed structurally. RHAC-008 | ✅ 12 Jul 2026 |

## GOAL B — No dark matter  *(split down the middle)*

| # | Item | Status |
|---|---|---|
| B-1 | a₀ = H₀cx₀/√3 derived (1.72% from empirical, zero free parameters); Step 1 now underwritten by the derived field equation | ✅ |
| B-2 | SPARC validation: 87 quality-1 galaxies, RMS = 15.6% (within BTFR scatter) | ✅ |
| B-3 | a₀ redshift constancy (H(z) cancels; 2.22×10⁻¹⁶) + parameter independence (3,600 H₀/Ωm combos; 8 datasets) | ✅ |
| B-4 | Ωm = x₀ — **downgraded (C2)** to a consistency relation, then **partially recovered (Door 2)**: principle P ("Ωm = R_H/r_p") inverts the circularity into the closed equation Ωm·I(Ωm) = 1 with a unique root — 0.3043 with zero inputs, **0.31408** with radiation (measured T_CMB, N_eff, h): **0.96% from Planck, 0.53σ inside its error bar**. Closure: a₀ improves 1.72% → **0.66%**; r_universe back-predicted to −1.07%; DESI DR2 χ²/N = **1.618** (fixed-ruler). Ledger after adopting P: measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ, x₀, r_universe, a₀}. **Conditional on P — Part B (deriving P from A1–A3) is open; three candidate routes, none attempted.** | 🔶 |
| B-4b | Gaztañaga comparison memo — the causal-universe scale (≈ 0.3176 H₀, via inflation) is adjacent to the bootstrap root. **Prerequisite for any Ωm novelty claim.** | ⬜ |
| B-4c | Part B: derive P from the axioms — **CLOSED 12 Jul 2026 (RHAC-009).** Obstruction: P's equality Ωm(a)=R_H/d_p holds only at a≈1 (crosses once; →1 vs ½ at high z), so no time-symmetric-axiom derivation exists. Route (ii)'s g_horizon=cH₀/2 IS that today-condition. Escapes closed (attractor breaks Phase-2; anthropic can't hit 0.3141). P = predictive postulate. | ✅ closed |
| B-5 | ρ_eddy = x₀ρ_crit from the homogeneous field equation (well-posed since Task 4) | ⬜ |
| B-6 | D_ESTIF(z) linear growth — ✅ derived 11 July 2026 (strict: decaying-only no-go, RHAC-005; A1′: D₊ = H·∫da/(aH)³ = GR growth, f(0.5)=0.76). No ESTIF-specific enhancement at linear order; C3 filter passed trivially (no boost) | ✅ |
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
| A-6 | **The distinguishing prediction.** FOUR exactness locks catalogued (RHAC-007): Ωm/Ω_k=0 (kill-shot), γ≈0.55 & slip=1 (Front 3), w=−1 (RHAC-004), c_gw=c (RHAC-008). Honest finding: none separates ESTIF-Core from FLAT ΛCDM — all separate it from GR's extra freedoms. Genuinely-distinct content lives OFF the linear sheet (nonlinear halos, N-body wall) + the exactness-as-law structure. | 🔶 catalogued |
| A-R | ~97,300-year recombination age — **retired with exponential shrinkage; EC-fork only; not a Path One prediction** | ❌ retired by design |

## Cross-cutting — publication track

| # | Item | Status |
|---|---|---|
| P-1 | v6.3 documentation suite + v6.3.1 errata + canonical Path One scripts (DESI parity; AIC/BIC) delivered | ✅ |
| P-2 | Apply the errata (C1–C6) across the repo; commit as v6.3.1 | ✅ 9 July 2026 |
| P-4 | Apply the bootstrap amendments across the repo; commit as v6.3.2 | ⬜ |
| P-3 | Gravity-letter final pass (C1/C6 wording + preferred-frame section) → submit (barrier 6) | ⬜ |

---

## Tally (35 tracked items)

| ✅ complete | 🔶 partial | ⬜ open | ❌ blocked/retired |
|---|---|---|---|
| **19** | **2** | **12** | **2** |

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
model E1, falsified against DESI DR2 (χ²/N = 3232). v6.4 (11 July 2026): axiom A1 amended to A1′ (even on average, local matter-sourced
dents) — see RHAC-006. Strict-A1 results survive as the exact-evenness limit.
"Frozen eddy" label deprecated (RHAC-004); full rename pending NAMING.md.

**Checklist version:** 1.0 · **Date:** 8 July 2026
