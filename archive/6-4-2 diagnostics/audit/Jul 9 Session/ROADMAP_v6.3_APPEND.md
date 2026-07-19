<!--
  APPEND THIS BLOCK TO THE END OF docs/guide/ROADMAP.md
  (it follows the existing "ROADMAP UPDATE — v6.1" pattern)
-->

---

## ROADMAP UPDATE — v6.3 "THE SPLIT" (8 July 2026)

The mission is unchanged ("replace ΛCDM, don't adjust it"), but the route forks.
Two derivation results since v6.2 restructure the roadmap. See
`MILESTONE_v6.3_THE_SPLIT.md`.

---

### Foundation upgrade — gravity is now DERIVED

The strong-field gravity floor was marked 100% complete in v6.2, but its
foundation was a Schwarzschild match (set n = ½). v6.3 replaces that: the field
equation `rho_eff = m′(r)/(4πr²)` (Poisson in integrated form) is now **forced**
by the flow axioms via Gauss–Codazzi, giving exact Schwarzschild in vacuum and
Newtonian ρ₀ for a source, with no Poisson postulate (Task 4). The foundation is
now derived, not borrowed.

```
Foundation (strong-field gravity):    ████████████  100%  ✅ Complete + DERIVED
Field equation from flow axioms:      ████████████  100%  ✅ NEW (Task 4)
Strong-field pressure/stress sector:  ███░░░░░░░░░   25%  🔄 full T_μν remaining
Cosmology (Path One = frozen eddy):   ████████████  100%  ✅ ties ΛCDM, derived
Cosmology (Path Two = thawing eddy):  █░░░░░░░░░░░   10%  🔬 hard, vorticity T_μν
Dark matter:                          ██░░░░░░░░░░   20%  🔄 analytical done, sim wall
```

---

### The fork replaces the old Phase 5–6 cosmology plan

The prior roadmap's Phase 5.4 (self-consistent Ω_tilt) and Phase 6 (CMB on top
of Ω_tilt) are **superseded**. Phase 5.4 was done (Task 5: 10.8 → 3.35) and
revealed the tilt shape is the problem; Task 6 then showed the frozen-eddy limit
(w = −1) ties ΛCDM and beats the tilt. So the cosmology roadmap becomes:

#### PATH ONE — ESTIF-Core (clean) ✅ recommended default

- **P1.1 Write the axioms into the theory.** Add A2 (universal speed c) and A3
  (empty space is not a source) to `ESTIF_CONCEPT.md`; the fidelity audit found
  they exist only in the derivation scripts.
- **P1.2 Resolve the A1 conflict.** Retire the shrinking-ruler narrative in
  favour of the flow (Painlevé–Gullstrand) picture. Relabel v_flow = cx₀ ≈ 0.31c
  as the *sideways component* of a total-c flow, not the full speed.
- **P1.3 Rewrite the cosmology sector** as "constant cosmic eddy → cosmological
  constant, χ²/N = 1.92 (ties ΛCDM)." Move Ω_tilt(z), N_MAX, B, the sign-flip and
  the z<2 cutoff to an "explored and set aside" appendix.
- **P1.4 Submit the gravity letter** (strengthened by the derived field equation).

#### PATH TWO — ESTIF-Extended (hard) 🔬 high-risk research

- **P2.1 Vorticity stress tensor.** Derive the leading correction to w = −1 from
  the full rotating-shear / vorticity T_μν (the cosmological half of the T_μν
  work). Target the mild DESI thawing (χ²/N ≈ 0.66). The naive reductions E1
  (stiff, χ²=3232) and E2 (tracker, χ²=754) are already falsified; the full
  off-diagonal tensor is required.
- **P2.2 CMB / ISW** only after P2.1 produces a well-behaved, DESI-consistent
  H(z). On Path One (pure Λ) the CMB check is the standard ΛCDM one.

---

### Immediate Next Steps (v6.3 priority order)

| Priority | Task | Path | Type | Estimate |
|---|---|---|---|---|
| 1 | Adopt Path One; write A2 + A3 into theory; retire shrinking-ruler | One | Writing | days |
| 2 | Rewrite cosmology sector around frozen-eddy = Λ | One | Writing | days |
| 3 | Submit gravity letter (now on derived field equation) | One | Submission | this week |
| 4 | Strong-field pressure/stress sector (full T_μν) | Both | Theory | unknown |
| 5 | Vorticity stress-tensor derivation of leading w(z) | Two | Theory | unknown (hard) |

---

**Roadmap Version:** 6.3 — "The Split"
**Last Updated:** 8 July 2026
