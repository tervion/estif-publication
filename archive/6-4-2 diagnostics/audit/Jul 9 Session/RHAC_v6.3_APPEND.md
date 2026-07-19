<!--
  APPEND THIS BLOCK TO THE END OF docs/plan/RHAC.md
  (it follows the existing "RHAC UPDATE — v6.1" pattern)
-->

---

## RHAC UPDATE — v6.3 "THE SPLIT" (8 July 2026)

**Major milestone. The project forks into two tracks.** See
`MILESTONE_v6.3_THE_SPLIT.md` in the repo root. Summary of what changed and the
scenarios it resolves.

---

### 🏁 MILESTONE: Project split — Path One (Core) and Path Two (Extended)

Applying this file's own decision filter (Step 2, filter 4: *"Does it improve
agreement with observations?"*), the tilt-based cosmology fails against its own
frozen-eddy limit and is marked 🔴 pivot. The project forks:

- **Path One — ESTIF-Core (clean) ✅ recommended default.** Gravity on a derived
  field equation; cosmology = frozen cosmic eddy → cosmological constant
  (χ²/N = 1.92, ties ΛCDM). Tilt apparatus retired to an appendix. Gravity
  letter stands and is strengthened.
- **Path Two — ESTIF-Extended (hard) 🔬 high-risk.** Derive the leading w(z)
  correction to w = −1 from the full rotating-shear/vorticity stress tensor,
  targeting the mild DESI thawing (χ²/N ≈ 0.66).

---

### ✅ Scenario N: Is gravity DERIVED or borrowed from Schwarzschild? — RESOLVED (DERIVED)

**What if:** the ESTIF force law only reproduces Newton because the flow profile
is matched to the Schwarzschild solution (n = ½ ⇒ β = √(1−x)) and then
differentiated — i.e. the answer is assumed, not derived.

**Resolution:** DERIVED (Task 4, July 2026). Starting only from the flow axioms
(A1 flat slices, A2 universal speed c, A3 empty space is not a source), the
Gauss–Codazzi engine FORCES `rho_eff = m′(r)/(4πr²)` — mass continuity, i.e.
Poisson in integrated form. Vacuum ⇒ v² = 2GM/r uniquely ⇒ exact Schwarzschild;
uniform ball ⇒ rho0 exactly. The old "D2 Poisson postulate" is now a theorem for
vacuum, Newton, and Schwarzschild. Scripts: `estif_task4_field_equation.py`
(5/5), `estif_flow_signature_dynamics.py` (18/18), `estif_converse_flow_law.py`
(Birkhoff in flow variables), `estif_tmunu_gauss_codazzi.py` (engine validated).

**Remaining:** strong-field pressure/stress sector (full off-diagonal T_μν).
Needed for none of {vacuum, Newton, Schwarzschild}.

**Status:** 🟢 Resolved — gravity field equation derived.

---

### ✅ Scenario O: Fidelity — are the flow axioms actually IN the theory? — RESOLVED (partially; writing task)

**What if:** the derivation in Scenario N proves a *neighbouring* principle, not
the one v6.2 actually states.

**Resolution:** audit run (`estif_fidelity_audit.py`, 71 files). Findings:
- A1 (flow, not stretch): **CONTESTED** — the flow picture is load-bearing, but
  `ESTIF_CONCEPT.md` also runs a competing "shrinking-ruler" narrative
  (lines ~49–57). Must be resolved in favour of flow.
- A2 (universal speed c): **absent from theory docs** — present only in the
  derivation scripts. The existing v_flow = cx₀ ≈ 0.31c is the *sideways
  component* of a total-c motion, not the full speed (a mislabel to fix).
- A3 (empty space is not a source): **absent from theory docs** — present only
  in scripts.

Author confirmed A2 and A3 match the intended physical picture. Therefore this
is a **writing/consolidation task**, not new physics: Path One must (i) write A2
and A3 into the theory documents, (ii) retire the shrinking-ruler narrative,
(iii) relabel 0.31c as the sideways component of a total-c flow.

**Status:** 🟡 Active (Path One writing task).

---

### ✅ Scenario E (re-update): DESI w(z) — circularity FIXED, tilt shape is the real problem

**Update to the v6.1 resolution.** Scenario E called for a self-consistent
Ω_tilt(z). Done (Task 5): `x(z) = x₀(1+z)H₀/H_ESTIF(z)` via fixed-point solve,
no ΛCDM. Result on real DESI DR2: χ²/N 10.80 → **3.35**. Large correctness
improvement, still short of ΛCDM (1.92). Task 5b localized the residual to the
tilt *shape* at low-to-mid z (not the ruler); ESTIF's self-consistent tilt w(z)
already thaws toward the DESI-preferred curve (within ~0.05), and 3.35 is where
DESI's own published w0,wa sit on the BAO-only subset (3.09).

**Status:** 🟢 Circularity resolved → 🔴 tilt shape superseded by frozen-eddy
(see Scenario P). Scripts: `estif_task5_desi_selfconsistent.py`,
`estif_task5b_cosmo_eos.py`.

---

### ✅ Scenario P: Can the eddy EoS be derived? — RESOLVED (naive routes fail; frozen-eddy reframe)

**What if:** the cosmic eddy's equation of state w(z) can be derived from the
rotating-hypersurface kinetic energy, reproducing the DESI-preferred thawing.

**Resolution (Task 6):** the two natural reductions FAIL badly against DESI DR2:
- E1 conserved-angular-momentum spin → w = +1 (stiff), χ²/N = 3232;
- E2 expansion tracker → thaws to ~0, χ²/N = 754.
Both dilute/blueshift the wrong way. **But** the comparison exposed the decisive
fact: the **frozen-eddy limit** (constant eddy density → de Sitter → w = −1,
which Task 4 derives for free) scores **χ²/N = 1.92, tying ΛCDM and beating the
tilt formula's 3.35**. On DESI the entire tilt apparatus is a net negative.

**Decision:** retire Ω_tilt(z) from the cosmology claim (Path One). The correct
first-principles route (Path Two) is to derive the *leading correction* to
w = −1 from the full vorticity stress tensor — a small perturbation, not the
strong evolution E1/E2 produce. Script: `estif_task6_eddy_eos.py`.

**Status:** 🟢 Resolved — frozen eddy is the honest cosmology; tilt retired;
Path Two is the remaining hard derivation.

---

### Scenario D (re-update): Ω_tilt high-z cutoff — moot under the reframe

The z < 2 hard cutoff was scaffolding for the tilt term. Under Path One the tilt
term is retired, so the cutoff is moot. Under Path Two the frozen w = −1 limit
does not diverge, so no cutoff is needed; only the leading correction must be
kept well-behaved. **Status:** 🟢 Superseded.

---

### Scenario H (re-update): Ωm = x₀ — now connected to the derived machinery

The Task 4 result derives the *local* field equation (rho_eff = m′/4πr²) from the
flow metric. The homogeneous version of the same calculation is the natural route
to test whether ρ_eddy = x₀ρ_crit emerges. This remains the most important
theoretical target for the dark-matter sector and is now well-posed rather than
abstract. **Status:** 🟡 Active, sharpened.

---

### Updated Summary Statistics (v6.3)

| Category | Count |
|---|---|
| Total scenarios documented | 24 |
| Resolved | 16 ✅ |
| Active | 6 🟡 |
| Budget wall (simulation) | 1 🔴 |
| Major pivots | 6 (H(t), dynamic n, Ω_tilt inversion, fluid→collisionless, cosmology→gravity letter, **tilt→frozen-eddy split**) |

---

**Document Version:** 6.3 (ESTIF v6.3 — "The Split")
**Last Updated:** 8 July 2026
