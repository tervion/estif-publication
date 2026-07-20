# ESTIF — Project Status, Checklist, Roadmap, Structure

**Version:** 6.4.2 · **Last updated:** 21 July 2026
**Status:** Gravity letter drafted and strengthened; send gated on the Bullet Cluster response (B-8). Path One (ESTIF-Core) is the default track; Path Two (ESTIF-Extended) is CLOSED — honorable null (11 July 2026, RHAC-004).
**Companions:** `docs/SCIENCE.md` (framework + evidence) · `docs/RHAC.md` (decision archive) · `CHANGELOG.md` (history) · `tests/docs/TEST_INDEX.md` (test-suite index)

This file consolidates and supersedes STATUS.md, ROADMAP.md, PATH_ONE_CHECKLIST.md, the Phase 2 declaration's status record, the document-update guide, and the dated progress reports. Their operative content is here; their history is in the changelog; the originals are archived.

---

## Executive summary

ESTIF is a geometric model deriving gravity, time, and galactic dynamics from the claim that 3D space is an even hypersurface (even on average — A1′) carried through a 4D bulk at the speed of light. As of v6.4:

- **Gravity is derived, not borrowed.** The flow axioms force the field equation (mass continuity = Poisson), giving exact Schwarzschild in vacuum and the Newtonian source in the weak field, with no Schwarzschild matching. The matter coupling (8πG) is adopted, not derived (precision C6).
- **Cosmology is honest.** The best result the axioms permit is an imported cosmological constant — w = −1 exactly, χ²/N = 1.92, tying ΛCDM on DESI DR2. Phase 2 proved the residual dark-energy sector is EMPTY (strong-form null): the axioms cannot manufacture dark energy. The evolving Ω_tilt(z) apparatus is retired.
- **A1 → A1′** (11 July 2026): slices even on average, matter dents locally, ⟨curvature⟩ ≡ 0 as law — the registered kill-shot falsifier (current measurement 0.0007 ± 0.0019). One amendment restored linear growth (f(z=0.5) = 0.76, DESI-consistent), the gravitational-wave sector (c_gw = c derived, GW170817-consistent), and a container for the ζ = 10⁻⁵ seed.
- **a₀ is a horizon quantity** (13 July 2026 doctrine, RHAC-010): the SCALE of a₀ (~c·H) is derived; the NUMBER (the O(1) prefactor) is not, and local flow provably cannot make it. This is the keystone open problem.
- **Ωm is not derived.** The bootstrap gives a unique root 0.31408 (0.53σ from Planck) conditional on postulate P, which is proven non-derivable (RHAC-009); and the causal-horizon → Ωm ≈ 0.3 result has prior art (Gaztañaga 2019–2023) — no novelty claim on Ωm is defensible.

---

## Three-sector status

### Gravity — ✅ solid, on a derived foundation

| Result | Value | Status |
|---|---|---|
| Field equation derived (mass continuity = Poisson) | ρ_eff = m′(r)/(4πr²) forced by the flow metric | ✅ (Task 4, 5/5) |
| Vacuum → exact Schwarzschild | v² = 2GM/r unique; full Einstein tensor = 0 | ✅ |
| Weak-field source → Newtonian ρ₀ | uniform ball returns ρ₀ exactly | ✅ |
| Signature + SR from Euclidean bulk + universal c | (−,+,+,+) and dτ/dt = √(1−v²/c²) | ✅ (18/18) |
| a₀ — horizon doctrine | scale ≈ c·H derived; working form H₀cx₀/√3 = 1.179×10⁻¹⁰ (1.72%; 0.66% with bootstrap Ωm, conditional); exact prefactor OPEN | 🔶 scale derived, number open (RHAC-010) |
| GW sector: c_gw = c derived (C-15) | single geometry; GW and light share the null cone; GW170817 passed structurally | ✅ (RHAC-008) |
| SPARC BTFR (87 quality-1 galaxies) | RMS = 15.6% | ✅ within observed scatter |
| Weak-lensing BTFR (Mistele+2024) | χ²/N ≈ 1.5 with the 0.1 dex M*/L systematic (5.6–8.0 stat-only); predictions ~5–13% low | ✅ executed 11 Jul; consistent within systematics (B-7) |
| Implied-a₀ inversion (correlated systematics) | +0.7σ to +2.0σ above derived (ETG worst) | ⚠️ mild tension, recorded |
| EHT M87* shadow / LISA delay | 0.00σ / 49.2σ | ⚠️ consistent; *deviation* claims conditional (C1) |
| a₀ redshift constancy · parameter independence | H(z) cancels exactly · 3,600 H₀/Ωm combos pass | ✅ |

> **Conditional (C1).** The ESTIF vacuum is exactly Schwarzschild, so any *deviation* from GR in shadows or GW propagation must be sourced by the non-vacuum eddy background — not yet derived. Observations are consistent; the deviation is what awaits derivation (single-speed flow forces p_r = −ρ, limitation L1).

**Open in gravity:** strong-field pressure/interior sector (TOV-level, full off-diagonal T_μν); the a₀ prefactor and its local mechanism (below).

### Cosmology — 🔄 reframed; imported Λ ties ΛCDM

| Test | Result | Status |
|---|---|---|
| Constant cosmic term = imported Λ | χ²/N = 1.92 (ties ΛCDM) | ✅ honest best; residual sector provably EMPTY (Phase 2) |
| Ω_tilt(z), circular | χ²/N = 10.80 | ❌ superseded |
| Ω_tilt(z), self-consistent (Task 5) | χ²/N = 3.35 | 🔄 circularity fixed, still short of ΛCDM — tilt retired |
| Fitted CPL bar (Task 5b / target lock) | χ²/N ≈ 0.66 at (w₀, wₐ) = (−0.85, −0.45) | ⚠️ under-marginalized (C4); pre-registered as the bar any outside-the-axioms proposal must beat |
| Derived eddies E1 / E2 | χ²/N = 3232 / 754 | ❌ falsified (Task 6) |
| Growth under A1′ | D₊ = H·∫da/(aH)³ = GR growth; f(0.5) = 0.76 | ✅ DESI RSD-consistent |
| Ωm bootstrap (conditional on P) | unique root 0.31408 (0.53σ); DESI χ²/N = 1.618 fixed-ruler | 🔶 conditional; P non-derivable (RHAC-009); no Ωm novelty (Gaztañaga) |
| Age of universe | 13.379 Gyr | ✅ |

### Dark matter — 🟡 analytical phase complete; simulation wall

Ωm = x₀ (0.12%, consistency relation C2 — final), Ωdm = x₀ − Ωb (0.10%), σ/v_esc = 0.5 exact, λ_Jeans = 2.565r. Halo internal structure (v_flat = 220 km/s, δ ~ 10⁵) requires N-body — off-machine, collaboration target (B-9, the budget wall).

---

## The two paths

**Path One — ESTIF-Core (clean) ✅ default.** Derived gravity + imported cosmological constant (ties ΛCDM). Publishable; no failing tests. The gravity letter rides on this path.

**Path Two — ESTIF-Extended — 🟢 CLOSED 11 July 2026 (honorable null; RHAC-004).** Phase 2 computed the *entire* residual stress the axioms permit — pre-registered, submitted unmodified — and found it empty: w = −1 exactly, Λ imported. Any future evolving-w must leave the single-speed flow class (a second metric function or explicit time dependence, outside A1–A3), and must beat the pre-registered χ²/N ≈ 0.66 bar under honest marginalization. Full record: `docs/SCIENCE.md` (Phase 2 section) and RHAC-004/005/006.

---

## Path One checklist (35 tracked items)

Maintenance rule (RHAC discipline): update status when a scenario closes; never delete items — mark ✅ / ❌ and date them. Legend: ✅ complete · 🔶 partial/downgraded · ⬜ open · ❌ externally blocked or retired.

### The six barriers (adversarial review)

| # | Barrier | Fix | Status |
|---|---|---|---|
| 1 | Exact-Schwarzschild vacuum vs EHT/LISA deviation claims | C1 applied (deviations conditional on the eddy-stress sector) | ✅ 9 Jul (P-2) |
| 2 | "Derived, not borrowed" overclaims the coupling | C6 applied (axioms select GR's constraint sector in PG gauge; coupling adopted) | ✅ 9 Jul (P-2) |
| 3 | Ωm = x₀ hides an Ωm-dependent horizon integral | C2 applied (consistency relation — final per RHAC-009) | ✅ |
| 4 | Pressure/interior sector (TOV-level) | new derivation | ⬜ (= C-12) |
| 5 | Preferred-frame nature of A2 | explicit simultaneity section | ⬜ (= C-11) |
| 6 | Peer review — the only external test yet faced | submit the letter | ⬜ (= C-14; gated on B-8) |

### Goal C — Gravity = Time

| # | Item | Status |
|---|---|---|
| C-1 | A2 stated as law: dw² + dσ² = c²dt²; proper time c·dτ = dw | ✅ |
| C-2 | Lorentzian signature (−,+,+,+) derived from Euclidean bulk + A2 | ✅ |
| C-3 | Exact SR recovered, zero free parameters | ✅ |
| C-4 | Static clock rate = local flow speed u(r); pull = −c²∇u | ✅ |
| C-5 | Field equation derived (Task 4, 5/5) | ✅ |
| C-6 | Vacuum (A3) ⇒ exact Schwarzschild (Birkhoff in flow variables) | ✅ |
| C-7 | Uniform ball ⇒ ρ₀ exactly | ✅ |
| C-8 | Axioms written into the theory; shrinking-ruler narrative retired | ✅ |
| C-9 | Apply C1 across the repo | ✅ 9 Jul 2026 (v6.3.1) |
| C-10 | Apply C6 across the repo | ✅ 9 Jul 2026 (v6.3.1) |
| C-11 | Preferred-frame / relativity-of-simultaneity section (barrier 5) | ⬜ pre-submission |
| C-12 | Pressure/interior sector — TOV-level star (barrier 4) | ⬜ |
| C-13 | Moving-observer bulk identity inside eddies (signature-suite caveat) | ⬜ |
| C-14 | Peer review of the gravity letter (barrier 6) | ⬜ |
| C-15 | GW sector under A1′: c_gw = c DERIVED; GW170817 structural pass | ✅ 12 Jul 2026 (RHAC-008) |

### Goal B — No dark matter

| # | Item | Status |
|---|---|---|
| B-1 | a₀ working form H₀cx₀/√3 (1.72%), Step 1 underwritten by the derived field equation — **doctrine note:** scale derived, prefactor open (RHAC-010) | ✅ with 🔶 doctrine caveat |
| B-2 | SPARC validation: 87 quality-1, RMS 15.6% | ✅ |
| B-3 | a₀ redshift constancy + parameter independence (3,600 combos; 8 datasets) | ✅ |
| B-4 | Ωm bootstrap: unique root 0.31408 (0.53σ), closure a₀ → 0.66%, DESI 1.618 fixed-ruler — conditional on P | 🔶 |
| B-4b | Gaztañaga comparison memo | 🔶 verdict delivered 13 Jul (🔴 no Ωm novelty — cite prominently, lead with a₀); item-4 numerical cross-check vs his Ω_Λ ≈ 0.70 still owed |
| B-4c | Part B: derive P from the axioms | ✅ CLOSED 12 Jul (RHAC-009) — non-derivable; P is a predictive postulate |
| B-5 | ρ_eddy = x₀ρ_crit from the homogeneous field equation | ⬜ |
| B-6 | D_ESTIF(z) linear growth | ✅ 11 Jul (no-go under strict A1; restored under A1′, f(0.5)=0.76) |
| B-7 | Weak-lensing BTFR vs derived a₀ | ✅ executed 11 Jul 2026 — consistent within the 0.1 dex M*/L systematic (χ²/N ≈ 1.5; ~9% low, degenerate with the calibration). Receipts: `tests/scripts/btfr_lensing.py`, `tests/scripts/a0_tension_corrected.py` |
| B-8 | Bullet Cluster response — hardest external test for a no-dark-matter framework | ⬜ **gates the letter send** |
| B-9 | Halo internal structure — N-body | ❌ blocked (compute/collaborator) |

### Goal A — Expansion is 4D infall

| # | Item | Status |
|---|---|---|
| A-1 | Mechanism as an equation: H = c·d(ln b)/dw (engine-verified Friedmann analog) | ✅ |
| A-2 | Constant-Λ limit ties ΛCDM on DESI DR2 (1.92; ΔAIC +0.59 with Ωm = x₀) | ✅ |
| A-3 | One-flow unification: v² = 2GM/r + H²r² (exact Schwarzschild–de Sitter) | ✅ |
| A-4 | Circularity removed (10.8 → 3.35); tilt cosmology retired (ΔAIC +18, decisive) | ✅ |
| A-5 | CMB anisotropy spectrum as an ESTIF prediction (Boltzmann code); background CMB inherited and passes | ⬜ |
| A-6 | Distinguishing prediction: FOUR exactness locks catalogued (Ω_k = 0 kill-shot · γ ≈ 0.55 & slip = 1 · w = −1 · c_gw = c). Honest: none separates ESTIF-Core from ΛCDM at Ω_k = 0 at linear order; all separate it from GR's extra freedoms. Distinct content lives off the linear sheet + in the exactness-as-law structure | 🔶 catalogued |
| A-R | ~97,300-yr recombination age (EC-fork only) | ❌ retired by design |

### Publication track

| # | Item | Status |
|---|---|---|
| P-1 | v6.3 doc suite + errata + canonical Path One scripts delivered | ✅ |
| P-2 | Errata C1–C6 applied repo-wide (v6.3.1) | ✅ 9 Jul 2026 |
| P-4 | Bootstrap amendments applied repo-wide (v6.3.2) | ⬜ |
| P-3 | Letter final pass (C1/C6 wording + preferred-frame section) → submit | ⬜ after C-11, gated on B-8 |

**Tally:** ✅ 20 · 🔶 3 · ⬜ 10 · ❌ 2 (of 35). The derivation core of Path One is finished; what remains is two referee-facing sections (C-11, C-12), the Bullet Cluster response (B-8), documentation hygiene (P-4), and the external gauntlet (C-14).

---

## Current priorities and open problems

1. **Bullet Cluster response (B-8)** — gates the letter send to McGaugh.
2. **a₀ keystone (RHAC-010):** (a) the O(1) prefactor — derive it from horizon geometry with zero fit freedom (the former √3 slot); (b) the LOCAL mechanism — what makes the horizon scale c·H act at galaxy radii? Every emergent-gravity approach (Verlinde, Padmanabhan) *assumes* this; deriving it from the 4D inward flow would be ESTIF's genuine contribution. `ripple_speed.py`'s `estif_eom` slot is the prepared instrument for the equation-of-motion insertion (parked).
3. **x_c = 0.272** — derive from black-hole-edge geometry; note the double-duty flag (crossover 0.272 vs data-preferred edge ≈ 0.221; possibly two conflated surfaces — see `docs/SCIENCE.md`).
4. **C-11 simultaneity section**, then **letter submission** (C-14).
5. **Interior/pressure sector** (C-12) — the remaining rigorous gravity step.

**Parked / queued (deliberately, not forgotten):** Bullet Cluster memo (item 1 above is the unblock); letter v6.4 held pending it; Zenodo DOI refresh for v6.4.2; PDF build of the letter; `ripple_speed.py` EoM insertion; Gaztañaga item-4 numerical cross-check; black-hole-seed front follow-up (front 2's direct-collapse channel).

---

## Repository structure (v6.4.2)

```
estif_publication/
├── README.md · CHANGELOG.md · LICENSE · CITATION.cff · setup.py · requirements.txt
├── src/                      estif_ec_gr_constants.py · estif_ec_gr_model.py ·
│                             estif_ec_gr_run_simulation.py  (21/21)
├── tests/
│   ├── scripts/              all 85 test/receipt scripts (pure .py)
│   │   └── phase2_a1prime/   Phase-2 doors + A1′ growth + re-census suite (kept whole)
│   ├── docs/                 TEST_INDEX.md + every cached data .txt
│   │                         (DESI dr1/dr2 family · pathone_* · dr2b_* ·
│   │                          pathtwo_target_lock_output.txt · SPARC VOTable cache)
│   └── plots/                generated figures (gitignored by default)
├── docs/                     PROJECT.md (this file) · SCIENCE.md · RHAC.md ·
│                             Letter/ · latex/  (+ superseded originals until archived)
├── data/ · results/          observational inputs · validated outputs
└── archive/                  6-2/ snapshot · "6-4-2 diagnostics"/ snapshot
                              (test_UKN.py, test_UKN2.py live HERE; audit/ lives
                               inside it — never modified)
```

Conventions: every CACHE_DIR targets `tests/docs` (scripts/ stays pure .py; caches never duplicate); all sys.path hops are `__file__`-relative and depth-aware, so scripts run from any working directory; `tests/docs/TEST_INDEX.md` is the single source of truth for the suite (the old `_index_view/` mirror is deleted). Note: the archive's plots and PDF are gitignored — they exist only on the author's machine unless force-added.

---

## Receipt map (where the key receipts live)

| Receipt | Path |
|---|---|
| Field equation, signature/SR, Birkhoff, ADM engine | `tests/scripts/estif_task4_field_equation.py` · `estif_flow_signature_dynamics.py` · `estif_converse_flow_law.py` · `estif_tmunu_gauss_codazzi.py` |
| Cosmology reframe (Tasks 5/5b/6) + target lock | `tests/scripts/estif_task5_desi_selfconsistent.py` · `estif_task5b_cosmo_eos.py` · `estif_task6_eddy_eos.py` · `estif_pathtwo_target_lock.py` |
| Phase-2 doors + A1′ growth + re-census (PASS 5/5) | `tests/scripts/phase2_a1prime/` (README inside lists verdicts and expected outputs) |
| Strict-A1 census + A1′ re-audit (UKN pair) | `archive/6-4-2 diagnostics/test_UKN.py` · `test_UKN2.py` |
| Fronts 1–3 (growth/σ₈/JWST · first holes · γ + slip) | `tests/scripts/estif_front1_growth_sigma8_jwst.py` · `estif_front2_first_hole_recipe.py` · `estif_front3_second_discriminator.py` |
| GW sector (C-15) · P non-derivability | `tests/scripts/estif_C15_gw_sector.py` · `estif_P_derivation_attempt.py` |
| a₀-horizon doctrine quartet | `tests/scripts/a0_horizon_test.py` · `a0_prefactor_derivation.py` · `estif_flow_sim.py` · `estif_horizon.py` |
| 11 Jul letter-session trio | `tests/scripts/btfr_lensing.py` · `a0_tension_corrected.py` · `mu_extraction.py` |
| Ωm bootstrap | `tests/scripts/estif_omega_bootstrap.py` · `estif_bootstrap_closure.py` |
| Fidelity audit · keystone instrument | `tests/scripts/estif_fidelity_audit.py` · `tests/scripts/ripple_speed.py` (prepared `estif_eom` slot) |

Full index: `tests/docs/TEST_INDEX.md`.

---

## Documentation map and update conventions

Four living documents (this consolidation replaces the old update guide):

1. **`docs/PROJECT.md`** — status, checklist, roadmap, structure. Update when an item's status changes; never delete checklist items; date every flip.
2. **`docs/SCIENCE.md`** — framework and evidence. Update when a claim's epistemic status changes (derived / conditional / retired), citing the RHAC entry and receipt.
3. **`docs/RHAC.md`** — the decision archive. Append-only in spirit: new RHAC-NNN entries for decisions; surgical status reconciliations only, dated.
4. **`CHANGELOG.md`** — reverse-chronological history; backdated entries marked as such.

Standing terminology rule (rename batch delivered v6.4.2, RHAC-011): spatial geometry is described as *even* / *even on average* — never the old term; behavioural ranges are *gradients* of one equation — never the old term. Every claim carries its status word (derived · adopted · conditional · consistency relation · retired) and its receipt path.

---

## Lineage

ESTIF-FD v1.0 (2024, exponential shrinkage; ruled out by SNe) → sibling branches: the EC fork (June 2025; shrinkage kept, variable H(t)) and the main line v3.0 → v6.x. v6.3 (8 Jul 2026): "The Split" — field equation derived, project forks into Path One/Two. v6.3.1 (9 Jul): adversarial errata C1–C6. v6.3.2 (9 Jul): Ωm bootstrap (conditional). v6.3.3 (11 Jul, backdated): letter measurement session (target lock, lensing BTFR, a₀ tension, μ extraction). v6.4.0 (11–12 Jul): Phase 2 honorable null; A1 → A1′; C-15 closed; P non-derivable. v6.4.1 (13 Jul): fronts 1–3, four-lock ledger, a₀-horizon doctrine. v6.4.2 (21 Jul): repository restructure + this documentation consolidation. EC's vortex mechanism corresponds to Task 6's model E1, falsified against DESI DR2 (χ²/N = 3232). Strict-A1 results survive as the exact-evenness limit of A1′.

---

**Document version:** 6.4.2 · 21 July 2026 · consolidates STATUS 6.4.1, ROADMAP 6.4.1, PATH_ONE_CHECKLIST 1.0 (+ reconciliations dated above)
