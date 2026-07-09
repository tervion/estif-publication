# ESTIF PROGRESS REPORT — the v6.3 / v6.3.1 / bootstrap campaign

**Report date:** 9 July 2026
**Work window:** July 2026 session series (Tasks 1–6, the Split, errata, bootstrap)
**Prepared for:** repo record + document-update planning
**Target location in repo:** `docs/report/PROGRESS_REPORT_2026-07-09.md`

---

## 1. Executive summary

The session began with three named bottlenecks — Ω_tilt(z) circularity, T_μν, and
N-body — and ended with: a **derived gravitational field equation** (the T_μν
energy sector closed), the **circularity fixed and then superseded** by the
frozen-eddy reframe, the **N-body wall correctly descoped** (population-level
tests are analytic; only halo interiors need simulation), the project **split
into Path One (ESTIF-Core) and Path Two (ESTIF-Extended)**, a six-item
**adversarial errata** (v6.3.1), a quantified **JWST early-structure
specification** with a σ8 hard filter, and — newest — the **Ωm bootstrap**:
conditional on one principle, the framework now computes Ωm, Ω_Λ, x₀,
r_universe, and a₀ from three directly measured inputs (H₀, T_CMB, N_eff).

## 2. The three original bottlenecks — before and after

| Bottleneck | Before | After |
|---|---|---|
| Ω_tilt(z) circularity | x(z) used H_ΛCDM as its own ruler; χ²/N = 10.8 vs DESI DR2 | Fixed (self-consistent solve → 3.35), then **superseded**: frozen eddy = Λ scores 1.92; tilt retired (ΔAIC +18) |
| T_μν | Postulated; force law matched to Schwarzschild (n = ½) | **Energy sector derived**: ρ_eff = m′(r)/(4πr²) forced by the flow metric; vacuum → exact Schwarzschild; ball → ρ₀. Open: pressure/stress sector, vorticity sector |
| N-body | Treated as blocking all structure claims | **Descoped**: halo-abundance/JWST tests are analytic (need only D_ESTIF(z)); the wall blocks only halo *interiors* (rotation curves) |

## 3. Results by workstream (all numbers are runtime output, reproduced on the author's machine where noted)

### A. Engine and the gravity derivation
- Gauss–Codazzi/ADM symbolic engine built and validated against exact GR: flat
  FRW → Friedmann + pressure; ADM Hamiltonian identity on flat and curved
  slices; global de Sitter (hyperboloid slicing) → G_μν = −Λg_μν.
  (`estif_tmunu_gauss_codazzi.py`; reproduced on author's Mac.)
- Signature + SR + dynamics suite, **18/18 PASS** (author-reproduced):
  Lorentzian (−,+,+,+) emerges from a Euclidean bulk + universal-c (minus sign
  produced by solving the constraint); exact SR; g_tt = −c² becomes a
  consequence; static clock rate = u(r); pull = −c²∇u; D1 "plughole" control
  fails (1/r⁵); Gullstrand–Painlevé flow metric = **exact Schwarzschild**
  (Einstein tensor ≡ 0), raindrop dτ = dt. (`estif_flow_signature_dynamics.py`)
- Converse theorem: vacuum + flat slices + universal-c **forces** v² = 2A/r
  uniquely (one constant; power-law scan → n = −1/2 unique) — Birkhoff in flow
  variables. Bonus: ρ_eff = const → v² = 2A/r + (λ/3)r² = exact
  Schwarzschild–de Sitter; one flow field carries gravity (small r) and Hubble
  flow (large r), crossover r³ = GM/H². (`estif_converse_flow_law.py`)
- **Task 4 — the field equation derived (5/5 PASS):** the engine forces
  ρ_eff = m′(r)/(4πr²) (mass continuity = Poisson in integrated form) from the
  flow metric; vacuum ⇒ exact Schwarzschild; uniform ball ⇒ ρ_eff = ρ₀ exactly.
  The old "D2 Poisson postulate" is a theorem for vacuum/Newton/Schwarzschild.
  (`estif_task4_field_equation.py`)
- Fidelity audit (71 files): A1 CONTESTED (flow vs shrinking-ruler narrative,
  187 support / 43 contradiction sentences); A2 and A3 **absent from theory
  docs** (present only in derivation scripts); v_flow = cx₀ ≈ 0.31c mislabeled
  as "the flow speed" (it is the in-slice component of a total-c motion).
  Author confirmed A2/A3 match intended physics → writing task, executed in the
  v6.3 concept rewrite. (`estif_fidelity_audit.py`)

### B. Cosmology reframe and the Split
- **Task 5 — de-circularization:** self-consistent x(z) via fixed point of
  ESTIF's own H(z). Real DESI DR2: 10.80 → **3.35** (11/13 bins within 2σ;
  worst bin z = 0.51 DM at −3.07σ); w_eff(0): −0.65 → −0.80.
  (`estif_task5_desi_selfconsistent.py`; author-reproduced)
- **Task 5b — what DESI wants:** best constant-w flow: w = −0.945, χ²/N = 0.884;
  best CPL: w0 = −0.85, wa = −0.45, χ²/N = **0.660**; DESI's published
  (w0, wa) = (−0.73, −0.66) scores 3.089 on this BAO-only subset — i.e. ESTIF's
  3.35 sits beside DESI's own combined-fit values here. ESTIF's derived w(z)
  tracks DESI's published curve within ~0.05 across 0 < z < 2.
  (`estif_task5b_cosmo_eos.py`; author-reproduced)
- **Task 6 — frozen-eddy reframe:** both naive derived eddy EoS models fail
  badly (E1 conserved-spin, stiff w = +1: χ²/N = 3232; E2 tracker: 754). The
  **frozen eddy** (constant density → de Sitter → w = −1, free from Task 4)
  scores **1.92 — ties ΛCDM and beats the tilt's 3.35**. Tilt cosmology
  retired. (`estif_task6_eddy_eos.py`)
- **The Split (v6.3):** Path One = ESTIF-Core (derived gravity + frozen eddy =
  Λ; publishable); Path Two = ESTIF-Extended (derive a small thawing on top of
  w = −1 from the vorticity stress tensor).

### C. Path One canonical tests
- Parity script: ESTIF-Core with Ωm = x₀ (imported r_u): χ²/N = 1.965 vs ΛCDM
  1.919. (`estif_pathone_cosmology.py`)
- AIC/BIC: ESTIF-Core vs ΛCDM **ΔAIC = +0.59 (indistinguishable)**; retired
  tilt vs frozen eddy **ΔAIC = +18 (k=0) / +22 (k=2) — decisively worse**;
  fitted CPL AIC = 12.54 → **Path Two bar:** a derived (k=0) thawing must reach
  χ²/N < 1.965 to beat frozen, **< 0.965** to beat the fitted CPL.
  (`estif_pathone_aic_bic.py`)

### D. Adversarial errata — v6.3.1 (CORRECTIONS_v6.3.1.md)
| # | Correction | Severity |
|---|---|---|
| C1 | Exact-Schwarzschild vacuum vs EHT/LISA deviation claims → deviations conditional on the un-derived eddy-stress sector | **Critical** |
| C2 | Ωm = x₀ → consistency relation (r_universe is the Ωm-dependent ΛCDM horizon) | **Critical** |
| C3 | JWST spec: σ8/S8 hard filter added (persistent +13% growth → σ8 = 0.917, ~18σ; boost must be transient/scale-dependent) | High |
| C4 | Path Two target (0.66) under-marginalized (rd, H₀, Ωm fixed) → re-derive before committing | High |
| C5 | Task 6 E1 explanatory text sloppy → clean replacement block provided | Low |
| C6 | "Derived, not borrowed" precision: axioms select GR's constraint sector in PG gauge; the matter coupling (8πG) is adopted | **Critical** |

### E. JWST early-structure specification
- ΛCDM baseline (colossus, Planck18, Sheth–Tormen): n(>10¹¹ M☉, z = 9.1) =
  2.3×10⁻⁵ Mpc⁻³; n(>3×10¹¹) = 8.0×10⁻⁷.
- Growth sensitivity: **+13% in D at z = 9.1 → 4.8× (>10¹¹) and 7.8× (>3×10¹¹)
  abundance boost**; 5× reservoir ← g = 1.13; 10× ← g = 1.20.
- σ8 hard filter (C3): the enhancement must be **two-sided** — ≥13% at z ≈ 9,
  → 1 by z ≲ 2 (or scale-dependent).
- Path assignment corrected: the boost cannot come from H(z) (matter-dominated
  at z ≈ 9 for both paths) — it is a **gravity-sector** (perturbation)
  prediction. Missing input: D_ESTIF(z), laptop-analytic; **not N-body
  blocked**. (`JWST_TEST_SPEC.md`, `estif_jwst_growth_spec.py`)

### F. The Ωm bootstrap and closure (Door 2 — newest)
- Principle P ("Ωm = R_H/r_p") inverts the C2 circularity into a closed
  equation **Ωm·I(Ωm) = 1** with a **unique** root:
  pure matter+Λ (zero inputs): **0.3043** (2.2% from Planck); with radiation
  (inputs = measured T_CMB, N_eff, h): **0.31408** — **0.96% from Planck,
  0.53σ inside Planck's error bar**. Back-predicts r_universe = 4.353×10²⁶ m
  (−1.07% vs the 4.4×10²⁶ import). Identity: P ⇔ mean-matter pull at the
  horizon = cH₀/2 — the same cH₀ that sets a₀ (ratio computes to 1.00000).
  (`estif_omega_bootstrap.py`)
- **Closure:** propagating the bootstrap Ωm: a₀ = 1.1920×10⁻¹⁰ — MOND agreement
  improves **1.72% → 0.66%** (SPARC insensitive: v_flat ×1.00269); DESI DR2
  with bootstrap Ωm: **χ²/N = 1.618** vs ΛCDM 1.919 — *within the fixed-(H₀,
  rd) test*; under C4-style marginalization the ordering could change. **Input
  ledger after adopting P: measured = {H₀, T_CMB, N_eff}; computed = {Ωm, Ω_Λ,
  x₀, r_universe, a₀}.** (`estif_bootstrap_closure.py`)
- Literature adjacency flagged: Gaztañaga's causal-universe scale
  (≈ 0.3176 H₀, via inflation) — comparison **required** before any novelty
  claim.

### G. Lineage and naming
- v6.3 is **not** the FD or EC fork: FD (2024, exponential shrinkage, ruled
  out) → siblings EC (June 2025; shrinkage + variable H(t); its vortex
  late-acceleration ≙ Task 6's falsified E1) and the main line → v6.3
  (flat-slice flow; shrinkage retired; field equation derived).
- Naming recommendation: lineage tag **ESTIF-FS** (drop the redundant second
  F); paths keep the published **ESTIF-Core / ESTIF-Extended**. Adoption =
  author's decision, pending.

## 4. Honest flags carried by the above (none hidden)

1. Bootstrap results are **conditional on principle P** (Part B: P not derived
   from A1–A3; three candidate routes, none attempted in earnest).
2. Gaztañaga comparison pending → no novelty claims on Ωm.
3. DESI 1.618 is a **fixed-ruler** result; C4 marginalization pending.
4. a₀'s empirical target carries ~10% scatter; 0.66% is pleasing, not decisive.
5. Strong-field deviation claims (EHT/LISA) are conditional (C1) until the
   eddy-stress sector is derived.
6. Preferred-frame nature of A2 needs its explicit section (checklist C-11).
7. Peer review remains the only external test not yet faced.

## 5. Script inventory (13 delivered; ensure all present in `tests/`)

| Script | One-line purpose |
|---|---|
| estif_tmunu_gauss_codazzi.py | ADM/Gauss–Codazzi engine, validated vs FRW + de Sitter |
| estif_flow_signature_dynamics.py | Signature + SR + Newton from Euclidean bulk + universal-c (18/18) |
| estif_converse_flow_law.py | Vacuum forces v² = 2A/r; SdS one-flow bonus |
| estif_fidelity_audit.py | Axiom-presence scan of the corpus (71 files) |
| estif_task4_field_equation.py | Field equation derived: mass continuity = Poisson (5/5) |
| estif_task5_desi_selfconsistent.py | De-circularized DESI test (10.8 → 3.35) |
| estif_task5b_cosmo_eos.py | DESI-preferred w(z); tilt tracks it within ~0.05 |
| estif_task6_eddy_eos.py | E1/E2 falsified; frozen-eddy reframe (1.92) |
| estif_pathone_cosmology.py | Canonical Path One parity test |
| estif_pathone_aic_bic.py | AIC/BIC parity + Path Two quantified bar |
| estif_jwst_growth_spec.py | JWST baseline + growth→abundance + σ8 filter |
| estif_omega_bootstrap.py | Ωm·I(Ωm) = 1: unique root, 0.53σ from Planck |
| estif_bootstrap_closure.py | Bootstrap Ωm propagated: a₀ 0.66%, ledger |

Housekeeping: two converse-law variants exist in `tests/`
(`estif_converse_flow_law.py`, `estif_converse_flow_law2.py`) — consolidate to
one canonical version.

---

## 6. DOCUMENT UPDATE MATRIX

### Tier A — already delivered, awaiting application to the repo (the P-2 backlog)

| Deliverable | Action | Target |
|---|---|---|
| MILESTONE_v6.3_THE_SPLIT.md | add (new) | repo root |
| README.md | replace | repo root |
| CITATION.cff | replace | repo root |
| STATUS.md | replace | docs/report/ |
| ESTIF_CONCEPT.md | replace | docs/report/ |
| VALIDATION_REPORT.md | replace | docs/report/ |
| SUMMARY_FOR_REVIEW.md | replace | docs/ |
| RHAC_v6.3_APPEND.md | append | docs/plan/RHAC.md |
| ROADMAP_v6.3_APPEND.md | append | docs/guide/ROADMAP.md |
| CHANGELOG_v6.3_PREPEND.md | prepend | CHANGELOG.md |
| CORRECTIONS_v6.3.1.md | add + **apply C1–C6 edits** across the files above | repo root |
| PATH_ONE_CHECKLIST.md | add (new) | docs/plan/ |
| JWST_TEST_SPEC.md | add (new) | docs/report/ |
| 13 scripts | ensure present | tests/ |

### Tier B — delivered docs now needing a v6.3.2 amendment (bootstrap + closure are not yet in any doc)

| Document | Exact amendment needed |
|---|---|
| PATH_ONE_CHECKLIST.md | B-4: 🔶 text → "conditional bootstrap prediction: Ωm·I(Ωm)=1 unique root 0.3141 (0.53σ Planck); closure: a₀ 0.66%, r_u −1.07%, DESI 1.618 (fixed-ruler); Part B (derive P) open — 3 candidate routes." Add new item B-4b: "Gaztañaga comparison memo" ⬜. |
| docs/plan/RHAC.md | New **Scenario Q — the Ωm bootstrap (Door 2)**: resolution of the C2 circularity into a conditional prediction; Scenario H updated: "sharpened — reduces to deriving principle P; routes: flow-budget amplitude, horizon-acceleration balance (P ⇔ g_horizon = cH₀/2), homogeneous field equation (needs vorticity attachment shared with Path Two)." |
| docs/report/STATUS.md | Cosmology table: add row "Ωm bootstrap (conditional on P): 0.3141, 0.53σ, DESI 1.618 fixed-ruler ✅NEW"; Dark-matter table: update Ωm = x₀ row per C2 + bootstrap; Priority actions: insert "Gaztañaga memo" and "attempt Part B". |
| docs/report/VALIDATION_REPORT.md | New subsection under Part 3 (or new Part 3.0): "The Ωm bootstrap" with the equation, both roots, sensitivity note, cH₀/2 identity, ledger, and all four honest flags. |
| docs/SUMMARY_FOR_REVIEW.md | "What is new" bullet for the bootstrap + ledger; add reviewer question #9: "Is principle P (Ωm = R_H/r_p) derivable from the flow axioms, and how does it relate to Gaztañaga's causal-boundary results?" |
| README.md | Key-results: add bootstrap line + input ledger; note a₀ refinement 1.72% → 0.66% (conditional). |
| CITATION.cff | Abstract: one clause for the conditional bootstrap; bump 6.3.0 → 6.3.2 on commit. |
| CHANGELOG.md | New **[6.3.2]** entry: "Ωm bootstrap (conditional on P): unique root 0.3141 (0.53σ); closure: a₀ 0.66%, r_universe output, DESI 1.618 fixed-ruler; ledger {H₀, T_CMB, N_eff}; Gaztañaga adjacency flagged." |
| estif_pathone_cosmology.py | Apply C2 wording (already specified in errata); optionally add a `--bootstrap` note pointing to the closure script as the successor test. |
| MILESTONE_v6.3_THE_SPLIT.md | **No change** — keep frozen as the split-day snapshot; the bootstrap enters via CHANGELOG + STATUS. |

### Tier C — new documents to create

| New document | Purpose | Priority |
|---|---|---|
| GAZTANAGA_COMPARISON.md | Same-chest-or-different memo vs the causal-universe papers; prerequisite for any Ωm novelty claim | High (before any Ωm claim ships) |
| Preferred-frame section (paper + concept doc) | Checklist C-11 / barrier 5: A2's absolute bulk time, simultaneity reinterpreted | High (pre-submission) |
| NAMING.md or README note | Adopt/decline ESTIF-FS lineage tag; record Core/Extended taxonomy | Low (author decision) |

### Recommended commit sequence
1. **v6.3** — apply Tier A (the Split suite).
2. **v6.3.1** — apply C1–C6 edits per CORRECTIONS_v6.3.1.md.
3. **v6.3.2** — apply Tier B amendments (bootstrap + closure) + add Tier C
   placeholders.
This ordering keeps each version's claims internally consistent at every commit.

---

**Report version:** 1.0 · 9 July 2026 · All figures traceable to the scripts in §5.
