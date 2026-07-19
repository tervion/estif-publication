<!--
  PREPEND THIS ENTRY TO THE TOP OF CHANGELOG.md
  (below the title/header, above the existing v6.2 entry)
-->

## [6.3.0] — 2026-07-08 — "The Split"

Major milestone. The project forks into **Path One (ESTIF-Core, clean)** and
**Path Two (ESTIF-Extended, hard)**. See `MILESTONE_v6.3_THE_SPLIT.md`.

### Added — gravity field equation is now DERIVED (Task 4)
- The field equation `rho_eff = m′(r)/(4πr²)` (mass continuity = Poisson in
  integrated form) is now **forced** by the flow axioms via the Gauss–Codazzi
  engine — no longer matched to the Schwarzschild solution.
  - Vacuum (rho_eff = 0) → v² = 2GM/r uniquely → **exact Schwarzschild**.
  - Uniform-density ball → rho_eff = rho0 **exactly** (correct Newtonian source).
  - The former "D2 Poisson postulate" is now a theorem for vacuum, Newton, and
    Schwarzschild.
- New scripts: `estif_task4_field_equation.py` (5/5),
  `estif_flow_signature_dynamics.py` (18/18: signature + SR + Newton from a
  Euclidean bulk with a universal speed-c constraint),
  `estif_converse_flow_law.py` (Birkhoff's theorem restated in flow variables),
  `estif_tmunu_gauss_codazzi.py` (ADM engine validated vs flat FRW and de Sitter).

### Changed — cosmology reframed; tilt apparatus retired (Tasks 5, 5b, 6)
- **Circularity fixed:** `x(z)` no longer uses ΛCDM as its own ruler. The
  self-consistent fixed-point solve (`estif_task5_desi_selfconsistent.py`) drops
  the DESI DR2 fit from χ²/N = 10.80 (circular) to **3.35** (self-consistent).
- **Frozen-eddy reframe:** the constant-eddy limit (w = −1, derived from Task 4)
  scores **χ²/N = 1.92, tying ΛCDM and beating the tilt formula's 3.35**. On DESI
  the entire Ω_tilt apparatus (N_MAX, B, sign-flip, z<2 cutoff) is a net negative
  and is **retired** from the cosmology claim under Path One.
  (`estif_task6_eddy_eos.py`.)
- **DESI-preferred w(z):** the self-consistent tilt already thaws toward the
  DESI-preferred curve (within ~0.05); best-fit evolving-w flow reaches χ²/N =
  0.66 (Path Two target). (`estif_task5b_cosmo_eos.py`.)

### Deprecated / retired
- `Ω_tilt(z)` as a cosmology claim (moved to an "explored and set aside"
  appendix under Path One). The tilt formula's local strong-field use (EHT, the
  n(x) dynamic exponent) is unaffected; only the *cosmological* dark-energy claim
  is retired.
- The z < 2 hard cutoff (scaffolding for the retired tilt term) is moot under
  both paths.

### Fixed / flagged (fidelity audit)
- `estif_fidelity_audit.py` found that axioms A2 (universal speed c) and A3
  (empty space is not a source) are present only in the derivation scripts, and
  that `ESTIF_CONCEPT.md` runs a competing "shrinking-ruler" narrative (A1
  conflict). Path One writing tasks: add A2 + A3 to the theory, retire the
  shrinking-ruler narrative, relabel v_flow = cx₀ ≈ 0.31c as the *sideways
  component* of a total-c flow (not the full flow speed).

### Portability
- The DESI test scripts now cache data beside the script
  (`os.path.dirname(__file__)`) instead of a hard-coded absolute path, so they
  run unmodified on any machine.

### Notes
- The gravity letter is unaffected and **strengthened**: the a₀/MOND derivation
  now rests on a derived field equation rather than a Schwarzschild match.
