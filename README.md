# ESTIF: Emergent Spacetime from Inward Flow

**Version 6.4.2** · 21 July 2026
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17261724-blue)](https://zenodo.org/records/17261724)

A geometric framework in which **3D space is an even hypersurface carried through a 4D bulk at the speed of light.** Gravity is the local tilt and slowing of that flow near mass — and its field equation is *derived* from the axioms, giving exact Schwarzschild in vacuum (the matter coupling, 8πG, is adopted — C6). Time is motion through the 4th dimension; time dilation is that motion partly diverted sideways. Cosmology is stated honestly: an imported cosmological constant that ties ΛCDM (χ²/N = 1.92), with the residual dark-energy sector *proven empty* (Phase 2, w = −1 exactly). The MOND scale a₀ is a horizon quantity — its scale (~c·H) derived, its O(1) prefactor the central open problem.

## Status at a glance

- **Gravity:** field equation derived (Task 4); SPARC BTFR RMS 15.6% (87 galaxies); c_gw = c derived (GW170817 structural pass). Strong-field interior sector open.
- **Cosmology:** Λ imported, ties ΛCDM; Ω_tilt retired; growth restored under A1′ (f(0.5) = 0.76, DESI-consistent); registered kill-shot: mean spatial curvature ≡ 0 at all epochs.
- **Dark matter:** analytics complete; Ωm = 0.31408 bootstrap conditional on a non-derivable postulate (no Ωm novelty — Gaztañaga 2019–2023 has priority); halo structure behind the N-body wall.
- **Publication:** gravity letter drafted; send gated on the Bullet Cluster response.

## Documentation (four living files)

| File | Contents |
|---|---|
| [`docs/PROJECT.md`](./docs/PROJECT.md) | Status, the 35-item checklist, priorities, repository map |
| [`docs/SCIENCE.md`](./docs/SCIENCE.md) | The framework and the evidence — axioms, derivations, receipts, limitations |
| [`docs/RHAC.md`](./docs/RHAC.md) | The decision archive (Scenarios A–Q, RHAC-001…011) |
| [`CHANGELOG.md`](./CHANGELOG.md) | Reverse-chronological history, v6.1.0 → v6.4.2 |

Test-suite index: [`tests/docs/TEST_INDEX.md`](./tests/docs/TEST_INDEX.md).

## Reproduce the key results

Every script runs on a laptop; the DESI/SPARC scripts fetch real data on first run and cache it in `tests/docs/`.

```bash
# Gravity: the field equation is derived, not postulated
python3 tests/scripts/estif_task4_field_equation.py       # 5/5  mass continuity = Poisson
python3 tests/scripts/estif_flow_signature_dynamics.py    # 18/18 signature + SR + Newton
python3 tests/scripts/estif_converse_flow_law.py          # vacuum forces v^2 = 2A/r
python3 tests/scripts/estif_tmunu_gauss_codazzi.py        # ADM engine vs FRW + de Sitter

# Cosmology: circularity fix, DESI parity, Phase-2 receipts
python3 tests/scripts/estif_task5_desi_selfconsistent.py  # DESI DR2: 10.8 -> 3.35
python3 tests/scripts/estif_task5b_cosmo_eos.py           # DESI-preferred w(z)
python3 tests/scripts/estif_task6_eddy_eos.py             # constant-term limit (1.92) beats tilt
python3 tests/scripts/phase2_a1prime/estif_a1prime_recensus.py   # A1' re-census PASS 5/5

# a0-horizon doctrine: the scale is derived, the number is not
python3 tests/scripts/a0_horizon_test.py
python3 tests/scripts/a0_prefactor_derivation.py
python3 tests/scripts/estif_flow_sim.py                   # local flow cannot make a0
python3 tests/scripts/estif_horizon.py                    # horizon background reproduces BTFR

# Housekeeping + prior analytical suite
python3 tests/scripts/estif_fidelity_audit.py
python3 src/estif_ec_gr_run_simulation.py                 # 21/21
```

## Citation

```bibtex
@software{angelov2026estif,
  author  = {Angelov, Peter},
  title   = {ESTIF: Emergent Spacetime from Inward Flow — v6.4.2},
  year    = {2026},
  version = {6.4.2},
  url     = {https://github.com/tervion/estif-publication},
  doi     = {10.5281/zenodo.17261724}
}
```

**Author:** Peter Angelov (Independent Researcher) · tervion@gmail.com
**License:** MIT
