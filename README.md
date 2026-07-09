# ESTIF: Emergent Spacetime from Inward Flow

**Version 6.3 — "The Split"** · 8 July 2026
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17261724-blue)](https://zenodo.org/records/17261724)

A geometric framework in which **3D space is a flat hypersurface carried through
a 4D bulk.** From that single picture:

- - **Gravity** is the local tilt/slowing of the flow near mass — and its field
  equation is now *derived* from the axioms rather than matched to the
  Schwarzschild solution (the matter coupling, 8πG, is still adopted — see C6).
- **Time** is motion through the 4th dimension at the speed of light; time
  dilation is that motion being partly diverted sideways by mass.
- **Cosmic expansion** is the projection of the 4D inward flow onto our 3D slice.
- **Dark matter** is the background eddy energy of the moving hypersurface.

---

## 🏁 v6.3 milestone: the project splits into two paths

A sequence of derivations (July 2026) restructured the project. Full record in
[`MILESTONE_v6.3_THE_SPLIT.md`](./MILESTONE_v6.3_THE_SPLIT.md).

### The headline result: gravity is now DERIVED

Earlier versions obtained Newtonian gravity by *matching* the flow profile to the
Schwarzschild solution (setting n = ½) and differentiating it — which presupposes
the answer. **v6.3 removes that.** Starting only from the flow axioms — flat
3-slices, everything moving through the bulk at speed *c*, and "empty space is not
a source" — the Gauss–Codazzi engine *forces* the field equation

```
rho_eff = m'(r) / (4 pi r^2)      (mass continuity = Poisson, integrated form)
```

with these engine-verified consequences:
- **vacuum → exact Schwarzschild** (v² = 2GM/r unique, full Einstein tensor = 0);
- **a uniform matter ball → ρ₀ exactly** (correct Newtonian source).

The old Poisson "postulate" is now a theorem for vacuum, Newton, and Schwarzschild.

> **Precision (C6).** What is *derived* is that the three flow axioms uniquely select
> the constraint (energy) sector of General Relativity in Painlevé–Gullstrand gauge —
> forcing mass continuity, hence exact Schwarzschild in vacuum and the Newtonian
> source in the weak field — **without matching to the Schwarzschild solution**. The
> gravitational *coupling* (geometric scalar ↔ 8πG × energy density) is adopted, not
> derived from below.

### The two paths

| | **Path One — ESTIF-Core (clean)** ✅ default | **Path Two — ESTIF-Extended (hard)** 🔬 |
|---|---|---|
| Gravity | Derived field equation (above) | same |
| Cosmology | Frozen cosmic eddy → **cosmological constant, ties ΛCDM** (χ²/N = 1.92, derived) | Derive a *small thawing* correction to w = −1 from the vorticity stress tensor (target χ²/N ≈ 0.66) |
| Tilt Ω_tilt(z) | Retired to appendix (net negative on DESI) | — |
| Risk | Low, publishable now | High, timeline unknown |

---

## Status at a glance

| Sector | Status |
|---|---|
| Gravity — field equation | ✅ **Derived** (Task 4); strong-field pressure sector open |
| Gravity — a₀ / MOND / SPARC | ✅ Solid (RMS 15.6%, 87 galaxies), now on a derived foundation |
| Gravity — EHT / Λ / LISA (local tilt) | ⚠️ Consistent (0.00σ / ratio 1.0000 / 49.2σ); *deviation* claims conditional — see note |
| Cosmology — honest best | ✅ Frozen eddy = Λ, ties ΛCDM (1.92) |
| Cosmology — tilt Ω_tilt(z) | 🔴 Retired (fits worse than the Λ limit under it) |
| Dark matter | 🟡 Analytical phase complete; N-body wall |

> ⚠️ **Conditional (C1).** The ESTIF vacuum is exactly Schwarzschild (Task 4), so any
> *deviation* from GR in shadows or GW propagation must be sourced by the non-vacuum
> eddy background — a sector not yet derived. The observations remain *consistent*
> with ESTIF; the deviation from GR is what awaits derivation. Planck Λ is a
> calibration match, not a vacuum deviation, and is unaffected.

---

## Reproduce the key results

Every script below runs on a laptop; the DESI scripts fetch real DESI DR2 data
(CobayaSampler) and cache it beside the script.

```bash
# Gravity: the field equation is derived, not postulated
python3 tests/estif_task4_field_equation.py         # 5/5 — mass continuity = Poisson
python3 tests/estif_flow_signature_dynamics.py      # 18/18 — signature + SR + Newton
python3 tests/estif_converse_flow_law.py            # vacuum forces v^2=2A/r (Birkhoff)
python3 tests/estif_tmunu_gauss_codazzi.py          # ADM engine, validated vs FRW + de Sitter

# Cosmology: circularity fix, DESI target, and the frozen-eddy reframe
python3 tests/estif_task5_desi_selfconsistent.py    # DESI DR2: 10.8 (circular) -> 3.35
python3 tests/estif_task5b_cosmo_eos.py             # DESI-preferred w(z); tilt tracks it
python3 tests/estif_task6_eddy_eos.py               # frozen eddy (1.92) beats tilt (3.35)

# Housekeeping: which axioms are actually written into the theory
python3 tests/estif_fidelity_audit.py

# Prior analytical suite (unchanged)
python3 src/estif_ec_gr_run_simulation.py           # 21/21
```

---

## The three project goals (unchanged in spirit, sharpened by v6.3)

1. **Gravity = Time = Eddies.** Now backed by a derived field equation:
   flow-speed gradients *are* gravitational acceleration, and the vacuum
   condition forces exactly Schwarzschild.
2. **Expansion = 4D inward fall.** Honest best form is the frozen-eddy
   cosmological constant (ties ΛCDM); the evolving-dark-energy version is Path Two.
3. 3. **No dark matter / dark energy.** Ωm = x₀ = (c/H₀)/r_universe holds to 0.12%, but
   this is a **consistency relation**, not an Ωm-independent prediction: r_universe is
   the ΛCDM particle horizon, which itself depends on Ωm (C2). a₀ = H₀cx₀/√3 (1.72%
   from empirical MOND). Halo structure needs N-body (documented wall).

---

## What is still open (honest)

- **Strong-field pressure/stress sector** of the gravity T_μν (a relativistic
  interior with pressure). Not needed for vacuum, Newton, or Schwarzschild.
- **Path Two cosmology:** deriving the leading w(z) correction from the full
  vorticity stress tensor. The two naive reductions (stiff spin, expansion
  tracker) are already falsified against DESI.
- ~~Writing tasks (Path One)~~ — **done in v6.3.** A2 and A3 are now stated in
  `ESTIF_CONCEPT.md`; the shrinking-ruler narrative is retired; v_flow = cx₀ is
  relabelled as the sideways component of a total-c flow.
- **Dark-matter halos:** N-body simulation (collaboration target).

---

## Citation

```bibtex
@software{angelov2026estif_v63,
  author  = {Angelov, Peter},
  title   = {ESTIF: Emergent Spacetime from Inward Flow — v6.3 "The Split"},
  year    = {2026},
  version = {6.3},
  url     = {https://github.com/tervion/estif-publication},
  doi     = {10.5281/zenodo.17261724}
}
```

**Author:** Peter Angelov (Independent Researcher) · tervion@gmail.com
**License:** MIT
