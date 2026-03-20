# ESTIF v6.0 — Work in Progress

This directory contains predictions that are either below current detection thresholds,
awaiting observational data, or blocked by computational resources not available on a
personal machine. Nothing here is speculative — all items below are analytically
motivated and have defined falsifiable criteria.

---

## Current Contents

### `galaxy_asymmetry_prediction.png`
**Status:** Below detection threshold — not currently testable
**Origin:** v2.0 friction-drag framework (October 2025)
**Note:** The galaxy asymmetry prediction has not been re-derived in v6.0.
In the current framework, the tilt formula is dormant at galactic scales
(x_local ≈ 10⁻⁷, correction ≈ 0 to 20 decimal places). A true v6.0 galactic
asymmetry prediction requires the N-body simulation output (Phase 7.2) — the
eddy halo structure will determine what, if any, asymmetry signature exists.

---

## What Will Go Here — Mac Mini Doable (Scripts Pending)

These three plots can be generated on a Mac Mini once the scripts are written.
They are the next analytical steps after v6.0.

### `desi_w_comparison.png`
**What it shows:** ESTIF w_eff = −1.08 plotted against published DESI DR2 w(z) data points
**Script to write:** `tests/test_desi_comparison.py`
**Expected result:** Consistent with DESI hints of w < −1 across z = 0.2–1.6
**Why it matters:** DESI DR2 (2024) reported hints of evolving dark energy. ESTIF predicted
w = −1.08 independently. Direct comparison would be a significant external validation.
**Blocked by:** Script not yet written (≤ 1 day of work)
**Phase:** 5.3

### `cmb_angle_estimate.png`
**What it shows:** ESTIF's modified H(z) effect on the CMB acoustic scale θ_s
**Script to write:** `tests/test_cmb_angle_estimate.py`
**Expected result:** Shift < 0.5% from ΛCDM (Ω_tilt capped at z=2, H(z) = ΛCDM at z > 2)
**Why it matters:** The CMB acoustic angle θ_s = r_s/D_A is measured to 0.03% precision
by Planck. If ESTIF shifts it by > 0.1%, that's a problem. If < 0.1%, it's a pass.
**Blocked by:** Script not yet written (≤ 1 day of work)
**Phase:** 6.1

### `isw_prediction.png`
**What it shows:** Integrated Sachs-Wolfe effect from Ω_tilt(z) evolution
**Script to write:** `tests/test_isw_prediction.py`
**Expected result:** Small ISW signal from the evolving dark energy (Ω_tilt peaks at z≈0.5)
**Why it matters:** The ISW effect is a direct signature of dark energy evolution.
ESTIF's bell-shaped Ω_tilt(z) produces a specific ISW pattern distinguishable from ΛCDM.
**Blocked by:** Script not yet written (≤ 1–2 days of work)
**Phase:** 6.2

---

## What Will Go Here — Requires N-body Simulation (Collaboration Target)

These plots require a university computing cluster or cloud HPC.
They cannot be computed on a personal machine. Each has a defined
falsifiable prediction — the simulation either confirms or rules out ESTIF.

### `flat_rotation_curves.png`
**What it shows:** v_flat vs galactic radius for simulated ESTIF halos
**Predicted result:** v_flat ≈ 220 km/s for MW-mass halos (currently gives ~2.7 km/s
without virialization — 82× below observed)
**Falsifiable criterion:** Does the ESTIF force law a = −∇(ω²/2) produce halos that
virialize to the correct circular velocity?
**Blocked by:** N-body simulation — requires tracking ~10⁷ particles over ~10,000 timesteps
**Phase:** 7.2

### `halo_concentration.png`
**What it shows:** Internal halo overdensity δ = ρ_halo/ρ_eddy vs halo mass
**Predicted result:** δ ~ 50,000–100,000 × background (derived analytically from
virial condition + free-fall time arguments)
**Falsifiable criterion:** If δ is in this range → flat rotation curves follow automatically.
If δ << 50,000 → ESTIF cannot explain v_flat without additional physics.
**Blocked by:** N-body simulation
**Phase:** 7.2

### `tully_fisher_simulation.png`
**What it shows:** v_flat ∝ M^α from ESTIF halos — does α = 1/4 emerge from simulation?
**Context:** Analytically, ESTIF gives M^(1/3) from geometry and M^(1/4) from the MOND
limit. The simulation will show which exponent the actual halo dynamics produce.
**Falsifiable criterion:** α = 1/4 (MOND) would confirm the ESTIF-MOND connection.
α = 1/3 (pure geometry) would mean additional physics is needed.
**Blocked by:** N-body simulation
**Phase:** 7.2

### `bullet_cluster_offset.png`
**What it shows:** Spatial separation between simulated eddy mass distribution and
simulated gas distribution during a cluster merger
**Context:** The Bullet Cluster shows ~0.7 Mpc offset between X-ray gas and lensing mass.
In ESTIF, the eddy background is collisionless (σ/v_esc = 0.5) so it should pass through
a merger like dark matter does — producing the offset.
**Falsifiable criterion:** Does the eddy background produce the correct offset magnitude
and direction? This is one of the strongest tests of any dark matter alternative.
**Blocked by:** N-body simulation + hydrodynamics (most demanding computation in the project)
**Phase:** 7.3 (future, after 7.2)

---

## Priority Order

| Priority | Plot | Blocker | Time estimate |
|---|---|---|---|
| 1 | `desi_w_comparison.png` | Script needed | ≤ 1 day |
| 2 | `cmb_angle_estimate.png` | Script needed | ≤ 1 day |
| 3 | `isw_prediction.png` | Script needed | 1–2 days |
| 4 | `flat_rotation_curves.png` | N-body simulation | Collaboration |
| 5 | `halo_concentration.png` | N-body simulation | Collaboration |
| 6 | `tully_fisher_simulation.png` | N-body simulation | Collaboration |
| 7 | `bullet_cluster_offset.png` | N-body + hydro | Future |

---

## How to Get the N-body Simulation Done

Three realistic paths, in order of likelihood:

**Path 1 — arXiv first, simulation second.** Publish the analytical results now.
Once the paper exists, computational cosmology groups have something to cite and test.
The paper explicitly invites simulation collaboration with co-authorship.

**Path 2 — Cloud HPC.** A small test run (10⁶ particles, simplified box) costs roughly
$200–500 on AWS/Google Cloud and takes a few days. Not a full cosmological simulation
but enough to check whether δ ~ 50,000 is achievable.

**Path 3 — University collaboration.** Approach a computational cosmology group with
the paper and the falsifiable prediction. Groups using Gadget-4, AREPO, or SWIFT can
run this in weeks on existing hardware.

---

**Document Version:** 6.0 | **Last Updated:** 17 March 2026
