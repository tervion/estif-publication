# ESTIF v6.0 — Work in Progress

This directory contains predictions that are either:
1. Below current detection thresholds
2. Awaiting observational data or computational resources
3. Analytically motivated but not yet confirmable

---

## Current Contents

### `galaxy_asymmetry_prediction.png`
**Status:** Below detection threshold — not currently testable  
**Origin:** v2.0 friction-drag framework (October 2025)  
**Note:** The galaxy asymmetry prediction has not been re-derived in v6.0.
In the current framework, the formula is dormant at galactic scales
(x_local ≈ 10⁻⁷, correction ≈ 0). A galactic-scale asymmetry prediction
from the eddy background requires Phase 7.2 N-body simulation.

### `jwst_ceers_comparison.png`
**Status:** Old framework — moved here from validated/  
**Origin:** v2.0 friction-drag framework comparison with CEERS data  
**Note:** Not a v6.0 result. Retained for historical reference only.
A v6.0 JWST comparison would require the ISW prediction and CMB work (Phase 6).

---

## What Will Go Here Next

### N-body Simulation Prediction (Phase 7.2)
**Predicted:** ESTIF halos should reach overdensity δ ~ 50,000–100,000 × background  
**Test:** N-body simulation with force law a = −∇(ω²/2) instead of −GM/r²  
**Resource:** University computing cluster or cloud HPC  
**Status:** Collaboration target — cannot be computed on a Mac Mini  
**Falsifiable:** If simulation gives δ in this range → v_flat = 220 km/s follows

### CMB Acoustic Scale Check (Phase 6.1)
**Test:** Does ESTIF's modified H(z) shift the CMB acoustic scale by < 0.5%?  
**Script:** `tests/test_cmb_angle_estimate.py` (pending)  
**Expected:** Small shift (Ω_tilt capped at z=2, H(z) identical to ΛCDM at z>2)

### DESI w(z) Comparison (Phase 5.3)
**Test:** ESTIF w = −1.08 vs published DESI DR2 w(z) measurements  
**Script:** `tests/test_desi_comparison.py` (pending)  
**Expected:** Consistent with DESI hints of w < −1

---

## Priority Order for Next Scripts (Mac Mini doable)

1. `test_desi_comparison.py` — ESTIF w=−1.08 vs DESI DR2 (Phase 5.3)
2. `test_cmb_angle_estimate.py` — CMB acoustic scale sanity check (Phase 6.1)
3. `test_isw_prediction.py` — ISW from Ω_tilt evolution (Phase 6.2)

---

**Document Version:** 6.0 | **Last Updated:** 17 March 2026

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2
