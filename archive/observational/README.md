# Observational Comparison Scripts

This directory contains three scripts that compare ESTIF-Gravity predictions with actual observational data from EHT, LIGO, and JWST.

---

## Overview

Each script:
1. Loads real observational data
2. Calculates ESTIF-Gravity predictions
3. Compares with General Relativity baseline
4. Generates publication-quality plots
5. Assesses detectability and statistical significance

**To run all comparisons:**
```bash
cd tests
python3 run_all_comparisons.py
```

This executes all three scripts and generates plots in the current directory.

---

## 1. `compare_eht_m87.py` - Event Horizon Telescope

### What It Does

Compares ESTIF-Gravity's black hole shadow prediction with EHT M87* observations.

### Key Physics

**ESTIF-Gravity lensing formula:**
```python
θ_ESTIF = θ_GR × (1 + β × R_s/(2r))
```

At the photon sphere (r = 1.5 R_s):
- **GR prediction:** 19.85 μas
- **ESTIF prediction:** 20.18 μas (1.67% larger)
- **EHT observation:** 42 ± 3 μas

### Why Both Models Are Wrong

Both GR and ESTIF predict ~20 μas, but EHT observes 42 μas (7σ tension). This is because:
- Simple Schwarzschild metric used (no spin)
- Black hole is actually rotating (Kerr metric needed)
- Accretion disk contributions not included
- Plasma effects near photon sphere

### Key Finding

The 0.33 μas difference between GR and ESTIF requires:
- Next-generation EHT precision (~0.3%)
- Expected timeline: 2030s
- **Status:** Marginally testable (5-10 years)

### Outputs

**Generated plots:**
- `eht_m87_comparison.png` - Shadow diameter comparison
- Saved to: `tests/eht_m87_comparison.png` (if run standalone)
- Or: `results/validated/eht_m87_comparison.png` (via main simulation)

**Console output:**
```
=== EHT M87* BLACK HOLE SHADOW COMPARISON ===

M87* Parameters:
  Mass: 6.5e+09 M_sun
  Distance: 16.8 Mpc
  Schwarzschild radius: 1.92e+13 m

Predictions at photon sphere (r = 1.5 R_s):
  GR (Schwarzschild): 19.85 μas
  ESTIF-Gravity: 20.18 μas
  Difference: 0.33 μas (1.67%)
  
EHT Observation: 42.0 ± 3.0 μas
  GR deviation: -22.2 μas (-7.38σ)
  ESTIF deviation: -21.8 μas (-7.27σ)
  
Status: Both models show large tension (likely due to spin/disk)
Detectability: Marginal (needs 0.3% precision, EHT+ by 2030)
```

---

## 2. `compare_ligo_gw.py` - LIGO/LISA Gravitational Waves

### What It Does

Calculates gravitational wave propagation delays predicted by ESTIF-Gravity and compares with detector sensitivities.

**This script produces the PRIMARY RESULT of ESTIF-Gravity.**

### Key Physics

**ESTIF-Gravity GW delay formula:**
```python
Δt = (G × M × β / c³) × (distance / Mpc)
```

For GW150914 (65 M☉ at 410 Mpc):
- **Predicted delay:** 32.03 μs
- **LIGO precision:** ~1 ms (insufficient)
- **LISA precision:** ~10 μs (sufficient)
- **Signal-to-noise:** 3.2σ

### Why This Is The Key Test

1. **Clear signal:** 32 μs is well above LISA threshold (10 μs)
2. **High significance:** 3.2σ detection confidence
3. **Definitive timeline:** LISA launch ~2034-2037
4. **Falsifiable:** Either detected or model ruled out

### Mass Dependence

The delay scales linearly with merger mass:
- **30 M☉:** 15 μs (1.5σ) - Marginal
- **65 M☉:** 32 μs (3.2σ) - **Detectable**
- **100 M☉:** 49 μs (4.9σ) - Strong detection

Heavier mergers provide higher significance.

### Outputs

**Generated plots:**
- `ligo_gw150914_comparison.png` - **Primary result figure**
- `gw_mass_dependence.png` - Scaling with binary mass

**Console output:**
```
=== LIGO/LISA GRAVITATIONAL WAVE COMPARISON ===

GW150914 Parameters:
  Component masses: 36 M_sun + 29 M_sun
  Total mass: 65 M_sun
  Distance: 410 Mpc
  
ESTIF-Gravity Prediction:
  Time delay: 3.203e-05 s (32.03 μs)
  
Detector Sensitivities:
  LIGO timing: ~1 ms (signal 31× below)
  LISA timing: ~10 μs (signal 3.2× above)
  
Signal-to-Noise Ratio:
  LIGO: 0.032 (undetectable)
  LISA: 3.2σ (DETECTABLE)
  
Status: TESTABLE when LISA launches (~2035)
This is the strongest prediction of ESTIF-Gravity.
```

---

## 3. `compare_jwst_galaxies.py` - JWST High-Redshift Galaxies

### What It Does

Predicts rotation curve asymmetries in high-redshift galaxies due to friction drag and compares with JWST capabilities.

### Key Physics

**ESTIF-Gravity asymmetry formula:**
```python
A_asym = β × (G × M_gal) / (r_half × c²) × 100%
```

For typical galaxy (10¹¹ M☉, r_half = 3 kpc) at z = 6:
- **Predicted asymmetry:** 0.0001%
- **JWST single-galaxy precision:** ~5%
- **JWST statistical precision (N=100):** ~0.5%

### Why This Test Fails

The predicted signal is **5,000× below** even the best statistical precision:
```
Signal:     0.0001%
Threshold:  0.5% (N=100 galaxies)
Ratio:      1/5000
```

Even with:
- Full CEERS survey (~10,000 galaxies)
- Next-generation TMT/ELT telescopes
- Perfect systematic control

The effect remains undetectable.

### Future Prospects

Possible (but difficult) alternatives:
- **Kinematic asymmetries:** IFU spectroscopy (different from morphology)
- **Statistical stacking:** Thousands of galaxies aligned by redshift
- **Indirect signatures:** Morphological correlations

**Current verdict:** Not a viable test.

### Outputs

**Generated plot:**
- `jwst_ceers_comparison.png` - Four-panel analysis

**Console output:**
```
=== JWST HIGH-REDSHIFT GALAXY COMPARISON ===

Typical Galaxy Parameters:
  Mass: 1.0e+11 M_sun
  Half-light radius: 3.0 kpc
  Redshift: z = 3
  
ESTIF-Gravity Prediction:
  Rotation asymmetry: 0.0001%
  
JWST Capabilities:
  Single galaxy precision: 5.0%
  Statistical (N=100): 0.5%
  Full CEERS (N=10000): 0.05%
  
Signal-to-Noise:
  Single galaxy: 2.0e-05 σ
  Statistical: 2.0e-04 σ
  Full survey: 2.0e-03 σ
  
Status: UNDETECTABLE (signal 5000× too weak)
```

---

## Running Individual Scripts

### Run single comparison:
```bash
cd tests/observational
python3 compare_eht_m87.py
python3 compare_ligo_gw.py
python3 compare_jwst_galaxies.py
```

### Run all comparisons:
```bash
cd tests
python3 run_all_comparisons.py
```

### Run full simulation (includes these + validation):
```bash
cd src
python3 estif_ec_gr_run_simulation.py
```

---

## Dependencies

All scripts require:
```python
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
```

And import from:
```python
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from estif_ec_gr_constants import *
from estif_ec_gr_model import *
```

Install via:
```bash
pip install numpy matplotlib
```

---

## Output Locations

**When run standalone:**
- Plots saved to: `tests/observational/[plot_name].png`

**When run via `run_all_comparisons.py`:**
- Plots saved to: `tests/[plot_name].png`

**When run via main simulation:**
- Plots saved to: `results/validated/[plot_name].png`

---

## Interpretation Guide

### EHT M87* Results

✅ **Both models consistent:** GR and ESTIF both predict ~20 μas  
⚠️ **Both have tension:** 7σ deviation from observation  
🔬 **Needs future data:** EHT+ with 0.3% precision (2030s)

**Interpretation:** Not a current test, but could distinguish models in ~10 years.

### LIGO/LISA Results (PRIMARY)

✅ **Clear prediction:** 32 μs delays  
✅ **High significance:** 3.2σ for LISA  
✅ **Definitive timeline:** LISA launch ~2035  
🎯 **This is the key test**

**Interpretation:** If LISA detects 32 μs delays → ESTIF validated. If null → ESTIF falsified.

### JWST Galaxy Results

❌ **Signal too weak:** 0.0001% vs 0.5% threshold  
❌ **Factor 5000× below:** Even statistical studies insufficient  
❌ **No foreseeable test:** Not viable with current technology

**Interpretation:** Interesting theoretical prediction, but not observationally accessible.

---

## Summary Table

| Script | Observable | Prediction | Status | Timeline |
|--------|-----------|-----------|--------|----------|
| `compare_eht_m87.py` | BH shadow | 1.67% larger | ⚠️ Marginal | 5-10 years |
| `compare_ligo_gw.py` | GW delays | **32 μs (3.2σ)** | ✅ **Testable** | **2034-2037** |
| `compare_jwst_galaxies.py` | Galaxy asymmetry | 0.0001% | ❌ Too weak | >2040 |

**Bottom line:** LISA gravitational wave timing is the definitive test of ESTIF-Gravity.

---

## For Developers

### Adding a New Comparison

Template for new observational test:

```python
"""
Compare ESTIF-Gravity prediction with [OBSERVATORY] data.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Import ESTIF-Gravity model
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from estif_ec_gr_constants import *
from estif_ec_gr_model import *

def main():
    print("=== [OBSERVATORY] COMPARISON ===\n")
    
    # 1. Define observational parameters
    observable_value = ...
    observable_uncertainty = ...
    
    # 2. Calculate ESTIF prediction
    estif_prediction = calculate_prediction(...)
    
    # 3. Calculate GR baseline
    gr_prediction = calculate_gr_baseline(...)
    
    # 4. Compare and assess significance
    deviation = estif_prediction - gr_prediction
    sigma = deviation / observable_uncertainty
    
    # 5. Print results
    print(f"Observable: {observable_value} ± {observable_uncertainty}")
    print(f"ESTIF: {estif_prediction}")
    print(f"GR: {gr_prediction}")
    print(f"Significance: {sigma}σ")
    
    # 6. Generate plot
    fig, ax = plt.subplots(figsize=(10, 6))
    # ... plotting code ...
    plt.savefig('new_comparison.png', dpi=150, bbox_inches='tight')
    print("\n✓ Plot saved: new_comparison.png")

if __name__ == "__main__":
    main()
```

### Updating run_all_comparisons.py

Add your script:
```python
scripts = [
    'observational/compare_eht_m87.py',
    'observational/compare_ligo_gw.py',
    'observational/compare_jwst_galaxies.py',
    'observational/compare_new_observatory.py',  # Add here
]
```

---

## Troubleshooting

### Import Errors

```
ModuleNotFoundError: No module named 'estif_ec_gr_constants'
```

**Fix:** Run from correct directory:
```bash
cd tests/observational  # For individual scripts
cd tests                # For run_all_comparisons.py
```

### Missing Plots

If plots don't appear, check:
1. Script completed without errors
2. Write permissions in output directory
3. `matplotlib` installed correctly

### Incorrect Values

If predictions don't match documentation:
1. Verify constants in `estif_ec_gr_constants.py`
2. Check β = 0.05 is set correctly
3. Ensure using `friction_drag_local()` not `friction_drag()`

---

## References

**EHT M87* Data:**
- Event Horizon Telescope Collaboration (2019)
- Astrophys. J. Lett. 875, L1

**LIGO GW150914 Data:**
- Abbott et al. (2016)
- Phys. Rev. Lett. 116, 061102

**JWST CEERS Data:**
- Finkelstein et al. (2022)
- Astrophys. J. Lett. 940, L55

**LISA Specifications:**
- LISA Science Requirements Document
- ESA/SRE(2018)1 revision 1.0

---

## Citation

If you use these comparison scripts in your research, please cite:

```bibtex
@software{angelov2025estif_gravity,
  author = {Angelov, Peter},
  title = {ESTIF-Gravity: Strong-Field Modifications to General Relativity},
  year = {2025},
  version = {2.0},
  url = {https://github.com/tervion/estif-publication}
}
```

---

**Last Updated:** October 16, 2025  
**Version:** 2.0 (ESTIF-Gravity Fork)  
**Status:** All three comparison scripts complete and validated


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2

