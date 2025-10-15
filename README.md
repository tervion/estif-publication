# ESTIF-Gravity: Strong-Field Modifications to General Relativity

[![Status](https://img.shields.io/badge/status-testable_prediction-green)]()
[![Version](https://img.shields.io/badge/version-2.0-blue)]()
[![Tests](https://img.shields.io/badge/tests-3%2F3_passing-success)]()

## 🎯 Key Result

**LISA-Detectable Gravitational Wave Signal:**
- **32 microsecond delays** in binary black hole mergers  
- **3.2σ detection significance**  
- **Testable when LISA launches (~2034-2037)**

This is currently the strongest testable prediction of ESTIF-Gravity.

---

## 📋 Project Overview

**ESTIF-Gravity** tests whether friction-like corrections to General Relativity produce observable effects in strong gravitational fields. Unlike the earlier ESTIF-FD version (which attempted to replace ΛCDM cosmology and was ruled out by supernova data), this approach:

- ✅ **Accepts standard ΛCDM cosmology**
- ✅ **Tests only strong-field gravity modifications**
- ✅ **Makes specific, falsifiable predictions**

### What Changed from ESTIF-FD?

| Aspect | ESTIF-FD (v1.0) | ESTIF-Gravity (v2.0) |
|--------|-----------------|----------------------|
| **Cosmology** | Custom exponential S(t) | Standard ΛCDM |
| **Status** | Ruled out (χ²=3.8× worse) | Uses validated cosmology |
| **Focus** | Universe expansion | Strong-field gravity only |
| **Key Prediction** | None testable | LISA GW delays (3.2σ) |
| **Archived at** | [Zenodo](https://zenodo.org/records/17261725) | Current version |

---

## 🔬 Scientific Approach

### Core Hypothesis

Friction-drag corrections modify light bending and gravitational wave propagation near massive objects via:

```
θ_ESTIF = θ_GR × (1 + β × R_s/(2r))
```

Where:
- `θ_GR` = General Relativity prediction
- `β` = friction coefficient (0.05)
- `R_s` = Schwarzschild radius
- `r` = distance from object

### Three Testable Predictions

| Observable | Prediction | Detector | Status |
|------------|-----------|----------|--------|
| **GW merger delays** | 32 μs | LISA | ✅ **3.2σ detectable** |
| **BH shadow size** | 1.67% larger | next-gen EHT | ⚠️ Marginal (needs ~0.3% precision) |
| **Galaxy asymmetry** | 0.0001% | JWST | ❌ Below threshold |

**Bottom line:** Gravitational wave timing with LISA is the key test.

---

## 📊 Results Summary

### Validation Tests (All Pass)

- ✅ **Weak-field GR compliance:** <1% deviation in solar system tests
- ✅ **Friction scaling:** Confirmed ∝ M, ∝ 1/r³
- ✅ **Mathematical consistency:** All equations self-consistent

### Observational Predictions

**1. LISA Gravitational Waves** (Strongest Prediction)
- Prediction: 32 microsecond delays in binary black hole mergers
- Detection significance: 3.2σ
- Timeline: LISA launch ~2034-2037
- **Status: TESTABLE** ✅

**2. EHT Black Hole Shadow**
- Prediction: 1.67% larger shadow than GR predicts for M87*
- Current status: Both GR and ESTIF show 7σ tension with observation
- Timeline: Next-generation EHT by 2030
- **Status: MARGINAL** ⚠️

**3. JWST Galaxy Morphology**
- Prediction: 0.0001% asymmetry at high redshift
- Required precision: 0.5% (with 100 galaxies)
- **Status: UNDETECTABLE** ❌

---

## 🗂️ Repository Structure

```
estif_publication/
├── README.md                    # This file
├── CHANGELOG.md                 # Version history
├── src/                         # Core implementation
│   ├── estif_ec_gr_constants.py
│   ├── estif_ec_gr_model.py
│   └── estif_ec_gr_run_simulation.py
├── tests/                       # Validation & predictions
│   ├── observational/           # EHT, LIGO, JWST comparisons
│   │   ├── compare_eht_m87.py
│   │   ├── compare_ligo_gw.py
│   │   └── compare_jwst_galaxies.py
│   ├── unit/                    # Basic functionality tests
│   └── run_all_comparisons.py   # Main test runner
├── results/                     # Generated plots
│   ├── validated/               # Publication-ready figures
│   └── work_in_progress/        # Future work
├── docs/                        # Documentation
│   ├── STATUS.md                # Technical status
│   ├── VALIDATION_REPORT.md     # Evidence summary
│   └── SUMMARY_FOR_REVIEW.md    # For expert review
├── data/                        # Observational data
└── archive/                     # Old ESTIF-FD version
```

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/tervion/estif-publication
cd estif_publication

# Install dependencies
pip install -r requirements.txt
```

### Run Tests

```bash
# Run all observational comparisons
cd tests
python3 run_all_comparisons.py

# Generated plots appear in tests/ directory
# - eht_m87_comparison.png
# - ligo_gw150914_comparison.png
# - gw_mass_dependence.png
# - jwst_ceers_comparison.png
```

### Run Main Simulation

```bash
# Full validation suite
cd src
python3 estif_ec_gr_run_simulation.py

# Generates plots in src/:
# - friction_scaling.png
# - lensing_comparison.png
# - predictions_summary.png
```

---

## 📈 Key Results

### Gravitational Wave Delays (LISA)

**Prediction:** Binary black hole mergers experience 32 microsecond delays due to friction drag.

```
For 65 M_sun merger:
- ESTIF delay: 3.2×10⁻⁵ s
- LIGO precision: 1 ms (too coarse)
- LISA precision: 10 μs
- LISA S/N: 3.2σ ← Detectable!
```

See: `results/validated/ligo_gw150914_comparison.png`

### Black Hole Shadow (EHT)

**Prediction:** M87* shadow 1.67% larger than GR predicts.

```
For M87* at photon sphere:
- GR: 19.85 μas diameter
- ESTIF: 20.18 μas diameter
- Observed: 42 μas (both models off by 7σ)
- ESTIF vs GR: 0.33 μas difference
- Needs: 0.3% precision (EHT goal: ~1% by 2030)
```

See: `results/validated/eht_m87_comparison.png`

---

## 📚 Documentation

### For Researchers

- **[STATUS.md](docs/STATUS.md)** - Current development status
- **[VALIDATION_REPORT.md](docs/VALIDATION_REPORT.md)** - Evidence summary
- **[CHANGELOG.md](CHANGELOG.md)** - Version history and major changes

### For Expert Review

- **[SUMMARY_FOR_REVIEW.md](docs/SUMMARY_FOR_REVIEW.md)** - One-page overview for experts
- **[Results Directory](results/validated/)** - All plots with descriptions

### Technical Details

- **[estif_ec_gr_model.py](src/estif_ec_gr_model.py)** - Implementation with inline documentation
- **[Comparison Scripts](tests/observational/)** - Full observational analysis code

---

## 🎓 Scientific Context

### Motivation

General Relativity has been extensively validated in weak fields (solar system, binary pulsars) but remains less constrained near black hole horizons. The Event Horizon Telescope's imaging of M87* and LIGO's detection of gravitational waves open new windows for testing gravity in extreme environments.

### Relation to Standard Physics

ESTIF-Gravity:
- ✅ **Preserves:** All weak-field GR predictions
- ✅ **Uses:** Standard ΛCDM cosmology
- 🔬 **Tests:** Whether friction-like corrections appear in strong fields
- ❌ **Not claiming:** To solve dark matter, dark energy, or replace quantum gravity

### Previous Work

An earlier version (ESTIF-FD) attempted to derive cosmic expansion from 4D flow dynamics. This was **ruled out by supernova data:**
- ESTIF-FD: χ² = 1428
- ΛCDM: χ² = 376
- Ratio: 3.8× worse

**Current approach** accepts standard cosmology and tests only strong-field modifications.

---

## 🔬 Current Status

### What's Validated

- ✅ Code works correctly
- ✅ Predictions are self-consistent
- ✅ LISA prediction is concrete and testable
- ✅ Mathematical framework is complete

### What's Not Validated

- ⚠️ Whether the friction mechanism is physically correct
- ⚠️ Whether β = 0.05 is the right value
- ⚠️ Whether GW and EM should use same scaling

### Honest Assessment

**Strengths:**
- Makes specific, testable prediction (LISA)
- Internal mathematical consistency
- Simpler than many alternative gravity theories

**Weaknesses:**
- Developed by non-expert (me)
- Previous cosmology version failed
- May have fundamental issues experts can identify

---

## 🤝 Contributing & Feedback

This is exploratory research by an independent researcher. **Expert feedback is crucial and welcome.**

### Questions for Experts

1. Is the lensing equation θ = θ_GR × (1 + β·R_s/2r) physically reasonable?
2. Should gravitational wave and electromagnetic predictions use the same scaling?
3. Are there obvious errors or fundamental problems?
4. Is this worth pursuing further?

### Contact

- **Repository:** https://github.com/tervion/estif-publication
- **Email:** [your email]
- **Zenodo (old version):** https://zenodo.org/records/17261725

---

## 📄 Citation

If you reference this work, please cite:

```bibtex
@software{angelov2025estif_gravity,
  author = {Angelov, Peter},
  title = {ESTIF-Gravity: Strong-Field Modifications to General Relativity},
  year = {2025},
  version = {2.0},
  url = {https://github.com/tervion/estif-publication}
}
```

Previous version (ESTIF-FD, ruled out):
```bibtex
@software{angelov2024estif_fd,
  author = {Angelov, Peter},
  title = {ESTIF: Emergent Spacetime from Inward Flow},
  year = {2024},
  doi = {10.5281/zenodo.17261725}
}
```

---

## 📜 License

MIT License - See [LICENSE](LICENSE) file for details.

---

## ⚠️ Disclaimer

This is independent research by a non-physicist. The framework makes testable predictions, but has not been peer-reviewed by experts in general relativity or gravitational physics. The code is provided for transparency and to enable expert evaluation.

**The strongest prediction (LISA gravitational wave delays) is falsifiable and will be tested when LISA launches in the 2030s.**

---

**Last Updated:** October 15, 2025  
**Version:** 2.0 (ESTIF-Gravity Fork)  
**Status:** Ready for expert review


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2

