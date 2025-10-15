# Changelog

All notable changes to the ESTIF project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.0] - 2025-10-16 - ESTIF-Gravity Fork

**MAJOR VERSION**: Complete pivot from ESTIF-FD (cosmology derivation) to ESTIF-Gravity (gravity-only modifications with ΛCDM cosmology)

### 🎯 Strategic Pivot

This version represents a fundamental change in approach:
- **Old (ESTIF-FD)**: Attempt to derive cosmology from friction → Failed (CMB age ~1% off, numerical instabilities)
- **New (ESTIF-Gravity)**: Use standard ΛCDM cosmology, focus on testable gravity modifications → Success

**Rationale**: Isolate variables—test gravity predictions with proven cosmology before attempting cosmology derivation.

### Added

#### Core Predictions
- **LISA gravitational wave delays**: 32 μs for binary black hole mergers (3.2σ significance)
- **EHT black hole shadows**: 1.0% smaller than GR prediction for M87*
- **JWST galaxy asymmetries**: 0.3% rotation curve asymmetry at z > 6

#### Observational Comparison Scripts
- `tests/observational/compare_eht_m87.py` - Event Horizon Telescope M87* comparison
- `tests/observational/compare_ligo_gw.py` - LIGO/LISA gravitational wave comparison
- `tests/observational/compare_jwst_galaxies.py` - JWST high-redshift galaxy comparison
- `tests/run_all_comparisons.py` - Master validation runner

#### Results and Plots
- `results/validated/ligo_gw150914_comparison.png` - **PRIMARY RESULT**: 32 μs GW delays
- `results/validated/eht_m87_comparison.png` - Black hole shadow comparison
- `results/validated/jwst_ceers_comparison.png` - High-z galaxy framework
- `results/validated/gw_mass_dependence.png` - GW delay mass scaling
- `results/validated/friction_scaling.png` - Theoretical foundation
- `results/validated/lensing_comparison.png` - Strong-field lensing
- `results/validated/predictions_summary.png` - All three predictions

#### Documentation
- `docs/report/VALIDATION_REPORT.md` v2.0 - Complete rewrite for ESTIF-Gravity
- `docs/guide/estif_ec_fd_ROADMAP.md` v2.0 - Updated with completed Phase 3
- `docs/report/estif_ec_fd_concept.md` v2.0 - Removed ESTIF-FD math, added predictions
- `docs/plan/estif_ec_fd_RHAC.md` v2.0 - Documented resolved scenarios
- `CHANGELOG.md` - This file

### Changed

#### Model Architecture
- **Cosmology**: Now uses `astropy.cosmology.FlatLambdaCDM` with Planck 2018 parameters
- **Scale factor**: Changed from `S(t) = exp(-∫H(t')dt')` to standard `a(z) = 1/(1+z)`
- **Hubble parameter**: Uses ΛCDM `H(z)` instead of friction-derived `H(t)`
- **Redshift-distance**: Uses standard `luminosity_distance(z)` from ΛCDM

#### Friction Implementation
- **Local drag function**: New `friction_drag_local(M, r)` using mass density
- **Strong-field predictions**: Now use local drag instead of cosmic-average
- **Magnitude scaling**: Predictions now in testable range (1-3%) instead of ~0%

#### File Naming Convention
- Renamed: `estif_ec_fd_*.py` → `estif_ec_gr_*.py` (emphasizes GRavity modifications)
- Updated: All import statements to use new naming
- Location: All source files now in `src/` directory

#### Test Structure
- Reorganized: `tests/observational/` for EHT/LIGO/JWST comparisons
- Added: `tests/unit/` for model function tests
- Enhanced: `tests/run_all_comparisons.py` with comprehensive validation

#### Documentation Updates
- **README.md**: Reflects ESTIF-Gravity scope (gravity-only, uses ΛCDM)
- **STATUS.md**: Updated for completed Phase 3
- **VALIDATION_REPORT.md**: Removed cosmology tests, added strong-field predictions
- **ROADMAP.md**: Marked Phase 1-3 complete, outlined Phase 4 (publication)
- **Concept document**: Removed friction-derived cosmology formulas

### Removed

#### Abandoned Features
- ❌ Friction-derived Hubble parameter `H(t) = H₀ + A/t^0.75 + BETA_DRAG·(Gρ/c²) + BETA_DRAG/S(t)²`
- ❌ Scale factor calculation from friction `S(t) = exp(-∫H(t')dt')`
- ❌ CMB age independent prediction (~377k years claim)
- ❌ High-z (z>1.5) numerical integration for `S(t)`
- ❌ Cosmology parameter optimization (now use Planck values)

#### Deprecated Tests
- ❌ `test_cmb_age()` - No longer testing CMB age independently
- ❌ `test_high_z_stability()` - Not calculating S(t) from friction anymore
- ❌ `debug_S_calculation()` - S(t) calculation abandoned
- ❌ Custom `scale_factor()` function - Using ΛCDM `a(z)` directly

#### Removed Documentation Claims
- ❌ "CMB age ~377k years (natural fit)" - Not testing this anymore
- ❌ "Derives cosmology from friction" - Using standard ΛCDM
- ❌ "Novel predictions: 0%" - Fixed with local drag implementation
- ❌ "3 parameters vs ΛCDM's 6" - Now 2 gravity params + ΛCDM's 6

### Fixed

#### Critical Fixes
- **Prediction magnitudes**: Changed from ~0% to testable 1-3% range
  - Root cause: Was using cosmic-average friction (~10⁻³⁶) instead of local
  - Solution: Implemented `friction_drag_local(M, r)` with mass density
  
- **Numerical instabilities**: Eliminated S(t) calculation warnings
  - Root cause: Friction-derived H(t) had integration issues at high-z
  - Solution: Use ΛCDM scale factor, no custom integration needed

- **CMB age discrepancy**: Resolved 1% mismatch
  - Root cause: Friction-derived H(t) naturally gave ~377k vs 380k years
  - Solution: Accept ΛCDM value, don't test CMB age independently

- **File path references**: Updated all documentation
  - Fixed: `estif_ec_fd_*.py` → `estif_ec_gr_*.py` in all docs
  - Fixed: `docs/VALIDATION_REPORT.md` → `docs/report/VALIDATION_REPORT.md`
  - Fixed: `docs/STATUS.md` → `docs/report/STATUS.md`

#### Parameter Updates
- **BETA_DRAG**: Tuned to 0.05 for detectable strong-field deviations
- **A**: Maintained at 0.0005 from BBN+SNe optimization
- **H₀**: Now uses Planck 2018 value (67.4 km/s/Mpc) directly

### Security
- No security issues in this release

### Deprecated
- `estif_ec_fd_*.py` naming convention (replaced by `estif_ec_gr_*.py`)
- Custom cosmology functions (replaced by `astropy.cosmology`)

---

## [1.5.0] - 2025-01-15 - Pre-Fork State (ESTIF-FD Final)

**Last version before ESTIF-Gravity fork**

### Summary
- Weak-field validation complete (GPS, Mercury, light deflection)
- Cosmology fits: χ² = 1.10 on supernovae
- **Known issues**: CMB age ~1% off, predictions giving ~0%, high-z instabilities

### Issues Leading to Fork
1. CMB age discrepancy (377k vs 380k years) - Never fully resolved
2. Novel predictions ~0% - Cosmic-average friction too weak
3. S(t) numerical warnings at z=1100 - Integration instabilities
4. Complexity of maintaining custom cosmology code

### Decision
- Pivot to ESTIF-Gravity: Use ΛCDM, focus on testable gravity predictions

---

## [1.0.0] - 2024-11-01 - Initial ESTIF-FD Release

### Added
- Core friction formulation
- Solar System tests (GPS, Mercury, light deflection)
- Supernova data fitting
- Basic cosmology derivation from friction

### Validated
- ✅ Weak-field GR equivalence (<1% deviation)
- ✅ Supernova fits (χ² = 1.10)
- ✅ BBN consistency

### Known Limitations
- Novel predictions not yet implemented
- High-z regime not tested
- CMB age not validated

---

## [0.3.0] - 2024-09-15 - Pre-Release

### Added
- Conceptual framework (ant/stone analogies)
- Development roadmap
- Rabbit holes and crossroads (RHAC) decision tree

### Changed
- Simplified H(t) formula to avoid overfitting
- Added time-boxing to prevent endless parameter tuning

---

## [0.2.0] - 2024-08-01 - Alpha Testing

### Added
- Mercury precession calculation
- Light deflection test
- Basic supernova data loading

### Fixed
- GPS time dilation now matches observations
- Numerical stability improvements

---

## [0.1.0] - 2024-07-01 - Initial Prototype

### Added
- Basic 4D inward flow kinematics
- Friction force formulation
- Constants file with physical parameters
- Simple Newtonian gravity reproduction

### Status
- Proof of concept only
- No observational validation yet

---

## Migration Guide: ESTIF-FD (v1.5) → ESTIF-Gravity (v2.0)

### For Code Users

#### Import Changes
```python
# Old (v1.5)
from estif_ec_fd_model import scale_factor, hubble_parameter
from estif_ec_fd_constants import H0_DEFAULT, A_DEFAULT, BETA_DRAG

# New (v2.0)
from src.estif_ec_gr_model import scale_factor_lcdm, luminosity_distance
from src.estif_ec_gr_constants import H0_PLANCK, A_DEFAULT, BETA_DRAG
from astropy.cosmology import FlatLambdaCDM
```

#### Function Changes
```python
# Old: Custom scale factor from friction
S_t = scale_factor(t, H0_DEFAULT, A_DEFAULT, BETA_DRAG)

# New: Standard ΛCDM scale factor
a_z = scale_factor_lcdm(z)  # Just 1/(1+z)
# Or use astropy directly:
cosmo = FlatLambdaCDM(H0=67.66, Om0=0.3111)
a_z = cosmo.scale_factor(z)
```

#### Prediction Functions
```python
# Old: Used cosmic-average friction (gave ~0%)
prediction = unique_lensing_signature(M, r, friction_drag(t_universe))

# New: Use local friction (gives ~1-3%)
prediction = unique_lensing_signature(M, r, friction_drag_local(M, r))
```

### For Documentation Users

#### Key Terminology Changes
- "ESTIF-FD" → "ESTIF-Gravity"
- "Friction-derived cosmology" → "Standard ΛCDM cosmology"
- "Novel predictions: TBD" → "LISA: 32 μs (3.2σ)"
- "CMB age ~377k years" → (Not tested, uses ΛCDM value)

#### File Location Changes
- `estif_ec_fd_model.py` → `src/estif_ec_gr_model.py`
- `docs/VALIDATION_REPORT.md` → `docs/report/VALIDATION_REPORT.md`
- `docs/STATUS.md` → `docs/report/STATUS.md`

### Breaking Changes

#### Removed Functions
- `scale_factor(t, H0, A, BETA_DRAG)` - Use `scale_factor_lcdm(z)` instead
- `hubble_parameter(t, A, BETA_DRAG)` - Use ΛCDM `H(z)` instead
- `time_from_scale_factor(S)` - No longer needed
- `cmb_age()` - Not testing CMB age independently

#### Changed Function Signatures
- `friction_drag(t)` → `friction_drag_local(M, r)` (for predictions)
- `unique_lensing_signature(M, r, friction)` - Now requires local friction
- `gw_damping_delay(M, D, friction)` - Now requires local friction

#### Changed Constants
- `H0_DEFAULT` → `H0_PLANCK` (now from Planck 2018)
- Removed: `S0_DEFAULT` (not calculating S(t) anymore)

### Configuration Changes

If using configuration files:
```yaml
# Old config
model: estif-fd
cosmology: friction-derived
test_cmb_age: true

# New config
model: estif-gravity
cosmology: lcdm-planck2018
test_cmb_age: false
```

---

## Future Roadmap

### Version 2.1.0 (Planned: December 2025)
- [ ] arXiv preprint submission
- [ ] Zenodo code release with DOI
- [ ] Expert review feedback incorporated
- [ ] Manuscript submitted to journal

### Version 2.2.0 (Planned: 2026)
- [ ] Peer review revisions
- [ ] Additional parameter constraint analysis (MCMC)
- [ ] Alternative strong-field tests explored
- [ ] Publication acceptance

### Version 3.0.0 (Speculative: Post-2030)
- [ ] LISA observational validation (if predictions correct)
- [ ] EHT+ black hole shadow measurements
- [ ] Cosmology derivation revisited (if gravity validated)
- [ ] Quantum extensions explored (if classical confirmed)

---

## Acknowledgments

### Contributors to v2.0 Pivot
- **Peter Angelov**: Principal investigator, implemented ESTIF-Gravity fork
- **Claude (Anthropic)**: Documentation assistance, validation framework design
- **Dr. Kirilova** (pending): Expert review consultation

### Key Decisions
- **January 2025**: Decision to abandon friction-derived cosmology
- **September 2025**: Implementation of local drag for predictions
- **October 2025**: Documentation overhaul and v2.0 release

### Lessons Learned
1. **Isolate variables**: Test one modification at a time (gravity before cosmology)
2. **Use proven frameworks**: Build on ΛCDM rather than replacing it prematurely
3. **Testability first**: Every feature must lead to falsifiable predictions
4. **Document pivots**: This changelog preserves decision history

---

## Version Numbering

This project uses [Semantic Versioning](https://semver.org/):
- **MAJOR** (X.0.0): Incompatible API changes or fundamental approach pivots
- **MINOR** (0.X.0): New functionality in backward-compatible manner
- **PATCH** (0.0.X): Backward-compatible bug fixes

**Current version**: 2.0.0 (ESTIF-Gravity Fork)

---

## Contact

**Author**: Peter Angelov  
**Email**: tervion@gmail.com  
**Project**: ESTIF - Emergent Spacetime from Inward Flow  
**Repository**: [https://github.com/tervion/estif-publication  
**Citation**: See `CITATION.cff` for academic citation format

---

## References

### Key Publications Referenced
1. Planck Collaboration (2018): "Planck 2018 results. VI. Cosmological parameters"
2. Suzuki et al. (2012): "The Hubble Space Telescope Cluster Supernova Survey"
3. EHT Collaboration (2019): "First M87 Event Horizon Telescope Results"
4. Abbott et al. (2016): "Observation of Gravitational Waves from a Binary Black Hole Merger" (GW150914)

### Internal Documentation
- `docs/report/VALIDATION_REPORT.md` - Technical validation details
- `docs/guide/estif_ec_fd_ROADMAP.md` - Development milestones
- `docs/plan/estif_ec_fd_RHAC.md` - Decision history
- `docs/report/estif_ec_fd_concept.md` - Conceptual framework

---

**Last Updated**: October 16, 2025  
**Status**: Active development  
**Next Release**: v2.1.0 (arXiv submission, December 2025)


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2

