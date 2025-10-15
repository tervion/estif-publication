# Data Directory

This directory contains observational data used for ESTIF model validation and fitting.

## Files

### sn_data.txt

**Type Ia Supernova Dataset**

- **Source**: Union2.1 compilation (Suzuki et al. 2012) and Pantheon dataset
- **Format**: Space-separated values
- **Columns**:
  1. Supernova name (string)
  2. Redshift (z) - dimensionless
  3. Distance modulus (μ) - magnitudes
  4. Distance modulus uncertainty (σ_μ) - magnitudes
  5. Third variance (systematic uncertainty) - magnitudes

- **Usage**: Used for cosmological distance-redshift relation validation
- **Sample size**: ~580 supernovae
- **Redshift range**: z ≈ 0.01 to 1.5

**Comments** (lines starting with #):
- Contains calibration parameters (alpha, beta, delta, M)
- Includes notes on Hubble constant assumptions (h=0.7)

## Data Provenance

The supernova data in this repository is derived from publicly available astronomical surveys:

1. **Union2.1 Compilation**: Suzuki et al. (2012), ApJ, 746, 85
   - DOI: 10.1088/0004-637X/746/1/85
   - URL: https://supernova.lbl.gov/Union/

2. **Pantheon Sample**: Scolnic et al. (2018), ApJ, 859, 101
   - DOI: 10.3847/1538-4357/aab9bb
   - URL: https://pantheonplussh0es.github.io/

## Usage in ESTIF

The supernova data is used to:
1. Validate distance modulus predictions: μ(z) = 5 log₁₀(d_L) + 25
2. Fit friction parameters (BETA_DRAG, A_DEFAULT) via χ² minimization
3. Compare ESTIF predictions to ΛCDM standard model
4. Test the model's cosmological distance-redshift relation

**Reference scripts**:
- `estif_ec_fd_run_simulation.py::test_and_plot_supernova_data()`
- `tests/test_cosmology.py`
- `tests/model_comparison.py`

## Citation

If you use this data in publications, please cite both:
1. The original supernova survey papers (listed above)
2. This ESTIF implementation (see main repository CITATION.cff)

## Notes

- Data quality flags are included in the fourth column
- Systematic uncertainties are accounted for in model fitting
- The dataset is filtered for z ≤ 1.5 in low-z validation tests
- Some entries may have been excluded based on quality cuts

---

*Last updated: 2025-01-15*
*Data format validated: ESTIF v0.3.0*


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2



