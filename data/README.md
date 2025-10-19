# Data Directory

This directory contains observational data used for model validation.

## Contents

### sn_data.txt
Type Ia Supernovae data for cosmological fitting and validation.

**Format:** Plain text file with columns for:
- Redshift (z)
- Distance modulus (μ)
- Uncertainty (σ)

**Source:** [Specify the source dataset, e.g., Pantheon, Union2.1, etc.]

**Usage:** 
This data is used in `src/estif_ec_gr_run_simulation.py` to validate the ESTIF model against standard cosmological observations.

## Adding New Data

When adding new observational datasets:
1. Document the source and format in this README
2. Include appropriate citation information
3. Update validation scripts in the `tests/` directory
4. Add entries to CITATION.cff if necessary


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2



