basic_functionality_test_estif_ec_gr_model.py

```

**Expected Output:**
```
======================================================================
ESTIF-GRAVITY MODEL - Basic Tests
======================================================================

1. Testing ΛCDM cosmology functions:
   z=0.5: d_L = 2726.7 Mpc, μ = 42.18 mag
   Expected: d_L ≈ 2700 Mpc, μ ≈ 42.2 mag
   ✓ ΛCDM functions working correctly

2. Testing friction drag scaling:
   ✓ Friction scaling test passed: 2.0M/1.0M = 2.000x drag

3. Comparing local vs cosmic drag:
   Cosmic drag (t=now): 1.025e-37
   Local drag (Sun surface): 2.456e+08
   Ratio (local/cosmic): 2.396e+45
   ✓ Local drag >> cosmic (as expected for predictions)

4. Testing lensing prediction:
   At r = 5 R_s from Sun:
   GR deflection: 8.494e-06 rad
   ESTIF deflection: 8.494e-06 rad
   Deviation: 0.0001%
   ⚠️  Note: Magnitude pending validation

5. Testing weak-field GR limit:
   At r = 1e+12 m:
   Passes <1% test: True

======================================================================
Basic tests complete!
======================================================================

📋 Summary:
   ✓ ΛCDM cosmology: Working
   ✓ Friction functions: Implemented
   ✓ Weak-field limit: GR compliance
   ⚠️  Strong-field predictions: Pending validation

Next steps:
   1. Run estif_ec_gr_run_simulation.py for full test suite
   2. Compare predictions to EHT/LIGO/JWST data
   3. Refine prediction magnitudes based on observations

   

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2


