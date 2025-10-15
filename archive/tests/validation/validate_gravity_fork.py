# validate_gravity_fork.py

"""
Quick validation that ESTIF-Gravity fork is working correctly.
"""

import estif_ec_gr_model as estif
import estif_ec_gr_constants as const
import numpy as np

def main():
    print("\n" + "="*70)
    print("ESTIF-GRAVITY FORK VALIDATION")
    print("="*70)
    
    checks_passed = 0
    checks_total = 5
    
    # Check 1: ΛCDM cosmology
    print("\n✓ Check 1: ΛCDM cosmology functions")
    try:
        z = 1.0
        d_L = estif.luminosity_distance(z)
        assert 4000 < d_L < 8000, f"d_L={d_L} Mpc out of range"
        print(f"   z=1.0 → d_L = {d_L:.0f} Mpc ✓")
        checks_passed += 1
    except Exception as e:
        print(f"   ❌ Failed: {e}")
    
    # Check 2: Local drag implemented
    print("\n✓ Check 2: Local drag function")
    try:
        drag = estif.friction_drag_local(const.M_sun, const.R_sun)
        assert drag > 0, "Drag should be positive"
        print(f"   Local drag at Sun surface: {drag:.3e} ✓")
        checks_passed += 1
    except Exception as e:
        print(f"   ❌ Failed: {e}")
    
    # Check 3: Deprecated functions commented out
    print("\n✓ Check 3: Deprecated functions not accessible")
    try:
        # These should raise NameError or be in commented block
        has_H_variable = hasattr(estif, 'H_variable')
        has_global_S = hasattr(estif, 'global_S')
        
        if not has_H_variable and not has_global_S:
            print("   Exponential cosmology functions properly deprecated ✓")
            checks_passed += 1
        else:
            print("   ⚠️  Warning: Some deprecated functions still accessible")
    except Exception as e:
        print(f"   ❌ Failed: {e}")
    
    # Check 4: Distance modulus uses ΛCDM
    print("\n✓ Check 4: Distance modulus using ΛCDM")
    try:
        mu = estif.distance_modulus_lcdm(0.5)
        assert 41 < mu < 43, f"μ={mu} out of range"
        print(f"   z=0.5 → μ = {mu:.2f} mag ✓")
        checks_passed += 1
    except Exception as e:
        print(f"   ❌ Failed: {e}")
    
    # Check 5: Lensing prediction exists
    print("\n✓ Check 5: Lensing prediction function")
    try:
        r = 5 * estif.schwarzschild_radius(const.M_sun)
        theta = estif.unique_lensing_signature(r, const.M_sun)
        assert theta > 0, "Deflection should be positive"
        print(f"   Lensing at 5 R_s: {theta:.3e} rad ✓")
        checks_passed += 1
    except Exception as e:
        print(f"   ❌ Failed: {e}")
    
    # Summary
    print("\n" + "="*70)
    print(f"VALIDATION RESULT: {checks_passed}/{checks_total} checks passed")
    print("="*70)
    
    if checks_passed == checks_total:
        print("\n✅ ESTIF-Gravity fork successfully implemented!")
        print("\nNext steps:")
        print("  1. Run full test suite: python estif_ec_gr_run_simulation.py")
        print("  2. Compare predictions to observational data")
        print("  3. Refine prediction magnitudes")
    else:
        print("\n⚠️  Some checks failed. Review output above.")
    
    print()

if __name__ == "__main__":
    main()

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-1

