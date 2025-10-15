# test_fixes_improved.py
"""
Improved diagnostic tool for ESTIF distance modulus calculations.
Addresses issues: multiple z values, error handling, consistency checks,
comparison to ΛCDM, validation of S(t) monotonicity.
"""

import numpy as np
import estif_ec_fd_model as estif
import estif_ec_fd_constants as const
import matplotlib.pyplot as plt
from scipy.integrate import quad

def distance_modulus_lcdm_simple(z, H0_kmsMpc=67.66, Omega_m=0.3111):
    """Quick ΛCDM distance modulus for comparison"""
    if z <= 0:
        return 0.0
    
    c_km_s = 299792.458
    Omega_Lambda = 1.0 - Omega_m
    
    def integrand(z_prime):
        E_z = np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_Lambda)
        return 1.0 / E_z
    
    try:
        integral, _ = quad(integrand, 0, z, limit=100)
        d_c = (c_km_s / H0_kmsMpc) * integral
        d_L = (1 + z) * d_c
        mu = 5 * np.log10(d_L) + 25
        return mu
    except:
        return np.nan

def test_distance_modulus_comprehensive():
    """
    Comprehensive test of distance modulus calculation across redshift range.
    Tests for consistency, monotonicity, and comparison to ΛCDM.
    """
    
    print("\n" + "="*80)
    print("COMPREHENSIVE DISTANCE MODULUS TEST")
    print("="*80)
    
    # Test parameters
    t_obs = 4.35e17  # Current universe age in seconds
    z_test_values = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]
    
    print(f"\nUniverse age: t_obs = {t_obs:.3e} s ({t_obs/(365.25*86400*1e9):.2f} Gyr)")
    print(f"Testing {len(z_test_values)} redshift values...")
    
    # Storage for results
    results = {
        'z': [],
        't_emit': [],
        'mu_estif': [],
        'mu_lcdm': [],
        'H_emit': [],
        'H_obs': [],
        'S_emit': [],
        'S_obs': [],
        'errors': []
    }
    
    # Test each redshift
    for z in z_test_values:
        print(f"\n{'─'*80}")
        print(f"Testing z = {z}")
        print(f"{'─'*80}")
        
        error_occurred = False
        
        try:
            # 1. Calculate emission time
            print("  1. Computing emission time...")
            t_emit = estif.t_from_z(z, t_obs)
            
            # Sanity check
            if t_emit >= t_obs:
                print(f"    ❌ ERROR: t_emit ({t_emit:.3e}) >= t_obs ({t_obs:.3e})")
                error_occurred = True
            else:
                print(f"    ✓ t_emit = {t_emit:.3e} s ({t_emit/(365.25*86400*1e9):.2f} Gyr)")
            
            # 2. Calculate scale factors
            print("  2. Computing scale factors...")
            S_emit = estif.global_S(t_emit)
            S_obs = estif.global_S(t_obs)
            
            # Check S(t) is decreasing (we flow inward)
            if S_emit <= S_obs:
                print(f"    ❌ ERROR: S not decreasing! S_emit ({S_emit:.6f}) <= S_obs ({S_obs:.6f})")
                error_occurred = True
            else:
                print(f"    ✓ S_emit = {S_emit:.6f}, S_obs = {S_obs:.6f}")
                
            # Check redshift consistency
            z_check = S_emit / S_obs - 1
            if not np.isclose(z_check, z, rtol=0.1):
                print(f"    ⚠️  WARNING: Redshift mismatch! Expected {z}, got {z_check:.4f}")
            
            # 3. Calculate H(t) values
            print("  3. Computing Hubble parameters...")
            H_emit = estif.H_variable(t_emit)
            H_obs = estif.H_variable(t_obs)
            
            print(f"    H(t_emit) = {H_emit:.6e} s⁻¹")
            print(f"    H(t_obs) = {H_obs:.6e} s⁻¹")
            print(f"    H₀ (const) = {const.H_0:.6e} s⁻¹")
            
            # 4. Calculate distance modulus (ESTIF)
            print("  4. Computing distance modulus (ESTIF)...")
            mu_estif = estif.distance_modulus_estif_numerical(z)
            
            if not np.isfinite(mu_estif):
                print(f"    ❌ ERROR: Invalid mu_estif = {mu_estif}")
                error_occurred = True
            else:
                print(f"    ✓ μ_ESTIF = {mu_estif:.2f} mag")
            
            # 5. Calculate distance modulus (ΛCDM for comparison)
            print("  5. Computing distance modulus (ΛCDM)...")
            mu_lcdm = distance_modulus_lcdm_simple(z)
            
            if np.isfinite(mu_lcdm):
                print(f"    ✓ μ_ΛCDM = {mu_lcdm:.2f} mag")
                diff = mu_estif - mu_lcdm
                print(f"    Δμ (ESTIF - ΛCDM) = {diff:.2f} mag ({100*diff/mu_lcdm:.1f}%)")
            
            # Store results
            results['z'].append(z)
            results['t_emit'].append(t_emit)
            results['mu_estif'].append(mu_estif)
            results['mu_lcdm'].append(mu_lcdm)
            results['H_emit'].append(H_emit)
            results['H_obs'].append(H_obs)
            results['S_emit'].append(S_emit)
            results['S_obs'].append(S_obs)
            results['errors'].append(error_occurred)
            
        except Exception as e:
            print(f"  ❌ EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            results['errors'].append(True)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    n_success = len([e for e in results['errors'] if not e])
    n_total = len(z_test_values)
    print(f"\nTests passed: {n_success}/{n_total}")
    
    if n_success > 0:
        print("\nDistance modulus comparison:")
        print(f"{'z':<8} {'μ_ESTIF':<12} {'μ_ΛCDM':<12} {'Δμ':<12} {'Δμ%':<10}")
        print("─"*60)
        for i in range(len(results['z'])):
            if not results['errors'][i]:
                z = results['z'][i]
                mu_e = results['mu_estif'][i]
                mu_l = results['mu_lcdm'][i]
                diff = mu_e - mu_l
                pct = 100 * diff / mu_l
                print(f"{z:<8.2f} {mu_e:<12.2f} {mu_l:<12.2f} {diff:<12.2f} {pct:<10.2f}%")
    
    # Check for systematic issues
    print("\nDiagnostic checks:")
    
    # Check S(t) monotonicity
    if len(results['S_emit']) > 1:
        S_decreasing = all(results['S_emit'][i] > results['S_emit'][i+1] 
                          for i in range(len(results['S_emit'])-1))
        if S_decreasing:
            print("  ✓ S(t) decreases monotonically with increasing z")
        else:
            print("  ❌ S(t) NOT monotonic - CRITICAL ISSUE")
    
    # Check mu monotonicity
    if len(results['mu_estif']) > 1:
        mu_increasing = all(results['mu_estif'][i] < results['mu_estif'][i+1] 
                           for i in range(len(results['mu_estif'])-1) 
                           if not results['errors'][i] and not results['errors'][i+1])
        if mu_increasing:
            print("  ✓ Distance modulus increases with z")
        else:
            print("  ⚠️  Distance modulus NOT monotonic")
    
    # Plot results if successful
    if n_success >= 3:
        plot_results(results)
    
    print("\n" + "="*80 + "\n")
    
    return results

def plot_results(results):
    """Create diagnostic plots"""
    
    # Filter out errors
    valid = [not e for e in results['errors']]
    z_valid = np.array([results['z'][i] for i, v in enumerate(valid) if v])
    mu_estif = np.array([results['mu_estif'][i] for i, v in enumerate(valid) if v])
    mu_lcdm = np.array([results['mu_lcdm'][i] for i, v in enumerate(valid) if v])
    
    if len(z_valid) < 2:
        return
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot 1: Distance modulus comparison
    axes[0].plot(z_valid, mu_estif, 'o-', label='ESTIF', linewidth=2)
    axes[0].plot(z_valid, mu_lcdm, 's--', label='ΛCDM', linewidth=2)
    axes[0].set_xlabel('Redshift z')
    axes[0].set_ylabel('Distance Modulus μ (mag)')
    axes[0].set_title('Distance Modulus: ESTIF vs ΛCDM')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    residuals = mu_estif - mu_lcdm
    axes[1].plot(z_valid, residuals, 'o-', linewidth=2, color='red')
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[1].set_xlabel('Redshift z')
    axes[1].set_ylabel('Δμ (ESTIF - ΛCDM) [mag]')
    axes[1].set_title('Residuals: ESTIF - ΛCDM')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    filename = 'distance_modulus_diagnostics.png'
    plt.savefig(filename, dpi=150)
    print(f"\nDiagnostic plots saved to '{filename}'")
    plt.close()

def test_s_monotonicity():
    """
    Explicit test of S(t) monotonicity to catch the "invalid target" issue
    """
    print("\n" + "="*80)
    print("S(t) MONOTONICITY TEST")
    print("="*80)
    
    # Sample times from Big Bang to now
    t_values = np.logspace(10, np.log10(4.35e17), 50)
    S_values = []
    
    print("\nComputing S(t) for 50 time values...")
    
    for t in t_values:
        try:
            S = estif.global_S(t)
            S_values.append(S)
        except Exception as e:
            print(f"ERROR at t={t:.3e}: {e}")
            S_values.append(np.nan)
    
    S_values = np.array(S_values)
    
    # Check monotonicity
    valid = np.isfinite(S_values)
    if np.sum(valid) < 2:
        print("❌ CRITICAL: Too few valid S(t) values")
        return
    
    dS = np.diff(S_values[valid])
    
    if np.all(dS < 0):
        print("✓ S(t) decreases monotonically (correct for inward flow)")
    elif np.all(dS > 0):
        print("❌ CRITICAL: S(t) INCREASES with time (wrong sign!)")
    else:
        print("❌ CRITICAL: S(t) is NOT monotonic")
        print(f"   Positive steps: {np.sum(dS > 0)}")
        print(f"   Negative steps: {np.sum(dS < 0)}")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(t_values[valid]/1e17, S_values[valid], 'o-')
    plt.xlabel('Time (×10¹⁷ s)')
    plt.ylabel('S(t)')
    plt.title('Scale Factor Evolution')
    plt.grid(True, alpha=0.3)
    plt.savefig('s_monotonicity_test.png', dpi=150)
    print("\nPlot saved to 's_monotonicity_test.png'")
    plt.close()
    
    print("="*80 + "\n")

if __name__ == "__main__":
    # Run comprehensive distance modulus test
    results = test_distance_modulus_comprehensive()
    
    # Run explicit S(t) monotonicity test
    test_s_monotonicity()
    
    print("All tests complete!")
