# test_fixes_low_z.py

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

def test_low_z_range():
    """Test ESTIF on valid z range (z≤1.5)"""
    print("\n" + "="*80)
    print("LOW-Z VALIDATION: ESTIF for z ≤ 1.5")
    print("="*80)
    
    t_obs = 4.35e17
    z_test_values = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 1.3, 1.5]
    
    results = {'z': [], 't_emit': [], 'mu_estif': [], 'mu_lcdm': [], 'errors': []}
    
    for z in z_test_values:
        print(f"\nTesting z = {z}")
        error_occurred = False
        
        try:
            t_emit = estif.t_from_z(z, t_obs)
            S_emit = estif.global_S(t_emit)
            S_obs = estif.global_S(t_obs)
            z_check = S_emit / S_obs - 1
            
            mu_estif = estif.distance_modulus_estif_numerical(z)
            mu_lcdm = distance_modulus_lcdm_simple(z)
            
            print(f"  t_emit = {t_emit:.3e} s ({t_emit/(365.25*86400*1e9):.2f} Gyr)")
            print(f"  S_emit = {S_emit:.6f}, S_obs = {S_obs:.6f}")
            print(f"  z_check = {z_check:.4f} (expected {z})")
            print(f"  μ_ESTIF = {mu_estif:.2f}, μ_ΛCDM = {mu_lcdm:.2f}")
            print(f"  Δμ = {mu_estif - mu_lcdm:.2f} mag")
            
            if abs(z_check - z) > 0.1:
                print(f"  ⚠️  WARNING: Redshift mismatch!")
                error_occurred = True
            
            results['z'].append(z)
            results['t_emit'].append(t_emit)
            results['mu_estif'].append(mu_estif)
            results['mu_lcdm'].append(mu_lcdm)
            results['errors'].append(error_occurred)
            
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            results['errors'].append(True)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    n_success = len([e for e in results['errors'] if not e])
    print(f"\nTests passed: {n_success}/{len(z_test_values)}")
    
    if n_success > 0:
        print("\nDistance modulus comparison (valid range):")
        print(f"{'z':<8} {'μ_ESTIF':<12} {'μ_ΛCDM':<12} {'Δμ':<12}")
        print("─"*50)
        for i in range(len(results['z'])):
            if not results['errors'][i]:
                z = results['z'][i]
                mu_e = results['mu_estif'][i]
                mu_l = results['mu_lcdm'][i]
                diff = mu_e - mu_l
                print(f"{z:<8.2f} {mu_e:<12.2f} {mu_l:<12.2f} {diff:<12.2f}")
    
    # Plot
    if n_success >= 3:
        valid = [not e for e in results['errors']]
        z_valid = np.array([results['z'][i] for i, v in enumerate(valid) if v])
        mu_estif = np.array([results['mu_estif'][i] for i, v in enumerate(valid) if v])
        mu_lcdm = np.array([results['mu_lcdm'][i] for i, v in enumerate(valid) if v])
        
        fig, axes = plt.subplots(2, 1, figsize=(10, 10))
        
        axes[0].plot(z_valid, mu_estif, 'o-', label='ESTIF (z≤1.5)', linewidth=2)
        axes[0].plot(z_valid, mu_lcdm, 's--', label='ΛCDM', linewidth=2)
        axes[0].set_xlabel('Redshift z')
        axes[0].set_ylabel('Distance Modulus μ (mag)')
        axes[0].set_title('ESTIF Low-z Model Performance')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        residuals = mu_estif - mu_lcdm
        axes[1].plot(z_valid, residuals, 'o-', linewidth=2, color='red')
        axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
        axes[1].set_xlabel('Redshift z')
        axes[1].set_ylabel('Δμ (ESTIF - ΛCDM) [mag]')
        axes[1].set_title('Residuals (ESTIF - ΛCDM)')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('estif_low_z_validation.png', dpi=150)
        print(f"\nPlot saved to 'estif_low_z_validation.png'")
        plt.close()
    
    print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    test_low_z_range()

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-1

