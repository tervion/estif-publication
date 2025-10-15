# model_comparison_fixed.py
"""
Fixed version: Compare ESTIF vs ΛCDM with proper calculations for both models.
Addresses issues: computes ΛCDM directly, adds validation, better error handling.
"""

import numpy as np
import estif_ec_fd_model as estif
from scipy.integrate import quad

# Physical constants
c_km_s = 299792.458  # Speed of light in km/s
H0_planck = 67.66    # Planck 2018 H0 in km/s/Mpc
Omega_m_planck = 0.3111  # Planck 2018 matter density

def calculate_aic_bic(chi_squared, n_data, n_params):
    """Calculate AIC and BIC for model comparison"""
    # Log-likelihood (for Gaussian errors)
    log_likelihood = -0.5 * chi_squared
    
    # Akaike Information Criterion
    aic = 2 * n_params - 2 * log_likelihood
    
    # Bayesian Information Criterion
    bic = n_params * np.log(n_data) - 2 * log_likelihood
    
    # Reduced chi-squared
    dof = n_data - n_params
    chi_squared_reduced = chi_squared / dof if dof > 0 else np.inf
    
    return {
        'chi2': chi_squared,
        'chi2_red': chi_squared_reduced,
        'aic': aic,
        'bic': bic,
        'n_params': n_params,
        'n_data': n_data,
        'dof': dof
    }

def distance_modulus_lcdm(z, H0, Omega_m):
    """
    Calculate ΛCDM distance modulus for given redshift.
    Uses flat ΛCDM: Omega_Lambda = 1 - Omega_m
    
    Args:
        z: Redshift
        H0: Hubble constant in km/s/Mpc
        Omega_m: Matter density parameter
    
    Returns:
        Distance modulus in magnitudes
    """
    if z <= 0:
        return 0.0
    
    Omega_Lambda = 1.0 - Omega_m
    
    # Integrand for comoving distance: c/H(z)
    def integrand(z_prime):
        E_z = np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_Lambda)
        return 1.0 / E_z
    
    try:
        # Integrate from 0 to z
        integral, _ = quad(integrand, 0, z, limit=100)
        
        # Comoving distance in Mpc
        d_c = (c_km_s / H0) * integral
        
        # Luminosity distance
        d_L = (1 + z) * d_c
        
        # Distance modulus
        mu = 5 * np.log10(d_L) + 25
        
        return mu
    except Exception as e:
        print(f"Warning: ΛCDM calculation failed for z={z}: {e}")
        return np.nan

def compare_models_improved():
    """
    Compare ESTIF and ΛCDM on supernova data with proper calculations.
    Both models computed on identical dataset.
    """
    
    print("\n" + "="*80)
    print("IMPROVED MODEL COMPARISON: ESTIF vs ΛCDM")
    print("="*80)
    
    # Load supernova data
    try:
        data = np.loadtxt('sn_data.txt', comments='#', usecols=(1, 2, 3))
        z_data = data[:, 0]
        mu_data = data[:, 1]
        sigma_mu = data[:, 2]
        
        # Filter valid data
        valid = np.isfinite(z_data) & np.isfinite(mu_data) & (z_data > 0) & (sigma_mu > 0)
        z_data = z_data[valid]
        mu_data = mu_data[valid]
        sigma_mu = sigma_mu[valid]
        
        n_data = len(z_data)
        print(f"\nDataset: {n_data} Type Ia Supernovae")
        print(f"Redshift range: z = {z_data.min():.4f} to {z_data.max():.4f}")
        
    except FileNotFoundError:
        print("\nERROR: 'sn_data.txt' not found!")
        print("Cannot perform comparison without data.")
        return
    except Exception as e:
        print(f"\nERROR loading data: {e}")
        return
    
    # Calculate ESTIF predictions
    print("\n--- Computing ESTIF predictions ---")
    try:
        mu_estif = estif.distance_modulus_estif_numerical(z_data)
        
        # Check for any failures
        invalid_estif = ~np.isfinite(mu_estif)
        if np.any(invalid_estif):
            print(f"WARNING: {np.sum(invalid_estif)} ESTIF predictions failed")
            # Remove failed points from both datasets
            valid_estif = np.isfinite(mu_estif)
            z_data = z_data[valid_estif]
            mu_data = mu_data[valid_estif]
            sigma_mu = sigma_mu[valid_estif]
            mu_estif = mu_estif[valid_estif]
            n_data = len(z_data)
        
        chi2_estif = np.sum(((mu_data - mu_estif) / sigma_mu)**2)
        print(f"ESTIF χ² = {chi2_estif:.2f}")
        
    except Exception as e:
        print(f"ERROR: ESTIF calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Calculate ΛCDM predictions
    print("\n--- Computing ΛCDM predictions ---")
    try:
        # Use Planck 2018 best-fit parameters
        mu_lcdm = np.array([distance_modulus_lcdm(z, H0_planck, Omega_m_planck) 
                            for z in z_data])
        
        # Check for any failures
        invalid_lcdm = ~np.isfinite(mu_lcdm)
        if np.any(invalid_lcdm):
            print(f"WARNING: {np.sum(invalid_lcdm)} ΛCDM predictions failed")
            valid_lcdm = np.isfinite(mu_lcdm)
            z_data = z_data[valid_lcdm]
            mu_data = mu_data[valid_lcdm]
            sigma_mu = sigma_mu[valid_lcdm]
            mu_estif = mu_estif[valid_lcdm]
            mu_lcdm = mu_lcdm[valid_lcdm]
            n_data = len(z_data)
        
        chi2_lcdm = np.sum(((mu_data - mu_lcdm) / sigma_mu)**2)
        print(f"ΛCDM χ² = {chi2_lcdm:.2f}")
        
    except Exception as e:
        print(f"ERROR: ΛCDM calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Calculate metrics
    # ESTIF: 3 parameters (H0, A, BETA_DRAG)
    estif_metrics = calculate_aic_bic(chi2_estif, n_data, 3)
    
    # ΛCDM: For distance-only fits, typically 3 params (H0, Omega_m, M_abs)
    # But using 6 for full cosmological model comparison per your docs
    lcdm_metrics_3param = calculate_aic_bic(chi2_lcdm, n_data, 3)
    lcdm_metrics_6param = calculate_aic_bic(chi2_lcdm, n_data, 6)
    
    # Print comparison tables
    print("\n" + "="*80)
    print("RESULTS: ESTIF (3 params) vs ΛCDM (3 params, distance-only)")
    print("="*80)
    print(f"\n{'Metric':<20} {'ESTIF':<20} {'ΛCDM':<20} {'Δ (ESTIF-ΛCDM)':<20}")
    print("-"*80)
    print(f"{'Parameters':<20} {estif_metrics['n_params']:<20} {lcdm_metrics_3param['n_params']:<20} {estif_metrics['n_params'] - lcdm_metrics_3param['n_params']:<20}")
    print(f"{'χ²':<20} {estif_metrics['chi2']:<20.2f} {lcdm_metrics_3param['chi2']:<20.2f} {estif_metrics['chi2'] - lcdm_metrics_3param['chi2']:<20.2f}")
    print(f"{'χ²_red':<20} {estif_metrics['chi2_red']:<20.3f} {lcdm_metrics_3param['chi2_red']:<20.3f} {estif_metrics['chi2_red'] - lcdm_metrics_3param['chi2_red']:<20.3f}")
    print(f"{'AIC':<20} {estif_metrics['aic']:<20.2f} {lcdm_metrics_3param['aic']:<20.2f} {estif_metrics['aic'] - lcdm_metrics_3param['aic']:<20.2f}")
    print(f"{'BIC':<20} {estif_metrics['bic']:<20.2f} {lcdm_metrics_3param['bic']:<20.2f} {estif_metrics['bic'] - lcdm_metrics_3param['bic']:<20.2f}")
    
    print("\n" + "="*80)
    print("RESULTS: ESTIF (3 params) vs ΛCDM (6 params, full cosmology)")
    print("="*80)
    print(f"\n{'Metric':<20} {'ESTIF':<20} {'ΛCDM':<20} {'Δ (ESTIF-ΛCDM)':<20}")
    print("-"*80)
    print(f"{'Parameters':<20} {estif_metrics['n_params']:<20} {lcdm_metrics_6param['n_params']:<20} {estif_metrics['n_params'] - lcdm_metrics_6param['n_params']:<20}")
    print(f"{'χ²':<20} {estif_metrics['chi2']:<20.2f} {lcdm_metrics_6param['chi2']:<20.2f} {estif_metrics['chi2'] - lcdm_metrics_6param['chi2']:<20.2f}")
    print(f"{'χ²_red':<20} {estif_metrics['chi2_red']:<20.3f} {lcdm_metrics_6param['chi2_red']:<20.3f} {estif_metrics['chi2_red'] - lcdm_metrics_6param['chi2_red']:<20.3f}")
    print(f"{'AIC':<20} {estif_metrics['aic']:<20.2f} {lcdm_metrics_6param['aic']:<20.2f} {estif_metrics['aic'] - lcdm_metrics_6param['aic']:<20.2f}")
    print(f"{'BIC':<20} {estif_metrics['bic']:<20.2f} {lcdm_metrics_6param['bic']:<20.2f} {estif_metrics['bic'] - lcdm_metrics_6param['bic']:<20.2f}")
    
    # Interpretation
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)
    
    delta_bic_3 = estif_metrics['bic'] - lcdm_metrics_3param['bic']
    delta_bic_6 = estif_metrics['bic'] - lcdm_metrics_6param['bic']
    
    print("\nGuidelines:")
    print("  • Negative Δ values favor ESTIF")
    print("  • |Δ| > 10 indicates strong evidence")
    print("  • |Δ| > 2 indicates moderate evidence")
    
    print(f"\nComparison to 3-parameter ΛCDM (distance-only fit):")
    if delta_bic_3 < -10:
        print(f"  → Strong evidence for ESTIF (ΔBIC = {delta_bic_3:.1f})")
    elif delta_bic_3 < -2:
        print(f"  → Moderate evidence for ESTIF (ΔBIC = {delta_bic_3:.1f})")
    elif delta_bic_3 > 10:
        print(f"  → Strong evidence for ΛCDM (ΔBIC = {delta_bic_3:.1f})")
    elif delta_bic_3 > 2:
        print(f"  → Moderate evidence for ΛCDM (ΔBIC = {delta_bic_3:.1f})")
    else:
        print(f"  → Models roughly equivalent (ΔBIC = {delta_bic_3:.1f})")
    
    print(f"\nComparison to 6-parameter ΛCDM (full cosmology):")
    if delta_bic_6 < -10:
        print(f"  → Strong evidence for ESTIF (ΔBIC = {delta_bic_6:.1f})")
    elif delta_bic_6 < -2:
        print(f"  → Moderate evidence for ESTIF (ΔBIC = {delta_bic_6:.1f})")
    elif delta_bic_6 > 10:
        print(f"  → Strong evidence for ΛCDM (ΔBIC = {delta_bic_6:.1f})")
    elif delta_bic_6 > 2:
        print(f"  → Moderate evidence for ΛCDM (ΔBIC = {delta_bic_6:.1f})")
    else:
        print(f"  → Models roughly equivalent (ΔBIC = {delta_bic_6:.1f})")
    
    print("\n" + "="*80)
    print("NOTES")
    print("="*80)
    print("• Both models computed on identical dataset")
    print("• ΛCDM uses Planck 2018 parameters (H₀=67.66 km/s/Mpc, Ω_m=0.3111)")
    print("• 3-param comparison: fair for distance-only fits")
    print("• 6-param comparison: accounts for full ΛCDM parameter space")
    print("• For publication, consider marginalizing over parameter uncertainties")
    print("="*80 + "\n")

if __name__ == "__main__":
    compare_models_improved()
