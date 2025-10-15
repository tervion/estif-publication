# model_comparison.py
"""
Compare ESTIF vs ΛCDM using information criteria (AIC/BIC) on supernova data.
Improved: Computes ΛCDM dynamically with astropy; optional param optimization;
plotting; consistent n_params for SNe fits; error handling.
"""

import numpy as np
import estif_ec_fd_model as estif
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def calculate_aic_bic(chi_squared, n_data, n_params):
    """Calculate AIC and BIC for model comparison"""
    log_likelihood = -0.5 * chi_squared
    aic = 2 * n_params - 2 * log_likelihood
    bic = n_params * np.log(n_data) - 2 * log_likelihood
    dof = n_data - n_params
    chi_squared_reduced = chi_squared / dof if dof > 0 else np.inf
    return {
        'chi2': chi_squared,
        'chi2_red': chi_squared_reduced,
        'aic': aic,
        'bic': bic,
        'n_params': n_params,
        'n_data': n_data
    }

def compute_lcdm_chi2(z_data, mu_data, sigma_mu, optimize=False):
    """Compute χ² for flat ΛCDM, with optional optimization over Om0 and M offset"""
    def chi2_func(params):
        Om0, M_offset = params
        lcdm = FlatLambdaCDM(H0=70, Om0=Om0)  # H0 arbitrary, absorbed in M
        mu_theory = 5 * np.log10(lcdm.luminosity_distance(z_data).to_value('Mpc')) + 25
        residuals = mu_data - mu_theory - M_offset
        return np.sum((residuals / sigma_mu)**2)

    if optimize:
        # Optimize Om0 (0.1-0.5) and M_offset (-1 to 1)
        res = minimize(chi2_func, [0.28, 0.0], bounds=[(0.1, 0.5), (-1, 1)])
        if res.success:
            return res.fun, res.x[0]  # min χ², best Om0
        else:
            raise ValueError("ΛCDM optimization failed")
    else:
        # Fixed params (from Suzuki et al. 2012 approx)
        Om0 = 0.28
        lcdm = FlatLambdaCDM(H0=70, Om0=Om0)
        mu_theory = 5 * np.log10(lcdm.luminosity_distance(z_data).to_value('Mpc')) + 25
        # Marginalize over M_offset analytically (mean shift)
        residuals = mu_data - mu_theory
        M_offset = np.average(residuals, weights=1/sigma_mu**2)
        adjusted_residuals = residuals - M_offset
        return np.sum((adjusted_residuals / sigma_mu)**2), Om0

def plot_residuals(z_data, mu_data, sigma_mu, mu_estif, mu_lcdm):
    """Plot residuals for both models"""
    plt.figure(figsize=(10, 6))
    plt.errorbar(z_data, mu_data - mu_estif, yerr=sigma_mu, fmt='b.', label='ESTIF residuals')
    plt.errorbar(z_data, mu_data - mu_lcdm, yerr=sigma_mu, fmt='r.', label='ΛCDM residuals')
    plt.axhline(0, color='k', linestyle='--')
    plt.xlabel('Redshift z')
    plt.ylabel('Δμ (data - model)')
    plt.legend()
    plt.title('Supernova Distance Modulus Residuals')
    plt.savefig('model_residuals.png')
    plt.close()
    print("Residual plot saved as 'model_residuals.png'")

def compare_models(optimize_lcdm=False, plot=False):
    """Compare ESTIF and ΛCDM on supernova data"""
    try:
        data = np.loadtxt('sn_data.txt', comments='#', usecols=(1, 2, 3))
        z_data = data[:, 0]
        mu_data = data[:, 1]
        sigma_mu = data[:, 2]
        valid = np.isfinite(z_data) & np.isfinite(mu_data) & (z_data > 0) & (sigma_mu > 0)
        z_data = z_data[valid]
        mu_data = mu_data[valid]
        sigma_mu = sigma_mu[valid]
        n_data = len(z_data)
        if n_data == 0:
            raise ValueError("No valid data points")
    except FileNotFoundError:
        print("Error: sn_data.txt not found!")
        return
    except Exception as e:
        print(f"Data loading error: {e}")
        return

    # ESTIF (fixed params, as in model; assumes no M offset needed, per project)
    mu_estif = estif.distance_modulus_estif_numerical(z_data)
    chi2_estif = np.sum(((mu_data - mu_estif) / sigma_mu)**2)
    estif_metrics = calculate_aic_bic(chi2_estif, n_data, 3)  # H0, A, BETA_DRAG

    # ΛCDM
    try:
        chi2_lcdm, best_om0 = compute_lcdm_chi2(z_data, mu_data, sigma_mu, optimize=optimize_lcdm)
        print(f"ΛCDM best Ω_m: {best_om0:.3f}")
    except Exception as e:
        print(f"ΛCDM computation error: {e}")
        return
    lcdm_metrics = calculate_aic_bic(chi2_lcdm, n_data, 2)  # Om0 + M (nuisance)

    # Print comparison
    print("\n" + "="*70)
    print("MODEL COMPARISON: ESTIF vs ΛCDM")
    print("="*70)
    print(f"\nDataset: {n_data} Type Ia Supernovae (optimized ΛCDM: {optimize_lcdm})\n")
    print(f"{'Metric':<20} {'ESTIF':<20} {'ΛCDM':<20} {'Δ (ESTIF-ΛCDM)':<20}")
    print("-"*70)
    for key in ['n_params', 'chi2', 'chi2_red', 'aic', 'bic']:
        est_val = estif_metrics[key]
        lcd_val = lcdm_metrics[key]
        delta = est_val - lcd_val
        print(f"{key.upper():<20} {est_val:<20.3f} {lcd_val:<20.3f} {delta:<20.3f}")
    print("="*70)

    # Interpretation
    delta_bic = estif_metrics['bic'] - lcdm_metrics['bic']
    print("\nInterpretation (based on BIC):")
    if delta_bic < -10:
        print(f"  • Strong evidence for ESTIF (ΔBIC = {delta_bic:.1f})")
    elif delta_bic < -2:
        print(f"  • Moderate evidence for ESTIF (ΔBIC = {delta_bic:.1f})")
    elif delta_bic > 10:
        print(f"  • Strong evidence for ΛCDM (ΔBIC = {delta_bic:.1f})")
    elif delta_bic > 2:
        print(f"  • Moderate evidence for ΛCDM (ΔBIC = {delta_bic:.1f})")
    else:
        print(f"  • Models roughly equivalent (ΔBIC = {delta_bic:.1f})")
    print("\nNote: ESTIF uses fixed params; ΛCDM χ² is minimized over Ω_m and M offset for fair comparison.")
    if plot:
        mu_lcdm = 5 * np.log10(FlatLambdaCDM(H0=70, Om0=best_om0).luminosity_distance(z_data).to_value('Mpc')) + 25
        plot_residuals(z_data, mu_data, sigma_mu, mu_estif, mu_lcdm)

if __name__ == "__main__":
    compare_models(optimize_lcdm=True, plot=True)  # Run with optimization and plot