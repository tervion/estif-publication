# estif_distance_tests.py (improved from test_fixes.py)
"""
Test ESTIF distance modulus for multiple z; compare to ΛCDM; assertions and plot.
"""

import numpy as np
import estif_ec_fd_model as estif
import estif_ec_fd_constants as const
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

def test_distances(z_tests=[0.1, 0.5, 1.0, 2.0], plot=False):
    """Test distance modulus for multiple z in ESTIF vs ΛCDM"""
    t_obs = 4.35e17  # ~13.8 Gyr in s; could derive from model if needed
    lcdm = FlatLambdaCDM(H0=const.H_0 * 3.08568e19 / 1e3, Om0=0.28)  # H0 in km/s/Mpc

    mu_estif_list = []
    mu_lcdm_list = []

    for z_test in z_tests:
        print(f"\n=== Testing at z={z_test} ===")
        try:
            t_emit = estif.t_from_z(z_test, t_obs=t_obs)
            print(f"t_emit = {t_emit:.3e} s ({t_emit / (365.25 * 86400 * 1e9):.2f} Gyr)")

            H_emit = estif.H_variable(t_emit)
            H_obs = estif.H_variable(t_obs)
            print(f"H_variable(t_emit) = {H_emit:.3e} s^-1")
            print(f"H_variable(t_obs) = {H_obs:.3e} s^-1")
            print(f"H0 = {const.H_0:.3e} s^-1")

            mu_estif = estif.distance_modulus_estif_numerical(z_test)
            mu_lcdm = 5 * np.log10(lcdm.luminosity_distance(z_test).to_value('Mpc')) + 25
            print(f"ESTIF μ = {mu_estif:.2f} mag")
            print(f"ΛCDM μ = {mu_lcdm:.2f} mag (fixed Om0=0.28)")
            print(f"Δμ (ESTIF - ΛCDM) = {mu_estif - mu_lcdm:.2f} mag")

            # Assertions
            assert 0 < mu_estif < 50, f"Unexpected μ for z={z_test}: {mu_estif}"
            assert abs(mu_estif - mu_lcdm) < 5, f"Large discrepancy for z={z_test}: Δ={mu_estif - mu_lcdm}"
            print("Assertions passed.")

            mu_estif_list.append(mu_estif)
            mu_lcdm_list.append(mu_lcdm)
        except Exception as e:
            print(f"Error at z={z_test}: {e}")

    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(z_tests, mu_estif_list, 'b-o', label='ESTIF')
        plt.plot(z_tests, mu_lcdm_list, 'r--o', label='ΛCDM')
        plt.xlabel('Redshift z')
        plt.ylabel('Distance Modulus μ (mag)')
        plt.legend()
        plt.title('Distance Modulus Comparison')
        plt.savefig('distance_comparison.png')
        plt.close()
        print("\nPlot saved as 'distance_comparison.png'")

if __name__ == "__main__":
    test_distances(plot=True)