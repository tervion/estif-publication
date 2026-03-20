"""
test_tilt_models.py

Tests two geometric suppression models for β against EHT M87* data.

BACKGROUND:
The full 4D Lorentz correction (β=1) is suppressed because 3D space sits
at an angle inside 4D space. Mass tilts the 3D hypersurface, hiding part
of the correction in the invisible 4th dimension.

The suppression factor = cos(θ) where θ is the tilt angle.

Two candidate tilt models:
  Model A: sin(θ) = √(Rs/r)  →  β(r) = √(1 - Rs/r)
  Model B: sin(θ) = Rs/r     →  β(r) = √(1 - (Rs/r)²)

EHT M87* data constrains which model (if either) is consistent.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../../estif_publication/src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# EHT M87* Observational Data
# ============================================================================
M_m87 = 6.5e9 * const.M_sun
shadow_observed_uas = 42.0
shadow_error_uas = 3.0
distance_m = 16.8 * 3.086e22

Rs = 2 * const.G * M_m87 / const.c**2

# GR shadow prediction (photon sphere geometry)
R_shadow_gr = np.sqrt(27) * Rs
theta_shadow_gr_uas = (R_shadow_gr / distance_m) * 206265 * 1e6

# ============================================================================
# Beta Functions: Two Tilt Models
# ============================================================================

def beta_model_A(r, Rs):
    """
    Model A: sin(θ) = √(Rs/r)
    Suppression: cos(θ) = √(1 - Rs/r)
    Stronger tilt — larger suppression near the black hole.
    """
    ratio = Rs / r
    if ratio >= 1.0:
        return 0.0  # Fully tilted at horizon
    return np.sqrt(1.0 - ratio)


def beta_model_B(r, Rs):
    """
    Model B: sin(θ) = Rs/r
    Suppression: cos(θ) = √(1 - (Rs/r)²)
    Weaker tilt — smaller suppression near the black hole.
    """
    ratio = Rs / r
    if ratio >= 1.0:
        return 0.0
    return np.sqrt(1.0 - ratio**2)


def lensing_with_beta(r, M, beta_func):
    """
    Modified lensing equation using geometric suppression:
    θ_ESTIF = θ_GR × (1 + β(r) × Rs/(2r))
    """
    Rs = 2 * const.G * M / const.c**2
    theta_gr = (4 * const.G * M) / (const.c**2 * r)
    beta = beta_func(r, Rs)
    correction = 1 + beta * (Rs / (2 * r))
    return theta_gr * correction


# ============================================================================
# Calculate Shadow Predictions for Both Models
# ============================================================================

r_photon = 1.5 * Rs  # Photon sphere

# Model A prediction
theta_gr_photon = (4 * const.G * M_m87) / (const.c**2 * r_photon)
theta_A = lensing_with_beta(r_photon, M_m87, beta_model_A)
deviation_A = (theta_A / theta_gr_photon - 1)
shadow_A_uas = theta_shadow_gr_uas * (1 + deviation_A)
sigma_A = (shadow_A_uas - shadow_observed_uas) / shadow_error_uas

# Model B prediction
theta_B = lensing_with_beta(r_photon, M_m87, beta_model_B)
deviation_B = (theta_B / theta_gr_photon - 1)
shadow_B_uas = theta_shadow_gr_uas * (1 + deviation_B)
sigma_B = (shadow_B_uas - shadow_observed_uas) / shadow_error_uas

# Original β=1 (no suppression)
theta_raw = theta_gr_photon * (1 + 1.0 * (Rs / (2 * r_photon)))
deviation_raw = (theta_raw / theta_gr_photon - 1)
shadow_raw_uas = theta_shadow_gr_uas * (1 + deviation_raw)
sigma_raw = (shadow_raw_uas - shadow_observed_uas) / shadow_error_uas

# ============================================================================
# Print Results
# ============================================================================

print("=" * 70)
print("GEOMETRIC SUPPRESSION: TWO TILT MODELS vs EHT M87*")
print("=" * 70)

print(f"\nSetup:")
print(f"   M87* mass:          {M_m87/const.M_sun:.2e} M_sun")
print(f"   Schwarzschild Rs:   {Rs:.3e} m")
print(f"   Photon sphere:      r = 1.5 Rs")
print(f"   EHT observed:       {shadow_observed_uas:.1f} ± {shadow_error_uas:.1f} μas")
print(f"   GR prediction:      {theta_shadow_gr_uas:.2f} μas")

print(f"\n{'='*70}")
print("β = 1.0 (No Suppression — Full Lorentz Correction)")
print(f"{'='*70}")
print(f"   β value at photon sphere: 1.000")
print(f"   Lensing deviation:        {deviation_raw*100:.2f}%")
print(f"   Shadow prediction:        {shadow_raw_uas:.2f} μas")
print(f"   EHT tension:              {sigma_raw:.2f}σ")
print(f"   Verdict:                  {'✅ Consistent' if abs(sigma_raw)<1 else ('⚠️  Tension' if abs(sigma_raw)<3 else '❌ Ruled out')}")

print(f"\n{'='*70}")
print("Model A: sin(θ) = √(Rs/r)  →  β(r) = √(1 - Rs/r)  [Stronger Tilt]")
print(f"{'='*70}")
beta_A_val = beta_model_A(r_photon, Rs)
print(f"   β value at photon sphere: {beta_A_val:.4f}")
print(f"   Lensing deviation:        {deviation_A*100:.2f}%")
print(f"   Shadow prediction:        {shadow_A_uas:.2f} μas")
print(f"   EHT tension:              {sigma_A:.2f}σ")
print(f"   Verdict:                  {'✅ Consistent' if abs(sigma_A)<1 else ('⚠️  Tension' if abs(sigma_A)<3 else '❌ Ruled out')}")

print(f"\n{'='*70}")
print("Model B: sin(θ) = Rs/r  →  β(r) = √(1 - (Rs/r)²)  [Weaker Tilt]")
print(f"{'='*70}")
beta_B_val = beta_model_B(r_photon, Rs)
print(f"   β value at photon sphere: {beta_B_val:.4f}")
print(f"   Lensing deviation:        {deviation_B*100:.2f}%")
print(f"   Shadow prediction:        {shadow_B_uas:.2f} μas")
print(f"   EHT tension:              {sigma_B:.2f}σ")
print(f"   Verdict:                  {'✅ Consistent' if abs(sigma_B)<1 else ('⚠️  Tension' if abs(sigma_B)<3 else '❌ Ruled out')}")

print(f"\n{'='*70}")
print("SUMMARY TABLE")
print(f"{'='*70}")
print(f"{'Model':<30} {'β at 1.5Rs':<12} {'Deviation':<12} {'Shadow(μas)':<14} {'EHT (σ)':<10} {'Status'}")
print("-"*90)
print(f"{'No suppression (β=1)':<30} {'1.000':<12} {f'{deviation_raw*100:.1f}%':<12} {f'{shadow_raw_uas:.2f}':<14} {f'{sigma_raw:.2f}':<10} {'❌ Ruled out'}")
print(f"{'Model A [√(1-Rs/r)]':<30} {f'{beta_A_val:.4f}':<12} {f'{deviation_A*100:.1f}%':<12} {f'{shadow_A_uas:.2f}':<14} {f'{sigma_A:.2f}':<10} {'✅ Consistent' if abs(sigma_A)<1 else ('⚠️ Tension' if abs(sigma_A)<3 else '❌ Ruled out')}")
print(f"{'Model B [√(1-(Rs/r)²)]':<30} {f'{beta_B_val:.4f}':<12} {f'{deviation_B*100:.1f}%':<12} {f'{shadow_B_uas:.2f}':<14} {f'{sigma_B:.2f}':<10} {'✅ Consistent' if abs(sigma_B)<1 else ('⚠️ Tension' if abs(sigma_B)<3 else '❌ Ruled out')}")
print(f"{'EHT Observed':<30} {'-':<12} {'-':<12} {f'{shadow_observed_uas:.1f}±{shadow_error_uas}':<14} {'0.00':<10} {'Reference'}")

# ============================================================================
# Visualization
# ============================================================================

r_range = np.linspace(1.5, 20, 500) * Rs

beta_A_curve = np.array([beta_model_A(r, Rs) for r in r_range])
beta_B_curve = np.array([beta_model_B(r, Rs) for r in r_range])
beta_raw_curve = np.ones_like(r_range)

dev_A_curve = beta_A_curve * (Rs / (2 * r_range)) * 100
dev_B_curve = beta_B_curve * (Rs / (2 * r_range)) * 100
dev_raw_curve = beta_raw_curve * (Rs / (2 * r_range)) * 100

fig, axes = plt.subplots(1, 3, figsize=(16, 6))
fig.suptitle('Geometric Suppression: Two Tilt Models vs EHT M87*',
             fontsize=14, fontweight='bold')

r_Rs = r_range / Rs

# Plot 1: β vs distance
axes[0].plot(r_Rs, beta_raw_curve, 'k--', linewidth=2, label='No suppression (β=1)')
axes[0].plot(r_Rs, beta_A_curve, 'b-', linewidth=2.5, label='Model A: √(1-Rs/r)')
axes[0].plot(r_Rs, beta_B_curve, 'r-', linewidth=2.5, label='Model B: √(1-(Rs/r)²)')
axes[0].axvline(1.5, color='gray', linestyle=':', alpha=0.7, label='Photon sphere (1.5Rs)')
axes[0].set_xlabel('Distance (r / Rs)', fontsize=12)
axes[0].set_ylabel('β(r) — Suppression Factor', fontsize=12)
axes[0].set_title('β vs Distance from Black Hole', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=10)
axes[0].grid(alpha=0.3)
axes[0].set_xlim(1.5, 20)
axes[0].set_ylim(0, 1.05)

# Plot 2: Lensing deviation vs distance
axes[1].plot(r_Rs, dev_raw_curve, 'k--', linewidth=2, label='No suppression (β=1)')
axes[1].plot(r_Rs, dev_A_curve, 'b-', linewidth=2.5, label='Model A')
axes[1].plot(r_Rs, dev_B_curve, 'r-', linewidth=2.5, label='Model B')
axes[1].axvline(1.5, color='gray', linestyle=':', alpha=0.7, label='Photon sphere')
axes[1].axhline(10, color='orange', linestyle='--', alpha=0.7, label='EHT upper bound (~10%)')
axes[1].axhline(3, color='green', linestyle='--', alpha=0.7, label='EHT 1σ safe zone (~3%)')
axes[1].set_xlabel('Distance (r / Rs)', fontsize=12)
axes[1].set_ylabel('Lensing Deviation from GR (%)', fontsize=12)
axes[1].set_title('Predicted Deviation vs Distance', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.3)
axes[1].set_xlim(1.5, 20)

# Plot 3: Shadow size comparison bar chart
labels = ['EHT\nObserved', 'GR', 'β=1\n(No supp.)', 'Model A\n[√(1-Rs/r)]', 'Model B\n[√(1-(Rs/r)²)]']
values = [shadow_observed_uas, theta_shadow_gr_uas, shadow_raw_uas, shadow_A_uas, shadow_B_uas]
colors = ['royalblue', 'green', 'black', 'blue', 'red']
alphas = [0.8, 0.7, 0.5, 0.8, 0.8]

bars = axes[2].bar(labels, values, color=colors, alpha=0.75, edgecolor='black', linewidth=1.5)
axes[2].errorbar([0], [shadow_observed_uas], yerr=[shadow_error_uas],
                 fmt='none', color='royalblue', capsize=8, linewidth=2)
axes[2].fill_between([-0.5, 4.5],
                     shadow_observed_uas - shadow_error_uas,
                     shadow_observed_uas + shadow_error_uas,
                     alpha=0.15, color='royalblue', label='EHT 1σ band')
axes[2].set_ylabel('Shadow Diameter (μas)', fontsize=12)
axes[2].set_title('Shadow Size: Models vs EHT', fontsize=12, fontweight='bold')
axes[2].legend(fontsize=10)
axes[2].grid(axis='y', alpha=0.3)
axes[2].set_ylim(15, 50)

for bar, val in zip(bars, values):
    axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{val:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('tilt_model_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: tilt_model_comparison.png")
print("=" * 70)
