"""
test_tilt_scan.py

Scans the tilt exponent n in sin(θ) = (Rs/r)^n to find which geometric
suppression model is consistent with EHT M87* observations.

BACKGROUND:
- Model A used n=0.5  → sin(θ) = √(Rs/r)   → 1.78σ tension
- Model B used n=1.0  → sin(θ) = Rs/r       → 2.52σ tension
- We scan n from 0.1 to 2.0 to find the EHT-consistent range.

The suppression factor is always: β(r) = cos(θ) = √(1 - (Rs/r)^(2n))
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# EHT M87* Observational Data
# ============================================================================
M_m87       = 6.5e9 * const.M_sun
shadow_obs  = 42.0   # μas
shadow_err  = 3.0    # μas
distance_m  = 16.8 * 3.086e22

Rs = 2 * const.G * M_m87 / const.c**2

# GR shadow diameter (corrected formula)
R_shadow_gr    = np.sqrt(27) * Rs
theta_gr_uas   = (R_shadow_gr / distance_m) * 206265 * 1e6
r_photon       = 1.5 * Rs

# ============================================================================
# General Tilt Model
# ============================================================================

def beta_tilt(r, Rs, n):
    """
    General tilt model: sin(θ) = (Rs/r)^n
    Suppression:        β(r)   = cos(θ) = √(1 - (Rs/r)^(2n))

    n=0.5  → Model A (stronger tilt)
    n=1.0  → Model B (weaker tilt)
    n→0    → no suppression (β→1)
    n→∞   → full suppression (β→0)
    """
    ratio = (Rs / r) ** (2 * n)
    if ratio >= 1.0:
        return 0.0
    return np.sqrt(1.0 - ratio)


def shadow_prediction(n):
    """
    Compute predicted shadow diameter (μas) for a given n.
    Returns (shadow_uas, deviation_percent, sigma_from_eht)
    """
    beta  = beta_tilt(r_photon, Rs, n)
    corr  = 1 + beta * (Rs / (2 * r_photon))
    shadow_uas = theta_gr_uas * corr
    deviation  = (corr - 1) * 100
    sigma      = (shadow_uas - shadow_obs) / shadow_err
    return shadow_uas, deviation, sigma


# ============================================================================
# Scan n from 0.05 to 3.0
# ============================================================================

n_values   = np.linspace(0.05, 3.0, 1000)
shadows    = []
deviations = []
sigmas     = []
betas      = []

for n in n_values:
    s, d, sig = shadow_prediction(n)
    shadows.append(s)
    deviations.append(d)
    sigmas.append(sig)
    betas.append(beta_tilt(r_photon, Rs, n))

shadows    = np.array(shadows)
deviations = np.array(deviations)
sigmas     = np.array(sigmas)
betas      = np.array(betas)

# ============================================================================
# Find Consistent Range
# ============================================================================

consistent_mask = np.abs(sigmas) < 1.0   # Within 1σ of EHT
tension_mask    = np.abs(sigmas) < 2.0   # Within 2σ of EHT

consistent_n = n_values[consistent_mask]
tension_n    = n_values[tension_mask]

# Find n where |sigma| is minimised (best fit)
best_idx = np.argmin(np.abs(sigmas))
best_n   = n_values[best_idx]
best_shadow, best_dev, best_sigma = shadow_prediction(best_n)
best_beta = beta_tilt(r_photon, Rs, best_n)

# ============================================================================
# Print Results
# ============================================================================

print("=" * 70)
print("TILT EXPONENT SCAN: sin(θ) = (Rs/r)^n")
print("Finding EHT-consistent geometric suppression")
print("=" * 70)

print(f"\nFixed Points:")
print(f"   EHT observed:    {shadow_obs:.1f} ± {shadow_err:.1f} μas")
print(f"   GR prediction:   {theta_gr_uas:.2f} μas (0.77σ from EHT)")
print(f"   Photon sphere:   r = 1.5 Rs")

print(f"\nPrevious Models:")
s_A, d_A, sig_A = shadow_prediction(0.5)
s_B, d_B, sig_B = shadow_prediction(1.0)
print(f"   Model A (n=0.5): {s_A:.2f} μas  |  {d_A:.1f}% deviation  |  {sig_A:.2f}σ  ⚠️")
print(f"   Model B (n=1.0): {s_B:.2f} μas  |  {d_B:.1f}% deviation  |  {sig_B:.2f}σ  ⚠️")

print(f"\n{'='*70}")
print(f"SCAN RESULTS")
print(f"{'='*70}")

print(f"\n✅ EHT-Consistent Range (within 1σ):")
if len(consistent_n) > 0:
    print(f"   n range:   {consistent_n[0]:.3f}  to  {consistent_n[-1]:.3f}")
    s_lo, _, _ = shadow_prediction(consistent_n[0])
    s_hi, _, _ = shadow_prediction(consistent_n[-1])
    print(f"   Shadow:    {s_lo:.2f} to {s_hi:.2f} μas")
else:
    print(f"   No values within 1σ found in scan range")

print(f"\n⚠️  EHT-Tension Range (within 2σ):")
if len(tension_n) > 0:
    print(f"   n range:   {tension_n[0]:.3f}  to  {tension_n[-1]:.3f}")
else:
    print(f"   No values within 2σ found in scan range")

print(f"\n🎯 Best-Fit Tilt Model:")
print(f"   Exponent n:         {best_n:.4f}")
print(f"   β at photon sphere: {best_beta:.4f}")
print(f"   Shadow prediction:  {best_shadow:.2f} μas")
print(f"   Deviation from GR:  {best_dev:.2f}%")
print(f"   EHT tension:        {best_sigma:.4f}σ")
print(f"   Tilt formula:       sin(θ) = (Rs/r)^{best_n:.3f}")
print(f"   Suppression:        β(r) = √(1 - (Rs/r)^{2*best_n:.3f})")

# Key n values table
print(f"\n{'='*70}")
print(f"KEY VALUES TABLE")
print(f"{'='*70}")
print(f"{'n':<8} {'β at 1.5Rs':<14} {'Shadow(μas)':<14} {'Deviation':<12} {'σ':<8} Status")
print("-" * 70)
for n_test in [0.1, 0.2, 0.3, best_n, 0.5, 0.75, 1.0, 1.5, 2.0]:
    s, d, sig = shadow_prediction(n_test)
    b = beta_tilt(r_photon, Rs, n_test)
    if abs(sig) < 1:
        status = "✅ Consistent"
    elif abs(sig) < 2:
        status = "⚠️  Tension"
    elif abs(sig) < 3:
        status = "⚠️  Moderate"
    else:
        status = "❌ Ruled out"
    marker = " ← BEST FIT" if abs(n_test - best_n) < 0.02 else ""
    print(f"{n_test:<8.3f} {b:<14.4f} {s:<14.2f} {d:<12.2f} {sig:<8.3f} {status}{marker}")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
fig.suptitle('Tilt Exponent Scan: Finding EHT-Consistent Suppression',
             fontsize=14, fontweight='bold')

# -- Plot 1: σ vs n --
axes[0].plot(n_values, sigmas, 'purple', linewidth=2.5)
axes[0].axhline(0,    color='black',  linewidth=1, linestyle='-')
axes[0].axhline(1,    color='green',  linewidth=1.5, linestyle='--', label='1σ boundary')
axes[0].axhline(-1,   color='green',  linewidth=1.5, linestyle='--')
axes[0].axhline(2,    color='orange', linewidth=1.5, linestyle='--', label='2σ boundary')
axes[0].axhline(-2,   color='orange', linewidth=1.5, linestyle='--')
axes[0].axhline(3,    color='red',    linewidth=1,   linestyle=':', label='3σ boundary')
axes[0].axhline(-3,   color='red',    linewidth=1,   linestyle=':')
axes[0].axvline(best_n, color='gold', linewidth=2,   linestyle='-', label=f'Best fit n={best_n:.3f}')

if len(consistent_n) > 0:
    axes[0].axvspan(consistent_n[0], consistent_n[-1],
                    alpha=0.2, color='green', label='1σ consistent range')

axes[0].scatter([0.5, 1.0], [sig_A, sig_B], color=['blue','red'],
                s=80, zorder=5, label='Model A (n=0.5) / B (n=1.0)')
axes[0].set_xlabel('Tilt Exponent n', fontsize=12)
axes[0].set_ylabel('EHT Tension (σ)', fontsize=12)
axes[0].set_title('EHT Tension vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)
axes[0].set_xlim(0.05, 3.0)
axes[0].set_ylim(-2, 6)

# -- Plot 2: Shadow size vs n --
axes[1].plot(n_values, shadows, 'purple', linewidth=2.5)
axes[1].axhline(shadow_obs, color='royalblue', linewidth=2,
                linestyle='-', label=f'EHT observed: {shadow_obs} μas')
axes[1].fill_between(n_values,
                     shadow_obs - shadow_err,
                     shadow_obs + shadow_err,
                     alpha=0.2, color='royalblue', label='EHT 1σ band')
axes[1].axhline(theta_gr_uas, color='green', linewidth=1.5,
                linestyle='--', label=f'GR: {theta_gr_uas:.1f} μas')
axes[1].axvline(best_n, color='gold', linewidth=2,
                linestyle='-', label=f'Best fit n={best_n:.3f}')
axes[1].scatter([0.5, 1.0],
                [shadow_prediction(0.5)[0], shadow_prediction(1.0)[0]],
                color=['blue','red'], s=80, zorder=5)
axes[1].set_xlabel('Tilt Exponent n', fontsize=12)
axes[1].set_ylabel('Predicted Shadow Diameter (μas)', fontsize=12)
axes[1].set_title('Shadow Size vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.3)
axes[1].set_xlim(0.05, 3.0)

# -- Plot 3: β at photon sphere vs n --
axes[2].plot(n_values, betas, 'purple', linewidth=2.5)
axes[2].axvline(best_n, color='gold', linewidth=2,
                linestyle='-', label=f'Best fit n={best_n:.3f}')
axes[2].axvline(0.5, color='blue', linewidth=1.5,
                linestyle='--', label='Model A (n=0.5)', alpha=0.7)
axes[2].axvline(1.0, color='red',  linewidth=1.5,
                linestyle='--', label='Model B (n=1.0)', alpha=0.7)
if len(consistent_n) > 0:
    axes[2].axvspan(consistent_n[0], consistent_n[-1],
                    alpha=0.2, color='green', label='1σ consistent range')
axes[2].set_xlabel('Tilt Exponent n', fontsize=12)
axes[2].set_ylabel('β(r) at Photon Sphere', fontsize=12)
axes[2].set_title('Suppression Factor β vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[2].legend(fontsize=9)
axes[2].grid(alpha=0.3)
axes[2].set_xlim(0.05, 3.0)

plt.tight_layout()
plt.savefig('tilt_scan_results.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: tilt_scan_results.png")
print("=" * 70)
