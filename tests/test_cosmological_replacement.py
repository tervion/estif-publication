"""
test_cosmological_replacement.py

Tests whether the ESTIF tilt geometry can modify the distance-redshift
relation in a way that fits or improves on supernova data.

THE IDEA:
At each redshift z, the universe was smaller:
    r_universe(z) = r_universe_0 / (1+z)

So the curvature ratio x changes with redshift:
    x(z) = R_H / r_universe(z) = x_0 × (1+z)

As z increases, x increases, n(x) decreases, β changes.
The tilt correction accumulated along the line of sight modifies
how bright supernovae appear — i.e. it modifies the distance modulus.

The tilt correction to distance modulus:
    Δμ(z) = -2.5 × log10(observable(x(z)) / observable(x_0))

where observable = √β from the combined formula.

TEST:
Does μ_ESTIF(z) = μ_ΛCDM(z) + Δμ(z) fit the supernova data
better, worse, or the same as pure ΛCDM?

If better → the tilt geometry contributes to explaining expansion
If same   → consistent but not distinguishable
If worse  → the correction moves in the wrong direction

DATA: data/sn_data.txt
Format: name  z  mu  sigma_mu  sigma_intrinsic
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Load supernova data
# ============================================================================

data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../data/sn_data.txt')

names, z_data, mu_data, sigma_mu, sigma_int = [], [], [], [], []

with open(data_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                names.append(parts[0])
                z_data.append(float(parts[1]))
                mu_data.append(float(parts[2]))
                sigma_mu.append(float(parts[3]))
                sigma_int.append(float(parts[4]))
            except ValueError:
                continue

z_data    = np.array(z_data)
mu_data   = np.array(mu_data)
sigma_mu  = np.array(sigma_mu)
sigma_int = np.array(sigma_int)
sigma_tot = np.sqrt(sigma_mu**2 + sigma_int**2)

print("=" * 70)
print("COSMOLOGICAL REPLACEMENT TEST")
print("=" * 70)
print(f"\n   Loaded {len(z_data)} supernovae")
print(f"   Redshift range: {z_data.min():.3f} – {z_data.max():.3f}")

# ============================================================================
# ΛCDM baseline
# ============================================================================

def mu_lcdm(z, M_offset=0.0):
    """Standard ΛCDM distance modulus with optional magnitude offset."""
    return estif.distance_modulus_lcdm(z) + M_offset

# Fit M_offset to data
def chi2_lcdm(M_offset):
    mu_pred = mu_lcdm(z_data, M_offset)
    return np.sum(((mu_data - mu_pred) / sigma_tot)**2)

result_lcdm = minimize_scalar(chi2_lcdm, bounds=(-2, 2), method='bounded')
M_lcdm      = result_lcdm.x
chi2_lcdm_best = chi2_lcdm(M_lcdm)
dof         = len(z_data) - 1
chi2_dof_lcdm = chi2_lcdm_best / dof

print(f"\n{'='*70}")
print("ΛCDM BASELINE")
print(f"{'='*70}")
print(f"\n   Best M_offset: {M_lcdm:.4f} mag")
print(f"   χ²:            {chi2_lcdm_best:.2f}")
print(f"   χ²/dof:        {chi2_dof_lcdm:.4f}  (target: ~1.0)")

# ============================================================================
# ESTIF tilt correction
# ============================================================================

R_H          = const.c / const.H_0
r_universe_0 = 4.4e26
x_0          = R_H / r_universe_0          # 0.311 today

def x_at_z(z):
    """Curvature ratio at redshift z — universe was smaller by (1+z)."""
    r_univ_z = r_universe_0 / (1.0 + z)
    return R_H / r_univ_z                  # = x_0 × (1+z)

def tilt_correction_mu(z):
    """
    Distance modulus correction from tilt geometry.

    As we look back in time (higher z), the universe was smaller,
    so x was larger, n was smaller, and the tilt observable √β was smaller.

    A smaller √β means less of the 4D correction projected into 3D —
    so light from the past was affected differently by the tilt geometry.

    Δμ = -2.5 × log10(observable(z) / observable(0))

    Positive Δμ means supernova appears dimmer (farther away).
    Negative Δμ means supernova appears brighter (closer).
    """
    obs_z   = estif.observable_combined(x_at_z(z))
    obs_now = estif.observable_combined(x_0)
    if obs_z <= 0 or obs_now <= 0:
        return 0.0
    return -2.5 * np.log10(obs_z / obs_now)

# Vectorise
tilt_corrections = np.array([tilt_correction_mu(z) for z in z_data])

print(f"\n{'='*70}")
print("TILT CORRECTION PROPERTIES")
print(f"{'='*70}")
print(f"\n   x_0 (today):           {x_0:.4f}")
print(f"   x at z=0.1:            {x_at_z(0.1):.4f}")
print(f"   x at z=0.5:            {x_at_z(0.5):.4f}")
print(f"   x at z=1.0:            {x_at_z(1.0):.4f}")
print(f"   x at z=1.5:            {x_at_z(1.5):.4f}")
print(f"\n   Δμ at z=0.1:           {tilt_correction_mu(0.1):+.4f} mag")
print(f"   Δμ at z=0.5:           {tilt_correction_mu(0.5):+.4f} mag")
print(f"   Δμ at z=1.0:           {tilt_correction_mu(1.0):+.4f} mag")
print(f"   Δμ at z=1.5:           {tilt_correction_mu(1.5):+.4f} mag")
print(f"\n   Correction range in data: "
      f"{tilt_corrections.min():+.4f} to {tilt_corrections.max():+.4f} mag")

# ============================================================================
# ESTIF fit — M_offset only
# ============================================================================

def mu_estif(z, M_offset=0.0):
    """ΛCDM + tilt correction."""
    return mu_lcdm(z, M_offset) + np.array([tilt_correction_mu(zi)
                                              for zi in np.atleast_1d(z)])

def chi2_estif(M_offset):
    mu_pred = mu_estif(z_data, M_offset)
    return np.sum(((mu_data - mu_pred) / sigma_tot)**2)

result_estif  = minimize_scalar(chi2_estif, bounds=(-2, 2), method='bounded')
M_estif       = result_estif.x
chi2_estif_best = chi2_estif(M_estif)
chi2_dof_estif  = chi2_estif_best / dof

print(f"\n{'='*70}")
print("ESTIF COMBINED FORMULA FIT")
print(f"{'='*70}")
print(f"\n   Best M_offset: {M_estif:.4f} mag")
print(f"   χ²:            {chi2_estif_best:.2f}")
print(f"   χ²/dof:        {chi2_dof_estif:.4f}")

# ============================================================================
# ESTIF fit with tilt amplitude free — scale factor s
# ============================================================================
# Allow the tilt correction to be scaled: Δμ → s × Δμ
# If best s ≈ 1 → formula is correct
# If best s > 1 → correction needs amplification
# If best s < 0 → correction has wrong sign

def chi2_estif_scaled(params):
    M_offset, s = params
    mu_pred = mu_lcdm(z_data, M_offset) + s * tilt_corrections
    return np.sum(((mu_data - mu_pred) / sigma_tot)**2)

result_scaled = minimize(chi2_estif_scaled, [M_lcdm, 1.0],
                         method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-8})
M_scaled, s_best = result_scaled.x
chi2_scaled = chi2_estif_scaled([M_scaled, s_best])
chi2_dof_scaled = chi2_scaled / (dof - 1)

print(f"\n{'='*70}")
print("ESTIF WITH FREE AMPLITUDE (s × Δμ)")
print(f"{'='*70}")
print(f"\n   Best M_offset: {M_scaled:.4f} mag")
print(f"   Best scale s:  {s_best:.4f}")
print(f"   χ²:            {chi2_scaled:.2f}")
print(f"   χ²/dof:        {chi2_dof_scaled:.4f}")
print(f"\n   s = {s_best:.4f} means:")
if abs(s_best) < 0.1:
    print(f"   Tilt correction is negligible — data prefers pure ΛCDM")
elif s_best < 0:
    print(f"   Correction has WRONG SIGN — tilt makes fit worse")
elif 0.5 < s_best < 2.0:
    print(f"   Correction has RIGHT SIGN and reasonable magnitude")
    print(f"   Tilt geometry contributes to explaining the data")
elif s_best > 2.0:
    print(f"   Correction has right sign but needs amplification")

# ============================================================================
# Summary comparison
# ============================================================================

print(f"\n{'='*70}")
print("COMPARISON SUMMARY")
print(f"{'='*70}")
print(f"\n   {'Model':<30} {'χ²':<12} {'χ²/dof':<12} {'Δχ²'}")
print("   " + "-"*55)
print(f"   {'ΛCDM (baseline)':<30} {chi2_lcdm_best:<12.2f} "
      f"{chi2_dof_lcdm:<12.4f} 0.00 (reference)")
print(f"   {'ESTIF (fixed formula)':<30} {chi2_estif_best:<12.2f} "
      f"{chi2_dof_estif:<12.4f} {chi2_estif_best-chi2_lcdm_best:+.2f}")
print(f"   {'ESTIF (free amplitude)':<30} {chi2_scaled:<12.2f} "
      f"{chi2_dof_scaled:<12.4f} {chi2_scaled-chi2_lcdm_best:+.2f}")

delta_chi2 = chi2_estif_best - chi2_lcdm_best

print(f"\n{'='*70}")
print("VERDICT")
print(f"{'='*70}")

if delta_chi2 < -2:
    print(f"\n   ✅ ESTIF IMPROVES THE FIT by Δχ² = {delta_chi2:.2f}")
    print(f"   The tilt geometry correction moves in the right direction.")
    print(f"   Expansion data is consistent with 4D inward flow.")
elif abs(delta_chi2) <= 2:
    print(f"\n   ⚠️  ESTIF IS CONSISTENT WITH DATA (Δχ² = {delta_chi2:.2f})")
    print(f"   Tilt correction is within noise — cannot distinguish from ΛCDM.")
    print(f"   Not ruled out, but not yet a detectable improvement.")
else:
    print(f"\n   ❌ ESTIF WORSENS THE FIT by Δχ² = {delta_chi2:.2f}")
    print(f"   The specific tilt correction formula moves in wrong direction.")
    print(f"   Formula needs adjustment before cosmological replacement.")

print(f"\n   Sign of correction: {'correct ✅' if s_best > 0 else 'wrong ❌'}")
print(f"   Magnitude: {abs(s_best):.3f}× (1.0 = perfect)")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Cosmological Replacement Test: ESTIF Tilt vs Supernova Data',
             fontsize=13, fontweight='bold')

z_smooth = np.linspace(0.01, z_data.max() * 1.05, 300)
mu_lcdm_smooth  = np.array([mu_lcdm(z, M_lcdm) for z in z_smooth])
mu_estif_smooth = np.array([mu_estif(z, M_estif)[0] for z in z_smooth])
corr_smooth     = np.array([tilt_correction_mu(z) for z in z_smooth])

# Plot 1: Hubble diagram
ax = axes[0]
ax.errorbar(z_data, mu_data, yerr=sigma_tot, fmt='.', color='gray',
            alpha=0.4, markersize=3, elinewidth=0.5, label='SN data')
ax.plot(z_smooth, mu_lcdm_smooth,  'blue', linewidth=2,   label='ΛCDM')
ax.plot(z_smooth, mu_estif_smooth, 'red',  linewidth=2,
        linestyle='--', label='ESTIF (tilt corrected)')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Distance Modulus μ', fontsize=11)
ax.set_title('Hubble Diagram', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 2: Residuals from ΛCDM
ax = axes[1]
res_lcdm  = mu_data - np.array([mu_lcdm(z, M_lcdm)  for z in z_data])
res_estif = mu_data - np.array([mu_estif(z, M_estif)[0] for z in z_data])

ax.scatter(z_data, res_lcdm,  color='blue', alpha=0.4, s=8,
           label=f'ΛCDM residuals (χ²/dof={chi2_dof_lcdm:.3f})')
ax.scatter(z_data, res_estif, color='red',  alpha=0.4, s=8,
           label=f'ESTIF residuals (χ²/dof={chi2_dof_estif:.3f})')
ax.axhline(0, color='black', linewidth=1)
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Residual (mag)', fontsize=11)
ax.set_title('Residuals from Best Fit', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_ylim(-1, 1)

# Plot 3: Tilt correction vs z
ax = axes[2]
ax.plot(z_smooth, corr_smooth * 1000, 'purple', linewidth=2.5)
ax.axhline(0, color='black', linewidth=1, linestyle='--')
ax.fill_between(z_smooth, 0, corr_smooth * 1000,
                alpha=0.2, color='purple',
                label='Tilt correction')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Tilt Correction Δμ (millimag)', fontsize=11)
ax.set_title('ESTIF Tilt Correction vs Redshift',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.text(0.05, 0.95,
        f'Max correction: {max(abs(corr_smooth))*1000:.2f} mmag\n'
        f'Best scale s = {s_best:.3f}\n'
        f'Δχ² = {delta_chi2:+.2f}',
        transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig('cosmological_replacement_test.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: cosmological_replacement_test.png")
print("=" * 70)
