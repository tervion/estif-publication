"""
test_cosmo_correction_scaling.py

The previous test showed:
  - Sign is correct (s > 0)
  - Magnitude is 16× too large (s = 0.061 ≈ 1/16)
  - Fixed formula: χ²/dof = 1.629 (worse than ΛCDM)
  - Free amplitude: χ²/dof = 0.370 (slightly better than ΛCDM)

THE PROBLEM:
The current formula uses x(z) = x_0 × (1+z) which grows from
0.311 today to 0.777 at z=1.5 — pushing into strong-field regime
where the combined formula drops steeply. The correction becomes
huge (2.3 mag at z=1.5) and overwhelms the data.

THIS SCRIPT:
Tests five different approaches to the cosmological distance correction
to find which one naturally gives s ≈ 1 against the supernova data.

APPROACH A: Slower x(z) growth — x scales as (1+z)^α for α < 1
APPROACH B: Logarithmic — Δμ ∝ ln(1+z) 
APPROACH C: Integrated line-of-sight correction
APPROACH D: Small-angle expansion — first-order tilt correction only
APPROACH E: n-based correction — use n(z) directly not β(z)
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

z_data, mu_data, sigma_mu, sigma_int = [], [], [], []
with open(data_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                z_data.append(float(parts[1]))
                mu_data.append(float(parts[2]))
                sigma_mu.append(float(parts[3]))
                sigma_int.append(float(parts[4]))
            except ValueError:
                continue

z_data    = np.array(z_data)
mu_data   = np.array(mu_data)
sigma_tot = np.sqrt(np.array(sigma_mu)**2 + np.array(sigma_int)**2)
dof       = len(z_data) - 1

R_H          = const.c / const.H_0
r_universe_0 = 4.4e26
x_0          = R_H / r_universe_0   # 0.311 today
N_MAX        = 33.265
B_EXP        = 15.429

# ============================================================================
# Base functions
# ============================================================================

def mu_lcdm(z, M=0.0):
    return estif.distance_modulus_lcdm(z) + M

def fit_chi2(corrections, M_init=0.0):
    """Fit M_offset and amplitude s to data given correction array."""
    def chi2(params):
        M, s = params
        mu_pred = np.array([mu_lcdm(z, M) for z in z_data]) + s * corrections
        return np.sum(((mu_data - mu_pred) / sigma_tot)**2)
    res = minimize(chi2, [M_init, 0.0], method='Nelder-Mead',
                   options={'xatol':1e-8, 'fatol':1e-8, 'maxiter':10000})
    M_best, s_best = res.x
    chi2_best = chi2([M_best, s_best])
    return M_best, s_best, chi2_best, chi2_best / (dof - 1)

# ΛCDM baseline
res_lcdm  = minimize_scalar(lambda M: np.sum(((mu_data - np.array(
    [mu_lcdm(z, M) for z in z_data])) / sigma_tot)**2),
    bounds=(-2,2), method='bounded')
M_lcdm    = res_lcdm.x
chi2_lcdm = res_lcdm.fun
chi2_dof_lcdm = chi2_lcdm / dof

print("=" * 70)
print("COSMOLOGICAL CORRECTION SCALING INVESTIGATION")
print("=" * 70)
print(f"\n   ΛCDM baseline: χ²={chi2_lcdm:.2f}  χ²/dof={chi2_dof_lcdm:.4f}")
print(f"   Target: find correction with s ≈ 1.0")

# ============================================================================
# APPROACH A: x(z) = x_0 × (1+z)^α  for different α
# ============================================================================

print(f"\n{'='*70}")
print("APPROACH A: x(z) = x_0 × (1+z)^α  [slower growth]")
print(f"{'='*70}")

alpha_vals = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 1.00]

print(f"\n   {'α':<8} {'s_best':<10} {'χ²/dof':<12} {'Δχ²':<10} {'verdict'}")
print("   " + "-"*55)

best_A = None
for alpha in alpha_vals:
    corrections = np.array([
        -2.5 * np.log10(
            max(estif.observable_combined(x_0*(1+z)**alpha), 1e-10) /
            max(estif.observable_combined(x_0), 1e-10)
        ) for z in z_data
    ])
    M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections)
    delta = chi2_b - chi2_lcdm
    verdict = "✅ s≈1" if 0.5 < s_b < 2.0 else ("← close" if 0.2 < s_b < 3 else "")
    print(f"   {alpha:<8.3f} {s_b:<10.4f} {chi2_dof_b:<12.4f} {delta:<10.2f} {verdict}")
    if best_A is None or abs(s_b - 1.0) < abs(best_A[1] - 1.0):
        best_A = (alpha, s_b, chi2_dof_b, corrections.copy())

print(f"\n   Best α = {best_A[0]:.3f}  (s = {best_A[1]:.4f})")

# ============================================================================
# APPROACH B: Logarithmic correction Δμ = A × ln(1+z)
# ============================================================================

print(f"\n{'='*70}")
print("APPROACH B: Δμ = A × ln(1+z)  [logarithmic]")
print(f"{'='*70}")

# This is the simplest possible correction — find best A
corrections_log = np.log(1 + z_data)
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_log)
delta_b = chi2_b - chi2_lcdm
print(f"\n   s_best = {s_b:.4f}  χ²/dof = {chi2_dof_b:.4f}  Δχ² = {delta_b:.2f}")
print(f"   Physical meaning: Δμ = {s_b:.4f} × ln(1+z)")
if 0.5 < abs(s_b) < 2.0:
    print(f"   ✅ Logarithmic correction has right magnitude")

# ============================================================================
# APPROACH C: Integrated line-of-sight correction
# ============================================================================

print(f"\n{'='*70}")
print("APPROACH C: Integrated correction ∫₀ᶻ n(x(z')) dz'")
print(f"{'='*70}")

# Instead of evaluating at endpoint z, integrate the n(x) change
# along the entire line of sight
def integrated_correction(z, n_steps=50):
    """Integral of n(x(z')) from 0 to z, normalized."""
    z_vals = np.linspace(0, z, n_steps)
    n_vals = np.array([N_MAX * np.exp(-B_EXP * x_0*(1+zp)) for zp in z_vals])
    integral = np.trapz(n_vals, z_vals)
    n_0 = N_MAX * np.exp(-B_EXP * x_0)
    return integral / (n_0 * z) if z > 0 else 1.0

corrections_int = np.array([
    -2.5 * np.log10(max(integrated_correction(z), 1e-10))
    for z in z_data
])
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_int)
delta_c = chi2_b - chi2_lcdm
print(f"\n   s_best = {s_b:.4f}  χ²/dof = {chi2_dof_b:.4f}  Δχ² = {delta_c:.2f}")
if 0.5 < abs(s_b) < 2.0:
    print(f"   ✅ Integrated correction has right magnitude")

# ============================================================================
# APPROACH D: First-order expansion (small correction limit)
# ============================================================================

print(f"\n{'='*70}")
print("APPROACH D: First-order tilt — Δμ ≈ C × (x(z) - x_0)")
print(f"{'='*70}")

# When the tilt correction is small, first-order:
# Δμ ≈ -2.5/ln(10) × δobs/obs  ≈ C × (x(z) - x_0)
corrections_first = np.array([x_0*(1+z) - x_0 for z in z_data])  # = x_0 × z
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_first)
delta_d = chi2_b - chi2_lcdm
print(f"\n   s_best = {s_b:.4f}  χ²/dof = {chi2_dof_b:.4f}  Δχ² = {delta_d:.2f}")
print(f"   (Equivalent to Δμ = {s_b*x_0:.4f} × z at low z)")
if 0.5 < abs(s_b) < 2.0:
    print(f"   ✅ First-order correction has right magnitude")

# ============================================================================
# APPROACH E: n(z) directly — use change in tilt exponent
# ============================================================================

print(f"\n{'='*70}")
print("APPROACH E: n-based correction — Δμ ∝ n(x_0) - n(x(z))")
print(f"{'='*70}")

n_0 = N_MAX * np.exp(-B_EXP * x_0)
corrections_n = np.array([
    N_MAX * np.exp(-B_EXP * x_0) - N_MAX * np.exp(-B_EXP * x_0*(1+z))
    for z in z_data
])
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_n)
delta_e = chi2_b - chi2_lcdm
print(f"\n   n(x_0) = {n_0:.4f}")
print(f"   s_best = {s_b:.4f}  χ²/dof = {chi2_dof_b:.4f}  Δχ² = {delta_e:.2f}")
if 0.5 < abs(s_b) < 2.0:
    print(f"   ✅ n-based correction has right magnitude")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

approaches = [
    ("ΛCDM baseline",               1,     chi2_dof_lcdm,   0),
    (f"A: x^α (α={best_A[0]:.2f})", best_A[1], best_A[2],  None),
    ("B: ln(1+z)",                   s_b,   chi2_dof_b,      delta_b),
    ("C: Integrated n",              s_b,   chi2_dof_b,      delta_c),
    ("D: First-order x",             s_b,   chi2_dof_b,      delta_d),
    ("E: Δn directly",               s_b,   chi2_dof_b,      delta_e),
]

# Recompute all properly
results = []

# A
corrections = best_A[3]
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections)
results.append(("A: x^α", best_A[0], s_b, chi2_dof_b, chi2_b - chi2_lcdm))

# B
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(np.log(1 + z_data))
results.append(("B: ln(1+z)", None, s_b, chi2_dof_b, chi2_b - chi2_lcdm))

# C
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_int)
results.append(("C: Integrated", None, s_b, chi2_dof_b, chi2_b - chi2_lcdm))

# D
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_first)
results.append(("D: First-order", None, s_b, chi2_dof_b, chi2_b - chi2_lcdm))

# E
M_b, s_b, chi2_b, chi2_dof_b = fit_chi2(corrections_n)
results.append(("E: Δn direct", None, s_b, chi2_dof_b, chi2_b - chi2_lcdm))

print(f"\n   {'Approach':<20} {'param':<10} {'s_best':<10} "
      f"{'χ²/dof':<12} {'Δχ²':<10} {'verdict'}")
print("   " + "-"*70)
print(f"   {'ΛCDM (baseline)':<20} {'-':<10} {'1.000':<10} "
      f"{chi2_dof_lcdm:<12.4f} {'0.00':<10} reference")

for name, param, s_b, chi2_dof_b, delta in results:
    param_str = f"{param:.3f}" if param is not None else "-"
    if abs(delta) < 2 and s_b > 0:
        verdict = "✅ consistent"
    elif s_b > 0 and delta < 0:
        verdict = "✅ improves"
    elif s_b > 0:
        verdict = "⚠️  wrong magnitude"
    else:
        verdict = "❌ wrong sign"
    print(f"   {name:<20} {param_str:<10} {s_b:<10.4f} "
          f"{chi2_dof_b:<12.4f} {delta:<10.2f} {verdict}")

# Find best approach
best_approach = min(results, key=lambda x: abs(x[2] - 1.0))
print(f"\n   Best approach (s closest to 1.0): {best_approach[0]}")
print(f"   s = {best_approach[2]:.4f}  χ²/dof = {best_approach[3]:.4f}")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Correction Scaling: Finding s ≈ 1',
             fontsize=13, fontweight='bold')

z_smooth = np.linspace(0.01, 1.5, 200)

# Plot 1: Correction shapes
ax = axes[0]
corr_A = [-2.5*np.log10(max(estif.observable_combined(x_0*(1+z)**best_A[0]),1e-10)/
           max(estif.observable_combined(x_0),1e-10)) for z in z_smooth]
corr_B = np.log(1 + np.array(z_smooth))
corr_D = x_0 * np.array(z_smooth)
corr_orig = [-2.5*np.log10(max(estif.observable_combined(x_0*(1+z)),1e-10)/
              max(estif.observable_combined(x_0),1e-10)) for z in z_smooth]

ax.plot(z_smooth, corr_orig,          'gray',       linewidth=1.5,
        linestyle='--', label='Original (too large)', alpha=0.7)
ax.plot(z_smooth, corr_A,             'blue',       linewidth=2,
        label=f'A: x^{best_A[0]:.2f}')
ax.plot(z_smooth, [0.05*c for c in corr_B], 'red', linewidth=2,
        label='B: ~0.05×ln(1+z)')
ax.plot(z_smooth, [0.1*c for c in corr_D],  'green', linewidth=2,
        label='D: ~0.1×x₀z')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Correction Δμ (mag)', fontsize=11)
ax.set_title('Correction Shapes', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_ylim(-0.1, 0.5)

# Plot 2: s values across approaches
ax = axes[1]
approach_names = [r[0] for r in results]
s_vals         = [r[2] for r in results]
colors = ['green' if 0.5 < s < 2.0 else
          'orange' if 0 < s else 'red' for s in s_vals]

bars = ax.bar(approach_names, s_vals, color=colors, alpha=0.8,
              edgecolor='black')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--',
           label='s = 1 (perfect)')
ax.axhspan(0.5, 2.0, alpha=0.1, color='green', label='Acceptable range')
ax.set_ylabel('Best-fit scale s', fontsize=11)
ax.set_title('Which Approach Gives s ≈ 1?', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(axis='y', alpha=0.3)
ax.set_ylim(-0.5, max(max(s_vals)*1.2, 2.5))

for bar, s in zip(bars, s_vals):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.05,
            f'{s:.3f}', ha='center', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('correction_scaling_results.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: correction_scaling_results.png")
print("=" * 70)
