"""
test_hubble_tilt_cosmology.py

EXPLORATORY: Tests whether substituting the Hubble radius for Rs in the
tilt formula produces a cosmological constant Λ consistent with ΛCDM.

THE QUESTION:
The tilt formula sin(θ) = (Rs/r)^n works locally near black holes.
Rs is the Schwarzschild radius — a curvature scale derived from local mass.

At cosmological scales there is no single Rs. But there is a cosmological
equivalent: the Hubble radius R_H = c/H_0 — the curvature scale of the
universe as a whole.

If we substitute:
    sin(θ)_cosmic = (R_H / r_universe)^n

And derive what Λ this implies, does it match the measured value?

    Λ_measured = 1.1056e-52 m⁻²  (Planck 2018)

If Λ_predicted / Λ_measured ≈ 1 for n in [0.05, 0.215] → bridge found.
If ratio is far from 1 → Hubble radius alone is insufficient.

Either result is scientifically useful.

NOTE: This is a new exploratory test. No existing files are modified.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Known Constants
# ============================================================================

# Hubble radius: the cosmological equivalent of Rs
# This is how far light can travel in one Hubble time — the curvature
# scale of the universe as a whole
R_H = const.c / const.H_0          # meters — ~1.37e26 m (~14.5 billion light years)

# Size of the observable universe (comoving radius)
# This plays the role of r — the distance from the "center" to the edge
r_universe = 4.4e26                 # meters (~46.5 billion light years)

# Measured cosmological constant (Planck 2018)
# This is what we need to reproduce
LAMBDA_MEASURED = 1.1056e-52        # m⁻²

# EHT-constrained tilt exponent range
N_LO   = 0.05
N_HI   = 0.215
N_BEST = 0.05

# ============================================================================
# Core Calculation
# ============================================================================

def beta_cosmic(n):
    """
    Cosmic suppression factor using Hubble radius substitution.

    Exactly the same structure as the local tilt formula:
        sin(θ) = (Rs/r)^n  →  (R_H / r_universe)^n
        β      = cos(θ)    = √(1 - sin²(θ))
    """
    ratio = (R_H / r_universe) ** (2 * n)
    if ratio >= 1.0:
        return 0.0
    return np.sqrt(1.0 - ratio)


def lambda_predicted(n):
    """
    Predicted cosmological constant from tilt geometry.

    The 4D flow energy density that projects into 3D as Λ:

        Λ = (3 / R_H²) × β_cosmic²

    Physical meaning:
    - 3 / R_H²  is the natural curvature scale of the universe
    - β_cosmic²  is how much of the 4D flow energy
                 projects into observable 3D space
    """
    beta = beta_cosmic(n)
    return (3.0 / R_H**2) * beta**2


def lambda_ratio(n):
    """Ratio of predicted to measured Λ. Goal: as close to 1.0 as possible."""
    return lambda_predicted(n) / LAMBDA_MEASURED


# ============================================================================
# Print Results
# ============================================================================

print("=" * 70)
print("HUBBLE RADIUS TILT: COSMOLOGICAL CONSTANT TEST")
print("=" * 70)

print(f"\nConstants:")
print(f"   Hubble radius R_H:         {R_H:.3e} m")
print(f"   Observable universe r:     {r_universe:.3e} m")
print(f"   R_H / r_universe:          {R_H/r_universe:.4f}")
print(f"   Λ measured (Planck 2018):  {LAMBDA_MEASURED:.4e} m⁻²")
print(f"   3/R_H² (base scale):       {3/R_H**2:.4e} m⁻²")

print(f"\nEHT-constrained n range: {N_LO} – {N_HI}")

print(f"\n{'='*70}")
print("PREDICTIONS ACROSS n RANGE")
print(f"{'='*70}")
print(f"\n{'n':<10} {'β_cosmic':<14} {'Λ_predicted':<18} {'Ratio Λp/Λm':<16} {'Status'}")
print("-" * 75)

best_n    = None
best_ratio = None
best_dist  = float('inf')

n_scan = np.linspace(0.01, 0.5, 1000)
ratios = []

for n in n_scan:
    r = lambda_ratio(n)
    ratios.append(r)
    dist = abs(np.log10(r)) if r > 0 else float('inf')
    if dist < best_dist:
        best_dist  = dist
        best_n     = n
        best_ratio = r

# Print key values
for n_test in [0.01, 0.05, 0.10, 0.215, 0.30, 0.50]:
    b  = beta_cosmic(n_test)
    lp = lambda_predicted(n_test)
    r  = lambda_ratio(n_test)

    if 0.1 <= r <= 10:
        status = "✅ Order of magnitude"
    elif 0.01 <= r <= 100:
        status = "⚠️  2 orders off"
    else:
        status = f"❌ {r:.1e}× off"

    in_eht = " ← EHT range" if N_LO <= n_test <= N_HI else ""
    print(f"{n_test:<10.3f} {b:<14.6f} {lp:<18.4e} {r:<16.4e} {status}{in_eht}")

print(f"\n{'='*70}")
print("BEST FIT")
print(f"{'='*70}")
print(f"\n   Best n:              {best_n:.4f}")
print(f"   β_cosmic at best n:  {beta_cosmic(best_n):.6f}")
print(f"   Λ_predicted:         {lambda_predicted(best_n):.4e} m⁻²")
print(f"   Λ_measured:          {LAMBDA_MEASURED:.4e} m⁻²")
print(f"   Ratio:               {best_ratio:.4e}")
print(f"   Log10(ratio):        {np.log10(best_ratio):.2f} orders of magnitude")

print(f"\n{'='*70}")
print("WITHIN EHT-CONSTRAINED RANGE (n = 0.05 – 0.215)")
print(f"{'='*70}")

ratio_lo   = lambda_ratio(N_LO)
ratio_hi   = lambda_ratio(N_HI)
ratio_best = lambda_ratio(N_BEST)

print(f"\n   n = {N_LO}:   Λ ratio = {ratio_lo:.4e}  ({np.log10(ratio_lo):.1f} orders from measured)")
print(f"   n = {N_HI}: Λ ratio = {ratio_hi:.4e}  ({np.log10(ratio_hi):.1f} orders from measured)")

print(f"\n{'='*70}")
print("VERDICT")
print(f"{'='*70}")

log_off = np.log10(ratio_lo)

if abs(log_off) < 1:
    print(f"\n   ✅ STRONG RESULT:")
    print(f"   Hubble radius substitution produces Λ within one order of magnitude")
    print(f"   of the measured value — using only EHT-constrained n.")
    print(f"   The same geometric parameter links black hole shadows and")
    print(f"   the expansion of the universe.")
elif abs(log_off) < 3:
    print(f"\n   ⚠️  PARTIAL RESULT:")
    print(f"   Hubble radius substitution is {abs(log_off):.1f} orders of magnitude off.")
    print(f"   The framework is directionally correct but needs a scaling factor.")
    print(f"   The bridge exists conceptually — the proportionality constant")
    print(f"   needs to be identified.")
else:
    print(f"\n   ❌ NEGATIVE RESULT:")
    print(f"   Hubble radius substitution is {abs(log_off):.1f} orders of magnitude off.")
    print(f"   Simple substitution is insufficient.")
    print(f"   A more fundamental derivation of the cosmological bridge is needed.")

print(f"\n   Note: Even a negative result is informative — it tells us how far")
print(f"   the simple substitution falls short and what scale of correction")
print(f"   is needed to bridge the gap.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
fig.suptitle('Hubble Radius Tilt: Cosmological Constant Test',
             fontsize=14, fontweight='bold')

n_plot  = np.linspace(0.01, 0.5, 500)
r_plot  = [lambda_ratio(n) for n in n_plot]
b_plot  = [beta_cosmic(n)  for n in n_plot]
lp_plot = [lambda_predicted(n) for n in n_plot]

eht_kw = dict(alpha=0.2, color='green', label='EHT-consistent range')

# -- Plot 1: Λ ratio vs n --
axes[0].semilogy(n_plot, r_plot, 'purple', linewidth=2.5)
axes[0].axhline(1.0,  color='blue',  linewidth=2, linestyle='--',
                label='Λ_predicted = Λ_measured (perfect)')
axes[0].axhline(10,   color='orange',linewidth=1.5, linestyle=':',
                label='10× off')
axes[0].axhline(0.1,  color='orange',linewidth=1.5, linestyle=':')
axes[0].axvspan(N_LO, N_HI, **eht_kw)
axes[0].axvline(best_n, color='gold', linewidth=2, linestyle='-',
                label=f'Best fit n={best_n:.3f}')
axes[0].set_xlabel('Tilt Exponent n', fontsize=12)
axes[0].set_ylabel('Λ_predicted / Λ_measured', fontsize=12)
axes[0].set_title('Λ Ratio vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)

# -- Plot 2: β_cosmic vs n --
axes[1].plot(n_plot, b_plot, 'purple', linewidth=2.5)
axes[1].axvspan(N_LO, N_HI, **eht_kw)
axes[1].axvline(best_n, color='gold', linewidth=2, linestyle='-',
                label=f'Best fit n={best_n:.3f}')
axes[1].set_xlabel('Tilt Exponent n', fontsize=12)
axes[1].set_ylabel('β_cosmic (suppression factor)', fontsize=12)
axes[1].set_title('Cosmic Suppression β vs n', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.3)
axes[1].set_ylim(0, 1.05)

# -- Plot 3: Λ_predicted vs Λ_measured bar --
labels = ['Λ measured\n(Planck 2018)',
          f'Λ predicted\n(n={N_LO})',
          f'Λ predicted\n(n={N_HI})',
          f'Λ predicted\n(best fit n={best_n:.3f})']
values = [LAMBDA_MEASURED,
          lambda_predicted(N_LO),
          lambda_predicted(N_HI),
          lambda_predicted(best_n)]
colors = ['royalblue', 'green', 'limegreen', 'gold']

bars = axes[2].bar(labels, values, color=colors, alpha=0.8,
                   edgecolor='black', linewidth=1.5)
axes[2].set_yscale('log')
axes[2].axhline(LAMBDA_MEASURED, color='royalblue', linewidth=2,
                linestyle='--', alpha=0.5, label='Measured Λ')
axes[2].set_ylabel('Cosmological Constant Λ (m⁻²)', fontsize=12)
axes[2].set_title('Predicted vs Measured Λ', fontsize=12, fontweight='bold')
axes[2].legend(fontsize=9)
axes[2].grid(axis='y', alpha=0.3)

for bar, val in zip(bars, values):
    axes[2].text(bar.get_x() + bar.get_width()/2,
                bar.get_height() * 1.5,
                f'{val:.1e}', ha='center', va='bottom',
                fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('hubble_tilt_cosmology.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: hubble_tilt_cosmology.png")
print("=" * 70)
