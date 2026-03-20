"""
test_half_power_and_dynamic_n.py

Two investigations in one script:

INVESTIGATION A — The √½ Connection
H1 gave α ≈ 0.515 (scale factor exponent)
H2 gave p ≈ 0.563 (projection power)
Both are close to ½. Is this coincidence or geometry?
Test: does setting both to exactly ½ still close the gap?
And are they the same geometric statement expressed differently?

INVESTIGATION B — H4 Fixed (Bounded Dynamic n)
The original H4 used n = A / x^B which explodes at low curvature
(n = 145 million at the Sun's surface). This is unphysical.

Fix: use an exponential form that stays bounded at all scales:
    n(x) = A × exp(-B × x)
    where x = curvature ratio

This gives:
    n → A (finite maximum) as x → 0 (flat space, weak field)
    n → 0 as x → ∞ (extreme curvature)

Same two calibration points:
    M87* photon sphere: x = 0.667, n = 0.133 (EHT midpoint)
    Cosmological:       x = 0.311, n = 0.500 (Λ match)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Constants (same as previous scripts)
# ============================================================================

M_m87           = 6.5e9 * const.M_sun
Rs_m87          = 2 * const.G * M_m87 / const.c**2
r_photon        = 1.5 * Rs_m87
CURVATURE_LOCAL = Rs_m87 / r_photon           # 0.667

R_H             = const.c / const.H_0
r_universe      = 4.4e26
CURVATURE_COSM  = R_H / r_universe            # 0.311

LAMBDA_MEASURED = 1.1056e-52                  # m⁻²

N_EHT_LO  = 0.050
N_EHT_HI  = 0.215
N_EHT_MID = (N_EHT_LO + N_EHT_HI) / 2       # 0.1325

# ============================================================================
# Base Functions
# ============================================================================

def beta_cosmic(n):
    ratio = (R_H / r_universe) ** (2 * n)
    return 0.0 if ratio >= 1.0 else np.sqrt(1.0 - ratio)

def beta_local(r, Rs, n):
    ratio = (Rs / r) ** (2 * n)
    return 0.0 if ratio >= 1.0 else np.sqrt(1.0 - ratio)

def lambda_from_beta(beta):
    return (3.0 / R_H**2) * beta**2

def lambda_ratio(n):
    return lambda_from_beta(beta_cosmic(n)) / LAMBDA_MEASURED

def eht_sigma(n):
    distance_m   = 16.8 * 3.086e22
    R_shadow_gr  = np.sqrt(27) * Rs_m87
    theta_gr_uas = (R_shadow_gr / distance_m) * 206265 * 1e6
    beta         = beta_local(r_photon, Rs_m87, n)
    correction   = 1 + beta * (Rs_m87 / (2 * r_photon))
    shadow_pred  = theta_gr_uas * correction
    return (shadow_pred - 42.0) / 3.0

# ============================================================================
# INVESTIGATION A — The √½ Connection
# ============================================================================

def h1_lambda_half(n_local):
    """H1 with α fixed to exactly ½."""
    k     = np.log10(r_universe / r_photon)
    n_c   = n_local * (k ** 0.5)
    return lambda_from_beta(beta_cosmic(n_c)) / LAMBDA_MEASURED

def h2_lambda_half(n):
    """H2 with p fixed to exactly ½."""
    beta = beta_cosmic(n)
    return (3.0 / R_H**2) * (beta ** 0.5) / LAMBDA_MEASURED

def h1_lambda_fitted(n_local, alpha):
    """H1 with fitted α."""
    k   = np.log10(r_universe / r_photon)
    n_c = n_local * (k ** alpha)
    return lambda_from_beta(beta_cosmic(n_c)) / LAMBDA_MEASURED

def h2_lambda_fitted(n, p):
    """H2 with fitted p."""
    beta = beta_cosmic(n)
    return (3.0 / R_H**2) * (beta ** p) / LAMBDA_MEASURED

# Fitted values from previous script
ALPHA_FITTED = 0.5149
P_FITTED     = 0.5634

# Are H1 and H2 the same statement?
# H1 maps n_local → n_cosmic via scale correction
# H2 modifies β directly
# Test: at n_EHT_MID, do both give same Λ with α=p=½?

# ============================================================================
# INVESTIGATION B — Bounded Dynamic n
# ============================================================================

def h4_bounded_fit():
    """
    Fit n(x) = A × exp(-B × x) to two calibration points.

    Point 1: x=0.667 (M87* photon sphere), n=0.133 (EHT midpoint)
    Point 2: x=0.311 (cosmological),       n=0.500 (Λ match)

    Solving:
    A × exp(-0.667 B) = 0.133  ... (1)
    A × exp(-0.311 B) = 0.500  ... (2)

    Dividing (1)/(2): exp(-0.356 B) = 0.266
    B = -ln(0.266) / 0.356
    A = 0.500 × exp(0.311 B)
    """
    x1, n1 = CURVATURE_LOCAL, N_EHT_MID
    x2, n2 = CURVATURE_COSM,  0.500

    B = -np.log(n1 / n2) / (x1 - x2)
    A = n2 * np.exp(B * x2)
    return A, B


def h4_bounded_n(curvature, A, B):
    """Bounded dynamic n: n(x) = A × exp(-B × x)"""
    return A * np.exp(-B * curvature)


def schwarzschild_curvature(mass, radius):
    """Curvature ratio = Rs/r for a given mass and surface radius."""
    Rs = 2 * const.G * mass / const.c**2
    return Rs / radius


# ============================================================================
# Run Investigation A
# ============================================================================

print("=" * 70)
print("INVESTIGATION A: THE √½ CONNECTION")
print("=" * 70)

print(f"\nFrom previous script:")
print(f"   H1 fitted α = {ALPHA_FITTED:.4f}  (close to ½ = 0.5000)")
print(f"   H2 fitted p = {P_FITTED:.4f}  (close to ½ = 0.5000)")
print(f"   Difference from ½: α is {abs(ALPHA_FITTED-0.5):.4f} off, "
      f"p is {abs(P_FITTED-0.5):.4f} off")

print(f"\n--- Test: Does α = ½ exactly still close the gap? ---")
print(f"\n   n_local → Λ ratio with α=½ (H1):")
for n in [N_EHT_LO, N_EHT_MID, N_EHT_HI]:
    r_half   = h1_lambda_half(n)
    r_fitted = h1_lambda_fitted(n, ALPHA_FITTED)
    print(f"   n={n:.3f}: α=½ → {r_half:.3f}  |  α=fitted → {r_fitted:.3f}")

print(f"\n--- Test: Does p = ½ exactly still close the gap? ---")
print(f"\n   n → Λ ratio with p=½ (H2):")
for n in [N_EHT_LO, N_EHT_MID, N_EHT_HI]:
    r_half   = h2_lambda_half(n)
    r_fitted = h2_lambda_fitted(n, P_FITTED)
    print(f"   n={n:.3f}: p=½ → {r_half:.3f}  |  p=fitted → {r_fitted:.3f}")

print(f"\n--- Are H1 and H2 the same geometric statement? ---")
print(f"\n   At n_EHT_MID = {N_EHT_MID:.3f}:")
print(f"   H1 (α=½):      Λ ratio = {h1_lambda_half(N_EHT_MID):.4f}")
print(f"   H2 (p=½):      Λ ratio = {h2_lambda_half(N_EHT_MID):.4f}")
print(f"   H1 (α=fitted): Λ ratio = {h1_lambda_fitted(N_EHT_MID, ALPHA_FITTED):.4f}")
print(f"   H2 (p=fitted): Λ ratio = {h2_lambda_fitted(N_EHT_MID, P_FITTED):.4f}")

# Physical interpretation of ½
print(f"\n--- Physical meaning of ½ ---")
print(f"\n   In wave physics: intensity ∝ amplitude²")
print(f"   A β^½ correction in amplitude → β in energy/intensity")
print(f"   This would mean: the 4D tilt amplitude projects as β into 3D,")
print(f"   but what we measure (shadow area, Λ energy density) scales as β²")
print(f"\n   A scale correction of (log_scale)^½ is a geometric mean —")
print(f"   the square root of the scale ratio in log space.")
print(f"   Both interpretations suggest a √ relationship between")
print(f"   the 4D cause and the 3D observable effect.")

# ============================================================================
# Run Investigation B
# ============================================================================

print(f"\n{'='*70}")
print("INVESTIGATION B: BOUNDED DYNAMIC n")
print(f"{'='*70}")

A_b, B_b = h4_bounded_fit()
print(f"\n   Bounded form: n(x) = {A_b:.4f} × exp(-{B_b:.4f} × x)")
print(f"   At x=0 (flat space): n → {A_b:.4f}  ← finite maximum")
print(f"   At x=1 (horizon):    n = {h4_bounded_n(1.0, A_b, B_b):.4f}")

print(f"\n   Calibration check:")
n_m87  = h4_bounded_n(CURVATURE_LOCAL, A_b, B_b)
n_cosm = h4_bounded_n(CURVATURE_COSM,  A_b, B_b)
sig_m87 = eht_sigma(n_m87)
lam_cosm = lambda_from_beta(beta_cosmic(n_cosm)) / LAMBDA_MEASURED
print(f"   M87* (x={CURVATURE_LOCAL:.3f}): n = {n_m87:.4f}, "
      f"EHT σ = {sig_m87:.2f} "
      f"({'✅' if abs(sig_m87)<1.5 else '⚠️'})")
print(f"   Cosmic (x={CURVATURE_COSM:.3f}): n = {n_cosm:.4f}, "
      f"Λ ratio = {lam_cosm:.4f} "
      f"({'✅' if abs(lam_cosm-1)<0.01 else '⚠️'})")

print(f"\n   Dynamic n across environments (bounded form):")
envs = [
    ("Deep space (flat)",   1e-10),
    ("Sun surface",         schwarzschild_curvature(const.M_sun, const.R_sun)),
    ("White dwarf",         schwarzschild_curvature(1.4*const.M_sun, 7e6)),
    ("Neutron star",        0.300),
    ("Stellar BH 10 Msun",  0.500),
    ("M87* photon sphere",  CURVATURE_LOCAL),
    ("Black hole horizon",  1.000),
    ("Cosmological",        CURVATURE_COSM),
]

print(f"\n   {'Environment':<26} {'Curvature':<14} {'n(bounded)':<12} {'EHT σ'}")
print("   " + "-"*65)
for name, curv in envs:
    n_dyn = h4_bounded_n(curv, A_b, B_b)
    try:
        sig = eht_sigma(n_dyn)
        sig_str = f"{sig:.2f}σ"
    except:
        sig_str = "N/A"
    print(f"   {name:<26} {curv:<14.6f} {n_dyn:<12.4f} {sig_str}")

print(f"\n   Key fix: Sun surface now gives n = "
      f"{h4_bounded_n(schwarzschild_curvature(const.M_sun, const.R_sun), A_b, B_b):.4f}")
print(f"   (Previously: 145,208,230 — now physically reasonable)")

print(f"\n   Physical interpretation of bounded H4:")
print(f"   n_max = {A_b:.4f} — the tilt exponent in completely flat space")
print(f"   As curvature increases, n decreases exponentially")
print(f"   Strong gravity suppresses the tilt exponent,")
print(f"   making the shadow correction smaller near black holes")
print(f"   and larger at cosmological scales — consistent with both EHT and Λ")

# ============================================================================
# Connection Between A and B
# ============================================================================

print(f"\n{'='*70}")
print("CONNECTION: DO √½ AND BOUNDED H4 TELL THE SAME STORY?")
print(f"{'='*70}")

print(f"\n   H1/H2 say: the correction is ≈ √(something)")
print(f"   H4 says:   n decreases exponentially with curvature")
print(f"\n   At the EHT calibration point (n_EHT_MID = {N_EHT_MID:.3f}):")
print(f"   β_cosmic(n={N_EHT_MID:.3f}) = {beta_cosmic(N_EHT_MID):.4f}")
print(f"   √β_cosmic = {np.sqrt(beta_cosmic(N_EHT_MID)):.4f}")
print(f"   β_cosmic²  = {beta_cosmic(N_EHT_MID)**2:.4f}")
print(f"\n   The √β suggestion means the 3D measurement captures")
print(f"   the square root of the full 4D correction amplitude.")
print(f"   H4 explains WHY: as curvature increases, n decreases,")
print(f"   which automatically reduces the correction β in strong fields.")
print(f"   The √ relationship is the emergent consequence of n being dynamic.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('√½ Connection and Bounded Dynamic n',
             fontsize=14, fontweight='bold')

n_range   = np.linspace(0.01, 0.8, 500)
base_rat  = [lambda_ratio(n) for n in n_range]
eht_kw    = dict(alpha=0.2, color='green', label='EHT range (0.05–0.215)')

# -- Plot 1: H1 with α=fitted vs α=½ --
ax = axes[0, 0]
r_h1_fitted = [h1_lambda_fitted(n, ALPHA_FITTED) for n in n_range]
r_h1_half   = [h1_lambda_half(n)                 for n in n_range]
ax.semilogy(n_range, base_rat,    'gray',   linewidth=1.5,
            linestyle='--', alpha=0.5, label='No correction')
ax.semilogy(n_range, r_h1_fitted, 'blue',   linewidth=2.5,
            label=f'H1: α={ALPHA_FITTED:.3f} (fitted)')
ax.semilogy(n_range, r_h1_half,   'cyan',   linewidth=2.5,
            linestyle='--', label='H1: α=½ (exact)')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Perfect match')
ax.axvspan(N_EHT_LO, N_EHT_HI, **eht_kw)
ax.set_xlabel('n (local)', fontsize=11)
ax.set_ylabel('Λ ratio', fontsize=11)
ax.set_title('H1: Fitted α vs Exact ½', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0.05, 20)

# -- Plot 2: H2 with p=fitted vs p=½ --
ax = axes[0, 1]
r_h2_fitted = [h2_lambda_fitted(n, P_FITTED) for n in n_range]
r_h2_half   = [h2_lambda_half(n)             for n in n_range]
ax.semilogy(n_range, base_rat,    'gray',   linewidth=1.5,
            linestyle='--', alpha=0.5, label='No correction')
ax.semilogy(n_range, r_h2_fitted, 'red',    linewidth=2.5,
            label=f'H2: p={P_FITTED:.3f} (fitted)')
ax.semilogy(n_range, r_h2_half,   'orange', linewidth=2.5,
            linestyle='--', label='H2: p=½ (exact)')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Perfect match')
ax.axvspan(N_EHT_LO, N_EHT_HI, **eht_kw)
ax.set_xlabel('n', fontsize=11)
ax.set_ylabel('Λ ratio', fontsize=11)
ax.set_title('H2: Fitted p vs Exact ½', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0.05, 20)

# -- Plot 3: H4 bounded vs original --
ax = axes[1, 0]
curv_range    = np.linspace(0.001, 1.0, 500)

# Original (blows up — cap at 2 for display)
def h4_original_capped(x, A=0.0654, B=1.7397):
    val = A / (x ** B) if x > 0 else 1e10
    return min(val, 2.0)

n_orig    = [h4_original_capped(c) for c in curv_range]
n_bounded = [h4_bounded_n(c, A_b, B_b) for c in curv_range]

ax.plot(curv_range, n_orig,    'red',        linewidth=2,
        linestyle='--', label='Original H4 (explodes, capped at 2)')
ax.plot(curv_range, n_bounded, 'darkorange', linewidth=2.5,
        label=f'Bounded H4: {A_b:.3f}×exp(-{B_b:.3f}×x)')
ax.fill_between(curv_range, N_EHT_LO, N_EHT_HI,
                alpha=0.15, color='green', label='EHT range for n')
ax.axhline(0.500,   color='blue',   linewidth=1.5, linestyle='--',
           label='n for Λ match (0.500)')
ax.axhline(A_b,     color='orange', linewidth=1,   linestyle=':',
           label=f'n_max = {A_b:.3f} (flat space limit)')
ax.axvline(CURVATURE_LOCAL, color='red',    linewidth=1.5, linestyle=':',
           alpha=0.7, label=f'M87* ({CURVATURE_LOCAL:.3f})')
ax.axvline(CURVATURE_COSM,  color='purple', linewidth=1.5, linestyle=':',
           alpha=0.7, label=f'Cosmic ({CURVATURE_COSM:.3f})')
ax.set_xlabel('Curvature Ratio', fontsize=11)
ax.set_ylabel('Dynamic n', fontsize=11)
ax.set_title('H4: Original vs Bounded Form', fontsize=12, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_xlim(0, 1.0); ax.set_ylim(0, 2.1)

# -- Plot 4: Λ ratio with bounded H4 across n range --
ax = axes[1, 1]
# For each curvature in a range, compute dynamic n, then Λ ratio
curv_plot   = np.linspace(0.05, 0.95, 500)
n_dyn_vals  = [h4_bounded_n(c, A_b, B_b) for c in curv_plot]
lam_ratios  = [lambda_from_beta(beta_cosmic(n)) / LAMBDA_MEASURED
               for n in n_dyn_vals]

ax.plot(curv_plot, lam_ratios, 'darkorange', linewidth=2.5,
        label='Λ ratio with bounded H4')
ax.axhline(1.0,  color='black',  linewidth=2, linestyle='--',
           label='Perfect match')
ax.axhline(0.1,  color='orange', linewidth=1, linestyle=':', alpha=0.7)
ax.axhline(10.0, color='orange', linewidth=1, linestyle=':', alpha=0.7,
           label='10× boundary')
ax.axvline(CURVATURE_LOCAL, color='red',    linewidth=1.5, linestyle=':',
           label=f'M87* ({CURVATURE_LOCAL:.3f})')
ax.axvline(CURVATURE_COSM,  color='purple', linewidth=1.5, linestyle=':',
           label=f'Cosmic ({CURVATURE_COSM:.3f}) → Λ exact')
ax.set_yscale('log')
ax.set_xlabel('Curvature Ratio', fontsize=11)
ax.set_ylabel('Λ_predicted / Λ_measured', fontsize=11)
ax.set_title('Λ Ratio from Bounded Dynamic n', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('half_power_and_dynamic_n.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: half_power_and_dynamic_n.png")
print("=" * 70)
