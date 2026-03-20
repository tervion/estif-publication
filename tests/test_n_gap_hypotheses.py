"""
test_n_gap_hypotheses.py

Tests four hypotheses for why the tilt exponent n differs between:
  - EHT M87* shadow measurement:    n = 0.05 – 0.215
  - Cosmological constant Λ match:  n ≈ 0.500

HYPOTHESES:
  H1 — Scale Factor:      A multiplicative correction based on scale ratio
  H2 — Projection Channel: EM vs gravitational projection differ by a power
  H3 — Dimensional:        1D path (lensing) vs 3D volume (Λ) introduces β³
  H4 — Dynamic n:          n is not constant but determined by local curvature
                           (the cannonball insight — g changes with environment)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Observational Anchors
# ============================================================================

M_m87           = 6.5e9 * const.M_sun
Rs_m87          = 2 * const.G * M_m87 / const.c**2
r_photon        = 1.5 * Rs_m87
CURVATURE_LOCAL = Rs_m87 / r_photon       # 0.667 — strong field

R_H             = const.c / const.H_0
r_universe      = 4.4e26
CURVATURE_COSM  = R_H / r_universe        # 0.311 — weak field

LAMBDA_MEASURED = 1.1056e-52              # m⁻² Planck 2018

N_EHT_LO  = 0.050
N_EHT_HI  = 0.215
N_EHT_MID = (N_EHT_LO + N_EHT_HI) / 2   # 0.1325

# ============================================================================
# Base Functions
# ============================================================================

def beta_local(r, Rs, n):
    ratio = (Rs / r) ** (2 * n)
    return 0.0 if ratio >= 1.0 else np.sqrt(1.0 - ratio)

def beta_cosmic(n):
    ratio = (R_H / r_universe) ** (2 * n)
    return 0.0 if ratio >= 1.0 else np.sqrt(1.0 - ratio)

def lambda_from_beta(beta):
    return (3.0 / R_H**2) * beta**2

def eht_shadow_sigma(n):
    distance_m   = 16.8 * 3.086e22
    R_shadow_gr  = np.sqrt(27) * Rs_m87
    theta_gr_uas = (R_shadow_gr / distance_m) * 206265 * 1e6
    beta         = beta_local(r_photon, Rs_m87, n)
    correction   = 1 + beta * (Rs_m87 / (2 * r_photon))
    shadow_pred  = theta_gr_uas * correction
    return (shadow_pred - 42.0) / 3.0

# ============================================================================
# H1: Scale Factor
# ============================================================================

def h1_effective_n_cosm(n_local, alpha):
    k = np.log10(r_universe / r_photon)
    return n_local * (k ** alpha)

def h1_find_alpha(n_local=N_EHT_MID, n_target=0.500):
    k = np.log10(r_universe / r_photon)
    return np.log(n_target / n_local) / np.log(k)

# ============================================================================
# H2: Projection Channel
# ============================================================================

def h2_lambda_with_power(n, p):
    beta = beta_cosmic(n)
    return (3.0 / R_H**2) * (beta ** p)

def h2_find_power(n=N_EHT_MID):
    beta = beta_cosmic(n)
    if beta <= 0:
        return None
    target = LAMBDA_MEASURED * R_H**2 / 3.0
    return np.log(target) / np.log(beta)

# ============================================================================
# H3: Dimensional Correction
# ============================================================================

def h3_lambda_cubic(n):
    beta = beta_cosmic(n)
    return (3.0 / R_H**2) * (beta ** 3)

def h3_find_n():
    n_scan  = np.linspace(0.01, 2.0, 10000)
    diffs   = [abs(h3_lambda_cubic(n) / LAMBDA_MEASURED - 1) for n in n_scan]
    best    = np.argmin(diffs)
    return n_scan[best], h3_lambda_cubic(n_scan[best])

# ============================================================================
# H4: Dynamic n
# ============================================================================

def h4_dynamic_n(curvature_ratio, A, B):
    return A / (curvature_ratio ** B)

def h4_fit():
    x1, n1 = CURVATURE_LOCAL, N_EHT_MID
    x2, n2 = CURVATURE_COSM,  0.500
    B = -np.log(n1 / n2) / (np.log(x1) - np.log(x2))
    A = n1 * (x1 ** B)
    return A, B

# ============================================================================
# Run All Hypotheses
# ============================================================================

print("=" * 70)
print("N-GAP HYPOTHESES: WHY DOES n_EHT ≠ n_Λ?")
print("=" * 70)
print(f"\n   EHT constrains:       n = {N_EHT_LO} – {N_EHT_HI}")
print(f"   Λ perfect match at:   n ≈ 0.500")
print(f"   EHT curvature ratio:  {CURVATURE_LOCAL:.4f}")
print(f"   Cosmic curvature:     {CURVATURE_COSM:.4f}")

# H1
print(f"\n{'='*70}")
print("H1 — SCALE FACTOR")
print(f"{'='*70}")
alpha  = h1_find_alpha(N_EHT_MID, 0.500)
print(f"   Scale ratio (log10): {np.log10(r_universe/r_photon):.1f} orders")
print(f"   Required alpha:      {alpha:.4f}")
print(f"\n   n_local → n_cosmic → Λ ratio:")
h1_results = []
for n_test in [N_EHT_LO, N_EHT_MID, N_EHT_HI]:
    n_c   = h1_effective_n_cosm(n_test, alpha)
    lam   = lambda_from_beta(beta_cosmic(n_c))
    ratio = lam / LAMBDA_MEASURED
    h1_results.append(ratio)
    print(f"   {n_test:.3f} → {n_c:.3f} → {ratio:.3f}")

# H2
print(f"\n{'='*70}")
print("H2 — PROJECTION CHANNEL (EM vs Gravitational)")
print(f"{'='*70}")
p = h2_find_power(N_EHT_MID)
print(f"   Required power p (at n={N_EHT_MID:.3f}): {p:.4f}")
print(f"\n   n → β^p → Λ ratio:")
h2_results = []
for n_test in [N_EHT_LO, N_EHT_MID, N_EHT_HI]:
    lam   = h2_lambda_with_power(n_test, p)
    ratio = lam / LAMBDA_MEASURED
    h2_results.append(ratio)
    print(f"   {n_test:.3f} → β^{p:.2f} → {ratio:.3f}")

# H3
print(f"\n{'='*70}")
print("H3 — DIMENSIONAL CORRECTION (β³)")
print(f"{'='*70}")
n_h3, lam_h3 = h3_find_n()
print(f"   Λ = (3/R_H²)×β³ — perfect match at n = {n_h3:.4f}")
inside = N_EHT_LO <= n_h3 <= N_EHT_HI
print(f"   {'✅ Inside' if inside else '❌ Outside'} EHT range ({N_EHT_LO}–{N_EHT_HI})")
print(f"\n   n → β³ → Λ ratio:")
h3_results = []
for n_test in [N_EHT_LO, N_EHT_MID, N_EHT_HI]:
    ratio = h3_lambda_cubic(n_test) / LAMBDA_MEASURED
    h3_results.append(ratio)
    print(f"   {n_test:.3f} → {ratio:.4f}")

# H4
print(f"\n{'='*70}")
print("H4 — DYNAMIC n (The Cannonball Insight)")
print(f"{'='*70}")
A, B = h4_fit()
print(f"   Fitted: n(x) = {A:.4f} / x^{B:.4f}")
print(f"   where x = curvature ratio at measurement scale")
n_local_pred = h4_dynamic_n(CURVATURE_LOCAL, A, B)
n_cosm_pred  = h4_dynamic_n(CURVATURE_COSM,  A, B)
eht_sigma    = eht_shadow_sigma(n_local_pred)
lam_h4       = lambda_from_beta(beta_cosmic(n_cosm_pred))
ratio_h4     = lam_h4 / LAMBDA_MEASURED
print(f"\n   At M87* (curvature={CURVATURE_LOCAL:.3f}):  n = {n_local_pred:.4f}")
print(f"   EHT shadow tension: {eht_sigma:.2f}σ "
      f"({'✅ consistent' if abs(eht_sigma)<1.5 else '⚠️  tension'})")
print(f"\n   At cosmic (curvature={CURVATURE_COSM:.3f}): n = {n_cosm_pred:.4f}")
print(f"   Λ ratio: {ratio_h4:.4f} "
      f"({'✅ exact' if abs(ratio_h4-1)<0.01 else f'{ratio_h4:.3f}×'})")
print(f"\n   Dynamic n across environments:")
envs = [
    ("Sun surface",        2*const.G*const.M_sun/const.c**2 / const.R_sun),
    ("Neutron star",       0.300),
    ("M87* photon sphere", CURVATURE_LOCAL),
    ("Stellar BH 10 Msun", 0.500),
    ("Cosmological",       CURVATURE_COSM),
]
print(f"   {'Environment':<26} {'Curvature':<12} {'n(dynamic)'}")
print("   " + "-"*50)
for name, curv in envs:
    n_dyn = h4_dynamic_n(curv, A, B)
    print(f"   {name:<26} {curv:<12.4f} {n_dyn:.4f}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

h1_ok = abs(np.mean(h1_results) - 1.0) < 0.5
h2_ok = abs(np.mean(h2_results) - 1.0) < 0.5
h3_ok = inside
h4_ok = abs(eht_sigma) < 1.5 and abs(ratio_h4 - 1) < 0.05

print(f"\n{'Hypothesis':<32} {'Closes gap?':<14} {'Key parameter'}")
print("-"*70)
print(f"{'H1: Scale factor':<32} {'✅ Yes' if h1_ok else '❌ No':<14} α={alpha:.3f}")
print(f"{'H2: Projection channel':<32} {'✅ Yes' if h2_ok else '❌ No':<14} p={p:.3f}")
print(f"{'H3: Dimensional (β³)':<32} {'✅ Yes' if h3_ok else '❌ No':<14} match at n={n_h3:.3f}")
print(f"{'H4: Dynamic n (cannonball)':<32} {'✅ Yes' if h4_ok else '⚠️  Partial':<14} A={A:.3f}, B={B:.3f}")

print(f"\n   H4 (dynamic n) is the most physically natural:")
print(f"   The tilt formula sin(θ)=(Rs/r)^n becomes universal —")
print(f"   n adjusts to the local gravitational environment,")
print(f"   just as g adjusts to the local planet in Newton's formula.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('N-Gap Hypotheses: Why Does n_EHT ≠ n_Λ?',
             fontsize=14, fontweight='bold')

n_range    = np.linspace(0.01, 0.8, 500)
base_ratios= [lambda_from_beta(beta_cosmic(n))/LAMBDA_MEASURED for n in n_range]
eht_kw     = dict(alpha=0.2, color='green', label='EHT range')

# H1
ax = axes[0, 0]
h1_ratios = []
for n in n_range:
    nc  = h1_effective_n_cosm(n, alpha)
    lam = lambda_from_beta(beta_cosmic(nc))
    h1_ratios.append(lam / LAMBDA_MEASURED)
ax.semilogy(n_range, base_ratios, 'gray', linewidth=1.5,
            linestyle='--', alpha=0.6, label='No correction')
ax.semilogy(n_range, h1_ratios, 'blue', linewidth=2.5, label='H1 corrected')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Perfect match')
ax.axvspan(N_EHT_LO, N_EHT_HI, **eht_kw)
ax.set_xlabel('n (local)', fontsize=11)
ax.set_ylabel('Λ ratio', fontsize=11)
ax.set_title(f'H1: Scale Factor (α={alpha:.3f})', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0.01, 100)

# H2
ax = axes[0, 1]
h2_ratios = [h2_lambda_with_power(n, p)/LAMBDA_MEASURED for n in n_range]
ax.semilogy(n_range, base_ratios, 'gray', linewidth=1.5,
            linestyle='--', alpha=0.6, label='No correction (β¹)')
ax.semilogy(n_range, h2_ratios, 'red', linewidth=2.5, label=f'H2: β^{p:.2f}')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Perfect match')
ax.axvspan(N_EHT_LO, N_EHT_HI, **eht_kw)
ax.set_xlabel('n', fontsize=11)
ax.set_ylabel('Λ ratio', fontsize=11)
ax.set_title(f'H2: Projection Channel (β^{p:.2f})', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0.01, 100)

# H3
ax = axes[1, 0]
h3_ratios = [h3_lambda_cubic(n)/LAMBDA_MEASURED for n in n_range]
ax.semilogy(n_range, base_ratios, 'gray', linewidth=1.5,
            linestyle='--', alpha=0.6, label='No correction (β¹)')
ax.semilogy(n_range, h3_ratios, 'purple', linewidth=2.5, label='H3: β³')
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Perfect match')
ax.axvline(n_h3, color='purple', linewidth=1.5, linestyle=':',
           label=f'H3 match n={n_h3:.3f}')
ax.axvspan(N_EHT_LO, N_EHT_HI, **eht_kw)
ax.set_xlabel('n', fontsize=11)
ax.set_ylabel('Λ ratio', fontsize=11)
ax.set_title('H3: Dimensional Correction (β³)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0.01, 100)

# H4
ax = axes[1, 1]
curv_range = np.linspace(0.05, 1.0, 500)
n_dyn_vals = [h4_dynamic_n(c, A, B) for c in curv_range]
ax.plot(curv_range, n_dyn_vals, 'darkorange', linewidth=2.5,
        label=f'n(x)={A:.3f}/x^{B:.3f}')
ax.fill_between(curv_range, N_EHT_LO, N_EHT_HI,
                alpha=0.15, color='green', label='EHT range')
ax.axhline(0.500, color='blue', linewidth=1.5, linestyle='--',
           label='n for Λ match')
ax.axvline(CURVATURE_LOCAL, color='red', linewidth=1.5, linestyle=':',
           label=f'M87* ({CURVATURE_LOCAL:.3f})')
ax.axvline(CURVATURE_COSM, color='purple', linewidth=1.5, linestyle=':',
           label=f'Cosmic ({CURVATURE_COSM:.3f})')
ax.set_xlabel('Curvature Ratio', fontsize=11)
ax.set_ylabel('Dynamic n', fontsize=11)
ax.set_title('H4: Dynamic n — Cannonball Insight', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_xlim(0.05, 1.0); ax.set_ylim(0, 1.0)

plt.tight_layout()
plt.savefig('n_gap_hypotheses.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: n_gap_hypotheses.png")
print("=" * 70)
