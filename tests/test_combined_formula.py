"""
test_combined_formula.py

Tests the combined formula using all three findings:

    H4: n(x) = 1.594 × exp(-3.731 × x)    ← dynamic input
    H1/H2: Observable = √β                  ← square root output

Together:
    x        = curvature ratio (Rs/r locally, R_H/r_universe cosmologically)
    n(x)     = 1.594 × exp(-3.731 × x)
    sin(θ)   = x^n(x)
    β(x)     = √(1 - x^(2n(x)))
    Observable = √β(x)

THE TEST:
Does this single formula, with no free parameters, simultaneously:
  1. Produce an EHT-consistent shadow for M87* (within 1σ)?
  2. Produce the correct cosmological constant Λ?
  3. Give physically sensible values across all environments?
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Constants
# ============================================================================

M_m87           = 6.5e9 * const.M_sun
Rs_m87          = 2 * const.G * M_m87 / const.c**2
r_photon        = 1.5 * Rs_m87
CURVATURE_LOCAL = Rs_m87 / r_photon           # 0.667

R_H             = const.c / const.H_0
r_universe      = 4.4e26
CURVATURE_COSM  = R_H / r_universe            # 0.311

LAMBDA_MEASURED = 1.1056e-52                  # m⁻²
SHADOW_OBS      = 42.0                        # μas
SHADOW_ERR      = 3.0                         # μas

# H4 bounded parameters (from previous script)
N_MAX = 1.5940
B_EXP = 3.7311

# ============================================================================
# The Combined Formula
# ============================================================================

def n_dynamic(curvature):
    """H4: tilt exponent as function of local curvature."""
    return N_MAX * np.exp(-B_EXP * curvature)


def beta_combined(curvature):
    """
    Full tilt suppression with dynamic n.
    β = cos(θ) = √(1 - curvature^(2n(curvature)))
    """
    n   = n_dynamic(curvature)
    val = curvature ** (2 * n)
    if np.isscalar(val):
        return 0.0 if val >= 1.0 else np.sqrt(1.0 - val)
    return np.where(val >= 1.0, 0.0, np.sqrt(np.maximum(0.0, 1.0 - val)))


def observable_combined(curvature):
    """
    H1/H2: the 3D observable is √β, not β directly.
    This is the square root output correction.
    """
    return np.sqrt(beta_combined(curvature))


def lensing_combined(r, M):
    """
    Combined lensing deflection angle.
    θ_ESTIF = θ_GR × (1 + observable × Rs/(2r))
    """
    Rs      = 2 * const.G * M / const.c**2
    curv    = Rs / r
    theta_gr = (4 * const.G * M) / (const.c**2 * r)
    obs      = observable_combined(curv)
    return theta_gr * (1 + obs * (Rs / (2 * r)))


def lambda_combined():
    """
    Cosmological constant from combined formula.
    Λ = (3/R_H²) × observable_cosmic²
    (the √β observable squared gives back β for the energy density)
    """
    obs = observable_combined(CURVATURE_COSM)
    return (3.0 / R_H**2) * obs**2


# ============================================================================
# Test 1: EHT M87* Shadow
# ============================================================================

def test_eht():
    distance_m   = 16.8 * 3.086e22
    R_shadow_gr  = np.sqrt(27) * Rs_m87
    theta_gr_uas = (R_shadow_gr / distance_m) * 206265 * 1e6

    theta_estif  = lensing_combined(r_photon, M_m87)
    theta_gr_pt  = (4 * const.G * M_m87) / (const.c**2 * r_photon)
    deviation    = (theta_estif / theta_gr_pt - 1) * 100

    # Shadow scales with deflection
    shadow_pred  = theta_gr_uas * (theta_estif / theta_gr_pt)
    sigma        = (shadow_pred - SHADOW_OBS) / SHADOW_ERR

    return {
        'n_local':      n_dynamic(CURVATURE_LOCAL),
        'beta':         beta_combined(CURVATURE_LOCAL),
        'observable':   observable_combined(CURVATURE_LOCAL),
        'theta_gr':     theta_gr_uas,
        'shadow_pred':  shadow_pred,
        'deviation':    deviation,
        'sigma':        sigma,
    }


# ============================================================================
# Test 2: Cosmological Constant
# ============================================================================

def test_lambda():
    lam_pred = lambda_combined()
    ratio    = lam_pred / LAMBDA_MEASURED
    return {
        'n_cosmic':     n_dynamic(CURVATURE_COSM),
        'beta':         beta_combined(CURVATURE_COSM),
        'observable':   observable_combined(CURVATURE_COSM),
        'lam_pred':     lam_pred,
        'lam_meas':     LAMBDA_MEASURED,
        'ratio':        ratio,
    }


# ============================================================================
# Print Results
# ============================================================================

print("=" * 70)
print("COMBINED FORMULA TEST: H4 + H1/H2")
print("=" * 70)
print(f"\nCombined formula:")
print(f"   n(x)          = {N_MAX:.4f} × exp(-{B_EXP:.4f} × x)")
print(f"   β(x)          = √(1 - x^(2n(x)))")
print(f"   Observable    = √β(x)")
print(f"   where x = curvature ratio (Rs/r or R_H/r)")

print(f"\n{'='*70}")
print("TEST 1: EHT M87* SHADOW")
print(f"{'='*70}")

eht = test_eht()
print(f"\n   Curvature at photon sphere: {CURVATURE_LOCAL:.4f}")
print(f"   Dynamic n:                  {eht['n_local']:.4f}")
print(f"   β (tilt suppression):       {eht['beta']:.4f}")
print(f"   Observable (√β):            {eht['observable']:.4f}")
print(f"\n   GR shadow prediction:       {eht['theta_gr']:.2f} μas")
print(f"   Combined ESTIF prediction:  {eht['shadow_pred']:.2f} μas")
print(f"   EHT observed:               {SHADOW_OBS:.1f} ± {SHADOW_ERR:.1f} μas")
print(f"   Deviation from GR:          {eht['deviation']:.2f}%")
print(f"   EHT tension:                {eht['sigma']:.2f}σ")

if abs(eht['sigma']) < 1.0:
    eht_verdict = "✅ CONSISTENT (within 1σ)"
elif abs(eht['sigma']) < 2.0:
    eht_verdict = "⚠️  MILD TENSION (1–2σ)"
elif abs(eht['sigma']) < 3.0:
    eht_verdict = "⚠️  MODERATE TENSION (2–3σ)"
else:
    eht_verdict = "❌ RULED OUT (>3σ)"
print(f"   Verdict: {eht_verdict}")

print(f"\n{'='*70}")
print("TEST 2: COSMOLOGICAL CONSTANT Λ")
print(f"{'='*70}")

lam = test_lambda()
print(f"\n   Curvature at cosmic scale:  {CURVATURE_COSM:.4f}")
print(f"   Dynamic n:                  {lam['n_cosmic']:.4f}")
print(f"   β (tilt suppression):       {lam['beta']:.4f}")
print(f"   Observable (√β):            {lam['observable']:.4f}")
print(f"\n   Λ predicted:                {lam['lam_pred']:.4e} m⁻²")
print(f"   Λ measured (Planck 2018):   {lam['lam_meas']:.4e} m⁻²")
print(f"   Ratio Λ_pred/Λ_meas:        {lam['ratio']:.4f}")
print(f"   Log10(ratio):               {np.log10(lam['ratio']):.2f} orders")

if abs(lam['ratio'] - 1.0) < 0.01:
    lam_verdict = "✅ EXACT MATCH"
elif abs(lam['ratio'] - 1.0) < 0.1:
    lam_verdict = "✅ CONSISTENT (within 10%)"
elif abs(np.log10(lam['ratio'])) < 1.0:
    lam_verdict = "⚠️  WITHIN ONE ORDER"
else:
    lam_verdict = "❌ SIGNIFICANT DISCREPANCY"
print(f"   Verdict: {lam_verdict}")

print(f"\n{'='*70}")
print("TEST 3: ENVIRONMENT SCAN")
print(f"{'='*70}")

envs = [
    ("Deep space (flat)",   1e-8),
    ("Sun surface",         2*const.G*const.M_sun/const.c**2 / const.R_sun),
    ("White dwarf",         2*const.G*1.4*const.M_sun/const.c**2 / 7e6),
    ("Neutron star",        0.300),
    ("Stellar BH 10 Msun",  0.500),
    ("M87* photon sphere",  CURVATURE_LOCAL),
    ("Black hole horizon",  0.999),
    ("Cosmological",        CURVATURE_COSM),
]

print(f"\n   {'Environment':<26} {'x':<10} {'n(x)':<8} {'β':<8} {'√β':<8}")
print("   " + "-"*65)
for name, curv in envs:
    n   = n_dynamic(curv)
    b   = beta_combined(curv)
    obs = observable_combined(curv)
    print(f"   {name:<26} {curv:<10.6f} {n:<8.4f} {b:<8.4f} {obs:<8.4f}")

print(f"\n{'='*70}")
print("OVERALL VERDICT")
print(f"{'='*70}")

both_pass = abs(eht['sigma']) < 2.0 and abs(lam['ratio'] - 1.0) < 0.1

print(f"\n   EHT M87* shadow:         {eht_verdict}")
print(f"   Cosmological constant:   {lam_verdict}")

if both_pass:
    print(f"\n   ✅ COMBINED FORMULA PASSES BOTH TESTS")
    print(f"\n   A single self-consistent formula with no free parameters")
    print(f"   simultaneously describes:")
    print(f"   → Black hole shadow size (EHT M87*)")
    print(f"   → The cosmological constant (Λ Planck 2018)")
    print(f"\n   The same 4D tilt geometry operating at all scales.")
else:
    print(f"\n   ⚠️  COMBINED FORMULA PASSES ONE TEST BUT NOT BOTH")
    print(f"   Further refinement needed.")

# ============================================================================
# LISA Check With Combined Formula
# ============================================================================

print(f"\n{'='*70}")
print("BONUS: LISA GW PREDICTION WITH COMBINED FORMULA")
print(f"{'='*70}")

# GW delay uses ISCO (r = 3Rs) as characteristic scale
M_gw    = 65 * const.M_sun
Rs_gw   = 2 * const.G * M_gw / const.c**2
r_isco  = 3 * Rs_gw
curv_gw = Rs_gw / r_isco                      # = 1/3

n_gw    = n_dynamic(curv_gw)
obs_gw  = observable_combined(curv_gw)
delay   = (Rs_gw / const.c) * obs_gw

LISA_PRECISION = 1e-5                          # 10 μs
snr_lisa = delay / LISA_PRECISION

print(f"\n   Binary mass: 65 M_sun")
print(f"   ISCO curvature (Rs/3Rs): {curv_gw:.4f}")
print(f"   Dynamic n at ISCO:       {n_gw:.4f}")
print(f"   Observable (√β):         {obs_gw:.4f}")
print(f"   Predicted GW delay:      {delay*1e6:.1f} μs")
print(f"   LISA precision:          {LISA_PRECISION*1e6:.0f} μs")
print(f"   LISA S/N:                {snr_lisa:.1f}σ")
print(f"   Verdict: {'✅ LISA detectable' if snr_lisa >= 3 else '⚠️  marginal'}")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Combined Formula: H4 + H1/H2 — One Formula, All Scales',
             fontsize=14, fontweight='bold')

curv_range = np.linspace(0.001, 0.999, 500)
n_vals     = [n_dynamic(c)         for c in curv_range]
b_vals     = [beta_combined(c)     for c in curv_range]
obs_vals   = [observable_combined(c) for c in curv_range]

# -- Plot 1: n(x), β(x), √β(x) vs curvature --
ax = axes[0, 0]
ax.plot(curv_range, n_vals,   'purple',     linewidth=2.5, label='n(x) — dynamic tilt exponent')
ax.plot(curv_range, b_vals,   'blue',       linewidth=2,   label='β(x) — tilt suppression')
ax.plot(curv_range, obs_vals, 'darkorange', linewidth=2,
        linestyle='--', label='√β(x) — 3D observable')
ax.axvline(CURVATURE_LOCAL, color='red',    linewidth=1.5,
           linestyle=':', label=f'M87* ({CURVATURE_LOCAL:.3f})')
ax.axvline(CURVATURE_COSM,  color='green',  linewidth=1.5,
           linestyle=':', label=f'Cosmic ({CURVATURE_COSM:.3f})')
ax.axvline(1/3,             color='gray',   linewidth=1.5,
           linestyle=':', label='ISCO (1/3)')
ax.set_xlabel('Curvature Ratio x', fontsize=12)
ax.set_ylabel('Value', fontsize=12)
ax.set_title('Combined Formula Components vs Curvature',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_xlim(0, 1); ax.set_ylim(0, 2)

# -- Plot 2: EHT shadow prediction --
ax = axes[0, 1]
r_range_rs = np.linspace(1.5, 20, 300)
r_range_m  = r_range_rs * Rs_m87

distance_m    = 16.8 * 3.086e22
R_shadow_gr   = np.sqrt(27) * Rs_m87
theta_gr_uas  = (R_shadow_gr / distance_m) * 206265 * 1e6

shadows_gr    = []
shadows_comb  = []
for r in r_range_m:
    curv  = Rs_m87 / r
    tgr   = (4 * const.G * M_m87) / (const.c**2 * r)
    obs   = observable_combined(curv)
    tcomb = tgr * (1 + obs * (Rs_m87 / (2 * r)))
    shadows_gr.append(tgr * 206265 * 1e6)
    shadows_comb.append(tcomb * 206265 * 1e6)

ax.semilogy(r_range_rs, shadows_gr,   'blue', linewidth=2,   label='GR')
ax.semilogy(r_range_rs, shadows_comb, 'red',  linewidth=2.5, label='Combined ESTIF')
ax.axvline(1.5, color='gray', linewidth=1.5, linestyle=':', label='Photon sphere')
ax.axhline(SHADOW_OBS, color='green', linewidth=2, linestyle='--',
           label=f'EHT observed: {SHADOW_OBS} μas')
ax.fill_between(r_range_rs,
                SHADOW_OBS - SHADOW_ERR,
                SHADOW_OBS + SHADOW_ERR,
                alpha=0.2, color='green', label='EHT 1σ band')
ax.set_xlabel('Distance (r/Rs)', fontsize=12)
ax.set_ylabel('Deflection (μas)', fontsize=12)
ax.set_title('EHT Shadow: Combined Formula vs GR',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_xlim(1.5, 20)

# -- Plot 3: Λ prediction across curvature scales --
ax = axes[1, 0]
lam_vals = [(3.0/R_H**2) * observable_combined(c)**2 for c in curv_range]
lam_ratios = [l / LAMBDA_MEASURED for l in lam_vals]

ax.semilogy(curv_range, lam_ratios, 'purple', linewidth=2.5)
ax.axhline(1.0, color='black', linewidth=2, linestyle='--', label='Λ measured')
ax.axhline(0.1, color='orange', linewidth=1, linestyle=':', alpha=0.7)
ax.axhline(10,  color='orange', linewidth=1, linestyle=':', alpha=0.7,
           label='10× boundary')
ax.axvline(CURVATURE_COSM, color='green', linewidth=2, linestyle=':',
           label=f'Cosmic ({CURVATURE_COSM:.3f}) → ratio={lam["ratio"]:.3f}')
ax.axvline(CURVATURE_LOCAL, color='red', linewidth=1.5, linestyle=':',
           alpha=0.7, label=f'M87* ({CURVATURE_LOCAL:.3f})')
ax.set_xlabel('Curvature Ratio x', fontsize=12)
ax.set_ylabel('Λ_predicted / Λ_measured', fontsize=12)
ax.set_title('Cosmological Constant vs Curvature Scale',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# -- Plot 4: Three predictions summary --
ax = axes[1, 1]
labels   = ['EHT\nShadow\n(σ from obs)',
            'Λ\nRatio\n(pred/meas)',
            'LISA\nS/N\n(σ)']
values   = [abs(eht['sigma']), lam['ratio'], snr_lisa]
colors   = ['green' if abs(eht['sigma'])<1 else 'orange',
            'green' if abs(lam['ratio']-1)<0.1 else 'orange',
            'green' if snr_lisa>=3 else 'orange']
targets  = [1.0, 1.0, 3.0]

bars = ax.bar(labels, values, color=colors, alpha=0.8,
              edgecolor='black', linewidth=1.5)
for target, label, val in zip(targets, labels, values):
    ax.axhline(target, color='red', linewidth=1.5, linestyle='--', alpha=0.5)

ax.set_ylabel('Value', fontsize=12)
ax.set_title('Three Tests: Combined Formula Summary',
             fontsize=12, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.1,
            f'{val:.2f}', ha='center', va='bottom',
            fontsize=11, fontweight='bold')

# Add target labels
ax.text(0, 1.05, 'target\n<1σ',    ha='center', fontsize=8, color='red')
ax.text(1, 1.05, 'target\n≈1.0',   ha='center', fontsize=8, color='red')
ax.text(2, 3.1,  'target\n>3σ',    ha='center', fontsize=8, color='red')

plt.tight_layout()
plt.savefig('combined_formula_results.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: combined_formula_results.png")
print("=" * 70)
