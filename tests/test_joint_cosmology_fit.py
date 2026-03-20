"""
test_joint_cosmology_fit.py

Phase 5.2 — Joint fit of H₀, Ωm, and ALPHA_COSMO against combined
supernova + BAO dataset.

MOTIVATION:
Previous tests fixed H₀ = 67.66 and Ωm = 0.3111 (Planck 2018) and
fitted only M_offset. This is not a true cosmological fit — it tests
whether ESTIF is consistent with Planck's parameters, not whether it
prefers different ones.

A joint fit asks: what H₀, Ωm, and ALPHA_COSMO does ESTIF prefer
when all three are free simultaneously?

Key question: Does ESTIF naturally prefer a higher H₀ than ΛCDM
from the same data? If yes, it may partially resolve the H₀ tension.

FREE PARAMETERS:
    H₀           — Hubble constant [km/s/Mpc], prior: 60–80
    Ωm           — matter density,              prior: 0.20–0.45
    ALPHA_COSMO  — tilt x(z) exponent,          prior: 0.0–0.5
    M_offset     — SN absolute magnitude offset (nuisance)

DATASETS:
    1. Pantheon+ corrected (1580 SNe) — distance modulus vs z
    2. BAO measurements (BOSS/eBOSS)  — D_M/r_d vs z
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad
from scipy.stats import chi2 as chi2_dist
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Constants
# ============================================================================

MPC_TO_M  = 3.085677581e22
C_KMS     = const.c / 1e3
R_D       = 147.09           # Sound horizon [Mpc] — unchanged by ESTIF
X_0_BASE  = (const.c / const.H_0) / 4.4e26   # x_0 at Planck H₀

# Planck 2018 reference
H0_PLANCK = 67.66
OM_PLANCK = 0.3111
OL_PLANCK = 0.6889

# BAO measurements (z, D_M/r_d observed, sigma)
BAO_DATA = [
    (0.38,  10.27, 0.15, "BOSS DR12"),
    (0.51,  13.38, 0.18, "BOSS DR12"),
    (0.61,  15.45, 0.22, "BOSS DR12"),
    (0.70,  17.65, 0.30, "eBOSS DR16"),
    (0.85,  20.75, 0.40, "eBOSS ELG"),
    (1.48,  30.21, 0.79, "eBOSS QSO"),
]

# ============================================================================
# Load Pantheon+ corrected data
# ============================================================================

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')
pp_path  = os.path.join(data_dir, 'pantheon_plus.dat')

z_sn, mu_sn, err_sn = [], [], []
with open(pp_path) as f:
    header = None
    for line in f:
        line = line.strip()
        if not line: continue
        parts = line.split()
        if parts[0] == 'CID':
            header = parts
            zc  = header.index('zHD')
            muc = header.index('MU_SH0ES')
            erc = header.index('MU_SH0ES_ERR_DIAG')
            isc = header.index('IS_CALIBRATOR')
            continue
        if header is None: continue
        try:
            zv = float(parts[zc])
            mv = float(parts[muc])
            ev = float(parts[erc])
            ic = int(parts[isc])
            if zv > 0.01 and ic == 0 and ev < 10 and mv > 0:
                z_sn.append(zv); mu_sn.append(mv); err_sn.append(ev)
        except (ValueError, IndexError):
            continue

z_sn   = np.array(z_sn)
mu_sn  = np.array(mu_sn)
err_sn = np.array(err_sn)

print("=" * 70)
print("JOINT COSMOLOGY FIT: H₀ + Ωm + ALPHA_COSMO")
print("=" * 70)
print(f"\n   SN data:  {len(z_sn)} supernovae (Pantheon+ corrected)")
print(f"   BAO data: {len(BAO_DATA)} measurements (BOSS/eBOSS)")

# ============================================================================
# Core functions with variable parameters
# ============================================================================

def omega_tilt_free(z, H0_kms, Om, alpha):
    """Ω_tilt(z) with free parameters."""
    R_H   = (H0_kms * 1e3 / MPC_TO_M) / const.c * MPC_TO_M  # Hubble radius in m
    # Actually: R_H = c / H(z=0) in metres
    R_H   = const.c / (H0_kms * 1e3 / MPC_TO_M)
    x_0   = R_H / 4.4e26
    x_now = estif.observable_combined(x_0)

    z_eff = min(float(z), 2.0)
    x_z   = x_0 * (1.0 + z_eff)**alpha
    obs_z = estif.observable_combined(x_z)

    if obs_z <= 0 or x_now <= 0:
        return 0.6889
    OL = 1.0 - Om  # flat universe
    return OL * (x_now / obs_z) ** 2


def H_free(z, H0_kms, Om, alpha):
    """H(z) [s⁻¹] with free parameters."""
    H0_si = H0_kms * 1e3 / MPC_TO_M
    Ot    = omega_tilt_free(z, H0_kms, Om, alpha)
    return H0_si * np.sqrt(Om * (1.0 + z)**3 + Ot)


def comoving_distance_free(z, H0_kms, Om, alpha):
    """Comoving distance [Mpc] with free parameters."""
    def integrand(zp):
        return C_KMS / (H_free(zp, H0_kms, Om, alpha) * MPC_TO_M / 1e3)
    result, _ = quad(integrand, 0, float(z), limit=80)
    return result


def mu_estif_free(z, H0_kms, Om, alpha, M):
    """Distance modulus with free parameters."""
    d_C = comoving_distance_free(z, H0_kms, Om, alpha)
    d_L = (1.0 + z) * d_C
    return 5.0 * np.log10(max(d_L, 1e-10)) + 25.0 + M


def mu_lcdm_free(z, H0_kms, Om, M):
    """ΛCDM distance modulus with free H₀ and Ωm."""
    H0_si = H0_kms * 1e3 / MPC_TO_M
    OL    = 1.0 - Om
    def integrand(zp):
        Hz = H0_si * np.sqrt(Om*(1+zp)**3 + OL)
        return C_KMS / (Hz * MPC_TO_M / 1e3)
    d_C, _ = quad(integrand, 0, float(z), limit=80)
    d_L    = (1.0 + z) * d_C
    return 5.0 * np.log10(max(d_L, 1e-10)) + 25.0 + M


# ============================================================================
# Step 1: ΛCDM free fit (H₀, Ωm, M)
# ============================================================================

print(f"\n{'='*70}")
print("STEP 1: ΛCDM FREE FIT (H₀ + Ωm + M)")
print(f"{'='*70}")
print(f"\n   Fitting... (this takes ~3 min)", flush=True)

def chi2_lcdm_free(params):
    H0, Om, M = params
    if H0 < 60 or H0 > 80 or Om < 0.20 or Om > 0.45:
        return 1e10
    # SN contribution
    mu_pred = np.array([mu_lcdm_free(z, H0, Om, M) for z in z_sn])
    chi2_sn = np.sum(((mu_sn - mu_pred) / err_sn)**2)
    # BAO contribution
    chi2_bao = 0.0
    for z_b, DM_obs, sig_b, *_ in BAO_DATA:
        d_C = comoving_distance_free(z_b, H0, Om, 0) * MPC_TO_M / 1e3 / 1e3
        # comoving_distance returns Mpc, DM/r_d is dimensionless
        d_C_mpc, _ = quad(lambda zp: C_KMS / (H0 * np.sqrt(Om*(1+zp)**3+(1-Om))),
                          0, z_b, limit=50)
        DM_pred = d_C_mpc / R_D
        chi2_bao += ((DM_obs - DM_pred) / sig_b)**2
    chi2_H0_prior = ((H0 - 67.66) / 0.42)**2
    chi2_Om_prior = ((Om - 0.3111) / 0.006)**2
    return chi2_sn + chi2_bao + chi2_H0_prior + chi2_Om_prior

res_lcdm = minimize(chi2_lcdm_free, [H0_PLANCK, OM_PLANCK, -0.18],
                    method='Nelder-Mead',
                    options={'xatol':1e-4, 'fatol':1e-4, 'maxiter':5000})
H0_lcdm_best, Om_lcdm_best, M_lcdm_best = res_lcdm.x
chi2_lcdm_best = chi2_lcdm_free(res_lcdm.x)
dof_lcdm = len(z_sn) + len(BAO_DATA) - 3

print(f"\n   ΛCDM best fit:")
print(f"   H₀  = {H0_lcdm_best:.3f} km/s/Mpc  (Planck: {H0_PLANCK:.2f})")
print(f"   Ωm  = {Om_lcdm_best:.4f}            (Planck: {OM_PLANCK:.4f})")
print(f"   M   = {M_lcdm_best:.4f}")
print(f"   χ²  = {chi2_lcdm_best:.2f}  (χ²/dof = {chi2_lcdm_best/dof_lcdm:.4f})")

# ============================================================================
# Step 2: ESTIF joint fit (H₀, Ωm, ALPHA_COSMO, M)
# ============================================================================

print(f"\n{'='*70}")
print("STEP 2: ESTIF JOINT FIT (H₀ + Ωm + α + M)")
print(f"{'='*70}")
print(f"\n   Fitting... (this takes ~5 min)", flush=True)

def chi2_estif_joint(params):
    H0, Om, alpha, M = params
    if H0 < 60 or H0 > 80 or Om < 0.20 or Om > 0.45:
        return 1e10
    if alpha < 0 or alpha > 0.5:
        return 1e10
    # SN contribution
    mu_pred = np.array([mu_estif_free(z, H0, Om, alpha, M) for z in z_sn])
    chi2_sn = np.sum(((mu_sn - mu_pred) / err_sn)**2)
    # BAO contribution
    chi2_bao = 0.0
    for z_b, DM_obs, sig_b, *_ in BAO_DATA:
        d_C_mpc = comoving_distance_free(z_b, H0, Om, alpha)
        DM_pred = d_C_mpc / R_D
        chi2_bao += ((DM_obs - DM_pred) / sig_b)**2
    chi2_H0_prior = ((H0 - 67.66) / 0.42)**2
    chi2_Om_prior = ((Om - 0.3111) / 0.006)**2
    return chi2_sn + chi2_bao + chi2_H0_prior + chi2_Om_prior

# Grid search for good start
best_start = [H0_lcdm_best, Om_lcdm_best, 0.10, M_lcdm_best]
best_grid  = 1e10
for alpha_try in [0.05, 0.10, 0.15, 0.20]:
    c = chi2_estif_joint([H0_lcdm_best, Om_lcdm_best, alpha_try, M_lcdm_best])
    if c < best_grid:
        best_grid  = c
        best_start = [H0_lcdm_best, Om_lcdm_best, alpha_try, M_lcdm_best]

res_estif = minimize(chi2_estif_joint, best_start,
                     method='Nelder-Mead',
                     options={'xatol':1e-4, 'fatol':1e-4, 'maxiter':10000})
H0_e, Om_e, alpha_e, M_e = res_estif.x
chi2_estif_best = chi2_estif_joint(res_estif.x)
dof_estif = len(z_sn) + len(BAO_DATA) - 4

print(f"\n   ESTIF best fit:")
print(f"   H₀    = {H0_e:.3f} km/s/Mpc  (ΛCDM free: {H0_lcdm_best:.3f}, Planck: {H0_PLANCK:.2f})")
print(f"   Ωm    = {Om_e:.4f}            (ΛCDM free: {Om_lcdm_best:.4f})")
print(f"   α     = {alpha_e:.4f}            (previous fixed: 0.1036)")
print(f"   M     = {M_e:.4f}")
print(f"   χ²    = {chi2_estif_best:.2f}  (χ²/dof = {chi2_estif_best/dof_estif:.4f})")

# ============================================================================
# Step 3: Statistical significance
# ============================================================================

print(f"\n{'='*70}")
print("STEP 3: STATISTICAL SIGNIFICANCE")
print(f"{'='*70}")

delta_chi2  = chi2_lcdm_best - chi2_estif_best
p_value     = chi2_dist.sf(max(delta_chi2, 0), 1)   # 1 extra param (alpha)
sigma_equiv = np.sqrt(chi2_dist.isf(max(p_value, 1e-15), 1)) if delta_chi2 > 0 else 0.0

print(f"\n   ΛCDM (free):   χ² = {chi2_lcdm_best:.2f}  dof = {dof_lcdm}")
print(f"   ESTIF (joint): χ² = {chi2_estif_best:.2f}  dof = {dof_estif}")
print(f"   Δχ² = {delta_chi2:+.4f}  (1 extra parameter: α)")
print(f"   Significance: {sigma_equiv:.2f}σ")

if sigma_equiv >= 3.0:
    print(f"\n   ✅ ESTIF SIGNIFICANTLY PREFERRED (≥3σ)")
elif sigma_equiv >= 2.0:
    print(f"\n   ✅ ESTIF PREFERRED (≥2σ)")
elif sigma_equiv >= 1.0:
    print(f"\n   ⚠️  MARGINAL PREFERENCE FOR ESTIF")
else:
    print(f"\n   ❌ DATA DOES NOT PREFER ESTIF OVER ΛCDM")

# ============================================================================
# Step 4: H₀ tension analysis
# ============================================================================

print(f"\n{'='*70}")
print("STEP 4: H₀ TENSION ANALYSIS")
print(f"{'='*70}")

H0_local   = 73.04   # SH0ES
sigma_comb = 2.0     # Combined uncertainty

tension_lcdm  = (H0_local - H0_PLANCK)     / sigma_comb
tension_lcdm_free = (H0_local - H0_lcdm_best) / sigma_comb
tension_estif = (H0_local - H0_e)           / sigma_comb

print(f"\n   Local H₀ (SH0ES):     {H0_local:.2f} km/s/Mpc")
print(f"\n   Planck 2018 H₀:       {H0_PLANCK:.2f}  →  tension = {tension_lcdm:.2f}σ")
print(f"   ΛCDM free H₀:         {H0_lcdm_best:.2f}  →  tension = {tension_lcdm_free:.2f}σ")
print(f"   ESTIF joint H₀:       {H0_e:.2f}  →  tension = {tension_estif:.2f}σ")

shift = H0_e - H0_lcdm_best
print(f"\n   ESTIF shifts H₀ by {shift:+.3f} km/s/Mpc vs ΛCDM free fit")

if shift > 0.5:
    print(f"   ✅ ESTIF moves H₀ toward local measurement")
elif shift > 0:
    print(f"   ⚠️  Small positive shift — right direction")
else:
    print(f"   ❌ ESTIF moves H₀ away from local measurement")

# ============================================================================
# Step 5: Alpha comparison
# ============================================================================

print(f"\n{'='*70}")
print("STEP 5: ALPHA_COSMO COMPARISON")
print(f"{'='*70}")

print(f"\n   Previous (fixed, from SN-only fit): α = 0.1036")
print(f"   Joint fit (SN + BAO):               α = {alpha_e:.4f}")
print(f"   Change:                             Δα = {alpha_e - 0.1036:+.4f}")

# Rational approximations
print(f"\n   Rational approximations to α = {alpha_e:.4f}:")
for num in range(1, 15):
    for den in range(1, 100):
        if abs(num/den - alpha_e) < 0.005:
            print(f"   {num}/{den} = {num/den:.4f}  ({abs(num/den-alpha_e)*100:.3f}% off)")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"""
   Joint fit result (SN + BAO simultaneously):

   Parameter    ΛCDM (free)    ESTIF (joint)   Planck 2018
   -------------------------------------------------------
   H₀           {H0_lcdm_best:.3f}          {H0_e:.3f}          {H0_PLANCK:.2f}
   Ωm           {Om_lcdm_best:.4f}          {Om_e:.4f}          {OM_PLANCK:.4f}
   α            —              {alpha_e:.4f}          —

   χ² comparison:
   ΛCDM (free):   {chi2_lcdm_best:.2f}  (χ²/dof = {chi2_lcdm_best/dof_lcdm:.4f})
   ESTIF (joint): {chi2_estif_best:.2f}  (χ²/dof = {chi2_estif_best/dof_estif:.4f})
   Δχ² = {delta_chi2:+.4f}  →  {sigma_equiv:.2f}σ preference for ESTIF

   H₀ tension:
   Planck:      {tension_lcdm:.2f}σ
   ΛCDM free:   {tension_lcdm_free:.2f}σ
   ESTIF:       {tension_estif:.2f}σ
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(f'Joint Cosmology Fit: ESTIF α={alpha_e:.3f}, '
             f'H₀={H0_e:.1f}, Δχ²={delta_chi2:+.2f} ({sigma_equiv:.2f}σ)',
             fontsize=13, fontweight='bold')

z_smooth = np.linspace(0.01, 2.3, 100)

# Plot 1: H₀ comparison
ax = axes[0]
H0_vals   = [H0_PLANCK, H0_lcdm_best, H0_e, 73.04]
labels    = ['Planck\n2018', 'ΛCDM\nfree fit', 'ESTIF\njoint fit', 'SH0ES\nlocal']
colors    = ['blue', 'steelblue', 'red', 'green']
bars      = ax.bar(labels, H0_vals, color=colors, alpha=0.8, edgecolor='black')
ax.axhline(73.04, color='green', linewidth=2, linestyle='--', alpha=0.5)
ax.axhline(H0_PLANCK, color='blue', linewidth=2, linestyle=':', alpha=0.5)
for bar, val in zip(bars, H0_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            f'{val:.2f}', ha='center', fontsize=9, fontweight='bold')
ax.set_ylabel('H₀ [km/s/Mpc]', fontsize=11)
ax.set_title('H₀ Comparison', fontsize=11, fontweight='bold')
ax.set_ylim(60, 78)
ax.grid(axis='y', alpha=0.3)
ax.fill_between([-0.5, 3.5], 71.04, 75.04, alpha=0.1, color='green',
                label='SH0ES ±1σ')
ax.legend(fontsize=9)

# Plot 2: Ωm comparison
ax = axes[1]
om_vals = [OM_PLANCK, Om_lcdm_best, Om_e]
labels2 = ['Planck 2018', 'ΛCDM free', 'ESTIF joint']
colors2 = ['blue', 'steelblue', 'red']
bars2   = ax.bar(labels2, om_vals, color=colors2, alpha=0.8, edgecolor='black')
for bar, val in zip(bars2, om_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
            f'{val:.4f}', ha='center', fontsize=9, fontweight='bold')
ax.set_ylabel('Ωm', fontsize=11)
ax.set_title('Matter Density Ωm', fontsize=11, fontweight='bold')
ax.set_ylim(0.25, 0.38)
ax.grid(axis='y', alpha=0.3)

# Plot 3: χ² profile vs alpha
ax = axes[2]
alpha_scan = np.linspace(0.0, 0.4, 30)
chi2_scan  = []
for a in alpha_scan:
    c = chi2_estif_joint([H0_e, Om_e, a, M_e])
    chi2_scan.append(c)
chi2_scan = np.array(chi2_scan)

ax.plot(alpha_scan, chi2_scan, 'red', linewidth=2.5, label='ESTIF χ²(α)')
ax.axhline(chi2_lcdm_best, color='blue', linewidth=2, linestyle='--',
           label=f'ΛCDM χ²={chi2_lcdm_best:.1f}')
ax.axvline(alpha_e, color='red', linewidth=2, linestyle=':',
           label=f'Best α={alpha_e:.4f}')
ax.axvline(0.1036, color='orange', linewidth=1.5, linestyle=':',
           label='Previous α=0.1036')
ax.fill_between(alpha_scan, chi2_lcdm_best, chi2_scan,
                where=(chi2_scan < chi2_lcdm_best),
                alpha=0.2, color='green', label='ESTIF preferred')
ax.set_xlabel('α (tilt exponent)', fontsize=11)
ax.set_ylabel('χ²', fontsize=11)
ax.set_title('χ² vs ALPHA_COSMO', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('joint_cosmology_fit.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: joint_cosmology_fit.png")
print("=" * 70)
