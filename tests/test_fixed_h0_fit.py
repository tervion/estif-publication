"""
test_fixed_h0_fit.py

Cleanest possible test: fix H₀ AND Ωm at Planck 2018 values.
Only one free cosmological parameter: ALPHA_COSMO.

WHY BOTH FIXED:
The Pantheon+ MU_SH0ES column is pre-calibrated assuming Planck
parameters. Freeing Ωm against these calibrated distance moduli
creates a degeneracy that drives Ωm to the boundary (0.45).
Fixing both H₀ and Ωm removes all parameter degeneracies.

THE QUESTION BEING ANSWERED:
Does the data prefer α > 0 (ESTIF tilt correction) over α = 0 (ΛCDM)
when both H₀ and Ωm are held at their known Planck values?

FREE PARAMETERS:
    ALPHA_COSMO  — tilt x(z) exponent  (0 = pure ΛCDM)
    M_offset     — SN absolute magnitude offset (nuisance, always free)

DATASETS:
    Pantheon+ corrected (1580 SNe) + BAO (6 points, BOSS/eBOSS)
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
# Fixed cosmological parameters (Planck 2018)
# ============================================================================

H0_FIXED  = 67.66    # km/s/Mpc
OM_FIXED  = 0.3111   # matter density
MPC_TO_M  = 3.085677581e22
C_KMS     = const.c / 1e3
R_D       = 147.09   # Sound horizon [Mpc]
H0_SI     = H0_FIXED * 1e3 / MPC_TO_M

# Pre-computed x_0 for ESTIF
_R_H  = const.c / H0_SI
_X_0  = _R_H / 4.4e26
_OBS0 = estif.observable_combined(_X_0)
_OL   = 1.0 - OM_FIXED

# Correct BAO data (BOSS DR12, eBOSS DR16 — D_M/r_d)
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

pp_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '../data/pantheon_plus.dat')

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
            zv = float(parts[zc]); mv = float(parts[muc])
            ev = float(parts[erc]); ic = int(parts[isc])
            if zv > 0.01 and ic == 0 and ev < 10 and mv > 0:
                z_sn.append(zv); mu_sn.append(mv); err_sn.append(ev)
        except (ValueError, IndexError):
            continue

z_sn = np.array(z_sn); mu_sn = np.array(mu_sn); err_sn = np.array(err_sn)

print("=" * 70)
print("FIXED H₀ + Ωm FIT: DOES DATA PREFER α > 0?")
print("=" * 70)
print(f"\n   H₀ = {H0_FIXED} km/s/Mpc  (Planck 2018, fixed)")
print(f"   Ωm = {OM_FIXED}            (Planck 2018, fixed)")
print(f"   Free parameters: α (ALPHA_COSMO), M_offset only")
print(f"   SN data: {len(z_sn)} SNe  |  BAO data: {len(BAO_DATA)} points")

# ============================================================================
# Distance functions — H₀ and Ωm fixed
# ============================================================================

def H_estif(z, alpha):
    """H(z) with fixed H₀, fixed Ωm, free α."""
    if alpha < 1e-6:
        Ot = _OL
    else:
        z_eff = min(float(z), 2.0)
        x_z   = _X_0 * (1.0 + z_eff)**alpha
        obs_z = estif.observable_combined(x_z)
        Ot    = _OL * (_OBS0 / obs_z)**2 if obs_z > 0 else _OL
    return H0_SI * np.sqrt(OM_FIXED*(1+z)**3 + Ot)


def comoving_dist(z, alpha):
    """Comoving distance [Mpc]."""
    def integrand(zp):
        return C_KMS / (H_estif(zp, alpha) * MPC_TO_M / 1e3)
    val, _ = quad(integrand, 0, float(z), limit=80)
    return val


def mu_model(z, alpha, M):
    """Distance modulus."""
    d_C = comoving_dist(z, alpha)
    d_L = (1.0 + z) * d_C
    return 5.0 * np.log10(max(d_L, 1e-10)) + 25.0 + M


# ============================================================================
# χ² function — only α and M free
# ============================================================================

def chi2_total(alpha, M):
    """Total χ² = SN + BAO with H₀ and Ωm fixed."""
    # SN
    mu_pred = np.array([mu_model(z, alpha, M) for z in z_sn])
    c2_sn   = np.sum(((mu_sn - mu_pred) / err_sn)**2)
    # BAO
    c2_bao = 0.0
    for z_b, DM_obs, sig_b, *_ in BAO_DATA:
        d_C     = comoving_dist(z_b, alpha)
        DM_pred = d_C / R_D
        c2_bao += ((DM_obs - DM_pred) / sig_b)**2
    return c2_sn + c2_bao


def chi2_wrapper(params, alpha_fixed=None):
    """Wrapper for scipy minimize."""
    if alpha_fixed is not None:
        M     = params[0]
        alpha = alpha_fixed
    else:
        alpha, M = params
        if alpha < 0 or alpha > 0.5: return 1e10
    return chi2_total(alpha, M)

# ============================================================================
# Step 1: Best M for ΛCDM (α = 0)
# ============================================================================

print(f"\n{'='*70}")
print("STEP 1: ΛCDM  (α = 0 fixed, only M free)")
print(f"{'='*70}")
print(f"   Fitting...", flush=True)

res_lcdm  = minimize(lambda p: chi2_wrapper(p, alpha_fixed=0.0),
                     [-0.18],
                     method='Nelder-Mead',
                     options={'xatol':1e-7, 'fatol':1e-7, 'maxiter':5000})
M_lcdm    = res_lcdm.x[0]
chi2_lcdm = chi2_wrapper(res_lcdm.x, alpha_fixed=0.0)
dof_lcdm  = len(z_sn) + len(BAO_DATA) - 1

print(f"\n   M       = {M_lcdm:.6f}")
print(f"   χ²      = {chi2_lcdm:.4f}")
print(f"   χ²/dof  = {chi2_lcdm/dof_lcdm:.6f}")

# ============================================================================
# Step 2: Best α and M for ESTIF
# ============================================================================

print(f"\n{'='*70}")
print("STEP 2: ESTIF  (α free, M free)")
print(f"{'='*70}")
print(f"   Fitting...", flush=True)

# Grid search for best starting α
best_start = [0.10, M_lcdm]
best_grid  = 1e10
for a_try in [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    c = chi2_wrapper([a_try, M_lcdm])
    if c < best_grid:
        best_grid  = c
        best_start = [a_try, M_lcdm]

res_estif  = minimize(chi2_wrapper, best_start,
                      method='Nelder-Mead',
                      options={'xatol':1e-8, 'fatol':1e-8, 'maxiter':10000})
alpha_e, M_e = res_estif.x
chi2_estif   = chi2_wrapper(res_estif.x)
dof_estif    = len(z_sn) + len(BAO_DATA) - 2

print(f"\n   α       = {alpha_e:.6f}  (0 = pure ΛCDM)")
print(f"   M       = {M_e:.6f}")
print(f"   χ²      = {chi2_estif:.4f}")
print(f"   χ²/dof  = {chi2_estif/dof_estif:.6f}")

# ============================================================================
# Step 3: Statistical significance
# ============================================================================

print(f"\n{'='*70}")
print("STEP 3: IS α > 0 PREFERRED?")
print(f"{'='*70}")

delta_chi2  = chi2_lcdm - chi2_estif
p_value     = chi2_dist.sf(max(delta_chi2, 0), 1)
sigma_equiv = np.sqrt(chi2_dist.isf(max(p_value, 1e-15), 1)) if delta_chi2 > 0 else 0.0

print(f"\n   ΛCDM (α=0):      χ² = {chi2_lcdm:.4f}  (dof={dof_lcdm})")
print(f"   ESTIF (α={alpha_e:.4f}): χ² = {chi2_estif:.4f}  (dof={dof_estif})")
print(f"   Δχ²     = {delta_chi2:+.4f}")
print(f"   p-value = {p_value:.6f}")
print(f"   Significance: {sigma_equiv:.2f}σ")

if sigma_equiv >= 3.0:
    print(f"\n   ✅ SIGNIFICANT: data prefers α > 0 at {sigma_equiv:.2f}σ")
elif sigma_equiv >= 2.0:
    print(f"\n   ✅ SUGGESTIVE: data prefers α > 0 at {sigma_equiv:.2f}σ")
elif sigma_equiv >= 1.0:
    print(f"\n   ⚠️  MARGINAL: {sigma_equiv:.2f}σ preference for α > 0")
else:
    print(f"\n   ❌ NO PREFERENCE: cannot distinguish α=0 from α>0")

# ============================================================================
# Step 4: χ² profile scan
# ============================================================================

print(f"\n{'='*70}")
print("STEP 4: χ² PROFILE vs α")
print(f"{'='*70}")
print(f"   Scanning...", flush=True)

alpha_scan = np.linspace(0.0, 0.4, 30)
chi2_scan  = []
for a in alpha_scan:
    res_a = minimize(lambda p: chi2_wrapper(p, alpha_fixed=a),
                     [M_lcdm],
                     method='Nelder-Mead',
                     options={'xatol':1e-7, 'fatol':1e-7})
    chi2_scan.append(chi2_wrapper(res_a.x, alpha_fixed=a))
chi2_scan = np.array(chi2_scan)

idx_min   = np.argmin(chi2_scan)
alpha_min = alpha_scan[idx_min]
chi2_min  = chi2_scan[idx_min]

below_1sig = alpha_scan[chi2_scan < chi2_min + 1.0]
below_2sig = alpha_scan[chi2_scan < chi2_min + 4.0]

print(f"\n   χ² minimum at α = {alpha_min:.3f}  (χ² = {chi2_min:.2f})")
if len(below_1sig) > 0:
    print(f"   1σ range:  α = {below_1sig[0]:.3f} – {below_1sig[-1]:.3f}")
if len(below_2sig) > 0:
    print(f"   2σ range:  α = {below_2sig[0]:.3f} – {below_2sig[-1]:.3f}")

delta_from_min = chi2_scan[0] - chi2_min
print(f"\n   Δχ²(α=0 vs minimum) = {delta_from_min:+.2f}")
if delta_from_min > 4.0:
    print(f"   α=0 is outside the 2σ region — ΛCDM disfavoured")
elif delta_from_min > 1.0:
    print(f"   α=0 is outside the 1σ region")
else:
    print(f"   α=0 is within the 1σ region — consistent with ΛCDM")

# ============================================================================
# Final verdict
# ============================================================================

print(f"\n{'='*70}")
print("FINAL VERDICT")
print(f"{'='*70}")
print(f"""
   Model: ESTIF Option A
   H₀ = {H0_FIXED} km/s/Mpc  (Planck, fixed)
   Ωm = {OM_FIXED}            (Planck, fixed)

   ΛCDM (α=0):          χ²/dof = {chi2_lcdm/dof_lcdm:.4f}
   ESTIF (α={alpha_e:.4f}):  χ²/dof = {chi2_estif/dof_estif:.4f}

   Δχ² = {delta_chi2:+.4f}  →  {sigma_equiv:.2f}σ preference for ESTIF
   Best-fit α = {alpha_e:.4f}  (previously fixed at 0.1036)
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    f'Fixed H₀={H0_FIXED}, Ωm={OM_FIXED} (Planck): '
    f'α={alpha_e:.4f}  [{sigma_equiv:.2f}σ]',
    fontsize=13, fontweight='bold')

# χ² profile
ax = axes[0]
ax.plot(alpha_scan, chi2_scan, 'red', linewidth=2.5, label='χ²(α) — SN+BAO')
ax.axhline(chi2_lcdm, color='blue', linewidth=2, linestyle='--',
           label=f'ΛCDM (α=0): χ²={chi2_lcdm:.1f}')
ax.axhline(chi2_min + 1.0, color='green', linewidth=1.5, linestyle=':',
           alpha=0.7, label='1σ above minimum')
ax.axhline(chi2_min + 4.0, color='orange', linewidth=1.5, linestyle=':',
           alpha=0.7, label='2σ above minimum')
ax.axvline(alpha_e, color='red', linewidth=2, linestyle=':',
           label=f'Best α={alpha_e:.4f}')
ax.axvline(0.1036, color='purple', linewidth=1.5, linestyle=':',
           alpha=0.7, label='Previous α=0.1036')
ax.fill_between(alpha_scan, chi2_min + 4.0, chi2_scan,
                where=(chi2_scan < chi2_min + 4.0),
                alpha=0.15, color='green', label='2σ preferred region')
ax.set_xlabel('α (ALPHA_COSMO)', fontsize=12)
ax.set_ylabel('χ² (SN + BAO)', fontsize=12)
ax.set_title(f'χ² Profile vs α\n(H₀={H0_FIXED}, Ωm={OM_FIXED} fixed)', fontsize=11)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Model comparison
ax = axes[1]
models  = [f'ΛCDM\n(α=0)', f'ESTIF\n(α={alpha_e:.4f})']
vals    = [chi2_lcdm/dof_lcdm, chi2_estif/dof_estif]
colors  = ['blue', 'green' if delta_chi2 > 0 else 'red']
bars    = ax.bar(models, vals, color=colors, alpha=0.8, edgecolor='black')
for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.001,
            f'{val:.4f}', ha='center', fontsize=11, fontweight='bold')
ax.text(0.98, 0.95,
        f'Δχ² = {delta_chi2:+.2f}\n{sigma_equiv:.2f}σ',
        transform=ax.transAxes, fontsize=13,
        ha='right', va='top', fontweight='bold',
        color='green' if sigma_equiv >= 2 else 'orange',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
ax.set_ylabel('χ²/dof', fontsize=12)
ax.set_title('Fit Quality: ΛCDM vs ESTIF\n(H₀ and Ωm both Planck-fixed)', fontsize=11)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('fixed_h0_fit.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: fixed_h0_fit.png")
print("=" * 70)
