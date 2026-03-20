"""
test_pantheon_plus_fit.py

Downloads the Pantheon+ dataset (1701 supernovae) and runs the
same ESTIF Approach A fit as test_cosmo_approach_a_fit.py.

Pantheon+ reference:
  Scolnic et al. 2022, ApJ 938, 113
  https://github.com/PantheonPlusSH0ES/DataRelease

The key prediction:
  If the 1.78σ signal from 580 SNe scales as √N,
  1701 SNe should give approximately:
  1.78 × √(1701/580) ≈ 3.0σ
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import chi2 as chi2_dist
import urllib.request
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Download Pantheon+ data
# ============================================================================

PANTHEON_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease"
    "/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
)
LOCAL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '../data/pantheon_plus.dat')

print("=" * 70)
print("PANTHEON+ FIT: ESTIF vs 1701 SUPERNOVAE")
print("=" * 70)

if not os.path.exists(LOCAL_PATH):
    print(f"\n   Downloading Pantheon+ data...")
    try:
        urllib.request.urlretrieve(PANTHEON_URL, LOCAL_PATH)
        print(f"   ✅ Downloaded to {LOCAL_PATH}")
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        print(f"\n   Manual download:")
        print(f"   1. Go to: https://github.com/PantheonPlusSH0ES/DataRelease")
        print(f"   2. Navigate to: Pantheon+_Data/4_DISTANCES_AND_COVAR/")
        print(f"   3. Download: Pantheon+SH0ES.dat")
        print(f"   4. Save to: {LOCAL_PATH}")
        sys.exit(1)
else:
    print(f"\n   Using cached Pantheon+ data: {LOCAL_PATH}")

# ============================================================================
# Parse Pantheon+ data
# ============================================================================
# Columns: CID IDSURVEY zHD zHDERR zCMB zCMBERR zHEL zHELERR
#          m_b_corr m_b_corr_err_DIAG MU_SH0ES MU_SH0ES_ERR_DIAG
#          CEPH_DIST IS_CALIBRATOR USED_IN_SH0ES_HF ...

z_data, mu_data, mu_err = [], [], []

with open(LOCAL_PATH) as f:
    header = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == 'CID':
            header = parts
            # Find column indices
            z_col  = header.index('zHD')
            mu_col = header.index('MU_SH0ES')
            er_col = header.index('MU_SH0ES_ERR_DIAG')
            calib_col = header.index('IS_CALIBRATOR')
            continue
        if header is None:
            continue
        try:
            z     = float(parts[z_col])
            mu    = float(parts[mu_col])
            err   = float(parts[er_col])
            is_calib = int(parts[calib_col])
            # Skip calibrators (Cepheid-anchored) and very low-z
            # Use only cosmological sample: z > 0.01, not calibrators
            if z > 0.01 and is_calib == 0 and err < 10 and mu > 0:
                z_data.append(z)
                mu_data.append(mu)
                mu_err.append(err)
        except (ValueError, IndexError):
            continue

z_data   = np.array(z_data)
mu_data  = np.array(mu_data)
sigma_tot= np.array(mu_err)
N_sn     = len(z_data)

print(f"\n   Loaded {N_sn} cosmological supernovae")
print(f"   Redshift range: {z_data.min():.4f} – {z_data.max():.4f}")
print(f"   Predicted significance (scaled from 580 SNe result):")
predicted_sigma = 1.78 * np.sqrt(N_sn / 580)
print(f"   1.78σ × √({N_sn}/580) ≈ {predicted_sigma:.2f}σ")

# ============================================================================
# Core functions (same as previous script)
# ============================================================================

R_H          = const.c / const.H_0
r_universe_0 = 4.4e26
x_0          = R_H / r_universe_0

def mu_lcdm(z, M=0.0):
    return estif.distance_modulus_lcdm(z) + M

def tilt_correction(z, alpha):
    x_z     = x_0 * (1.0 + z)**alpha
    obs_z   = estif.observable_combined(x_z)
    obs_now = estif.observable_combined(x_0)
    if obs_z <= 0 or obs_now <= 0:
        return 0.0
    return -2.5 * np.log10(obs_z / obs_now)

def mu_lcdm_vec(z_arr, M):
    return np.array([mu_lcdm(z, M) for z in z_arr])

def tilt_corr_vec(z_arr, alpha):
    return np.array([tilt_correction(z, alpha) for z in z_arr])

# ============================================================================
# Step 1: ΛCDM baseline
# ============================================================================

print(f"\n   Fitting ΛCDM baseline (this may take a moment)...")

def chi2_lcdm(M):
    return np.sum(((mu_data - mu_lcdm_vec(z_data, M)) / sigma_tot)**2)

res_lcdm       = minimize_scalar(chi2_lcdm, bounds=(-2, 2), method='bounded')
M_lcdm         = res_lcdm.x
chi2_lcdm_best = chi2_lcdm(M_lcdm)
dof_lcdm       = N_sn - 1

print(f"\n{'='*70}")
print("STEP 1: ΛCDM BASELINE")
print(f"{'='*70}")
print(f"\n   M_offset   = {M_lcdm:.6f}")
print(f"   χ²         = {chi2_lcdm_best:.4f}")
print(f"   χ²/dof     = {chi2_lcdm_best/dof_lcdm:.6f}")

# ============================================================================
# Step 2: ESTIF fixed α = 0.10
# ============================================================================

print(f"\n   Fitting ESTIF (α=0.10 fixed)...")

def chi2_estif_fixed(M):
    return np.sum(((mu_data - mu_lcdm_vec(z_data, M)
                    - tilt_corr_vec(z_data, 0.10)) / sigma_tot)**2)

res_fixed       = minimize_scalar(chi2_estif_fixed, bounds=(-2, 2),
                                   method='bounded')
M_fixed         = res_fixed.x
chi2_fixed      = chi2_estif_fixed(M_fixed)
dof_fixed       = N_sn - 1

print(f"\n{'='*70}")
print("STEP 2: ESTIF FIXED α = 0.10")
print(f"{'='*70}")
print(f"\n   M_offset   = {M_fixed:.6f}")
print(f"   χ²         = {chi2_fixed:.4f}")
print(f"   χ²/dof     = {chi2_fixed/dof_fixed:.6f}")
print(f"   Δχ²        = {chi2_fixed - chi2_lcdm_best:+.4f}")

# ============================================================================
# Step 3: Joint fit
# ============================================================================

print(f"\n   Fitting ESTIF (joint M + α)...")

def chi2_joint(params):
    M, alpha = params
    if alpha < 0 or alpha > 1:
        return 1e10
    mu_pred = mu_lcdm_vec(z_data, M) + tilt_corr_vec(z_data, alpha)
    return np.sum(((mu_data - mu_pred) / sigma_tot)**2)

res_joint        = minimize(chi2_joint, [M_lcdm, 0.10],
                             method='Nelder-Mead',
                             options={'xatol':1e-8, 'fatol':1e-8,
                                      'maxiter':50000})
M_joint, alpha_joint = res_joint.x
chi2_joint_best  = chi2_joint([M_joint, alpha_joint])
dof_joint        = N_sn - 2

print(f"\n{'='*70}")
print("STEP 3: JOINT FIT (M + α free)")
print(f"{'='*70}")
print(f"\n   M_offset   = {M_joint:.6f}")
print(f"   Best α     = {alpha_joint:.6f}")
print(f"   χ²         = {chi2_joint_best:.4f}")
print(f"   χ²/dof     = {chi2_joint_best/dof_joint:.6f}")
print(f"   Δχ²        = {chi2_joint_best - chi2_lcdm_best:+.4f}")

# ============================================================================
# Step 4: Statistical significance
# ============================================================================

delta_chi2  = chi2_lcdm_best - chi2_joint_best
p_value     = chi2_dist.sf(delta_chi2, 1)
sigma_equiv = np.sqrt(chi2_dist.isf(p_value, 1)) if p_value > 0 else np.inf

print(f"\n{'='*70}")
print("STEP 4: STATISTICAL SIGNIFICANCE")
print(f"{'='*70}")
print(f"\n   Δχ² = {delta_chi2:.4f}  (Δdof = 1)")
print(f"   p-value = {p_value:.6f}")
print(f"   Significance: {sigma_equiv:.2f}σ")
print(f"   Predicted (from scaling): {predicted_sigma:.2f}σ")

if p_value < 0.003:
    verdict = "✅ STATISTICALLY SIGNIFICANT (>3σ)"
elif p_value < 0.05:
    verdict = "✅ STATISTICALLY SIGNIFICANT (>2σ)"
elif p_value < 0.10:
    verdict = "⚠️  MARGINAL (1-2σ)"
else:
    verdict = "❌ NOT SIGNIFICANT"
print(f"   Verdict: {verdict}")

# ============================================================================
# Step 5: Residual analysis
# ============================================================================

res_lcdm_arr  = mu_data - mu_lcdm_vec(z_data, M_lcdm)
res_estif_arr = (mu_data
                 - mu_lcdm_vec(z_data, M_joint)
                 - tilt_corr_vec(z_data, alpha_joint))

z_bins = np.percentile(z_data, [0, 25, 50, 75, 100])
bin_labels = ['low-z', 'mid-low-z', 'mid-high-z', 'high-z']

print(f"\n{'='*70}")
print("STEP 5: RESIDUALS BY REDSHIFT BIN")
print(f"{'='*70}")
print(f"\n   {'Bin':<14} {'z range':<16} {'ΛCDM mean':>12} "
      f"{'ESTIF mean':>12} {'better?'}")
print("   " + "-"*60)

improvements = 0
for i in range(4):
    mask = (z_data >= z_bins[i]) & (z_data < z_bins[i+1])
    if i == 3:
        mask = (z_data >= z_bins[i]) & (z_data <= z_bins[i+1])
    mean_l = np.mean(res_lcdm_arr[mask])
    mean_e = np.mean(res_estif_arr[mask])
    better = abs(mean_e) < abs(mean_l)
    if better:
        improvements += 1
    flag = "✅" if better else "⚠️"
    print(f"   {bin_labels[i]:<14} {z_bins[i]:.3f}–{z_bins[i+1]:.3f}     "
          f"{mean_l:>+10.4f}   {mean_e:>+10.4f}   {flag}")

# ============================================================================
# Final verdict
# ============================================================================

print(f"\n{'='*70}")
print("FINAL VERDICT")
print(f"{'='*70}")
print(f"""
   Dataset:        Pantheon+ ({N_sn} SNe, z = {z_data.min():.3f}–{z_data.max():.3f})
   ΛCDM:           χ²/dof = {chi2_lcdm_best/dof_lcdm:.6f}
   ESTIF (α=0.10): χ²/dof = {chi2_fixed/dof_fixed:.6f}
   ESTIF (joint):  χ²/dof = {chi2_joint_best/dof_joint:.6f}

   Best-fit α:     {alpha_joint:.4f}  (previous: 0.1036 from 580 SNe)
   Δχ²:            {chi2_joint_best - chi2_lcdm_best:+.4f}
   Significance:   {sigma_equiv:.2f}σ  (predicted: {predicted_sigma:.2f}σ)
   Residuals improved in {improvements}/4 bins

   {verdict}
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle(f'Pantheon+ ({N_sn} SNe): ESTIF α={alpha_joint:.3f} vs ΛCDM '
             f'[{sigma_equiv:.2f}σ improvement]',
             fontsize=13, fontweight='bold')

z_smooth   = np.linspace(0.01, z_data.max()*1.02, 300)
mu_lcdm_s  = np.array([mu_lcdm(z, M_lcdm) for z in z_smooth])
mu_estif_s = np.array([mu_lcdm(z, M_joint) + tilt_correction(z, alpha_joint)
                        for z in z_smooth])
corr_s     = np.array([tilt_correction(z, alpha_joint) for z in z_smooth])

# Plot 1: Hubble diagram
ax = axes[0, 0]
ax.errorbar(z_data, mu_data, yerr=sigma_tot, fmt='.', color='gray',
            alpha=0.2, markersize=2, elinewidth=0.3, zorder=1)
ax.plot(z_smooth, mu_lcdm_s,  'blue', linewidth=2.5, label='ΛCDM', zorder=3)
ax.plot(z_smooth, mu_estif_s, 'red',  linewidth=2.5, linestyle='--',
        label=f'ESTIF (α={alpha_joint:.3f})', zorder=3)
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Distance Modulus μ', fontsize=12)
ax.set_title(f'Hubble Diagram — Pantheon+ ({N_sn} SNe)', fontsize=12,
             fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# Plot 2: Residuals
ax = axes[0, 1]
ax.scatter(z_data, res_lcdm_arr,  color='blue', alpha=0.2, s=4,
           label=f'ΛCDM (χ²/dof={chi2_lcdm_best/dof_lcdm:.4f})')
ax.scatter(z_data, res_estif_arr, color='red',  alpha=0.2, s=4,
           label=f'ESTIF (χ²/dof={chi2_joint_best/dof_joint:.4f})')
ax.axhline(0, color='black', linewidth=1.5)
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Residual (mag)', fontsize=12)
ax.set_title('Residuals', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_ylim(-1, 1)

# Plot 3: Tilt correction
ax = axes[1, 0]
ax.plot(z_smooth, corr_s * 1000, 'purple', linewidth=2.5)
ax.fill_between(z_smooth, 0, corr_s * 1000, alpha=0.2, color='purple')
ax.axhline(0, color='black', linewidth=1, linestyle='--')
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Tilt Correction Δμ (millimag)', fontsize=12)
ax.set_title('ESTIF Correction vs Redshift', fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)
ax.text(0.05, 0.95,
        f'α = {alpha_joint:.4f}\n'
        f'Δχ² = {chi2_joint_best-chi2_lcdm_best:+.2f}\n'
        f'Significance: {sigma_equiv:.2f}σ\n'
        f'N_SN = {N_sn}',
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Plot 4: Significance comparison
ax = axes[1, 1]
datasets   = ['Previous\n(580 SNe)', 'Pantheon+\n({} SNe)'.format(N_sn),
              'Predicted\n(3σ threshold)']
sigmas     = [1.78, sigma_equiv, 3.0]
colors     = ['orange',
              'green' if sigma_equiv >= 3 else
              ('orange' if sigma_equiv >= 2 else 'red'),
              'gray']
bars       = ax.bar(datasets, sigmas, color=colors, alpha=0.8,
                    edgecolor='black', linewidth=1.5)
ax.axhline(3.0, color='green',  linewidth=2, linestyle='--',
           label='3σ significance')
ax.axhline(2.0, color='orange', linewidth=1.5, linestyle=':',
           label='2σ significance')
ax.axhline(1.0, color='red',    linewidth=1, linestyle=':',
           label='1σ significance')
for bar, s in zip(bars, sigmas):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.05,
            f'{s:.2f}σ', ha='center', fontsize=11, fontweight='bold')
ax.set_ylabel('Statistical Significance (σ)', fontsize=12)
ax.set_title('Signal Strength vs Dataset Size', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0, max(sigmas) * 1.3)

plt.tight_layout()
plt.savefig('pantheon_plus_fit.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: pantheon_plus_fit.png")
print("=" * 70)
