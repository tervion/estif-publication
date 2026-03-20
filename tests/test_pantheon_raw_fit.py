"""
test_pantheon_raw_fit.py

Tests ESTIF against RAW Pantheon+ magnitudes — before bias corrections.

THE PREVIOUS PROBLEM:
MU_SH0ES is bias-corrected by the Pantheon+ team. Their corrections
are tuned to minimise residuals against ΛCDM. A smooth redshift-dependent
correction like ESTIF produces gets absorbed into their pipeline.
Result: best-fit α = 0.000, signal = 0.00σ.

THIS APPROACH:
Use the raw light curve parameters directly:
  mB  = raw observed peak magnitude
  x1  = light curve width (stretch)
  c   = colour

Apply the Tripp formula ourselves:
  μ_raw = mB - M_B + α_salt × x1 - β_salt × c

Then run the ESTIF fit on μ_raw. If the signal was buried in the
Pantheon+ corrections, it should reappear here.

TRIPP PARAMETERS (Pantheon+ best fit values):
  α_salt = 0.1545  (stretch correction)
  β_salt = 3.1     (colour correction)
  M_B    = -19.3   (absolute magnitude, fitted as free parameter)

IMPORTANT CAVEATS:
Without bias corrections the data will be noisier — selection effects
and survey-to-survey calibration offsets inflate scatter. We add a
conservative intrinsic scatter σ_int = 0.15 mag in quadrature.
The signal test is cleaner but the absolute χ² values are not
directly comparable to the corrected-data results.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import chi2 as chi2_dist
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Tripp parameters
# ============================================================================

ALPHA_SALT = 0.1545   # stretch correction (Pantheon+ best fit)
BETA_SALT  = 3.1      # colour correction  (Pantheon+ best fit)
SIGMA_INT  = 0.15     # intrinsic scatter (mag) — conservative

# ============================================================================
# Load raw Pantheon+ data
# ============================================================================

LOCAL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '../data/pantheon_plus.dat')

if not os.path.exists(LOCAL_PATH):
    print(f"ERROR: {LOCAL_PATH} not found.")
    print("Run test_pantheon_plus_fit.py first to download the data.")
    sys.exit(1)

print("=" * 70)
print("PANTHEON+ RAW FIT: TRIPP FORMULA (no bias corrections)")
print("=" * 70)

z_arr, mB_arr, x1_arr, c_arr, mB_err_arr, x1_err_arr, c_err_arr = \
    [], [], [], [], [], [], []

with open(LOCAL_PATH) as f:
    header = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == 'CID':
            header = parts
            z_col    = header.index('zHD')
            mB_col   = header.index('mB')
            mBe_col  = header.index('mBERR')
            x1_col   = header.index('x1')
            x1e_col  = header.index('x1ERR')
            c_col    = header.index('c')
            ce_col   = header.index('cERR')
            calib_col= header.index('IS_CALIBRATOR')
            ndof_col = header.index('NDOF')
            prob_col = header.index('FITPROB')
            continue
        if header is None:
            continue
        try:
            z       = float(parts[z_col])
            mB      = float(parts[mB_col])
            mB_err  = float(parts[mBe_col])
            x1      = float(parts[x1_col])
            x1_err  = float(parts[x1e_col])
            c       = float(parts[c_col])
            c_err   = float(parts[ce_col])
            is_calib= int(parts[calib_col])
            ndof    = int(parts[ndof_col])
            fitprob = float(parts[prob_col])

            # Quality cuts:
            # - Not a Cepheid calibrator
            # - z > 0.01 (avoid very local peculiar velocity dominated)
            # - Good light curve fit (fitprob > 0.001, ndof > 0)
            # - Reasonable colour and stretch ranges
            # - No bad values
            if (is_calib == 0 and z > 0.01 and
                    ndof > 0 and fitprob > 0.001 and
                    abs(x1) < 3.0 and abs(c) < 0.3 and
                    mB_err < 0.5 and mB > 0):
                z_arr.append(z)
                mB_arr.append(mB)
                mB_err_arr.append(mB_err)
                x1_arr.append(x1)
                x1_err_arr.append(x1_err)
                c_arr.append(c)
                c_err_arr.append(c_err)
        except (ValueError, IndexError):
            continue

z_arr    = np.array(z_arr)
mB_arr   = np.array(mB_arr)
x1_arr   = np.array(x1_arr)
c_arr    = np.array(c_arr)
mB_err_arr = np.array(mB_err_arr)

print(f"\n   Loaded {len(z_arr)} supernovae (after quality cuts)")
print(f"   Redshift range: {z_arr.min():.4f} – {z_arr.max():.4f}")
print(f"\n   Tripp formula: μ = mB - M_B + {ALPHA_SALT}×x1 - {BETA_SALT}×c")
print(f"   Intrinsic scatter: {SIGMA_INT} mag added in quadrature")

# ============================================================================
# Compute raw distance modulus
# ============================================================================

R_H          = const.c / const.H_0
r_universe_0 = 4.4e26
x_0          = R_H / r_universe_0

def tripp_mu(mB, x1, c, M_B):
    return mB - M_B + ALPHA_SALT * x1 - BETA_SALT * c

def tripp_sigma(mB_err):
    """Propagated uncertainty — simplified (dominant term is mB_err)."""
    return np.sqrt(mB_err**2 + SIGMA_INT**2)

sigma_arr = tripp_sigma(mB_err_arr)

def mu_lcdm(z, M=0.0):
    return estif.distance_modulus_lcdm(z) + M

def tilt_correction(z, alpha):
    x_z     = x_0 * (1.0 + z)**alpha
    obs_z   = estif.observable_combined(x_z)
    obs_now = estif.observable_combined(x_0)
    if obs_z <= 0 or obs_now <= 0:
        return 0.0
    return -2.5 * np.log10(obs_z / obs_now)

def tilt_corr_vec(z_a, alpha):
    return np.array([tilt_correction(z, alpha) for z in z_a])

# ============================================================================
# Step 1: ΛCDM fit (M_B free — absorbs absolute magnitude)
# ============================================================================

print(f"\n{'='*70}")
print("STEP 1: ΛCDM FIT (M_B free)")
print(f"{'='*70}")

def chi2_lcdm(M_B):
    mu_obs  = tripp_mu(mB_arr, x1_arr, c_arr, M_B)
    mu_pred = np.array([mu_lcdm(z) for z in z_arr])
    return np.sum(((mu_obs - mu_pred) / sigma_arr)**2)

res_lcdm       = minimize_scalar(chi2_lcdm, bounds=(-21, -17), method='bounded')
M_B_lcdm       = res_lcdm.x
chi2_lcdm_best = chi2_lcdm(M_B_lcdm)
dof_lcdm       = len(z_arr) - 1

print(f"\n   Best M_B       = {M_B_lcdm:.6f}")
print(f"   χ²             = {chi2_lcdm_best:.4f}")
print(f"   χ²/dof         = {chi2_lcdm_best/dof_lcdm:.6f}")

# ============================================================================
# Step 2: ESTIF joint fit (M_B and α free)
# ============================================================================

print(f"\n{'='*70}")
print("STEP 2: ESTIF JOINT FIT (M_B + α free)")
print(f"{'='*70}")

print(f"\n   Fitting... (this takes ~1 minute)")

def chi2_estif(params):
    M_B, alpha = params
    if alpha < -0.5 or alpha > 1.0:
        return 1e10
    mu_obs  = tripp_mu(mB_arr, x1_arr, c_arr, M_B)
    mu_pred = (np.array([mu_lcdm(z) for z in z_arr])
               + tilt_corr_vec(z_arr, alpha))
    return np.sum(((mu_obs - mu_pred) / sigma_arr)**2)

# Grid search for starting point
best_start = [M_B_lcdm, 0.10]
best_grid  = 1e10
for alpha_try in [0.0, 0.05, 0.10, 0.15, 0.20]:
    for M_try in [M_B_lcdm - 0.1, M_B_lcdm, M_B_lcdm + 0.1]:
        c = chi2_estif([M_try, alpha_try])
        if c < best_grid:
            best_grid  = c
            best_start = [M_try, alpha_try]

res_joint   = minimize(chi2_estif, best_start,
                       method='Nelder-Mead',
                       options={'xatol':1e-8, 'fatol':1e-8,
                                'maxiter':50000})
M_B_joint, alpha_joint = res_joint.x
chi2_joint_best = chi2_estif([M_B_joint, alpha_joint])
dof_joint   = len(z_arr) - 2

print(f"\n   Best M_B       = {M_B_joint:.6f}")
print(f"   Best α         = {alpha_joint:.6f}")
print(f"   χ²             = {chi2_joint_best:.4f}")
print(f"   χ²/dof         = {chi2_joint_best/dof_joint:.6f}")
print(f"   Δχ²            = {chi2_joint_best - chi2_lcdm_best:+.4f}")

# ============================================================================
# Step 3: Statistical significance
# ============================================================================

delta_chi2  = chi2_lcdm_best - chi2_joint_best
p_value     = chi2_dist.sf(delta_chi2, 1)
sigma_equiv = np.sqrt(chi2_dist.isf(max(p_value, 1e-15), 1))

print(f"\n{'='*70}")
print("STEP 3: STATISTICAL SIGNIFICANCE")
print(f"{'='*70}")
print(f"\n   Δχ² = {delta_chi2:.4f}")
print(f"   p-value = {p_value:.6f}")
print(f"   Significance: {sigma_equiv:.2f}σ")

if p_value < 0.003:
    verdict = "✅ SIGNIFICANT (>3σ)"
elif p_value < 0.05:
    verdict = f"✅ SIGNIFICANT ({sigma_equiv:.1f}σ)"
elif p_value < 0.10:
    verdict = f"⚠️  MARGINAL ({sigma_equiv:.1f}σ)"
else:
    verdict = "❌ NOT SIGNIFICANT"
print(f"   {verdict}")

# ============================================================================
# Step 4: Residual analysis
# ============================================================================

mu_obs_best = tripp_mu(mB_arr, x1_arr, c_arr, M_B_lcdm)
mu_obs_estif= tripp_mu(mB_arr, x1_arr, c_arr, M_B_joint)
mu_pred_lcdm_arr  = np.array([mu_lcdm(z) for z in z_arr])
mu_pred_estif_arr = mu_pred_lcdm_arr + tilt_corr_vec(z_arr, alpha_joint)

res_lcdm_arr  = mu_obs_best  - mu_pred_lcdm_arr
res_estif_arr = mu_obs_estif - mu_pred_estif_arr

z_bins      = np.percentile(z_arr, [0, 25, 50, 75, 100])
bin_labels  = ['low-z', 'mid-low', 'mid-high', 'high-z']

print(f"\n{'='*70}")
print("STEP 4: RESIDUALS BY REDSHIFT BIN")
print(f"{'='*70}")
print(f"\n   {'Bin':<12} {'z range':<16} {'ΛCDM':<10} {'ESTIF':<10} {'better?'}")
print("   " + "-"*55)
improvements = 0
for i in range(4):
    mask = (z_arr >= z_bins[i]) & (z_arr <= z_bins[i+1])
    ml = np.mean(res_lcdm_arr[mask])
    me = np.mean(res_estif_arr[mask])
    better = abs(me) < abs(ml)
    if better: improvements += 1
    print(f"   {bin_labels[i]:<12} {z_bins[i]:.3f}–{z_bins[i+1]:.3f}    "
          f"{ml:>+8.4f}  {me:>+8.4f}  {'✅' if better else '⚠️'}")

# ============================================================================
# Verdict
# ============================================================================

print(f"\n{'='*70}")
print("OVERALL VERDICT")
print(f"{'='*70}")
print(f"""
   Raw Pantheon+ ({len(z_arr)} SNe after quality cuts)
   Using Tripp formula — NO bias corrections from Pantheon+ pipeline

   ΛCDM:   χ²/dof = {chi2_lcdm_best/dof_lcdm:.4f}
   ESTIF:  χ²/dof = {chi2_joint_best/dof_joint:.4f}
   Δχ² = {delta_chi2:+.4f}
   α_best = {alpha_joint:.4f}
   Significance: {sigma_equiv:.2f}σ

   Residuals improved in {improvements}/4 bins

   {verdict}

   INTERPRETATION:
   {'Signal emerged from raw data — Pantheon+ corrections were suppressing it.' if delta_chi2 > 2 else
    'Signal still absent in raw data — correction pipeline not the issue.' if delta_chi2 < 0.5 else
    'Marginal signal in raw data — needs further investigation.'}
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle(f'Raw Pantheon+ Tripp Fit: ESTIF α={alpha_joint:.3f} '
             f'[{sigma_equiv:.2f}σ]  — No Bias Corrections',
             fontsize=13, fontweight='bold')

z_smooth   = np.linspace(0.01, z_arr.max()*1.02, 300)
mu_lcdm_s  = np.array([mu_lcdm(z, 0) for z in z_smooth])
mu_estif_s = mu_lcdm_s + np.array([tilt_correction(z, alpha_joint)
                                     for z in z_smooth])
corr_s     = np.array([tilt_correction(z, alpha_joint) for z in z_smooth])

# Plot 1: Hubble diagram
ax = axes[0, 0]
ax.errorbar(z_arr, mu_obs_best, yerr=sigma_arr, fmt='.', color='gray',
            alpha=0.15, markersize=2, elinewidth=0.3, zorder=1)
ax.plot(z_smooth, mu_lcdm_s + M_B_lcdm,  'blue', linewidth=2.5,
        label='ΛCDM', zorder=3)
ax.plot(z_smooth, mu_estif_s + M_B_joint, 'red',  linewidth=2.5,
        linestyle='--', label=f'ESTIF (α={alpha_joint:.3f})', zorder=3)
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Distance Modulus μ', fontsize=12)
ax.set_title('Hubble Diagram (raw Tripp μ)', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# Plot 2: Residuals
ax = axes[0, 1]
ax.scatter(z_arr, res_lcdm_arr,  color='blue', alpha=0.15, s=3,
           label=f'ΛCDM (χ²/dof={chi2_lcdm_best/dof_lcdm:.3f})')
ax.scatter(z_arr, res_estif_arr, color='red',  alpha=0.15, s=3,
           label=f'ESTIF (χ²/dof={chi2_joint_best/dof_joint:.3f})')
ax.axhline(0, color='black', linewidth=1.5)
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Residual (mag)', fontsize=12)
ax.set_title('Residuals', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_ylim(-1.5, 1.5)

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
        f'Δχ² = {delta_chi2:+.2f}\n'
        f'Significance: {sigma_equiv:.2f}σ\n'
        f'N = {len(z_arr)} SNe (raw)',
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Plot 4: χ² scan vs α
ax = axes[1, 1]
alpha_scan = np.linspace(-0.1, 0.4, 40)
chi2_scan  = []
for a in alpha_scan:
    def f(M): return chi2_estif([M, a])
    r = minimize_scalar(f, bounds=(-21, -17), method='bounded')
    chi2_scan.append(r.fun)
chi2_scan = np.array(chi2_scan)

ax.plot(alpha_scan, chi2_scan, 'purple', linewidth=2.5)
ax.axhline(chi2_lcdm_best, color='blue', linewidth=2, linestyle='--',
           label=f'ΛCDM χ²={chi2_lcdm_best:.1f}')
ax.axvline(alpha_joint, color='red', linewidth=2, linestyle=':',
           label=f'Best α={alpha_joint:.4f}')
ax.fill_between(alpha_scan, chi2_lcdm_best, chi2_scan,
                where=(chi2_scan < chi2_lcdm_best),
                alpha=0.2, color='green', label='Improvement over ΛCDM')
ax.set_xlabel('Tilt exponent α', fontsize=12)
ax.set_ylabel('χ²', fontsize=12)
ax.set_title('χ² vs α (raw data)', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('pantheon_raw_fit.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: pantheon_raw_fit.png")
print("=" * 70)
