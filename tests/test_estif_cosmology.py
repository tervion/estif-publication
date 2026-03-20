"""
test_estif_cosmology.py

Full comparison of ESTIF Option A cosmology vs ΛCDM against
both supernova datasets.

ESTIF Option A:
    H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
    Ω_tilt(z) = Ω_Λ × (obs(x(z)) / obs(x_0))²

This is now a STANDALONE cosmological model — not a perturbative
correction to ΛCDM, but a replacement of the dark energy sector
with geometry derived from the tilt formula.

THREE TESTS:
1. Original 580 SNe dataset
2. Pantheon+ bias-corrected (1580 SNe)
3. Pantheon+ raw Tripp magnitudes (1443 SNe)

For each: fit M_offset only (same degrees of freedom as ΛCDM)
then compare χ², residuals, and significance.
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
# Helper: load original 580 SNe
# ============================================================================

def load_580(path):
    z, mu, s_mu, s_int = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    z.append(float(parts[1]))
                    mu.append(float(parts[2]))
                    s_mu.append(float(parts[3]))
                    s_int.append(float(parts[4]))
                except ValueError:
                    continue
    z  = np.array(z);  mu = np.array(mu)
    st = np.sqrt(np.array(s_mu)**2 + np.array(s_int)**2)
    return z, mu, st

# ============================================================================
# Helper: load Pantheon+ corrected
# ============================================================================

def load_pantheon_corrected(path):
    z, mu, err = [], [], []
    with open(path) as f:
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
                    z.append(zv); mu.append(mv); err.append(ev)
            except (ValueError, IndexError):
                continue
    return np.array(z), np.array(mu), np.array(err)

# ============================================================================
# Helper: load Pantheon+ raw Tripp
# ============================================================================

ALPHA_SALT = 0.1545
BETA_SALT  = 3.1
SIGMA_INT  = 0.15

def load_pantheon_raw(path):
    z, mu, err = [], [], []
    with open(path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split()
            if parts[0] == 'CID':
                header = parts
                zc   = header.index('zHD')
                mBc  = header.index('mB')
                mBec = header.index('mBERR')
                x1c  = header.index('x1')
                cc   = header.index('c')
                isc  = header.index('IS_CALIBRATOR')
                ndc  = header.index('NDOF')
                prc  = header.index('FITPROB')
                continue
            if header is None: continue
            try:
                zv   = float(parts[zc])
                mBv  = float(parts[mBc])
                mBev = float(parts[mBec])
                x1v  = float(parts[x1c])
                cv   = float(parts[cc])
                icv  = int(parts[isc])
                ndv  = int(parts[ndc])
                prv  = float(parts[prc])
                if (icv == 0 and zv > 0.01 and ndv > 0 and prv > 0.001
                        and abs(x1v) < 3 and abs(cv) < 0.3
                        and mBev < 0.5 and mBv > 0):
                    mu_raw = mBv + ALPHA_SALT*x1v - BETA_SALT*cv
                    sigma  = np.sqrt(mBev**2 + SIGMA_INT**2)
                    z.append(zv); mu.append(mu_raw); err.append(sigma)
            except (ValueError, IndexError):
                continue
    return np.array(z), np.array(mu), np.array(err)

# ============================================================================
# Core fit functions
# ============================================================================

def mu_lcdm_vec(z_arr, M):
    return np.array([estif.distance_modulus_lcdm(z) + M for z in z_arr])

def mu_estif_vec(z_arr, M):
    return estif.distance_modulus_estif(np.asarray(z_arr)) + M

def fit_both(z_arr, mu_obs, sigma):
    """Fit M_offset for both ΛCDM and ESTIF. Returns results dict."""
    dof = len(z_arr) - 1

    # Detect if raw Tripp data (μ values around 10-15 without M_B offset)
    mu_median = np.median(mu_obs)
    bounds = (-21, -17) if mu_median < 30 else (-3, 3)

    # ΛCDM
    def c2_lcdm(M):
        return np.sum(((mu_obs - mu_lcdm_vec(z_arr, M)) / sigma)**2)
    rl     = minimize_scalar(c2_lcdm,  bounds=bounds, method='bounded')
    chi2_l = c2_lcdm(rl.x)

    # ESTIF
    def c2_estif(M):
        return np.sum(((mu_obs - mu_estif_vec(z_arr, M)) / sigma)**2)
    re     = minimize_scalar(c2_estif, bounds=bounds, method='bounded')
    chi2_e = c2_estif(re.x)

    delta = chi2_l - chi2_e
    p_val = chi2_dist.sf(max(delta, 0), 1)
    sig   = np.sqrt(chi2_dist.isf(max(p_val, 1e-15), 1)) if delta > 0 else 0.0

    # Compute residual arrays using best-fit M values
    res_lcdm_arr  = mu_obs - mu_lcdm_vec(z_arr,  rl.x)
    res_estif_arr = mu_obs - mu_estif_vec(z_arr,  re.x)

    # Residuals by z bin
    z_bins  = np.percentile(z_arr, [0, 25, 50, 75, 100])
    bin_imp = 0
    for i in range(4):
        lo, hi = z_bins[i], z_bins[i+1]
        mask = (z_arr >= lo) & (z_arr <= hi)
        if abs(np.mean(res_estif_arr[mask])) < abs(np.mean(res_lcdm_arr[mask])):
            bin_imp += 1

    return {
        'M_lcdm':    rl.x,        'chi2_lcdm':  chi2_l,
        'M_estif':   re.x,        'chi2_estif': chi2_e,
        'dof':       dof,         'delta_chi2': delta,
        'p_value':   p_val,       'sigma':      sig,
        'res_lcdm':  res_lcdm_arr,'res_estif':  res_estif_arr,
        'bin_imp':   bin_imp,
    }

# ============================================================================
# Run all three tests
# ============================================================================

print("=" * 70)
print("ESTIF OPTION A COSMOLOGY — FULL COMPARISON")
print("=" * 70)
print(f"\n   H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]")
print(f"   Ω_tilt(z) = Ω_Λ × (obs(x(z)) / obs(x_0))²")
print(f"   x(z) = x_0 × (1+z)^{estif.ALPHA_COSMO}")

data_dir  = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')
sn_path   = os.path.join(data_dir, 'sn_data.txt')
pp_path   = os.path.join(data_dir, 'pantheon_plus.dat')

datasets  = []

# Dataset 1 — 580 SNe
if os.path.exists(sn_path):
    z1, mu1, s1 = load_580(sn_path)
    datasets.append(('Original 580 SNe', z1, mu1, s1))
    print(f"\n   ✅ Loaded original dataset: {len(z1)} SNe")
else:
    print(f"\n   ⚠️  {sn_path} not found — skipping")

# Dataset 2 — Pantheon+ corrected
if os.path.exists(pp_path):
    z2, mu2, s2 = load_pantheon_corrected(pp_path)
    datasets.append(('Pantheon+ corrected', z2, mu2, s2))
    print(f"   ✅ Loaded Pantheon+ corrected: {len(z2)} SNe")

    # Dataset 3 — Pantheon+ raw
    z3, mu3, s3 = load_pantheon_raw(pp_path)
    datasets.append(('Pantheon+ raw Tripp', z3, mu3, s3))
    print(f"   ✅ Loaded Pantheon+ raw: {len(z3)} SNe")
else:
    print(f"\n   ⚠️  Pantheon+ data not found.")
    print(f"   Run test_pantheon_plus_fit.py first to download it.")

print(f"\n   Fitting all datasets (ESTIF comoving distance is slow — ~2 min)...")

results = {}
for name, z_d, mu_d, sig_d in datasets:
    print(f"\n   [{name}] fitting...", flush=True)
    r = fit_both(z_d, mu_d, sig_d)
    results[name] = (z_d, mu_d, sig_d, r)

# ============================================================================
# Print results table
# ============================================================================

print(f"\n{'='*70}")
print("RESULTS SUMMARY")
print(f"{'='*70}")

print(f"\n   {'Dataset':<24} {'N':>6} {'ΛCDM χ²/dof':>13} "
      f"{'ESTIF χ²/dof':>13} {'Δχ²':>8} {'σ':>6} {'bins':>5}")
print("   " + "-"*75)

for name, (z_d, mu_d, sig_d, r) in results.items():
    chi2_dof_l = r['chi2_lcdm']  / r['dof']
    chi2_dof_e = r['chi2_estif'] / r['dof']
    print(f"   {name:<24} {len(z_d):>6} {chi2_dof_l:>13.4f} "
          f"{chi2_dof_e:>13.4f} {r['delta_chi2']:>8.2f} "
          f"{r['sigma']:>6.2f} {r['bin_imp']:>4}/4")

# ============================================================================
# Detailed verdict per dataset
# ============================================================================

print(f"\n{'='*70}")
print("DETAILED VERDICT")
print(f"{'='*70}")

for name, (z_d, mu_d, sig_d, r) in results.items():
    print(f"\n   {name}:")
    print(f"   Δχ² = {r['delta_chi2']:+.2f}  "
          f"significance = {r['sigma']:.2f}σ  "
          f"residuals {r['bin_imp']}/4 bins improved")
    if r['sigma'] >= 3.0:
        print(f"   ✅ STATISTICALLY SIGNIFICANT (≥3σ)")
    elif r['sigma'] >= 2.0:
        print(f"   ✅ STRONG HINT (≥2σ)")
    elif r['sigma'] >= 1.0:
        print(f"   ⚠️  MARGINAL SIGNAL (≥1σ)")
    else:
        print(f"   ❌ NO SIGNAL")
    if name == 'Original 580 SNe':
        print(f"   ℹ️  Note: 0.118 mag calibration offset vs Pantheon+")
        print(f"   (known systematic between SN compilations — not a model failure)")

# ============================================================================
# Ω_tilt evolution table
# ============================================================================

print(f"\n{'='*70}")
print("Ω_TILT EVOLUTION vs REDSHIFT")
print(f"{'='*70}")
print(f"\n   {'z':<8} {'Ω_tilt(z)':<14} {'Ω_Λ':<14} {'ratio':<10}")
print("   " + "-"*44)
for z_val in [0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0]:
    ot = estif.omega_tilt(z_val)
    print(f"   {z_val:<8.1f} {ot:<14.6f} {estif.OMEGA_LAMBDA:<14.6f} "
          f"{ot/estif.OMEGA_LAMBDA:.6f}")

# ============================================================================
# Visualization
# ============================================================================

n_plots = len(results)
fig, axes = plt.subplots(n_plots, 3, figsize=(18, 5*n_plots))
if n_plots == 1:
    axes = axes[None, :]

fig.suptitle('ESTIF Option A Cosmology: Full Comparison vs ΛCDM',
             fontsize=14, fontweight='bold')

z_smooth = np.linspace(0.01, 2.3, 200)
mu_lcdm_s = np.array([estif.distance_modulus_lcdm(z) for z in z_smooth])
mu_estif_s = estif.distance_modulus_estif(z_smooth)
omega_tilt_s = np.array([estif.omega_tilt(z) for z in z_smooth])

for idx, (name, (z_d, mu_d, sig_d, r)) in enumerate(results.items()):
    ax_hub, ax_res, ax_om = axes[idx]

    # Hubble diagram
    ax_hub.errorbar(z_d, mu_d, yerr=sig_d, fmt='.', color='gray',
                    alpha=0.15, markersize=2, elinewidth=0.3)
    ax_hub.plot(z_smooth, mu_lcdm_s  + r['M_lcdm'],  'blue', linewidth=2,
                label='ΛCDM')
    ax_hub.plot(z_smooth, mu_estif_s + r['M_estif'], 'red',  linewidth=2,
                linestyle='--', label=f'ESTIF (Δχ²={r["delta_chi2"]:+.1f})')
    ax_hub.set_xlabel('Redshift z', fontsize=10)
    ax_hub.set_ylabel('μ', fontsize=10)
    ax_hub.set_title(f'{name} — Hubble Diagram', fontsize=10, fontweight='bold')
    ax_hub.legend(fontsize=9); ax_hub.grid(alpha=0.3)

    # Residuals
    ax_res.scatter(z_d, r['res_lcdm'],  color='blue', alpha=0.15, s=3,
                   label=f'ΛCDM  χ²/dof={r["chi2_lcdm"]/r["dof"]:.3f}')
    ax_res.scatter(z_d, r['res_estif'], color='red',  alpha=0.15, s=3,
                   label=f'ESTIF χ²/dof={r["chi2_estif"]/r["dof"]:.3f}')
    ax_res.axhline(0, color='black', linewidth=1)
    ax_res.set_xlabel('Redshift z', fontsize=10)
    ax_res.set_ylabel('Residual (mag)', fontsize=10)
    ax_res.set_title(f'{name} — Residuals ({r["sigma"]:.2f}σ)',
                     fontsize=10, fontweight='bold')
    ax_res.legend(fontsize=8); ax_res.grid(alpha=0.3)
    ax_res.set_ylim(-1, 1)

    # Ω_tilt evolution (same for all datasets)
    ax_om.plot(z_smooth, omega_tilt_s, 'purple', linewidth=2.5,
               label='Ω_tilt(z)')
    ax_om.axhline(estif.OMEGA_LAMBDA, color='blue', linewidth=2,
                  linestyle='--', label=f'Ω_Λ = {estif.OMEGA_LAMBDA:.4f}')
    ax_om.fill_between(z_smooth, estif.OMEGA_LAMBDA, omega_tilt_s,
                       alpha=0.2, color='purple',
                       label='Departure from Λ')
    ax_om.set_xlabel('Redshift z', fontsize=10)
    ax_om.set_ylabel('Dark energy density', fontsize=10)
    ax_om.set_title('Ω_tilt(z) Evolution', fontsize=10, fontweight='bold')
    ax_om.legend(fontsize=9); ax_om.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('estif_cosmology_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: estif_cosmology_comparison.png")
print("=" * 70)
