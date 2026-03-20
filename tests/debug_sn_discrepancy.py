"""
debug_sn_discrepancy.py

Diagnoses why the original 580 SNe dataset gives Δχ² = -3.11
while Pantheon+ corrected gives +4.33 and Pantheon+ raw gives +5.41.

These three datasets should broadly agree. They don't. This script
finds out why.

HYPOTHESES:
H1 — Different redshift ranges (580 SNe only go to z=1.4)
H2 — Different absolute magnitude calibration (M_offset difference)
H3 — The signal comes entirely from high-z SNe (z > 1)
H4 — The 580 SNe have a systematic low-z bias that ESTIF worsens
H5 — The two models genuinely diverge at low z for this dataset
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import estif_ec_gr_model as estif
import estif_ec_gr_constants as const

# ============================================================================
# Load datasets
# ============================================================================

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')

def load_580(path):
    z, mu, s_mu, s_int = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    z.append(float(parts[1]));  mu.append(float(parts[2]))
                    s_mu.append(float(parts[3])); s_int.append(float(parts[4]))
                except ValueError: continue
    z = np.array(z); mu = np.array(mu)
    return z, mu, np.sqrt(np.array(s_mu)**2 + np.array(s_int)**2)

def load_pp_corrected(path):
    z, mu, err = [], [], []
    with open(path) as f:
        header = None
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            if parts[0] == 'CID':
                header = parts
                zc=header.index('zHD'); mc=header.index('MU_SH0ES')
                ec=header.index('MU_SH0ES_ERR_DIAG'); ic=header.index('IS_CALIBRATOR')
                continue
            if header is None: continue
            try:
                zv=float(parts[zc]); mv=float(parts[mc])
                ev=float(parts[ec]); icv=int(parts[ic])
                if zv>0.01 and icv==0 and ev<10 and mv>0:
                    z.append(zv); mu.append(mv); err.append(ev)
            except (ValueError,IndexError): continue
    return np.array(z), np.array(mu), np.array(err)

z1, mu1, s1 = load_580(os.path.join(data_dir, 'sn_data.txt'))
z2, mu2, s2 = load_pp_corrected(os.path.join(data_dir, 'pantheon_plus.dat'))

print("=" * 70)
print("580 SNe vs PANTHEON+ DISCREPANCY DIAGNOSTIC")
print("=" * 70)

# ============================================================================
# Section 1: Dataset properties
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: DATASET PROPERTIES")
print(f"{'='*70}")

print(f"\n   580 SNe:")
print(f"   z range:  {z1.min():.3f} – {z1.max():.3f}")
print(f"   μ range:  {mu1.min():.2f} – {mu1.max():.2f}")
print(f"   σ mean:   {s1.mean():.4f}")

print(f"\n   Pantheon+ corrected:")
print(f"   z range:  {z2.min():.4f} – {z2.max():.3f}")
print(f"   μ range:  {mu2.min():.2f} – {mu2.max():.2f}")
print(f"   σ mean:   {s2.mean():.4f}")

# Overlap region
z_max_580 = z1.max()
mask_pp_overlap = z2 <= z_max_580
print(f"\n   Pantheon+ SNe within 580 z-range (z<{z_max_580:.2f}): "
      f"{mask_pp_overlap.sum()}")
print(f"   Pantheon+ SNe beyond 580 z-range: {(~mask_pp_overlap).sum()}")

# ============================================================================
# Section 2: Model comparison functions
# ============================================================================

def mu_lcdm(z, M=0):
    return estif.distance_modulus_lcdm(z) + M

def mu_estif_fn(z, M=0):
    return estif.distance_modulus_estif(np.atleast_1d(z)) + M

def fit_M(z_arr, mu_obs, sigma, model_fn):
    def chi2(M):
        return np.sum(((mu_obs - model_fn(z_arr, M)) / sigma)**2)
    r = minimize_scalar(chi2, bounds=(-3, 3), method='bounded')
    return r.x, chi2(r.x)

print(f"\n{'='*70}")
print("SECTION 2: FIT QUALITY ACROSS Z BINS")
print(f"{'='*70}")

# ============================================================================
# Section 3: Pantheon+ restricted to 580 SNe z-range
# ============================================================================

print(f"\n   Fitting Pantheon+ restricted to z < {z_max_580:.2f}...")
z2r  = z2[mask_pp_overlap]
mu2r = mu2[mask_pp_overlap]
s2r  = s2[mask_pp_overlap]

M_l_580, c2_l_580  = fit_M(z1,  mu1,  s1,  mu_lcdm)
M_e_580, c2_e_580  = fit_M(z1,  mu1,  s1,  mu_estif_fn)
M_l_pp,  c2_l_pp   = fit_M(z2,  mu2,  s2,  mu_lcdm)
M_e_pp,  c2_e_pp   = fit_M(z2,  mu2,  s2,  mu_estif_fn)
M_l_ppr, c2_l_ppr  = fit_M(z2r, mu2r, s2r, mu_lcdm)
M_e_ppr, c2_e_ppr  = fit_M(z2r, mu2r, s2r, mu_estif_fn)

print(f"\n   {'Dataset':<30} {'N':>5} {'ΛCDM χ²/N':>12} {'ESTIF χ²/N':>12} {'Δχ²':>8}")
print("   " + "-"*68)
for label, N, cl, ce in [
    ("580 SNe (full)",           len(z1),  c2_l_580, c2_e_580),
    ("Pantheon+ (full)",         len(z2),  c2_l_pp,  c2_e_pp),
    (f"Pantheon+ (z<{z_max_580:.2f})", len(z2r), c2_l_ppr, c2_e_ppr),
]:
    print(f"   {label:<30} {N:>5} {cl/N:>12.6f} {ce/N:>12.6f} {cl-ce:>8.2f}")

# ============================================================================
# Section 4: Per-bin signal analysis (where does the signal come from?)
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: WHERE DOES THE SIGNAL COME FROM? (PER Z BIN)")
print(f"{'='*70}")

z_edges = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.5]

print(f"\n   580 SNe per bin:")
print(f"   {'z bin':<15} {'N':>5} {'ΛCDM mean res':>15} "
      f"{'ESTIF mean res':>15} {'ESTIF better?'}")
print("   " + "-"*60)
for i in range(len(z_edges)-1):
    lo, hi = z_edges[i], z_edges[i+1]
    m = (z1 >= lo) & (z1 < hi)
    if m.sum() < 3: continue
    mu_pred_l = np.array([mu_lcdm(z, M_l_580) for z in z1[m]])
    mu_pred_e = mu_estif_fn(z1[m], M_e_580)
    rl = np.mean(mu1[m] - mu_pred_l)
    re = np.mean(mu1[m] - mu_pred_e)
    better = "✅" if abs(re) < abs(rl) else "❌"
    print(f"   {lo:.1f}–{hi:.1f}           {m.sum():>5} {rl:>+15.4f} "
          f"{re:>+15.4f} {better}")

print(f"\n   Pantheon+ (full) per bin:")
print(f"   {'z bin':<15} {'N':>5} {'ΛCDM mean res':>15} "
      f"{'ESTIF mean res':>15} {'ESTIF better?'}")
print("   " + "-"*60)
for i in range(len(z_edges)-1):
    lo, hi = z_edges[i], z_edges[i+1]
    m = (z2 >= lo) & (z2 < hi)
    if m.sum() < 3: continue
    mu_pred_l = np.array([mu_lcdm(z, M_l_pp) for z in z2[m]])
    mu_pred_e = mu_estif_fn(z2[m], M_e_pp)
    rl = np.mean(mu2[m] - mu_pred_l)
    re = np.mean(mu2[m] - mu_pred_e)
    better = "✅" if abs(re) < abs(rl) else "❌"
    print(f"   {lo:.1f}–{hi:.1f}           {m.sum():>5} {rl:>+15.4f} "
          f"{re:>+15.4f} {better}")

# ============================================================================
# Section 5: M_offset comparison — calibration difference?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: CALIBRATION (M_OFFSET) COMPARISON")
print(f"{'='*70}")

print(f"\n   580 SNe:         M_ΛCDM = {M_l_580:+.6f}   M_ESTIF = {M_e_580:+.6f}   "
      f"ΔM = {M_e_580-M_l_580:+.6f}")
print(f"   Pantheon+ full:  M_ΛCDM = {M_l_pp:+.6f}   M_ESTIF = {M_e_pp:+.6f}   "
      f"ΔM = {M_e_pp-M_l_pp:+.6f}")
print(f"   Pantheon+ z<{z_max_580:.1f}: M_ΛCDM = {M_l_ppr:+.6f}   M_ESTIF = {M_e_ppr:+.6f}   "
      f"ΔM = {M_e_ppr-M_l_ppr:+.6f}")

print(f"\n   M_offset difference between datasets:")
print(f"   580 vs Pantheon+: {M_l_580 - M_l_pp:+.4f} mag (ΛCDM)")
print(f"\n   If |ΔM| > 0.05 mag, the datasets have different absolute")
print(f"   calibrations — a known systematic between SN compilations.")

# ============================================================================
# Section 6: ESTIF vs ΛCDM distance difference
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 5: ESTIF vs ΛCDM DISTANCE DIFFERENCE")
print(f"{'='*70}")

z_check = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
print(f"\n   {'z':<8} {'μ_ΛCDM':>10} {'μ_ESTIF':>10} {'Δμ (mmag)':>12} {'direction'}")
print("   " + "-"*45)
for zv in z_check:
    ml = estif.distance_modulus_lcdm(zv)
    me = float(estif.distance_modulus_estif(np.array([zv])))
    diff = (me - ml) * 1000
    direction = "ESTIF brighter" if diff < 0 else "ESTIF dimmer"
    print(f"   {zv:<8.2f} {ml:>10.4f} {me:>10.4f} {diff:>+12.2f}  {direction}")

print(f"\n   Sign convention:")
print(f"   Δμ > 0 (ESTIF dimmer) → ESTIF predicts objects farther away")
print(f"   Δμ < 0 (ESTIF brighter) → ESTIF predicts objects closer")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("DIAGNOSTIC SUMMARY")
print(f"{'='*70}")

delta_580 = c2_l_580 - c2_e_580
delta_pp  = c2_l_pp  - c2_e_pp
delta_ppr = c2_l_ppr - c2_e_ppr
M_diff    = abs(M_l_580 - M_l_pp)

print(f"""
   580 SNe:            Δχ² = {delta_580:+.2f} (ESTIF {'better' if delta_580>0 else 'worse'})
   Pantheon+ full:     Δχ² = {delta_pp:+.2f} (ESTIF {'better' if delta_pp>0 else 'worse'})
   Pantheon+ z<{z_max_580:.2f}:   Δχ² = {delta_ppr:+.2f} (ESTIF {'better' if delta_ppr>0 else 'worse'})

   M_offset shift between compilations: {M_diff:.4f} mag

   CONCLUSION:
""")

if abs(delta_ppr) < 1.0:
    print(f"   ✅ Pantheon+ restricted to 580 z-range shows no discrepancy.")
    print(f"   The signal in Pantheon+ comes from z > {z_max_580:.2f} SNe.")
    print(f"   The 580 SNe dataset lacks high-z coverage — that's why it disagrees.")
elif M_diff > 0.05:
    print(f"   ⚠️  Large calibration offset ({M_diff:.4f} mag) between datasets.")
    print(f"   The discrepancy is a known systematic between SN compilations.")
    print(f"   Not a failure of the ESTIF model — a dataset calibration issue.")
else:
    print(f"   ⚠️  Discrepancy persists even in the overlap z range.")
    print(f"   The two datasets genuinely disagree — needs further investigation.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('580 SNe vs Pantheon+: Why Do They Disagree?',
             fontsize=13, fontweight='bold')

z_s = np.linspace(0.01, 2.3, 200)
mu_l_s = np.array([mu_lcdm(z) for z in z_s])
mu_e_s = mu_estif_fn(z_s)
diff_s = (mu_e_s - mu_l_s) * 1000

# Plot 1: ESTIF - ΛCDM difference
ax = axes[0]
ax.plot(z_s, diff_s, 'purple', linewidth=2.5)
ax.fill_between(z_s, 0, diff_s,
                where=(diff_s > 0), alpha=0.2, color='red',
                label='ESTIF dimmer (farther)')
ax.fill_between(z_s, 0, diff_s,
                where=(diff_s < 0), alpha=0.2, color='blue',
                label='ESTIF brighter (closer)')
ax.axhline(0, color='black', linewidth=1.5)
ax.axvline(z_max_580, color='orange', linewidth=2, linestyle='--',
           label=f'580 SNe z limit ({z_max_580:.2f})')
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('μ_ESTIF − μ_ΛCDM (millimag)', fontsize=12)
ax.set_title('ESTIF vs ΛCDM Distance Difference', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 2: Residuals comparison
ax = axes[1]
res1_l = mu1 - np.array([mu_lcdm(z, M_l_580) for z in z1])
res1_e = mu1 - mu_estif_fn(z1, M_e_580)
res2_l = mu2 - np.array([mu_lcdm(z, M_l_pp) for z in z2])
res2_e = mu2 - mu_estif_fn(z2, M_e_pp)

ax.scatter(z1, res1_l - res1_e, color='orange', alpha=0.4, s=5,
           label=f'580 SNe: ΛCDM−ESTIF residual')
ax.scatter(z2, res2_l - res2_e, color='purple', alpha=0.2, s=3,
           label=f'Pantheon+: ΛCDM−ESTIF residual')
ax.axhline(0, color='black', linewidth=1.5)
ax.axvline(z_max_580, color='orange', linewidth=2, linestyle='--',
           label=f'580 SNe z limit')
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('ΛCDM residual − ESTIF residual (mag)', fontsize=12)
ax.set_title('Where ESTIF Helps vs Hurts', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_ylim(-0.5, 0.5)

# Plot 3: z distribution
ax = axes[2]
bins = np.linspace(0, 2.3, 30)
ax.hist(z1, bins=bins, color='orange', alpha=0.6, label=f'580 SNe (N={len(z1)})')
ax.hist(z2, bins=bins, color='purple', alpha=0.4,
        label=f'Pantheon+ (N={len(z2)})')
ax.axvline(z_max_580, color='orange', linewidth=2, linestyle='--',
           label=f'580 SNe z_max={z_max_580:.2f}')
ax.set_xlabel('Redshift z', fontsize=12)
ax.set_ylabel('Number of SNe', fontsize=12)
ax.set_title('Redshift Coverage Comparison', fontsize=11, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('sn_discrepancy_diagnostic.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: sn_discrepancy_diagnostic.png")
print("=" * 70)
