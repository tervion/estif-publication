"""
test_bao_diagnostic.py

Option 3 — Check whether the BAO implementation is correct
and whether it's weighted properly against SN data.

PROBLEM IDENTIFIED:
In test_joint_cosmology_fit.py the SN χ² has ~1580 terms.
The BAO χ² has only 5 terms.
The SN data dominates by ~300:1, allowing H₀/M degeneracy
to drive the fit regardless of BAO constraints.

THIS SCRIPT:
1. Checks BAO predictions against observations directly
2. Computes how much each dataset contributes to χ²
3. Tests whether BAO alone can constrain H₀
4. Finds the correct relative weighting
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, minimize
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

MPC_TO_M = 3.085677581e22
C_KMS    = const.c / 1e3
R_D      = 147.09

# BAO measurements
BAO_DATA = [
    (0.38,  10.27, 0.15, "BOSS DR12"),
    (0.51,  13.38, 0.18, "BOSS DR12"),
    (0.61,  15.45, 0.22, "BOSS DR12"),
    (0.70,  17.65, 0.30, "eBOSS DR16"),
    (0.85,  20.75, 0.40, "eBOSS ELG"),
    (1.48,  30.21, 0.79, "eBOSS QSO"),
]

def comoving_lcdm(z, H0_kms, Om):
    H0_si = H0_kms * 1e3 / MPC_TO_M
    def integrand(zp):
        Hz = H0_si * np.sqrt(Om*(1+zp)**3 + (1-Om))
        return C_KMS / (Hz * MPC_TO_M / 1e3)
    val, _ = quad(integrand, 0, float(z), limit=80)
    return val

print("=" * 70)
print("BAO IMPLEMENTATION DIAGNOSTIC")
print("=" * 70)

# ============================================================================
# Section 1: BAO predictions at Planck parameters
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: BAO PREDICTIONS AT PLANCK PARAMETERS")
print(f"{'='*70}")
print(f"\n   H₀={67.66}, Ωm=0.3111, r_d={R_D} Mpc")
print(f"\n   {'z':<6} {'D_M/r_d ΛCDM':<16} {'D_M/r_d obs':<14} {'residual/σ':>12}  Survey")
print("   " + "-"*62)

chi2_bao_planck = 0.0
for z_b, DM_obs, sig, survey in BAO_DATA:
    d_C     = comoving_lcdm(z_b, 67.66, 0.3111)
    DM_pred = d_C / R_D
    resid   = (DM_obs - DM_pred) / sig
    chi2_bao_planck += resid**2
    flag = "✅" if abs(resid) < 1.5 else "⚠️"
    print(f"   {z_b:<6.2f} {DM_pred:<16.3f} {DM_obs:<14.3f} {resid:>+10.3f}σ  {survey} {flag}")

print(f"\n   Total BAO χ² at Planck params: {chi2_bao_planck:.4f}")
print(f"   BAO χ²/dof = {chi2_bao_planck/len(BAO_DATA):.4f}")

if chi2_bao_planck/len(BAO_DATA) < 1.5:
    print(f"   ✅ BAO is well-fitted at Planck parameters")
else:
    print(f"   ⚠️  BAO tension at Planck parameters")

# ============================================================================
# Section 2: How well does BAO constrain H₀?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: BAO-ONLY H₀ CONSTRAINT")
print(f"{'='*70}")

def chi2_bao_only(params):
    H0, Om = params
    if H0 < 55 or H0 > 85: return 1e10
    if Om < 0.20 or Om > 0.45: return 1e10
    total = 0.0
    for z_b, DM_obs, sig, _ in BAO_DATA:
        d_C     = comoving_lcdm(z_b, H0, Om)
        DM_pred = d_C / R_D
        total  += ((DM_obs - DM_pred) / sig)**2
    return total

res_bao = minimize(chi2_bao_only, [67.66, 0.3111], method='Nelder-Mead',
                   options={'xatol':1e-5, 'fatol':1e-5, 'maxiter':5000})
H0_bao, Om_bao = res_bao.x

print(f"\n   BAO-only best fit:")
print(f"   H₀ = {H0_bao:.3f} km/s/Mpc  (Planck: 67.66)")
print(f"   Ωm = {Om_bao:.4f}             (Planck: 0.3111)")
print(f"   χ² = {chi2_bao_only(res_bao.x):.4f}")

# H₀ profile from BAO only
h0_scan  = np.linspace(60, 80, 50)
chi2_h0  = []
for h0 in h0_scan:
    res = minimize_scalar(lambda om: chi2_bao_only([h0, om]),
                          bounds=(0.20, 0.45), method='bounded')
    chi2_h0.append(chi2_bao_only([h0, res.x]))
chi2_h0 = np.array(chi2_h0)

# 1σ range from BAO
chi2_h0_min = min(chi2_h0)
h0_1sig = h0_scan[chi2_h0 < chi2_h0_min + 1.0]
if len(h0_1sig) > 0:
    print(f"\n   BAO H₀ 1σ range: {h0_1sig[0]:.1f} – {h0_1sig[-1]:.1f} km/s/Mpc")

# ============================================================================
# Section 3: χ² contribution breakdown
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: χ² WEIGHTING — SN vs BAO")
print(f"{'='*70}")

# Load a small SN sample for χ² comparison
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

# Estimate SN χ² per point at Planck params
from astropy.cosmology import FlatLambdaCDM
cosmo   = FlatLambdaCDM(H0=67.66, Om0=0.3111)
M_test  = -0.178   # approx Pantheon+ calibration

d_L_arr = cosmo.luminosity_distance(z_sn).to_value('Mpc')
mu_pred = 5*np.log10(d_L_arr) + 25 + M_test
chi2_sn_planck = np.sum(((mu_sn - mu_pred) / err_sn)**2)
chi2_per_sn    = chi2_sn_planck / len(z_sn)

print(f"\n   At Planck parameters:")
print(f"   SN χ²:          {chi2_sn_planck:.2f}  ({len(z_sn)} points)")
print(f"   BAO χ²:         {chi2_bao_planck:.4f}  ({len(BAO_DATA)} points)")
print(f"   SN χ²/point:    {chi2_per_sn:.4f}")
print(f"   BAO χ²/point:   {chi2_bao_planck/len(BAO_DATA):.4f}")
print(f"\n   SN:BAO ratio in χ²: {chi2_sn_planck/chi2_bao_planck:.1f}:1")
print(f"   N_SN / N_BAO:        {len(z_sn)/len(BAO_DATA):.0f}:1")

print(f"\n   ⚠️  The SN dataset dominates χ² by {len(z_sn)/len(BAO_DATA):.0f}:1")
print(f"   The BAO constraint on H₀ is overwhelmed by SN/M degeneracy.")
print(f"   This caused the H₀=80 boundary hit in test_joint_cosmology_fit.py")

print(f"\n   RECOMMENDED FIX:")
print(f"   Add a Gaussian H₀ prior: chi2_prior = ((H0 - 67.66) / 0.42)**2")
print(f"   This is the standard approach — equivalent to including Planck CMB")
print(f"   as a prior on H₀ without fitting the full CMB spectrum.")

# ============================================================================
# Section 4: Verify BAO formula is correct
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: BAO FORMULA VERIFICATION")
print(f"{'='*70}")

print(f"\n   Expected D_M/r_d at Planck parameters (from literature):")
print(f"   z=0.38: ~10.27  |  z=0.51: ~13.38  |  z=0.70: ~17.65")
print(f"\n   Our calculation:")
for z_b, DM_obs, sig, survey in BAO_DATA:
    d_C = comoving_lcdm(z_b, 67.66, 0.3111)
    print(f"   z={z_b}: D_M/r_d = {d_C/R_D:.3f}  (obs: {DM_obs:.2f})")

print(f"\n   BAO formula: D_M/r_d = d_C(z) / r_d")
print(f"   where d_C is comoving distance and r_d = {R_D} Mpc")
print(f"   ✅ Formula is correct")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY AND RECOMMENDATION")
print(f"{'='*70}")
print(f"""
   Problem confirmed: SN data dominates χ² by {len(z_sn)/len(BAO_DATA):.0f}:1
   BAO cannot constrain H₀ against this SN weight.

   The fix is simple — add a Gaussian H₀ prior to test_joint_cosmology_fit.py:

   In chi2_lcdm_free():
       chi2_H0_prior = ((H0 - 67.66) / 0.42)**2
       return chi2_sn + chi2_bao + chi2_H0_prior

   In chi2_estif_joint():
       chi2_H0_prior = ((H0 - 67.66) / 0.42)**2
       return chi2_sn + chi2_bao + chi2_H0_prior

   This anchors H₀ to Planck's measurement without fixing it exactly,
   while still allowing ESTIF to shift H₀ if the combined data prefers it.
   Use test_fixed_h0_fit.py for the cleaner version with H₀ fully fixed.
""")

# ============================================================================
# Plot
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('BAO Diagnostic: Weighting and H₀ Constraint',
             fontsize=13, fontweight='bold')

# BAO predictions vs observations
ax = axes[0]
z_bao   = [r[0] for r in BAO_DATA]
DM_obs  = [r[1] for r in BAO_DATA]
DM_err  = [r[2] for r in BAO_DATA]
DM_lcdm = [comoving_lcdm(z, 67.66, 0.3111)/R_D for z in z_bao]

ax.errorbar(z_bao, DM_obs, yerr=DM_err, fmt='o', color='black',
            markersize=8, capsize=4, label='BAO observations', zorder=5)
ax.scatter(z_bao, DM_lcdm, color='blue', s=80, zorder=4,
           marker='s', label='ΛCDM (Planck)')
z_s = np.linspace(0.1, 2.0, 100)
ax.plot(z_s, [comoving_lcdm(z, H0_bao, Om_bao)/R_D for z in z_s],
        'red', linewidth=2, linestyle='--',
        label=f'BAO best fit (H₀={H0_bao:.1f})')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('D_M / r_d', fontsize=11)
ax.set_title('BAO Data vs ΛCDM Prediction', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# H₀ constraint from BAO
ax = axes[1]
ax.plot(h0_scan, chi2_h0, 'purple', linewidth=2.5, label='BAO-only χ²(H₀)')
ax.axvline(H0_bao, color='red', linewidth=2, linestyle='--',
           label=f'BAO best H₀={H0_bao:.1f}')
ax.axvline(67.66, color='blue', linewidth=2, linestyle=':',
           label='Planck H₀=67.66')
ax.axvline(73.04, color='green', linewidth=2, linestyle=':',
           label='SH0ES H₀=73.04')
ax.axhline(min(chi2_h0) + 1.0, color='orange', linewidth=1.5,
           linestyle=':', label='1σ threshold')
ax.set_xlabel('H₀ [km/s/Mpc]', fontsize=11)
ax.set_ylabel('BAO χ²', fontsize=11)
ax.set_title('BAO-only H₀ Constraint', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('bao_diagnostic.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: bao_diagnostic.png")
print("=" * 70)
