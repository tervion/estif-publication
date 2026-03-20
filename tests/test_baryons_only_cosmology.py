"""
test_baryons_only_cosmology.py

THE QUESTION:
Can the ESTIF 4D tilt geometry carry ALL non-baryonic energy density,
replacing both dark energy AND dark matter from the Friedmann equation?

CURRENT MODEL (ESTIF v4.0):
    H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
    Ωm = 0.3111  ← includes dark matter (0.262) + baryons (0.049)

THIS TEST:
    H²(z) = H₀² × [Ωb(1+z)³ + Ω_tilt_total(z)]
    Ωb = 0.049   ← baryons only, dark matter removed
    Ω_tilt_total carries EVERYTHING else = 0.951 today

KEY DIFFERENCE FROM CURRENT MODEL:
In the current model: Ω_tilt(z=0) = 0.6889 (dark energy only)
In this test:         Ω_tilt(z=0) = 0.9510 (dark energy + dark matter)

The tilt formula is the same — only the normalization and background change.

USING THE EXACT GEOMETRIC FORMULA (from test_alpha_from_geometry.py):
    x(z) = x_0 × (1+z) × H_0 / H_background(z)

where H_background uses baryons-only to avoid circular dependency.

TESTS:
1. Age of universe — must be ≥ 13.5 Gyr
2. Supernova distances — Δχ² vs ΛCDM
3. BAO scale — D_M/r_d comparison
4. Ω_tilt shape — does it look physically reasonable?
5. Compare to current ESTIF v4.0 (Ωm = 0.3111)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

MPC_TO_M  = 3.085677581e22
GYR_TO_SEC= 3.15576e16
C_KMS     = const.c / 1e3
H0_SI     = const.H_0
R_D       = 147.09   # Sound horizon [Mpc] — unchanged

# ============================================================================
# Cosmological parameters
# ============================================================================

# Standard ΛCDM
OM_LCDM  = 0.3111
OL_LCDM  = 0.6889

# Baryons only — this test
OB        = 0.049    # Ωb — measured from BBN + CMB acoustic peaks
OT_0      = 1.0 - OB  # Ω_tilt today must carry everything else = 0.951

# Pre-computed x_0 (same as current model — universe size hasn't changed)
X_0       = (const.c / H0_SI) / 4.4e26
OBS_0     = estif.observable_combined(X_0)

print("=" * 70)
print("BARYONS-ONLY COSMOLOGY — 4D TILT CARRIES ALL NON-BARYONIC ENERGY")
print("=" * 70)
print(f"\n   Current ESTIF:    Ωm = {OM_LCDM}  (dark matter + baryons)")
print(f"   This test:        Ωb = {OB}   (baryons only — BBN measured)")
print(f"   Dark matter:      removed — tilt geometry carries it")
print(f"   Ω_tilt today:     {OT_0:.4f}  (vs {OL_LCDM:.4f} in current model)")
print(f"   x_0:              {X_0:.6f}")
print(f"   obs_now:          {OBS_0:.6f}")

# ============================================================================
# Core functions — baryons-only background
# ============================================================================

def H_background_baryons(z):
    """
    Background H(z) using baryons only — used as the ruler for x(z).
    Avoids circular dependency.
    H_bg²(z) = H₀² × [Ωb(1+z)³ + Ω_tilt_0]
    (flat universe: Ωb + Ω_tilt_0 = 1)
    """
    z = np.asarray(z, dtype=float)
    return H0_SI * np.sqrt(OB*(1+z)**3 + OT_0)


def x_baryons(z):
    """
    Exact geometric x(z) using baryons-only background as ruler.
    x(z) = x_0 × (1+z) × H_0 / H_background(z)
    """
    z    = np.asarray(z, dtype=float)
    H_bg = H_background_baryons(z)
    return X_0 * (1.0 + z) * H0_SI / H_bg


def omega_tilt_baryons(z):
    """
    Ω_tilt for baryons-only universe.
    Normalized so Ω_tilt(0) = OT_0 = 0.951.
    """
    z    = np.asarray(z, dtype=float)
    x_z  = x_baryons(z)
    obs_z = estif.observable_combined(x_z)

    if np.isscalar(obs_z):
        if obs_z <= 0: return OT_0
    else:
        obs_z = np.where(obs_z <= 0, OBS_0, obs_z)

    return OT_0 * (OBS_0 / obs_z) ** 2


def H_baryons(z):
    """
    Full ESTIF H(z) with baryons only + tilt.
    H²(z) = H₀² × [Ωb(1+z)³ + Ω_tilt_baryons(z)]
    """
    z  = np.asarray(z, dtype=float)
    return H0_SI * np.sqrt(OB*(1+z)**3 + omega_tilt_baryons(z))


def H_lcdm(z):
    """Standard ΛCDM for comparison."""
    z = np.asarray(z, dtype=float)
    return H0_SI * np.sqrt(OM_LCDM*(1+z)**3 + OL_LCDM)


def comoving_dist(z, H_func):
    """Comoving distance [Mpc] for any H(z)."""
    def integrand(zp):
        return C_KMS / (H_func(zp) * MPC_TO_M / 1e3)
    val, _ = quad(integrand, 0, float(z), limit=100)
    return val


def mu_model(z, H_func, M):
    """Distance modulus."""
    d_C = comoving_dist(z, H_func)
    d_L = (1.0 + z) * d_C
    return 5.0 * np.log10(max(d_L, 1e-10)) + 25.0 + M

# ============================================================================
# Section 1: Ω_tilt shape
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: Ω_tilt EVOLUTION")
print(f"{'='*70}")

z_vals = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 5.0]
print(f"\n   {'z':<8} {'x(z)':<12} {'Ω_tilt':<14} {'ratio/Ω_tilt_0':<16} {'current ESTIF Ω_tilt'}")
print("   " + "-"*65)
for z in z_vals:
    ot  = omega_tilt_baryons(z)
    xz  = x_baryons(z)
    ot_current = estif.omega_tilt(z)
    print(f"   {z:<8.1f} {xz:<12.6f} {ot:<14.6f} {ot/OT_0:<16.6f} {ot_current:.6f}")

# ============================================================================
# Section 2: Age of universe
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: AGE OF UNIVERSE")
print(f"{'='*70}")

def age_gyr(H_func):
    def integrand(z):
        return 1.0 / ((1.0 + z) * H_func(z))
    result, _ = quad(integrand, 0, 1000, limit=200)
    return result / GYR_TO_SEC

t_lcdm    = age_gyr(H_lcdm)
t_current = age_gyr(estif.H_estif)
t_baryons = age_gyr(H_baryons)

print(f"\n   ΛCDM:              {t_lcdm:.3f} Gyr")
print(f"   Current ESTIF:     {t_current:.3f} Gyr  (Ωm=0.3111, exact formula)")
print(f"   Baryons-only ESTIF:{t_baryons:.3f} Gyr  (Ωb=0.049)")
print(f"   Hard limit:        13.500 Gyr (oldest observed stars)")

if t_baryons >= 13.5:
    print(f"\n   ✅ PASS — baryons-only universe old enough")
elif t_baryons >= 13.0:
    print(f"\n   ⚠️  MARGINAL — {13.5 - t_baryons:.3f} Gyr short")
else:
    print(f"\n   ❌ FAIL — {13.5 - t_baryons:.3f} Gyr too young")

# ============================================================================
# Section 3: BAO scale
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: BAO ANGULAR SCALE")
print(f"{'='*70}")

BAO_DATA = [
    (0.38,  10.27, 0.15, "BOSS DR12"),
    (0.51,  13.38, 0.18, "BOSS DR12"),
    (0.61,  15.45, 0.22, "BOSS DR12"),
    (0.70,  17.65, 0.30, "eBOSS DR16"),
    (0.85,  20.75, 0.40, "eBOSS ELG"),
    (1.48,  30.21, 0.79, "eBOSS QSO"),
]

print(f"\n   r_d = {R_D} Mpc (unchanged — early universe BBN physics)")
print(f"\n   {'z':<6} {'ΛCDM':<10} {'Current':<12} {'Baryons':<12} {'Obs':<10} {'σ_bary':>8}  Survey")
print("   " + "-"*70)

n_better_current = 0
n_better_baryons = 0
chi2_bao_baryons = 0.0
chi2_bao_lcdm    = 0.0

for z_b, DM_obs, sig, survey in BAO_DATA:
    dm_lcdm = comoving_dist(z_b, H_lcdm) / R_D
    dm_cur  = comoving_dist(z_b, estif.H_estif) / R_D
    dm_bar  = comoving_dist(z_b, H_baryons) / R_D

    resid_bary = (DM_obs - dm_bar) / sig
    resid_lcdm = (DM_obs - dm_lcdm) / sig
    chi2_bao_baryons += resid_bary**2
    chi2_bao_lcdm    += resid_lcdm**2

    if abs(DM_obs - dm_bar) < abs(DM_obs - dm_lcdm):
        n_better_baryons += 1
    if abs(DM_obs - dm_cur) < abs(DM_obs - dm_lcdm):
        n_better_current += 1

    flag = "✅" if abs(resid_bary) < 2 else "⚠️"
    print(f"   {z_b:<6.2f} {dm_lcdm:<10.3f} {dm_cur:<12.3f} {dm_bar:<12.3f} "
          f"{DM_obs:<10.2f} {resid_bary:>+7.2f}σ  {survey} {flag}")

print(f"\n   BAO χ²: ΛCDM={chi2_bao_lcdm:.3f}  Baryons-only={chi2_bao_baryons:.3f}")
print(f"   ΛCDM better at {len(BAO_DATA)-n_better_baryons}/{len(BAO_DATA)} redshifts")
print(f"   Baryons better at {n_better_baryons}/{len(BAO_DATA)} redshifts")

# ============================================================================
# Section 4: Supernova distances
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: SUPERNOVA DISTANCES")
print(f"{'='*70}")
print(f"   Loading Pantheon+ data...", flush=True)

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
print(f"   Loaded {len(z_sn)} SNe. Fitting M_offset for each model...", flush=True)

def best_M_and_chi2(H_func):
    """Find optimal M_offset and return chi2."""
    def chi2_fn(M):
        mu_pred = np.array([mu_model(z, H_func, M[0]) for z in z_sn])
        return np.sum(((mu_sn - mu_pred) / err_sn)**2)
    res = minimize(chi2_fn, [-0.18], method='Nelder-Mead',
                   options={'xatol':1e-6, 'fatol':1e-6, 'maxiter':3000})
    return res.x[0], res.fun

print(f"   Fitting ΛCDM...", flush=True)
M_lcdm, chi2_lcdm_sn = best_M_and_chi2(H_lcdm)
print(f"   Fitting current ESTIF...", flush=True)
M_cur, chi2_cur_sn   = best_M_and_chi2(estif.H_estif)
print(f"   Fitting baryons-only...", flush=True)
M_bar, chi2_bar_sn   = best_M_and_chi2(H_baryons)

dof = len(z_sn) - 1

print(f"\n   {'Model':<22} {'χ²':<12} {'χ²/dof':<12} {'Δχ² vs ΛCDM':<14} {'σ'}")
print("   " + "-"*65)

for name, chi2_val, M_val in [
    ("ΛCDM", chi2_lcdm_sn, M_lcdm),
    ("Current ESTIF", chi2_cur_sn, M_cur),
    ("Baryons-only ESTIF", chi2_bar_sn, M_bar),
]:
    dchi2 = chi2_lcdm_sn - chi2_val
    p     = chi2_dist.sf(max(dchi2, 0), 1) if dchi2 > 0 else 1.0
    sig   = np.sqrt(chi2_dist.isf(max(p, 1e-15), 1)) if dchi2 > 0 else 0.0
    print(f"   {name:<22} {chi2_val:<12.2f} {chi2_val/dof:<12.4f} "
          f"{dchi2:>+10.2f}     {sig:.2f}σ")

# ============================================================================
# Section 5: Effective equation of state
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 5: DARK ENERGY BEHAVIOUR")
print(f"{'='*70}")

print(f"\n   {'z':<8} {'Ω_tilt':<14} {'Current Ω_tilt':<18} {'Ratio (bar/cur)'}")
print("   " + "-"*55)
for z in [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0]:
    ot_b = omega_tilt_baryons(z)
    ot_c = estif.omega_tilt(z)
    print(f"   {z:<8.1f} {ot_b:<14.6f} {ot_c:<18.6f} {ot_b/ot_c:.4f}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY: CAN THE 4D SHADOW REPLACE DARK MATTER?")
print(f"{'='*70}")

delta_sn = chi2_lcdm_sn - chi2_bar_sn
sn_sig   = 0.0
if delta_sn > 0:
    p = chi2_dist.sf(delta_sn, 1)
    sn_sig = np.sqrt(chi2_dist.isf(max(p, 1e-15), 1))

print(f"""
   Baryons-only ESTIF (Ωb=0.049, Ω_tilt carries rest):

   Age of universe:    {t_baryons:.3f} Gyr  ({'✅ PASS' if t_baryons >= 13.5 else '❌ FAIL — too young'})
   BAO scale:          {n_better_baryons}/{len(BAO_DATA)} redshifts better than ΛCDM
   SN improvement:     Δχ² = {delta_sn:+.2f}  ({sn_sig:.2f}σ)

   For comparison — current ESTIF (Ωm=0.3111):
   Age:  {t_current:.3f} Gyr
   SN:   Δχ² = {chi2_lcdm_sn - chi2_cur_sn:+.2f}

   INTERPRETATION:
""")

if t_baryons >= 13.5 and delta_sn > 0:
    print("   ✅ The 4D tilt CAN geometrically replace dark matter in the")
    print("   Friedmann equation. The universe is old enough and SN data")
    print("   is better fit. This is a major result — proceed to RHAC.")
elif t_baryons >= 13.5 and delta_sn <= 0:
    print("   ⚠️  Age passes but SN fit is worse. The tilt shape with")
    print("   baryons-only background doesn't match distance-redshift data.")
    print("   Dark matter is needed for the correct expansion history.")
elif t_baryons < 13.5 and delta_sn > 0:
    print("   ⚠️  SN improves but universe is too young. The tilt is too")
    print("   strong at low-z. Needs a lower Ω_tilt_0 or a modified shape.")
else:
    print("   ❌ Both tests fail. The 4D shadow in its current form cannot")
    print("   replace dark matter in the cosmological Friedmann equation.")
    print("   The galactic tilt mechanism (Phase 7) is the correct approach.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Baryons-Only ESTIF: Can 4D Tilt Replace Dark Matter?',
             fontsize=14, fontweight='bold')

z_plot = np.linspace(0.001, 2.5, 150)

# Plot 1: H(z) comparison
ax = axes[0, 0]
H_lcdm_arr = np.array([H_lcdm(z) * MPC_TO_M / 1e3 for z in z_plot])
H_cur_arr  = np.array([estif.H_estif(z) * MPC_TO_M / 1e3 for z in z_plot])
H_bar_arr  = np.array([H_baryons(z) * MPC_TO_M / 1e3 for z in z_plot])
ax.plot(z_plot, H_lcdm_arr, 'blue',  linewidth=2.5, label='ΛCDM')
ax.plot(z_plot, H_cur_arr,  'red',   linewidth=2,   linestyle='--',
        label='Current ESTIF (Ωm=0.311)')
ax.plot(z_plot, H_bar_arr,  'green', linewidth=2,   linestyle=':',
        label='Baryons-only (Ωb=0.049)')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('H(z) [km/s/Mpc]', fontsize=11)
ax.set_title('Hubble Parameter', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 2: Ω_tilt evolution
ax = axes[0, 1]
OT_cur_arr = np.array([estif.omega_tilt(z) for z in z_plot])
OT_bar_arr = np.array([omega_tilt_baryons(z) for z in z_plot])
ax.plot(z_plot, OT_cur_arr, 'red',   linewidth=2.5, label=f'Current ESTIF (×Ω_Λ={OL_LCDM})')
ax.plot(z_plot, OT_bar_arr, 'green', linewidth=2.5, linestyle='--',
        label=f'Baryons-only (×Ω_tilt={OT_0})')
ax.axhline(OL_LCDM, color='blue', linewidth=1.5, linestyle=':', label='Ω_Λ ΛCDM')
ax.axhline(OT_0,    color='gray', linewidth=1.5, linestyle=':', label='Ω_tilt_0=0.951')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Dark Energy Density', fontsize=11)
ax.set_title('Ω_tilt Evolution', fontsize=11, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Plot 3: BAO comparison
ax = axes[1, 0]
z_bao   = [r[0] for r in BAO_DATA]
DM_obs  = [r[1] for r in BAO_DATA]
DM_err  = [r[2] for r in BAO_DATA]
DM_l    = [comoving_dist(z, H_lcdm) / R_D for z in z_bao]
DM_c    = [comoving_dist(z, estif.H_estif) / R_D for z in z_bao]
DM_b    = [comoving_dist(z, H_baryons) / R_D for z in z_bao]
ax.errorbar(z_bao, DM_obs, yerr=DM_err, fmt='ko', markersize=8,
            capsize=4, label='BAO observations', zorder=5)
z_s = np.linspace(0.1, 2.0, 100)
ax.plot(z_s, [comoving_dist(z, H_lcdm) / R_D for z in z_s],
        'blue', linewidth=2.5, label='ΛCDM')
ax.plot(z_s, [comoving_dist(z, estif.H_estif) / R_D for z in z_s],
        'red', linewidth=2, linestyle='--', label='Current ESTIF')
ax.plot(z_s, [comoving_dist(z, H_baryons) / R_D for z in z_s],
        'green', linewidth=2, linestyle=':', label='Baryons-only')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('D_M / r_d', fontsize=11)
ax.set_title('BAO Distance Ratio', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 4: Results summary
ax = axes[1, 1]
ax.axis('off')
results_text = (
    f"RESULTS SUMMARY\n"
    f"{'─'*35}\n\n"
    f"{'Metric':<22} {'ΛCDM':<10} {'Current':<10} {'Baryons'}\n"
    f"{'─'*55}\n"
    f"{'Age (Gyr)':<22} {t_lcdm:<10.3f} {t_current:<10.3f} {t_baryons:.3f}\n"
    f"{'Age pass (≥13.5)':<22} {'✅':<10} {'{'+'⚠️' if t_current < 13.5 else '✅'+'}':<10} {'✅' if t_baryons >= 13.5 else '❌'}\n"
    f"{'SN χ²/dof':<22} {chi2_lcdm_sn/dof:<10.4f} {chi2_cur_sn/dof:<10.4f} {chi2_bar_sn/dof:.4f}\n"
    f"{'SN Δχ²':<22} {'0':<10} {chi2_lcdm_sn-chi2_cur_sn:>+9.2f} {delta_sn:>+9.2f}\n"
    f"{'BAO better':<22} {'—':<10} {n_better_current}/6{'':5} {n_better_baryons}/6\n"
    f"{'─'*55}\n\n"
    f"Ωm = 0.3111 (current) vs Ωb = 0.049 (this test)\n"
    f"Difference: dark matter removed from Friedmann eq."
)
ax.text(0.05, 0.95, results_text, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig('baryons_only_cosmology.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: baryons_only_cosmology.png")
print("=" * 70)
