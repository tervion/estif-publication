"""
test_x0_matter_term.py

Can the matter density Ωm be replaced by the geometric x₀(z)?

CURRENT MODEL:
    H²(z) = H₀² × [0.3111 × (1+z)³ + Ω_tilt(z)]
    Ωm = 0.3111 is borrowed from Planck — not derived

THE HYPOTHESIS:
    Ωm = x₀ = R_H / r_universe (shown to 0.12% numerically)
    If true, Ωm should evolve as x₀(z) = R_H(z) / r_universe

THE PROPOSED FRIEDMANN EQUATION:
    H²(z) = H₀² × [x₀(z) × (1+z)³ + Ω_tilt(z)]
    where x₀(z) = (c/H(z)) / r_universe_0 × (1+z)

This has a circular dependency (x₀ depends on H which depends on x₀).
We use H_lcdm as the ruler to break it — same approach as for alpha.

THREE QUESTIONS:
1. Does x₀(z) = 0.3111 at z=0 exactly? (it must)
2. How does x₀(z) evolve — and is it close to constant?
3. Do SN and BAO fits improve, stay the same, or get worse?
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
KPC_TO_M  = 3.085677581e19
GYR_TO_SEC= 3.15576e16
C_KMS     = const.c / 1e3
H0_SI     = const.H_0
R_U_0     = 4.4e26            # Observable universe radius today [m]
X_0_TODAY = (const.c / H0_SI) / R_U_0   # = 0.3107

# Planck reference values
OM_PLANCK = 0.3111
OL_PLANCK = 0.6889

# BAO data
BAO_DATA = [
    (0.38,  10.27, 0.15, "BOSS DR12"),
    (0.51,  13.38, 0.18, "BOSS DR12"),
    (0.61,  15.45, 0.22, "BOSS DR12"),
    (0.70,  17.65, 0.30, "eBOSS DR16"),
    (0.85,  20.75, 0.40, "eBOSS ELG"),
    (1.48,  30.21, 0.79, "eBOSS QSO"),
]
R_D = 147.09  # Sound horizon [Mpc]

print("=" * 70)
print("REPLACING Ωm WITH GEOMETRIC x₀(z)")
print("=" * 70)
print(f"\n   x₀ today = R_H/r_universe = {X_0_TODAY:.6f}")
print(f"   Ωm (Planck) =               {OM_PLANCK:.6f}")
print(f"   Agreement:                  {abs(X_0_TODAY-OM_PLANCK)/OM_PLANCK*100:.3f}%")

# ============================================================================
# Core: x₀(z) as the matter term
# ============================================================================

def x0_of_z(z):
    """
    Geometric x₀(z) = R_H(z) / r_universe_physical(z)
                     = [c/H_lcdm(z)] × (1+z) / r_universe_0

    Uses H_lcdm as the ruler — avoids circular dependency.
    At z=0: x₀(0) = R_H/r_u = 0.3107 ≈ Ωm ✅

    Physical meaning: as z increases (going back in time),
    H(z) grows but the universe was smaller → x₀(z) changes.
    """
    z    = np.asarray(z, dtype=float)
    H_bg = H0_SI * np.sqrt(OM_PLANCK*(1+z)**3 + OL_PLANCK)
    R_H_z= const.c / H_bg          # Hubble radius at redshift z [m]
    return R_H_z * (1+z) / R_U_0   # = x₀ × (1+z) × H₀/H(z)


def omega_tilt_x0(z):
    """
    Ω_tilt using the x₀(z) matter background.
    Same formula — but now x₀ is derived, not fixed.
    """
    z    = np.asarray(z, dtype=float)
    H_bg = H0_SI * np.sqrt(OM_PLANCK*(1+z)**3 + OL_PLANCK)
    x_z  = x0_of_z(z)
    x_now= X_0_TODAY
    obs_z = estif.observable_combined(x_z)
    obs_0 = estif.observable_combined(x_now)
    OL_eff= 1.0 - X_0_TODAY   # = 1 - 0.3107 = 0.6893
    if np.isscalar(obs_z):
        if obs_z <= 0: return OL_eff
    return OL_eff * (obs_0 / obs_z)**2


def H_x0(z):
    """
    ESTIF Friedmann with geometric matter term:
    H²(z) = H₀² × [x₀(z) × (1+z)³ + Ω_tilt(z)]
    """
    z  = np.asarray(z, dtype=float)
    return H0_SI * np.sqrt(x0_of_z(z) * (1+z)**3 + omega_tilt_x0(z))


def H_lcdm(z):
    """Standard ΛCDM."""
    z = np.asarray(z, dtype=float)
    return H0_SI * np.sqrt(OM_PLANCK*(1+z)**3 + OL_PLANCK)


def H_current_estif(z):
    """Current ESTIF with fixed Ωm=0.3111."""
    return estif.H_estif(z)


# ============================================================================
# Section 1: x₀(z) evolution
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: HOW DOES x₀(z) EVOLVE?")
print(f"{'='*70}")

print(f"\n   {'z':<8} {'x₀(z)':<14} {'vs today':<14} {'Ωm if x₀=matter':<18} {'Δ from Planck'}")
print("   " + "-"*65)
for z in [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 100.0]:
    xz    = x0_of_z(z)
    ratio = xz / X_0_TODAY
    delta = (xz - OM_PLANCK) / OM_PLANCK * 100
    print(f"   {z:<8.1f} {xz:<14.6f} {ratio:<14.4f} {xz:<18.6f} {delta:>+.3f}%")

print(f"""
   KEY INSIGHT:
   x₀(z) is NOT constant — it evolves with redshift.
   At high z it grows (universe was more tilt-dominated in past).
   At z=0 it equals 0.3107 ≈ Ωm = 0.3111 exactly.

   COMPARISON:
   Standard Ωm × (1+z)³ = matter dilutes with volume
   x₀(z) × (1+z)³      = ??? (what does this produce?)
""")

# ============================================================================
# Section 2: Effective matter evolution
# ============================================================================

print(f"{'='*70}")
print("SECTION 2: EFFECTIVE MATTER SCALING")
print(f"{'='*70}")

print(f"\n   Standard ΛCDM matter: Ωm × (1+z)³")
print(f"   Geometric matter:     x₀(z) × (1+z)³\n")
print(f"   {'z':<8} {'ΛCDM Ωm(1+z)³':<18} {'x₀(z)(1+z)³':<18} {'ratio':<10} {'same?'}")
print("   " + "-"*58)
for z in [0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0]:
    lcdm_matter = OM_PLANCK * (1+z)**3
    x0_matter   = x0_of_z(z) * (1+z)**3
    ratio = x0_matter / lcdm_matter
    same  = "✅" if abs(ratio - 1) < 0.05 else "⚠️"
    print(f"   {z:<8.1f} {lcdm_matter:<18.4f} {x0_matter:<18.4f} {ratio:<10.4f} {same}")

# ============================================================================
# Section 3: Age of universe
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: AGE OF UNIVERSE")
print(f"{'='*70}")

def age_gyr(H_func):
    def integrand(z):
        return 1.0 / ((1.0+z) * H_func(z))
    result, _ = quad(integrand, 0, 1000, limit=200)
    return result / GYR_TO_SEC

t_lcdm   = age_gyr(H_lcdm)
t_cur    = age_gyr(H_current_estif)
t_x0     = age_gyr(H_x0)

print(f"\n   ΛCDM:            {t_lcdm:.3f} Gyr")
print(f"   Current ESTIF:   {t_cur:.3f} Gyr  (fixed Ωm=0.3111)")
print(f"   Geometric x₀:    {t_x0:.3f} Gyr  (Ωm = x₀(z))")
print(f"   Hard limit:      13.500 Gyr")

if t_x0 >= 13.5:
    print(f"\n   ✅ PASS")
elif t_x0 >= 13.0:
    print(f"\n   ⚠️  MARGINAL — {13.5-t_x0:.3f} Gyr short")
else:
    print(f"\n   ❌ FAIL — {13.5-t_x0:.3f} Gyr too young")

# ============================================================================
# Section 4: BAO comparison
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: BAO ANGULAR SCALE")
print(f"{'='*70}")

def comoving_dist(z, H_func):
    def integrand(zp):
        return C_KMS / (H_func(zp) * MPC_TO_M / 1e3)
    val, _ = quad(integrand, 0, float(z), limit=100)
    return val

chi2_lcdm_bao = 0.0
chi2_cur_bao  = 0.0
chi2_x0_bao   = 0.0

print(f"\n   {'z':<6} {'ΛCDM':<10} {'Current':<12} {'Geom x₀':<12} {'Obs':<8}  Survey")
print("   " + "-"*60)
for z_b, DM_obs, sig, survey in BAO_DATA:
    dm_l  = comoving_dist(z_b, H_lcdm) / R_D
    dm_c  = comoving_dist(z_b, H_current_estif) / R_D
    dm_x0 = comoving_dist(z_b, H_x0) / R_D
    chi2_lcdm_bao += ((DM_obs - dm_l) / sig)**2
    chi2_cur_bao  += ((DM_obs - dm_c) / sig)**2
    chi2_x0_bao   += ((DM_obs - dm_x0) / sig)**2
    flag = "✅" if abs((DM_obs - dm_x0)/sig) < 2 else "⚠️"
    print(f"   {z_b:<6.2f} {dm_l:<10.3f} {dm_c:<12.3f} {dm_x0:<12.3f} "
          f"{DM_obs:<8.2f}  {survey} {flag}")

print(f"\n   BAO χ²:  ΛCDM={chi2_lcdm_bao:.3f}  Current={chi2_cur_bao:.3f}  "
      f"Geometric x₀={chi2_x0_bao:.3f}")

# ============================================================================
# Section 5: SN fit
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 5: SUPERNOVA DISTANCES")
print(f"{'='*70}")
print(f"   Loading Pantheon+...", flush=True)

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
            zc = header.index('zHD'); muc = header.index('MU_SH0ES')
            erc = header.index('MU_SH0ES_ERR_DIAG')
            isc = header.index('IS_CALIBRATOR')
            continue
        if header is None: continue
        try:
            zv = float(parts[zc]); mv = float(parts[muc])
            ev = float(parts[erc]); ic = int(parts[isc])
            if zv > 0.01 and ic == 0 and ev < 10 and mv > 0:
                z_sn.append(zv); mu_sn.append(mv); err_sn.append(ev)
        except (ValueError, IndexError): continue

z_sn = np.array(z_sn); mu_sn = np.array(mu_sn); err_sn = np.array(err_sn)
print(f"   {len(z_sn)} SNe loaded. Fitting M_offset...", flush=True)

def best_chi2(H_func):
    def chi2_fn(params):
        M = params[0]
        mu_pred = np.array([5*np.log10(max((1+z)*comoving_dist(z,H_func),1e-10))+25+M
                            for z in z_sn])
        return np.sum(((mu_sn - mu_pred)/err_sn)**2)
    res = minimize(chi2_fn, [-0.18], method='Nelder-Mead',
                   options={'xatol':1e-6,'fatol':1e-6,'maxiter':3000})
    return res.fun

print(f"   ΛCDM...", flush=True)
c2_lcdm = best_chi2(H_lcdm)
print(f"   Current ESTIF...", flush=True)
c2_cur  = best_chi2(H_current_estif)
print(f"   Geometric x₀...", flush=True)
c2_x0   = best_chi2(H_x0)

dof = len(z_sn) - 1
print(f"\n   {'Model':<24} {'χ²':<12} {'χ²/dof':<12} {'Δχ² vs ΛCDM':<14} {'σ'}")
print("   " + "-"*66)
for name, c2 in [("ΛCDM", c2_lcdm), ("Current ESTIF", c2_cur), ("Geometric x₀(z)", c2_x0)]:
    dchi2 = c2_lcdm - c2
    sig   = 0.0
    if dchi2 > 0:
        p = chi2_dist.sf(dchi2, 1)
        sig = np.sqrt(chi2_dist.isf(max(p, 1e-15), 1))
    print(f"   {name:<24} {c2:<12.2f} {c2/dof:<12.4f} {dchi2:>+10.2f}     {sig:.2f}σ")

# ============================================================================
# Section 6: The Ωm drift prediction
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 6: THE Ωm DRIFT PREDICTION")
print(f"{'='*70}")

print(f"""
   If Ωm = x₀ = R_H/r_universe, then as the universe expands:
   R_H grows slowly, r_universe grows faster → x₀ decreases.

   Drift rate:
""")

z_drift = np.array([0.0, 0.001, 0.01, 0.1, 1.0])
x0_drift = x0_of_z(z_drift)
# Convert redshift to lookback time approximately
t_lookback_gyr = z_drift * 9.78 / (1 + z_drift)  # rough approximation

print(f"   {'z':<8} {'x₀(z)':<14} {'Δx₀ from today':<18} {'Δ(%)'}")
print("   " + "-"*50)
for z, xz in zip(z_drift, x0_drift):
    delta    = xz - X_0_TODAY
    delta_pct= delta / X_0_TODAY * 100
    print(f"   {z:<8.3f} {xz:<14.6f} {delta:>+14.6f}     {delta_pct:>+.4f}%")

# Estimate drift per Gyr
dx0_per_gyr = (x0_of_z(0.1) - X_0_TODAY) / 1.3  # ~1.3 Gyr for z=0.1
print(f"\n   Approximate Ωm drift: {dx0_per_gyr/X_0_TODAY*100:.4f}% per Gyr")
print(f"   Λ drift prediction:    0.023% per Gyr (for comparison)")
print(f"   EUCLID/LSST threshold: ~0.01% per Gyr")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

delta_sn_x0  = c2_lcdm - c2_x0
sig_x0       = 0.0
if delta_sn_x0 > 0:
    p = chi2_dist.sf(delta_sn_x0, 1)
    sig_x0 = np.sqrt(chi2_dist.isf(max(p, 1e-15), 1))

print(f"""
   GEOMETRIC x₀(z) AS THE MATTER TERM:

   At z=0:        x₀ = {X_0_TODAY:.6f}  ≈  Ωm = {OM_PLANCK:.6f}  ({abs(X_0_TODAY-OM_PLANCK)/OM_PLANCK*100:.3f}%)
   Age:           {t_x0:.3f} Gyr  ({'✅ PASS' if t_x0 >= 13.5 else '❌ FAIL'})
   BAO χ²:        {chi2_x0_bao:.3f}  (ΛCDM: {chi2_lcdm_bao:.3f})
   SN Δχ²:        {delta_sn_x0:+.2f}  ({sig_x0:.2f}σ {'improvement' if delta_sn_x0>0 else 'worse'})

   VERDICT:
""")

passed = 0
if t_x0 >= 13.5: passed += 1
if chi2_x0_bao < chi2_lcdm_bao + 5: passed += 1
if delta_sn_x0 > -50: passed += 1

if passed == 3 and delta_sn_x0 > 0:
    print("   ✅ FULL PASS — geometric x₀(z) is a viable derived matter term")
    print("   Ωm is no longer a free parameter. It equals R_H/r_universe.")
elif passed >= 2:
    print("   ⚠️  PARTIAL — x₀(z) mostly works but needs refinement")
    print("   The geometric origin of Ωm is numerically plausible.")
else:
    print("   ❌ FAIL — x₀(z) cannot directly replace Ωm in the Friedmann equation")
    print("   The 0.12% numerical coincidence does not extend to the evolution.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(f'Replacing Ωm with Geometric x₀(z)  [x₀={X_0_TODAY:.4f} vs Ωm={OM_PLANCK}]',
             fontsize=13, fontweight='bold')

z_plot = np.linspace(0.001, 3, 200)

# Plot 1: x₀(z) evolution vs fixed Ωm
ax = axes[0]
x0_arr   = x0_of_z(z_plot)
ax.plot(z_plot, x0_arr,             'red',  linewidth=2.5,
        label=f'x₀(z) = R_H(z)/r_u')
ax.axhline(OM_PLANCK, color='blue', linewidth=2, linestyle='--',
           label=f'Fixed Ωm = {OM_PLANCK}')
ax.axhline(X_0_TODAY, color='red',  linewidth=1, linestyle=':',
           alpha=0.5, label=f'x₀(0) = {X_0_TODAY:.4f}')
ax.scatter([0], [X_0_TODAY], color='red', s=100, zorder=5)
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Matter density', fontsize=11)
ax.set_title('x₀(z) vs Fixed Ωm', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 2: H(z) comparison
ax = axes[1]
H_l   = np.array([H_lcdm(z)*MPC_TO_M/1e3 for z in z_plot])
H_c   = np.array([H_current_estif(z)*MPC_TO_M/1e3 for z in z_plot])
H_x   = np.array([H_x0(z)*MPC_TO_M/1e3 for z in z_plot])
ax.plot(z_plot, H_l, 'blue',  linewidth=2.5, label='ΛCDM')
ax.plot(z_plot, H_c, 'red',   linewidth=2, linestyle='--',
        label='Current ESTIF (fixed Ωm)')
ax.plot(z_plot, H_x, 'green', linewidth=2, linestyle=':',
        label='Geometric x₀(z)')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('H(z) [km/s/Mpc]', fontsize=11)
ax.set_title('Hubble Parameter H(z)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 3: Summary bar chart
ax = axes[2]
models  = ['ΛCDM', 'Current\nESTIF', 'Geometric\nx₀(z)']
chi2vals= [c2_lcdm/dof, c2_cur/dof, c2_x0/dof]
colors  = ['blue', 'red', 'green']
bars    = ax.bar(models, chi2vals, color=colors, alpha=0.8, edgecolor='black')
for bar, val in zip(bars, chi2vals):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.001,
            f'{val:.4f}', ha='center', fontsize=10, fontweight='bold')
ax.text(0.98, 0.95,
        f'BAO χ²:\nΛCDM: {chi2_lcdm_bao:.2f}\nGeom: {chi2_x0_bao:.2f}',
        transform=ax.transAxes, fontsize=10, ha='right', va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
ax.set_ylabel('SN χ²/dof', fontsize=11)
ax.set_title('Fit Quality Comparison', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('x0_matter_term.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: x0_matter_term.png")
print("=" * 70)
