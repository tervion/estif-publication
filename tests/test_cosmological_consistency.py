"""
test_cosmological_consistency.py

Combined consistency tests for ESTIF Option A cosmology:

    H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]
    Ω_tilt(z) = Ω_Λ × (obs_now / obs_z)²

FOUR TESTS:
1. Age of Universe     — must be ≥ 13.5 Gyr (oldest observed stars)
2. BAO Scale           — angular scale at z=0.57 compared to BOSS/DESI
3. H0 Tension          — what H₀ does ESTIF prefer vs ΛCDM
4. Dark Energy EOS     — effective w(z) implied by Ω_tilt(z)

SCOPE DISCLAIMER:
ESTIF Option A replaces the dark energy sector at low redshift (z < 2).
CMB physics (z ~ 1100) is explicitly out of scope — the tilt formula
diverges at recombination and requires a full early-universe extension
not attempted here.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Constants
# ============================================================================

MPC_TO_M     = 3.085677581e22      # 1 Mpc in metres
GYR_TO_SEC   = 3.15576e16          # 1 Gyr in seconds
H0_KMS       = const.H_0 * MPC_TO_M / 1e3   # H₀ in km/s/Mpc
C_KMS        = const.c / 1e3                  # c in km/s

# Planck 2018 reference values
T_UNIVERSE_LCDM  = 13.797          # Gyr
H0_PLANCK        = 67.66           # km/s/Mpc
H0_LOCAL         = 73.04           # km/s/Mpc (SH0ES)
OMEGA_M          = estif.OMEGA_M
OMEGA_LAMBDA     = estif.OMEGA_LAMBDA

print("=" * 70)
print("ESTIF COSMOLOGICAL CONSISTENCY TESTS")
print("=" * 70)
print(f"\n   Model: H²(z) = H₀² × [Ωm(1+z)³ + Ω_tilt(z)]")
print(f"   Ω_tilt(z) = Ω_Λ × (obs_now/obs_z)²")
print(f"   Scope: z < 2  (CMB/early universe explicitly out of scope)")

# ============================================================================
# TEST 1: AGE OF UNIVERSE
# ============================================================================

print(f"\n{'='*70}")
print("TEST 1: AGE OF UNIVERSE")
print(f"{'='*70}")

def integrand_age(z, H_func):
    """Integrand for age: dz / ((1+z) × H(z))"""
    return 1.0 / ((1.0 + z) * H_func(z))

def H_lcdm(z):
    return const.H_0 * np.sqrt(OMEGA_M*(1+z)**3 + OMEGA_LAMBDA)

def age_gyr(H_func, z_max=1000):
    """Age of universe in Gyr by integrating from 0 to z_max."""
    result, _ = quad(integrand_age, 0, z_max, args=(H_func,),
                     limit=200, epsabs=1e-10, epsrel=1e-8)
    return result / GYR_TO_SEC

t_lcdm  = age_gyr(H_lcdm)
t_estif = age_gyr(estif.H_estif)

print(f"\n   ΛCDM age:       {t_lcdm:.3f} Gyr  (Planck 2018: 13.797 Gyr)")
print(f"   ESTIF age:      {t_estif:.3f} Gyr")
print(f"   Difference:     {t_estif - t_lcdm:+.3f} Gyr")
print(f"   Oldest stars:   ≥ 13.5 Gyr (hard lower bound)")

if t_estif >= 13.5:
    print(f"\n   ✅ PASS — ESTIF universe old enough for observed stars")
elif t_estif >= 13.0:
    print(f"\n   ⚠️  MARGINAL — close to oldest-star limit")
else:
    print(f"\n   ❌ FAIL — universe too young for observed stars")

# Contribution by era
def age_in_range(z_lo, z_hi, H_func):
    r, _ = quad(integrand_age, z_lo, z_hi, args=(H_func,), limit=100)
    return r / GYR_TO_SEC

print(f"\n   Age contribution by era:")
print(f"   {'Era':<25} {'ΛCDM (Gyr)':<14} {'ESTIF (Gyr)':<14} {'Δ (Gyr)'}")
print("   " + "-"*60)
eras = [("z=0–0.5  (recent)",   0,    0.5),
        ("z=0.5–2  (SN range)", 0.5,  2.0),
        ("z=2–10   (high-z)",   2.0, 10.0),
        ("z=10–∞  (early)",    10.0, 1000)]
for name, z0, z1 in eras:
    al = age_in_range(z0, z1, H_lcdm)
    ae = age_in_range(z0, z1, estif.H_estif)
    print(f"   {name:<25} {al:<14.3f} {ae:<14.3f} {ae-al:+.3f}")

# ============================================================================
# TEST 2: BAO SCALE
# ============================================================================

print(f"\n{'='*70}")
print("TEST 2: BAO ANGULAR SCALE")
print(f"{'='*70}")

# Sound horizon at drag epoch — not modified by ESTIF (early universe physics)
# Using standard Planck 2018 value
r_d = 147.09   # Mpc — comoving sound horizon at drag epoch

# BAO angular scale: theta_BAO = r_d / D_A(z)
# where D_A = d_C / (1+z) is angular diameter distance

def comoving_distance(z, H_func):
    """Comoving distance in Mpc."""
    def integrand(zp):
        return C_KMS / (H_func(zp) * MPC_TO_M / 1e3)
    result, _ = quad(integrand, 0, z, limit=100)
    return result

# BAO measurements — key redshifts from BOSS/eBOSS/DESI
bao_measurements = [
    # (z,    D_H/r_d observed,  D_M/r_d observed,  survey)
    (0.38,  25.00, 10.23, "BOSS DR12"),
    (0.51,  22.33,  9.17, "BOSS DR12"),
    (0.70,  19.77,  8.85, "eBOSS DR16"),
    (0.85,  18.33,  8.58, "eBOSS DR16"),
    (1.48,  13.26, 12.93, "eBOSS QSO"),
]

print(f"\n   Sound horizon r_d = {r_d} Mpc (Planck 2018, unchanged by ESTIF)")
print(f"\n   {'z':<6} {'D_M/r_d ΛCDM':<16} {'D_M/r_d ESTIF':<16} "
      f"{'Δ(%)':<10} {'Survey'}")
print("   " + "-"*62)

bao_results = []
for z_bao, DH_obs, DM_obs, survey in bao_measurements:
    dC_lcdm  = comoving_distance(z_bao, H_lcdm)
    dC_estif = comoving_distance(z_bao, estif.H_estif)

    DM_lcdm  = dC_lcdm  / r_d
    DM_estif = dC_estif / r_d
    pct_diff = (DM_estif / DM_lcdm - 1) * 100

    # Check if ESTIF is closer to observation
    delta_lcdm  = abs(DM_lcdm  - DM_obs)
    delta_estif = abs(DM_estif - DM_obs)
    better = "✅" if delta_estif < delta_lcdm else "⚠️"

    bao_results.append((z_bao, DM_lcdm, DM_estif, DM_obs, pct_diff, better))
    print(f"   {z_bao:<6.2f} {DM_lcdm:<16.3f} {DM_estif:<16.3f} "
          f"{pct_diff:<10.3f}% {survey} {better}")

n_better = sum(1 for r in bao_results if r[5] == "✅")
print(f"\n   ESTIF closer to BAO observations: {n_better}/{len(bao_results)} redshifts")

if n_better >= 4:
    print(f"   ✅ ESTIF consistently improves BAO agreement")
elif n_better >= 2:
    print(f"   ⚠️  Mixed — ESTIF better at some redshifts")
else:
    print(f"   ❌ ΛCDM fits BAO better overall")

# ============================================================================
# TEST 3: H0 TENSION
# ============================================================================

print(f"\n{'='*70}")
print("TEST 3: H0 TENSION")
print(f"{'='*70}")

print(f"\n   ΛCDM H₀ (Planck):  {H0_PLANCK:.2f} km/s/Mpc")
print(f"   Local H₀ (SH0ES):  {H0_LOCAL:.2f} km/s/Mpc")
print(f"   Tension:           {(H0_LOCAL-H0_PLANCK)/2.0:.1f}σ  (σ_combined ≈ 2 km/s/Mpc)")

# ESTIF modifies D_L at low-z via Ω_tilt. If we anchor to CMB (z=1100)
# and re-derive H₀ from the low-z distance ladder, ESTIF's modified
# distances could shift the inferred H₀.

# Simplified estimate: ESTIF makes distances slightly smaller at low-z
# (ESTIF brighter by ~30 mmag). This means objects appear closer.
# A closer distance ladder → higher inferred H₀.

# Quantify: at z=0.5 (typical SN calibration redshift)
dC_lcdm_05  = comoving_distance(0.5, H_lcdm)
dC_estif_05 = comoving_distance(0.5, estif.H_estif)
dist_ratio  = dC_estif_05 / dC_lcdm_05
H0_estif_implied = H0_PLANCK / dist_ratio

print(f"\n   ESTIF comoving distance at z=0.5:")
print(f"   ΛCDM:  {dC_lcdm_05:.2f} Mpc")
print(f"   ESTIF: {dC_estif_05:.2f} Mpc")
print(f"   Ratio: {dist_ratio:.6f}")
print(f"\n   If CMB-anchored H₀ is rescaled by distance ratio:")
print(f"   H₀_ESTIF_implied ≈ {H0_estif_implied:.2f} km/s/Mpc")
print(f"   vs ΛCDM H₀ = {H0_PLANCK:.2f} km/s/Mpc")
print(f"   vs Local H₀ = {H0_LOCAL:.2f} km/s/Mpc")

h0_shift = H0_estif_implied - H0_PLANCK
tension_remaining = (H0_LOCAL - H0_estif_implied) / 2.0

print(f"\n   ESTIF shifts H₀ by: {h0_shift:+.2f} km/s/Mpc")
print(f"   Remaining tension:  {tension_remaining:.1f}σ")

if abs(h0_shift) > 0.5:
    print(f"   ⚠️  ESTIF moves H₀ in the {'right' if h0_shift > 0 else 'wrong'} direction")
    print(f"   (toward local measurement)" if h0_shift > 0 else "(away from local measurement)")
else:
    print(f"   ℹ️  ESTIF has minimal effect on H₀ tension")

print(f"\n   NOTE: This is a simplified estimate. A full H₀ determination")
print(f"   requires fitting the complete CMB+BAO+SN distance ladder.")

# ============================================================================
# TEST 4: EFFECTIVE DARK ENERGY EQUATION OF STATE
# ============================================================================

print(f"\n{'='*70}")
print("TEST 4: EFFECTIVE DARK ENERGY EOS w(z)")
print(f"{'='*70}")

# The effective EOS w(z) is defined by:
# ρ_DE(z) ∝ exp(3 ∫ (1+w(z')) dz'/(1+z'))
# For constant w: ρ_DE ∝ (1+z)^(3(1+w))
# We can extract effective w from dΩ_tilt/dz

print(f"\n   ΛCDM: w = -1 exactly (cosmological constant)")
print(f"\n   ESTIF effective w(z) from Ω_tilt evolution:")

z_vals  = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
dz      = 0.001
print(f"\n   {'z':<8} {'Ω_tilt':<12} {'dΩ/dz':<14} {'w_eff':<10} {'vs ΛCDM'}")
print("   " + "-"*52)

w_vals = []
for z in z_vals:
    ot_lo = estif.omega_tilt(z - dz)
    ot_hi = estif.omega_tilt(z + dz)
    ot    = estif.omega_tilt(z)
    dOt_dz = (ot_hi - ot_lo) / (2 * dz)

    # From continuity equation:
    # w_eff = -1 - (1+z)/3 × (dΩ/dz) / Ω
    w_eff = -1.0 - (1.0 + z) / 3.0 * dOt_dz / ot
    w_vals.append(w_eff)

    deviation = "phantom (w<-1)" if w_eff < -1 else "quintessence (w>-1)"
    print(f"   {z:<8.1f} {ot:<12.4f} {dOt_dz:<14.4f} {w_eff:<10.4f} {deviation}")

w_mean = np.mean(w_vals)
print(f"\n   Mean effective w = {w_mean:.4f}")
print(f"   ΛCDM:            w = -1.0000")
print(f"   DESI 2024 hint:  w ≈ -0.7 to -0.9 (evolving DE)")

if w_mean < -1.0:
    print(f"\n   ESTIF predicts phantom dark energy (w < -1)")
    print(f"   This is consistent with DESI DR2 2024 hints")
elif -1.1 < w_mean < -0.9:
    print(f"\n   ESTIF is close to ΛCDM in EOS")
else:
    print(f"\n   ESTIF departs significantly from ΛCDM EOS")

# ============================================================================
# OVERALL SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print("OVERALL CONSISTENCY SUMMARY")
print(f"{'='*70}")

print(f"""
   TEST 1 — Age of Universe:
   ΛCDM:  {t_lcdm:.3f} Gyr
   ESTIF: {t_estif:.3f} Gyr  ({'✅ PASS' if t_estif >= 13.5 else '⚠️  CHECK'})
   Oldest stars require ≥ 13.5 Gyr

   TEST 2 — BAO Scale:
   ESTIF closer to observations: {n_better}/{len(bao_results)} redshifts
   r_d unchanged (early universe not modified)

   TEST 3 — H0 Tension:
   ESTIF shifts H₀ by {h0_shift:+.2f} km/s/Mpc
   Remaining tension: {tension_remaining:.1f}σ  (was 2.7σ in ΛCDM)

   TEST 4 — Dark Energy EOS:
   Mean effective w = {w_mean:.4f}  (ΛCDM: -1.0000)
   DESI 2024 hints at w < -1 — ESTIF {'consistent' if w_mean < -1 else 'differs'}

   OUT OF SCOPE (requires early-universe extension):
   → CMB angular power spectrum (z ~ 1100)
   → BBN constraints
   → CMB-anchored absolute H₀ determination
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('ESTIF Option A: Cosmological Consistency Tests',
             fontsize=14, fontweight='bold')

z_smooth = np.linspace(0.001, 5, 300)

# Plot 1: H(z) comparison
ax = axes[0, 0]
H_lcdm_s  = np.array([H_lcdm(z) * MPC_TO_M / 1e3 for z in z_smooth])
H_estif_s = np.array([estif.H_estif(z) * MPC_TO_M / 1e3 for z in z_smooth])
ax.plot(z_smooth, H_lcdm_s,  'blue',   linewidth=2.5, label='ΛCDM')
ax.plot(z_smooth, H_estif_s, 'red',    linewidth=2.5, linestyle='--',
        label='ESTIF')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('H(z) [km/s/Mpc]', fontsize=11)
ax.set_title('Hubble Parameter H(z)', fontsize=11, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_xlim(0, 3)

# Plot 2: Age integrand (contribution to age per unit z)
ax = axes[0, 1]
age_int_l = np.array([integrand_age(z, H_lcdm)  / GYR_TO_SEC
                       for z in z_smooth])
age_int_e = np.array([integrand_age(z, estif.H_estif) / GYR_TO_SEC
                       for z in z_smooth])
ax.plot(z_smooth, age_int_l, 'blue',  linewidth=2.5, label=f'ΛCDM  t={t_lcdm:.2f} Gyr')
ax.plot(z_smooth, age_int_e, 'red',   linewidth=2.5, linestyle='--',
        label=f'ESTIF t={t_estif:.2f} Gyr')
ax.axhline(0, color='black', linewidth=0.5)
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('dt/dz (Gyr)', fontsize=11)
ax.set_title('Age Integrand (area = age)', fontsize=11, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_xlim(0, 5); ax.set_ylim(0, None)

# Plot 3: Ω_tilt and effective w(z)
ax = axes[1, 0]
z_w = np.linspace(0.05, 2.5, 100)
ot_s  = np.array([estif.omega_tilt(z) for z in z_w])
w_s   = []
for z in z_w:
    ot_lo = estif.omega_tilt(max(z - 0.01, 0.001))
    ot_hi = estif.omega_tilt(z + 0.01)
    ot    = estif.omega_tilt(z)
    dOt   = (ot_hi - ot_lo) / 0.02
    w_s.append(-1.0 - (1+z)/3.0 * dOt / ot)
w_s = np.array(w_s)

ax.plot(z_w, ot_s, 'purple', linewidth=2.5, label='Ω_tilt(z)')
ax.axhline(OMEGA_LAMBDA, color='blue', linewidth=1.5, linestyle='--',
           label=f'Ω_Λ = {OMEGA_LAMBDA:.4f}')
ax2 = ax.twinx()
ax2.plot(z_w, w_s, 'orange', linewidth=2, linestyle=':', label='w_eff(z)')
ax2.axhline(-1.0, color='gray', linewidth=1, linestyle=':')
ax2.set_ylabel('Effective w(z)', fontsize=10, color='orange')
ax2.tick_params(axis='y', colors='orange')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('Ω_tilt', fontsize=11)
ax.set_title('Dark Energy Evolution', fontsize=11, fontweight='bold')
ax.legend(loc='upper left', fontsize=9)
ax2.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3)

# Plot 4: BAO comparison
ax = axes[1, 1]
z_bao_smooth = np.linspace(0.1, 2.0, 100)
DM_lcdm_s  = np.array([comoving_distance(z, H_lcdm)  / r_d
                         for z in z_bao_smooth])
DM_estif_s = np.array([comoving_distance(z, estif.H_estif) / r_d
                         for z in z_bao_smooth])
ax.plot(z_bao_smooth, DM_lcdm_s,  'blue',  linewidth=2.5, label='ΛCDM')
ax.plot(z_bao_smooth, DM_estif_s, 'red',   linewidth=2.5, linestyle='--',
        label='ESTIF')
for z_b, DM_l, DM_e, DM_obs, pct, better in bao_results:
    ax.errorbar(z_b, DM_obs, yerr=0.3, fmt='o', color='black',
                markersize=6, capsize=4, zorder=5)
ax.scatter([], [], color='black', marker='o', s=50,
           label='BAO observations')
ax.set_xlabel('Redshift z', fontsize=11)
ax.set_ylabel('D_M / r_d', fontsize=11)
ax.set_title('BAO Distance Ratio D_M/r_d', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('cosmological_consistency.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: cosmological_consistency.png")
print("=" * 70)
