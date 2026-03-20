"""
test_virialized_eddy.py

Option 2: What happens when we use the virialized halo density?

The previous test used ρ_eddy = x₀ × ρ_crit (background density).
This gave σ(10 kpc) = 0.2 km/s — 1000× too small.

The fix: inside a collapsed halo, the density is not the background.
It is ρ_halo = δ_virial × ρ_eddy where δ_virial ≈ 200.

This factor of 200 comes from spherical collapse theory:
when a region collapses and virializes, it reaches an average
overdensity of ~200× the background. This is the same δ_virial
used in ALL dark matter models including ΛCDM.

What changes with δ_virial = 200:
    σ(r) → σ(r) × √200 ≈ 14×
    v_flat → v_flat × √200 ≈ 14×
    t_ff → t_ff / √200 ≈ 14× faster

We also test a range of δ_virial to see where the model breaks or succeeds.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

MPC_TO_M  = 3.085677581e22
KPC_TO_M  = 3.085677581e19
GYR_TO_SEC= 3.15576e16

H0        = const.H_0
X_0       = (const.c / H0) / 4.4e26
RHO_CRIT  = 3 * H0**2 / (8 * np.pi * const.G)
RHO_EDDY_BG = X_0 * RHO_CRIT     # background eddy density
OBS_NOW   = estif.observable_combined(X_0)

print("=" * 70)
print("VIRIALIZED EDDY HALO — δ_virial CORRECTION")
print("=" * 70)
print(f"\n   Background: ρ_eddy = x₀ × ρ_crit = {RHO_EDDY_BG:.4e} kg/m³")
print(f"   Virialized: ρ_halo = δ_virial × ρ_eddy")
print(f"\n   Previous result (δ=1): σ(10 kpc) ≈ 0.2 km/s  ← too small by 1000×")
print(f"   This test:   δ=200:  σ(10 kpc) ≈ {0.2*np.sqrt(200):.1f} km/s  ← getting closer")

# ============================================================================
# Section 1: How δ_virial changes everything
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: THE EFFECT OF δ_virial ON KEY QUANTITIES")
print(f"{'='*70}")

print(f"""
   All physical quantities scale with √δ_virial (from ρ → σ ∝ √ρ):

   σ(r)    ∝ √ρ  →  σ_halo   = σ_bg × √δ
   v_flat  ∝ √ρ  →  v_flat   = v_bg × √δ
   t_ff    ∝ 1/√ρ → t_ff_halo = t_ff_bg / √δ
   λ_Jeans ∝ 1/√ρ → λ_halo   = λ_bg / √δ  (smaller → more concentrated)
""")

delta_vals = [1, 10, 50, 100, 200, 500, 1000]
print(f"   {'δ_virial':<12} {'σ(10kpc) [km/s]':<20} {'v_flat [km/s]':<18} "
      f"{'t_ff [Gyr]':<14} {'λ_J(10kpc) [kpc]'}")
print("   " + "-"*78)

sigma_bg_10kpc = 10 * KPC_TO_M * np.sqrt(2*np.pi*const.G*RHO_EDDY_BG/3) / 1e3
vflat_bg = 220  # km/s target — we'll compute properly below
tff_bg = np.sqrt(3*np.pi / (32*const.G*RHO_EDDY_BG)) / GYR_TO_SEC
lambda_bg_10kpc = 10 * np.sqrt(2*np.pi**2/3)  # kpc, self-similar

for delta in delta_vals:
    sigma_halo = sigma_bg_10kpc * np.sqrt(delta)
    tff_halo   = tff_bg / np.sqrt(delta)
    lambda_halo= lambda_bg_10kpc / np.sqrt(delta)

    flag = ""
    if 50 < sigma_halo < 500:  flag = "← galactic range ✅"
    elif sigma_halo > 500:      flag = "← too hot ⚠️"

    print(f"   {delta:<12} {sigma_halo:<20.1f} {'?':<18} "
          f"{tff_halo:<14.3f} {lambda_halo:.2f}  {flag}")

# ============================================================================
# Section 2: Full calculation at δ = 200
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: FULL CALCULATION AT δ_virial = 200")
print(f"{'='*70}")

DELTA_VIRIAL = 200
RHO_HALO = DELTA_VIRIAL * RHO_EDDY_BG

print(f"\n   ρ_halo = {DELTA_VIRIAL} × ρ_eddy = {RHO_HALO:.4e} kg/m³")

def sigma_halo(r_m):
    return r_m * np.sqrt(2 * np.pi * const.G * RHO_HALO / 3)

def v_circular_halo(r_m):
    return r_m * np.sqrt(4 * np.pi * const.G * RHO_HALO / 3)

print(f"\n   {'Scale r':<25} {'σ_bg [km/s]':<16} {'σ_halo [km/s]':<18} {'Observed σ'}")
print("   " + "-"*72)

scales = [
    (0.1,  "100 pc (GMC)",              "~5-15 km/s"),
    (1,    "1 kpc (bulge)",             "~50-100 km/s"),
    (5,    "5 kpc (inner MW)",          "~150-200 km/s"),
    (8,    "8 kpc (Sun's orbit)",       "~220 km/s ← target"),
    (15,   "15 kpc (MW disk edge)",     "~200 km/s"),
    (50,   "50 kpc (MW halo)",          "~150-200 km/s"),
    (200,  "200 kpc (MW virial)",       "~100-200 km/s"),
    (1000, "1 Mpc (group)",             "~300-1000 km/s"),
]

for r_kpc, name, obs in scales:
    r_m      = r_kpc * KPC_TO_M
    sig_bg   = r_m * np.sqrt(2*np.pi*const.G*RHO_EDDY_BG/3) / 1e3
    sig_halo = sigma_halo(r_m) / 1e3
    flag     = "✅" if abs(sig_halo - 220)/220 < 0.5 else ""
    print(f"   {name:<25} {sig_bg:<16.1f} {sig_halo:<18.1f} {obs} {flag}")

# ============================================================================
# Section 3: Flat rotation velocity at δ = 200
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: FLAT ROTATION VELOCITY")
print(f"{'='*70}")

print(f"""
   For an isothermal halo with ρ ∝ 1/r²:
   v_flat² = 4πG × ρ_0 × r_0²

   The scale radius r_0 sets where the halo transitions from
   rising (inner baryonic dominated) to flat (halo dominated).

   Instead of r_0 = Rs/x₀ (which was too small), we use:
   r_0 = r_virial × x₀  (tidal truncation from the eddy background)

   This is physically motivated: x₀ is the ratio of Hubble radius
   to universe size — it sets the fraction of the virial radius
   that the eddy can "see" from the galaxy's perspective.
""")

M_gals = {
    "Dwarf galaxy (10⁸ M☉)":   1e8  * const.M_sun,
    "Milky Way (10¹² M☉)":     1e12 * const.M_sun,
    "Massive galaxy (10¹³ M☉)":1e13 * const.M_sun,
    "Galaxy cluster (10¹⁴ M☉)":1e14 * const.M_sun,
}

print(f"   {'Galaxy type':<30} {'r_virial [kpc]':<18} {'r₀ [kpc]':<12} "
      f"{'v_flat [km/s]':<16} {'Observed'}")
print("   " + "-"*82)

observed_vflat = {
    "Dwarf galaxy (10⁸ M☉)":    "~20-50 km/s",
    "Milky Way (10¹² M☉)":      "~220 km/s",
    "Massive galaxy (10¹³ M☉)": "~300-400 km/s",
    "Galaxy cluster (10¹⁴ M☉)": "~800-1200 km/s",
}

for name, M_gal in M_gals.items():
    # Virial radius: M = (4π/3) × δ_virial × ρ_crit × r_virial³
    r_virial = (M_gal / ((4*np.pi/3) * DELTA_VIRIAL * RHO_CRIT))**(1/3)
    r_0      = r_virial * X_0   # scale radius from eddy curvature
    v_flat   = np.sqrt(4 * np.pi * const.G * RHO_HALO * r_0**2) / 1e3
    obs      = observed_vflat[name]
    flag     = ""
    try:
        obs_mid = float(obs.split("~")[1].split(" ")[0].split("-")[0])
        if abs(v_flat - obs_mid)/obs_mid < 0.5: flag = "✅"
        elif abs(v_flat - obs_mid)/obs_mid < 2:  flag = "⚠️"
        else:                                     flag = "❌"
    except: pass
    print(f"   {name:<30} {r_virial/KPC_TO_M:<18.1f} {r_0/KPC_TO_M:<12.2f} "
          f"{v_flat:<16.1f} {obs} {flag}")

# ============================================================================
# Section 4: The missing factor — where does it come from?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: FINDING THE MISSING FACTOR")
print(f"{'='*70}")

# Compute what v_flat the Milky Way model gives
M_MW     = 1e12 * const.M_sun
r_vir_MW = (M_MW / ((4*np.pi/3) * DELTA_VIRIAL * RHO_CRIT))**(1/3)
r_0_MW   = r_vir_MW * X_0
v_model  = np.sqrt(4 * np.pi * const.G * RHO_HALO * r_0_MW**2) / 1e3
v_target = 220  # km/s

ratio_needed = v_target / v_model if v_model > 0 else np.inf
factor_needed= ratio_needed**2   # since v ∝ √(ρ r²)

print(f"""
   Milky Way:
   r_virial = {r_vir_MW/KPC_TO_M:.1f} kpc
   r₀ = r_virial × x₀ = {r_0_MW/KPC_TO_M:.2f} kpc
   v_flat (model) = {v_model:.1f} km/s
   v_flat (target) = {v_target} km/s
   Missing factor in v²: {factor_needed:.1f}×

   This factor needs to come from somewhere in the geometry.
   Options:
   A) Different scale radius: r₀ = r_virial × n(x₀) = {r_vir_MW/KPC_TO_M * estif.n_dynamic(X_0):.2f} kpc
   B) Concentration parameter: c = r_virial / r_s where r_s is NFW scale
   C) The eddy density inside the halo needs a correction from
      the local tilt at r_virial: x_local = Rs_galaxy/r_virial

   Testing option C — local tilt correction:
""")

Rs_MW    = 2 * const.G * M_MW / const.c**2
x_local  = Rs_MW / r_vir_MW
n_local  = estif.n_dynamic(x_local)
obs_local= estif.observable_combined(x_local)
print(f"   At r_virial of Milky Way:")
print(f"   x_local = Rs/r_virial = {x_local:.6e}")
print(f"   n(x_local) = {n_local:.4f}")
print(f"   obs(x_local) = {obs_local:.6f}")
print(f"   1/obs(x_local) = {1/obs_local:.4f}  ← tilt amplification")

# The eddy density is amplified by 1/obs at the virial radius
# because the tilt geometry concentrates flow toward mass
rho_local_corrected = RHO_HALO / obs_local  # eddy amplified by tilt
r_0_corrected = r_vir_MW * X_0
v_corrected = np.sqrt(4 * np.pi * const.G * rho_local_corrected * r_0_corrected**2) / 1e3
print(f"\n   With tilt amplification (ρ_halo / obs(x_local)):")
print(f"   ρ_corrected = {rho_local_corrected:.4e} kg/m³")
print(f"   v_flat = {v_corrected:.1f} km/s  (target: {v_target} km/s)")
print(f"   Ratio: {v_corrected/v_target:.3f}  "
      f"{'✅ Close!' if 0.5 < v_corrected/v_target < 2 else '⚠️  still off'}")

# ============================================================================
# Section 5: Mass-velocity relation (Tully-Fisher)
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 5: TULLY-FISHER RELATION — v_flat vs M_galaxy")
print(f"{'='*70}")

print(f"""
   The Tully-Fisher relation is observed: v_flat ∝ M^(1/4)
   (discovered empirically, not derived in ΛCDM — it just works)

   In ESTIF, can we derive it from geometry?

   v_flat² = 4πG × ρ_halo × r₀²
   r₀ = x₀ × r_virial = x₀ × (M / ((4π/3) × δ × ρ_crit))^(1/3)

   v_flat² = 4πG × (δ × x₀ × ρ_crit) × x₀² × (M/(4π/3×δ×ρ_crit))^(2/3)

   v_flat² ∝ M^(2/3)  →  v_flat ∝ M^(1/3)

   The ESTIF model gives v_flat ∝ M^(1/3).
   The observed Tully-Fisher gives v_flat ∝ M^(1/4).

   Exponent comparison:
   ESTIF (geometric): 1/3 = 0.333
   Observed (T-F):    1/4 = 0.250
   MOND prediction:   1/4 = 0.250

   ESTIF gets the scaling direction right but the exponent is off.
   The difference (1/3 vs 1/4) might come from the tilt correction
   to ρ_halo which adds an additional M-dependent factor.
""")

m_range = np.logspace(8, 15, 50) * const.M_sun
v_estif = np.zeros(len(m_range))
v_tf    = np.zeros(len(m_range))

# Calibrate Tully-Fisher to MW
M_MW_calib = 1e12 * const.M_sun
r_vir_calib= (M_MW_calib / ((4*np.pi/3) * DELTA_VIRIAL * RHO_CRIT))**(1/3)
r_0_calib  = r_vir_calib * X_0
v_calib    = np.sqrt(4*np.pi*const.G * RHO_HALO * r_0_calib**2)
# T-F calibrated to 220 km/s at 10^12 Msun
tf_norm = (220e3) / (M_MW_calib**(1/4))

for i, M in enumerate(m_range):
    r_vir = (M / ((4*np.pi/3) * DELTA_VIRIAL * RHO_CRIT))**(1/3)
    r_0   = r_vir * X_0
    v_estif[i] = np.sqrt(4*np.pi*const.G * RHO_HALO * r_0**2) / 1e3
    v_tf[i]    = tf_norm * M**(1/4) / 1e3

# Fit power law to ESTIF
log_m = np.log10(m_range/const.M_sun)
log_v = np.log10(np.maximum(v_estif, 1e-10))
coeffs = np.polyfit(log_m, log_v, 1)
print(f"   ESTIF power law fit: v_flat ∝ M^{coeffs[0]:.3f}")
print(f"   Observed (T-F):      v_flat ∝ M^0.250")

# ============================================================================
# Section 6: Free-fall time at δ=200
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 6: COLLAPSE TIMESCALE AT δ = 200")
print(f"{'='*70}")

t_ff_halo = np.sqrt(3*np.pi / (32*const.G*RHO_HALO)) / GYR_TO_SEC
print(f"\n   t_ff(background, z=0): {np.sqrt(3*np.pi/(32*const.G*RHO_EDDY_BG))/GYR_TO_SEC:.2f} Gyr")
print(f"   t_ff(δ=200 halo, z=0): {t_ff_halo:.3f} Gyr")
print(f"   t_ff(δ=200 halo, z=5): {t_ff_halo/np.sqrt((1+5)**3):.3f} Gyr")
print(f"   t_ff(δ=200 halo, z=10):{t_ff_halo/np.sqrt((1+10)**3):.4f} Gyr")
print(f"\n   At z=5: t_ff = {t_ff_halo/np.sqrt(6**3):.3f} Gyr  "
      f"{'✅ < galaxy formation timescale' if t_ff_halo/np.sqrt(6**3) < 1 else '⚠️'}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY: WHAT δ_virial = 200 CHANGES")
print(f"{'='*70}")

print(f"""
   BEFORE (background density only):
   σ(10 kpc) = 0.2 km/s     → 1000× too small
   v_flat    = ~0 km/s       → completely off
   t_ff(z=0) = 40.7 Gyr     → too slow

   AFTER (δ_virial = 200):
   σ(10 kpc) = {sigma_halo(10*KPC_TO_M)/1e3:.1f} km/s    → still {220/max(sigma_halo(10*KPC_TO_M)/1e3,0.01):.0f}× below observed 220 km/s
   v_flat    = {v_model:.1f} km/s      → {v_target/max(v_model,0.01):.0f}× below 220 km/s
   t_ff(z=5) = {t_ff_halo/np.sqrt(6**3):.3f} Gyr    → {'✅ reasonable' if t_ff_halo/np.sqrt(6**3) < 2 else '⚠️'}

   THE REMAINING GAP:
   The {ratio_needed:.0f}× gap in v_flat suggests the effective overdensity
   inside the optical disk of a galaxy is not 200 but closer to
   {200 * factor_needed:.0f} × ρ_eddy.

   This is physically reasonable: the core of a dark matter halo
   (where most of the gravitational effect is observed) has overdensity
   much larger than the average virial overdensity of 200.
   NFW profiles have central densities 10,000-100,000× ρ_crit.

   HONEST CONCLUSION:
   δ_virial = 200 closes the gap from 1000× to {220/max(sigma_halo(10*KPC_TO_M)/1e3,0.01):.0f}×.
   The remaining gap requires knowing the concentration parameter c
   of the eddy halo — which is a simulation output, not analytical.
   This is where the analytical path ends and simulation begins.
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(f'Virialized Eddy Halo (δ={DELTA_VIRIAL}): Closing The Gap',
             fontsize=13, fontweight='bold')

r_range_kpc = np.logspace(-1, 4, 300)
r_range_m   = r_range_kpc * KPC_TO_M

# Plot 1: σ(r) — background vs virialized
ax = axes[0]
sig_bg_arr   = np.array([r * np.sqrt(2*np.pi*const.G*RHO_EDDY_BG/3) for r in r_range_m]) / 1e3
sig_halo_arr = np.array([sigma_halo(r) for r in r_range_m]) / 1e3

ax.loglog(r_range_kpc, sig_bg_arr,   'blue',  linewidth=2,
          linestyle='--', label=f'Background (δ=1): σ ∝ r', alpha=0.6)
ax.loglog(r_range_kpc, sig_halo_arr, 'red',   linewidth=2.5,
          label=f'Virialized (δ={DELTA_VIRIAL}): σ × √{DELTA_VIRIAL}')

# Observed dispersions
obs_r   = [0.1,  1,    8,   50,   200,  1000]
obs_sig = [10,   80,  220,  175,  200,   700]
obs_err = [5,    30,   30,   30,   50,   300]
ax.errorbar(obs_r, obs_sig, yerr=obs_err, fmt='k*', markersize=10,
            capsize=4, label='Observed', zorder=5)

ax.axhline(220, color='gray', linewidth=1, linestyle=':', alpha=0.5)
ax.text(0.15, 250, '220 km/s (MW flat)', fontsize=8)
ax.set_xlabel('Radius [kpc]', fontsize=11)
ax.set_ylabel('Velocity dispersion σ [km/s]', fontsize=11)
ax.set_title('σ(r): Background vs Virialized\n(δ=200 closes gap by √200 = 14×)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')

# Plot 2: Tully-Fisher — v_flat vs galaxy mass
ax = axes[1]
m_solar = m_range / const.M_sun
ax.loglog(m_solar, v_estif, 'red',   linewidth=2.5, label=f'ESTIF (δ={DELTA_VIRIAL}): ∝M^{coeffs[0]:.2f}')
ax.loglog(m_solar, v_tf,    'blue',  linewidth=2, linestyle='--', label='Observed T-F: ∝M^0.25')
ax.scatter([1e8, 1e12, 1e13, 1e14],
           [30, 220, 350, 900], color='black', s=100, zorder=5,
           marker='*', label='Reference galaxies')
ax.set_xlabel('Galaxy mass [M☉]', fontsize=11)
ax.set_ylabel('Flat rotation velocity [km/s]', fontsize=11)
ax.set_title('Tully-Fisher Relation\nESTIF ∝M^(1/3) vs Observed ∝M^(1/4)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3, which='both')

# Plot 3: Progress summary — gap closing
deltas = np.logspace(0, 5, 100)
v_flat_arr = np.zeros(len(deltas))
for i, d in enumerate(deltas):
    rho_h = d * RHO_EDDY_BG
    r_v   = (M_MW / ((4*np.pi/3) * d * RHO_CRIT))**(1/3)
    r_0_d = r_v * X_0
    v_flat_arr[i] = np.sqrt(4*np.pi*const.G * rho_h * r_0_d**2) / 1e3

ax = axes[2]
ax.loglog(deltas, v_flat_arr, 'red', linewidth=2.5, label='v_flat vs δ_virial')
ax.axhline(220, color='blue', linewidth=2, linestyle='--', label='Target: 220 km/s')
ax.axvline(200, color='gray', linewidth=1.5, linestyle=':', label='δ=200 (standard)')
# Find where it crosses 220
crossing = deltas[np.argmin(np.abs(v_flat_arr - 220))]
ax.axvline(crossing, color='green', linewidth=2, linestyle='--',
           label=f'220 km/s at δ≈{crossing:.0f}')
ax.fill_between([200, crossing], [1, 1], [1000, 1000],
                alpha=0.15, color='orange', label='Simulation range')
ax.set_xlabel('Overdensity δ_virial', fontsize=11)
ax.set_ylabel('v_flat Milky Way [km/s]', fontsize=11)
ax.set_title('Required δ for v_flat = 220 km/s\n(where analytical path ends)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('virialized_eddy.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"✓ Plot saved: virialized_eddy.png")
print("=" * 70)
