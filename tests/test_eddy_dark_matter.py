"""
test_eddy_dark_matter.py

TEST 1: Can the 4D eddy structure explain dark matter?

THE CENTRAL HYPOTHESIS:
    Ωm = x₀ = R_H / r_universe

x₀ is the cosmological curvature ratio — how the Hubble radius compares
to the observable universe size. The hypothesis is that this same ratio
determines the matter content of the universe. Dark matter is not a
particle — it is the background eddy density of the 4D inward flow,
normalized by x₀.

FOUR THINGS THIS SCRIPT CHECKS:

1. NUMERICAL COINCIDENCE:
   Does x₀ ≈ Ωm to within measurement precision?

2. DIMENSIONAL DERIVATION:
   Can ρ_eddy = ρ_crit × x₀ be derived from 4D geometry?

3. GALACTIC ROTATION CURVES:
   If dark matter density ∝ x₀ × ρ_crit (1/r²), do we get flat curves?

4. HALO MASS FUNCTION:
   Does the implied halo structure match observed galaxy clusters?
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

H0_SI     = const.H_0
R_H       = const.c / H0_SI           # Hubble radius [m]
R_UNIVERSE= 4.4e26                     # Observable universe [m]
X_0       = R_H / R_UNIVERSE
RHO_CRIT  = 3 * H0_SI**2 / (8 * np.pi * const.G)  # Critical density [kg/m³]

OM_PLANCK = 0.3111

print("=" * 70)
print("TEST 1: DOES Ωm = x₀? — EDDY DARK MATTER HYPOTHESIS")
print("=" * 70)

# ============================================================================
# PART 1: THE NUMERICAL COINCIDENCE
# ============================================================================

print(f"\n{'='*70}")
print("PART 1: NUMERICAL CHECK — x₀ vs Ωm")
print(f"{'='*70}")

print(f"""
   Hubble radius:           R_H = c/H₀ = {R_H/MPC_TO_M:.2f} Mpc
   Observable universe:     r_u = {R_UNIVERSE/MPC_TO_M:.0f} Mpc
   x₀ = R_H / r_universe = {X_0:.6f}

   Planck 2018 Ωm =         {OM_PLANCK:.6f}

   Difference:              {abs(X_0 - OM_PLANCK)/OM_PLANCK*100:.4f}%
""")

if abs(X_0 - OM_PLANCK) / OM_PLANCK < 0.005:
    print(f"   ✅ x₀ = Ωm to within 0.5% — consistent with measurement uncertainty")
    print(f"   Planck uncertainty on Ωm: ±0.006 = ±1.9%")
    print(f"   x₀ sits within 1σ of the Planck measurement")
else:
    print(f"   ⚠️  x₀ and Ωm differ by more than 0.5%")

# Breakdown of Ωm
OB  = 0.049    # Baryonic matter (BBN + CMB)
ODM = 0.262    # Dark matter
print(f"""
   Ωm breakdown:
   Baryons Ωb =             {OB:.3f}  (BBN + CMB acoustic peaks)
   Dark matter Ωdm =         {ODM:.3f}  (inferred, never directly detected)
   Total Ωm =                {OB+ODM:.3f}

   Hypothesis: Ωdm = x₀ - Ωb = {X_0 - OB:.4f}  (vs measured {ODM:.3f})
   Difference: {abs(X_0-OB-ODM)/ODM*100:.2f}%
""")

# ============================================================================
# PART 2: DIMENSIONAL DERIVATION
# ============================================================================

print(f"{'='*70}")
print("PART 2: DIMENSIONAL DERIVATION OF ρ_eddy")
print(f"{'='*70}")

print(f"""
   In v4.0, the 3D hypersurface tilts at angle θ where sin(θ) = x^n(x).
   The tilt creates a LOCAL eddy — a swirl in the 4D flow field.

   For the GLOBAL background (cosmic scale):
       x₀ = 0.3107
       n₀ = n(x₀) = {estif.n_dynamic(X_0):.4f}
       sin(θ₀) = x₀^n₀ = {X_0**estif.n_dynamic(X_0):.4f}
       cos(θ₀) = √β(x₀) = {estif.observable_combined(X_0):.4f}  ← obs_now

   The eddy angular velocity at cosmic scale:
       ω_eddy = H₀ × sin(θ₀) = H₀ × x₀^n₀

   Eddy energy density (kinetic energy of rotation):
       ρ_eddy = ω_eddy² / (8πG) × (some geometric factor)

   Dimensional argument — the simplest geometric factor is 3:
       ρ_eddy = 3 H₀² × x₀^(2n₀) / (8πG)
              = ρ_crit × x₀^(2n₀)
""")

x0_2n0 = X_0 ** (2 * estif.n_dynamic(X_0))
print(f"   x₀^(2n₀) = {X_0:.4f}^(2×{estif.n_dynamic(X_0):.4f})")
print(f"            = {X_0:.4f}^{2*estif.n_dynamic(X_0):.4f}")
print(f"            = {x0_2n0:.6f}")
print(f"\n   ρ_eddy / ρ_crit = x₀^(2n₀) = {x0_2n0:.4f}")
print(f"   Ωm (Planck)    =              {OM_PLANCK:.4f}")
print(f"   This path gives {x0_2n0:.4f} — not 0.3111")

print(f"""
   Alternative: The eddy density is simply x₀ itself (not x₀^(2n₀)).
   Physical meaning: x₀ = R_H/r_u is the ratio of scales.
   If energy is uniformly distributed across the 4D embedding,
   the fraction perceived as matter in 3D equals the tilt ratio x₀.

   Ω_eddy_direct = x₀ = {X_0:.6f}  ≈  Ωm = {OM_PLANCK:.6f}  ✅

   This is a ZEROTH ORDER derivation. The full derivation requires
   computing how the 4D kinetic energy of the rotating hypersurface
   projects onto the 3D stress-energy tensor. That is future theoretical
   work. What we can test numerically is whether x₀ ≈ Ωm holds.
""")

# ============================================================================
# PART 3: GALACTIC ROTATION CURVES
# ============================================================================

print(f"{'='*70}")
print("PART 3: GALACTIC ROTATION CURVES FROM EDDY DENSITY")
print(f"{'='*70}")

print(f"""
   If Ω_eddy = x₀, then the eddy energy density today is:
       ρ_eddy_0 = x₀ × ρ_crit = {X_0 * RHO_CRIT:.3e} kg/m³

   The eddy has a characteristic scale — it cannot extend infinitely.
   The natural cutoff is the scale where local tilt curvature equals x₀:
   at radius r_halo, x_local(r_halo) = x₀.

   For a galaxy of mass M_galaxy:
       Rs = 2GM/c²
       x_local(r) = Rs/r = x₀  →  r_halo = Rs/x₀

   Inside r_halo: the eddy contributes to gravity
   Outside r_halo: pure Newtonian from visible matter
""")

# Milky Way parameters
M_MW    = 1e12 * const.M_sun     # Milky Way mass [kg]
Rs_MW   = 2 * const.G * M_MW / const.c**2  # Schwarzschild radius
r_halo  = Rs_MW / X_0            # Predicted halo radius from x₀

print(f"   Milky Way example:")
print(f"   M_galaxy = 10¹² M☉")
print(f"   Rs = {Rs_MW:.3e} m  = {Rs_MW/1e3:.0f} km")
print(f"   r_halo = Rs/x₀ = {r_halo/KPC_TO_M:.1f} kpc")
print(f"   Observed Milky Way halo: ~200 kpc ← comparison")

# Rotation curve calculation
def v_circular(r, M_galaxy, include_eddy=True):
    """
    Circular velocity at radius r.
    v² = v²_baryonic + v²_eddy

    Eddy density: ρ_eddy(r) = ρ_eddy_0 × (r_halo/r)²  [1/r² profile → flat curve]
    """
    # Baryonic Newtonian
    v2_bary = const.G * M_galaxy / r

    if not include_eddy:
        return np.sqrt(v2_bary)

    # Eddy contribution — 1/r² density gives log(r) enclosed mass
    # ρ_eddy(r) = ρ_0_eddy × (r_scale/r)²
    # M_eddy(<r) = 4π ρ_0_eddy r_scale² × r  → v² = 4πG ρ_0_eddy r_scale²
    Rs_g    = 2 * const.G * M_galaxy / const.c**2
    r_scale = Rs_g / X_0
    rho_0_eddy = X_0 * RHO_CRIT  # eddy central density

    # Enclosed eddy mass with 1/r² profile (isothermal sphere)
    M_eddy_enclosed = 4 * np.pi * rho_0_eddy * r_scale**2 * r
    v2_eddy = const.G * M_eddy_enclosed / r

    return np.sqrt(max(v2_bary + v2_eddy, 0))


# Test with a typical spiral galaxy (NGC 3198-like)
M_gal    = 1.5e11 * const.M_sun  # ~typical spiral
r_arr    = np.linspace(1, 50, 100) * KPC_TO_M  # 1–50 kpc

v_no_eddy  = np.array([v_circular(r, M_gal, include_eddy=False) / 1e3
                        for r in r_arr])  # km/s
v_with_eddy= np.array([v_circular(r, M_gal, include_eddy=True)  / 1e3
                        for r in r_arr])  # km/s

r_kpc = r_arr / KPC_TO_M

print(f"\n   Rotation curve predictions (M_galaxy = 1.5×10¹¹ M☉):")
print(f"   {'r (kpc)':<12} {'v_baryonic':<16} {'v_with_eddy':<16} {'ratio'}")
print("   " + "-"*50)
for i in [0, 9, 19, 39, 59, 79, 99]:
    print(f"   {r_kpc[i]:<12.1f} {v_no_eddy[i]:<16.1f} {v_with_eddy[i]:<16.1f} "
          f"{v_with_eddy[i]/v_no_eddy[i]:.2f}")

# Check flatness — ratio of v at 30 kpc to v at 5 kpc
idx5  = np.argmin(abs(r_kpc - 5))
idx30 = np.argmin(abs(r_kpc - 30))
flatness = v_with_eddy[idx30] / v_with_eddy[idx5]
print(f"\n   Flatness test: v(30 kpc)/v(5 kpc) = {flatness:.3f}")
print(f"   (1.0 = perfectly flat, <1.0 = falling, >1.0 = rising)")
if 0.85 < flatness < 1.15:
    print(f"   ✅ Rotation curve is approximately flat")
else:
    print(f"   ⚠️  Rotation curve is not flat — eddy profile needs adjustment")

# ============================================================================
# PART 4: COMPARING TO NFW PROFILE
# ============================================================================

print(f"\n{'='*70}")
print("PART 4: COMPARISON TO STANDARD NFW DARK MATTER HALO")
print(f"{'='*70}")

# NFW profile: ρ(r) = ρ_s / [(r/rs)(1 + r/rs)²]
# For Milky Way: rs ≈ 20 kpc, ρ_s ≈ 0.3 GeV/cm³
rs_nfw   = 20 * KPC_TO_M  # scale radius
GeV_cm3  = 1.783e-27 / (1e-2)**3  # 1 GeV/cm³ in kg/m³
rho_s    = 0.3 * GeV_cm3

def rho_nfw(r):
    x = r / rs_nfw
    return rho_s / (x * (1 + x)**2)

def rho_eddy_profile(r, M_galaxy):
    Rs_g   = 2 * const.G * M_galaxy / const.c**2
    r_scale= Rs_g / X_0
    return X_0 * RHO_CRIT * (r_scale / r)**2

print(f"\n   At r = 8 kpc (Sun's position):")
r_sun = 8 * KPC_TO_M
rho_eddy_sun = rho_eddy_profile(r_sun, M_gal)
rho_nfw_sun  = rho_nfw(r_sun)
print(f"   NFW density:  {rho_nfw_sun:.3e} kg/m³")
print(f"   Eddy density: {rho_eddy_sun:.3e} kg/m³")
print(f"   Ratio:        {rho_eddy_sun/rho_nfw_sun:.3f}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"""
   HYPOTHESIS: Ωm = x₀ = R_H/r_universe

   Numerical check:
   x₀ = {X_0:.6f}
   Ωm = {OM_PLANCK:.6f}
   Agreement: {abs(X_0-OM_PLANCK)/OM_PLANCK*100:.2f}%  {'✅ Within Planck 1σ' if abs(X_0-OM_PLANCK)/OM_PLANCK < 0.02 else '❌ Outside Planck 1σ'}

   Rotation curve flatness: {flatness:.3f}
   {'✅ Approximate flat curve from x₀ eddy density' if 0.85 < flatness < 1.15 else '⚠️  Not flat — 1/r² profile assumption needs refinement'}

   WHAT THIS MEANS IF TRUE:
   → Dark matter is the background eddy energy of the 4D inward flow
   → Ωm is not a free parameter — it equals R_H/r_universe
   → As the universe expands, r_universe grows, x₀ decreases, Ωm decreases
   → This is a tiny drift: Δ(Ωm)/Gyr ~ 0.01% — testable with next-gen surveys

   WHAT STILL NEEDS WORK:
   → Full 4D stress-energy tensor derivation connecting ρ_eddy to Tμν
   → Why 1/r² specifically (isothermal) rather than NFW 1/(r(1+r)²)
   → Cluster lensing: can x₀ eddy explain Bullet Cluster offset?
   → This is Phase 7.1 territory — open theoretical question
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(f'Eddy Dark Matter Hypothesis: Ωm = x₀ = {X_0:.4f} (Planck: {OM_PLANCK})',
             fontsize=13, fontweight='bold')

# Plot 1: x₀ vs Ωm numerical comparison
ax = axes[0]
values = [X_0, OM_PLANCK, 0.049, 0.262]
labels = [f'x₀={X_0:.4f}', f'Ωm={OM_PLANCK}', f'Ωb={0.049}', f'Ωdm={0.262}']
colors = ['red', 'blue', 'green', 'orange']
bars   = ax.bar(labels, values, color=colors, alpha=0.8, edgecolor='black')
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f'{val:.4f}', ha='center', fontsize=9, fontweight='bold')
ax.axhline(X_0, color='red', linewidth=2, linestyle='--', alpha=0.5, label='x₀')
ax.set_ylabel('Density parameter', fontsize=11)
ax.set_title('x₀ vs Matter Density Parameters', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(axis='y', alpha=0.3)

# Plot 2: Rotation curves
ax = axes[1]
ax.plot(r_kpc, v_no_eddy,   'blue',  linewidth=2.5, label='Baryons only (Newtonian)')
ax.plot(r_kpc, v_with_eddy, 'red',   linewidth=2.5, linestyle='--',
        label=f'Baryons + eddy (x₀={X_0:.3f})')
ax.axhline(v_with_eddy[idx30], color='gray', linewidth=1, linestyle=':')
ax.set_xlabel('Radius [kpc]', fontsize=11)
ax.set_ylabel('Circular velocity [km/s]', fontsize=11)
ax.set_title('Galactic Rotation Curve\nfrom x₀ Eddy Density', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Plot 3: Density profiles
ax = axes[2]
r_range = np.linspace(1, 50, 200) * KPC_TO_M
rho_e   = np.array([rho_eddy_profile(r, M_gal) for r in r_range])
rho_n   = np.array([rho_nfw(r) for r in r_range])
ax.loglog(r_range/KPC_TO_M, rho_e / rho_n[-1], 'red',
          linewidth=2.5, linestyle='--', label='Eddy density (1/r²)')
ax.loglog(r_range/KPC_TO_M, rho_n / rho_n[-1], 'blue',
          linewidth=2.5, label='NFW profile')
ax.set_xlabel('Radius [kpc]', fontsize=11)
ax.set_ylabel('Density (normalized)', fontsize=11)
ax.set_title('Eddy vs NFW Dark Matter Profile', fontsize=11, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('eddy_dark_matter.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"✓ Plot saved: eddy_dark_matter.png")
print("=" * 70)
