"""
test_jeans_length_eddy.py

Does the eddy background density collapse at galactic scales?

THE SETUP:
The background eddy density is ρ_eddy = x₀ × ρ_crit = 0.3107 × ρ_crit.
This is present everywhere in the universe as a uniform background.

A uniform background with small perturbations is unstable to collapse
when the perturbation scale exceeds the Jeans length:

    λ_Jeans = c_s × √(π / G × ρ)

Perturbations larger than λ_Jeans collapse into structure.
Perturbations smaller than λ_Jeans oscillate away as waves.

THE KEY QUESTION:
What is the effective "sound speed" c_s of the eddy background?

In ordinary gas: c_s = √(kT/m) — thermal pressure resists collapse
In the eddy:     c_s = ? — what resists the eddy from collapsing?

THE CANDIDATE SOUND SPEEDS (three hypotheses):
1. c_s = c × √(x₀)          — tilt geometry sets the pressure
2. c_s = c × obs_now         — observable (√β at x₀) sets it
3. c_s = H₀ × R_H / √(3)    — Hubble flow pressure (cosmological)

For each, compute λ_Jeans and compare to:
- Galactic scale:       1–100 kpc   ← where halos form
- Galaxy cluster scale: 1–10 Mpc    ← large scale structure
- Hubble scale:         ~4000 Mpc   ← entire observable universe

If λ_Jeans falls in the galactic range → Option B works geometrically.
If λ_Jeans is too large → eddy doesn't clump at galactic scales.
If λ_Jeans is too small → eddy clumps too easily (over-collapses).

ALSO CHECKS:
- Jeans mass: what mass of eddy material collapses together?
- Free-fall time: how long does collapse take?
- Compare to observed halo formation timescales (~1–5 Gyr)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ============================================================================
# Constants and setup
# ============================================================================

MPC_TO_M  = 3.085677581e22
KPC_TO_M  = 3.085677581e19
PC_TO_M   = 3.085677581e16
GYR_TO_SEC= 3.15576e16

H0_SI     = const.H_0
R_U_0     = 4.4e26
X_0       = (const.c / H0_SI) / R_U_0      # = 0.3107
OBS_NOW   = estif.observable_combined(X_0)  # = 0.8300
RHO_CRIT  = 3 * H0_SI**2 / (8 * np.pi * const.G)
RHO_EDDY  = X_0 * RHO_CRIT                  # Background eddy density

print("=" * 70)
print("JEANS LENGTH OF THE EDDY DARK MATTER BACKGROUND")
print("=" * 70)
print(f"\n   x₀         = {X_0:.6f}")
print(f"   obs_now    = {OBS_NOW:.6f}")
print(f"   ρ_crit     = {RHO_CRIT:.4e} kg/m³")
print(f"   ρ_eddy     = x₀ × ρ_crit = {RHO_EDDY:.4e} kg/m³")
print(f"   ρ_eddy/ρ_b = {RHO_EDDY/(0.049*RHO_CRIT):.2f}  (ratio to baryon density)")

# ============================================================================
# Section 1: The Jeans Length Formula
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: JEANS LENGTH FOR THREE SOUND SPEED HYPOTHESES")
print(f"{'='*70}")

print(f"""
   λ_Jeans = c_s × √(π / G × ρ_eddy)

   This is the critical scale — perturbations LARGER than λ_Jeans collapse.
   Perturbations SMALLER than λ_Jeans oscillate and disperse.

   The only unknown is c_s — the effective sound speed of the eddy.
""")

# Three candidate sound speeds
cs_candidates = [
    ("H1: c_s = c × √x₀",        const.c * np.sqrt(X_0),         "tilt geometry pressure"),
    ("H2: c_s = c × obs_now",     const.c * OBS_NOW,              "observable sets pressure"),
    ("H3: c_s = H₀ × R_H / √3",  H0_SI * (const.c/H0_SI) / np.sqrt(3), "Hubble flow"),
    ("H4: c_s = c × x₀",         const.c * X_0,                  "direct curvature"),
    ("H5: c_s = c × x₀² ",       const.c * X_0**2,               "curvature squared"),
]

print(f"   {'Hypothesis':<28} {'c_s [m/s]':<16} {'c_s/c':<10} {'Description'}")
print("   " + "-"*70)

jeans_results = []
for name, cs, desc in cs_candidates:
    ratio = cs / const.c
    lambda_j = cs * np.sqrt(np.pi / (const.G * RHO_EDDY))
    lambda_kpc = lambda_j / KPC_TO_M
    lambda_mpc = lambda_j / MPC_TO_M
    jeans_results.append((name, cs, ratio, lambda_j, lambda_kpc, lambda_mpc, desc))
    print(f"   {name:<28} {cs:<16.4e} {ratio:<10.6f} {desc}")

print(f"\n   {'Hypothesis':<28} {'λ_Jeans [kpc]':<18} {'λ_Jeans [Mpc]':<16} {'Scale regime'}")
print("   " + "-"*68)

for name, cs, ratio, lambda_j, lambda_kpc, lambda_mpc, desc in jeans_results:
    if lambda_kpc < 100:
        regime = "✅ Galactic (1–100 kpc)"
    elif lambda_kpc < 1000:
        regime = "⚠️  Sub-cluster (100 kpc–1 Mpc)"
    elif lambda_kpc < 10000:
        regime = "⚠️  Cluster (1–10 Mpc)"
    else:
        regime = "❌ Cosmic (>10 Mpc — too large)"
    print(f"   {name:<28} {lambda_kpc:<18.2f} {lambda_mpc:<16.4f} {regime}")

# ============================================================================
# Section 2: Jeans Mass
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: JEANS MASS — WHAT MASS COLLAPSES TOGETHER?")
print(f"{'='*70}")

print(f"""
   M_Jeans = (4π/3) × ρ_eddy × (λ_Jeans/2)³
   This is the characteristic mass that collapses into one halo.

   Reference masses:
   Milky Way:          ~10¹² M☉
   Dwarf galaxy:       ~10⁸ M☉
   Galaxy cluster:     ~10¹⁴ M☉
   Observable universe: ~10²³ M☉
""")

print(f"   {'Hypothesis':<28} {'M_Jeans [M☉]':<18} {'log₁₀(M)':<12} {'Structure match'}")
print("   " + "-"*72)

for name, cs, ratio, lambda_j, lambda_kpc, lambda_mpc, desc in jeans_results:
    M_jeans = (4*np.pi/3) * RHO_EDDY * (lambda_j/2)**3
    M_solar = M_jeans / const.M_sun
    log_M   = np.log10(max(M_solar, 1e-10))

    if 8 < log_M < 10:
        match = "✅ Dwarf galaxy"
    elif 10 < log_M < 13:
        match = "✅ Milky Way-scale"
    elif 13 < log_M < 15:
        match = "✅ Galaxy cluster"
    elif log_M < 8:
        match = "⚠️  Sub-galactic"
    else:
        match = "❌ Too massive"

    print(f"   {name:<28} {M_solar:<18.3e} {log_M:<12.2f} {match}")

# ============================================================================
# Section 3: Free-fall Time
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: FREE-FALL TIME — HOW LONG DOES COLLAPSE TAKE?")
print(f"{'='*70}")

t_ff = np.sqrt(3*np.pi / (32 * const.G * RHO_EDDY))
t_ff_gyr = t_ff / GYR_TO_SEC

print(f"""
   t_freefall = √(3π / 32Gρ)   [depends only on density, not sound speed]

   t_freefall = {t_ff:.4e} s
              = {t_ff_gyr:.3f} Gyr

   Reference timescales:
   First galaxies formed:  ~0.5 Gyr after Big Bang
   Milky Way formed:       ~1–2 Gyr after Big Bang
   Universe age:            13.8 Gyr

""")

if t_ff_gyr < 1.0:
    print(f"   ✅ Free-fall time {t_ff_gyr:.2f} Gyr — consistent with early galaxy formation")
elif t_ff_gyr < 5.0:
    print(f"   ✅ Free-fall time {t_ff_gyr:.2f} Gyr — consistent with galaxy formation era")
elif t_ff_gyr < 13.0:
    print(f"   ⚠️  Free-fall time {t_ff_gyr:.2f} Gyr — slow but possible")
else:
    print(f"   ❌ Free-fall time {t_ff_gyr:.2f} Gyr — longer than universe age, too slow")

# ============================================================================
# Section 4: The NFW comparison — what profile does eddy collapse produce?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: COLLAPSE PROFILE — WHAT SHAPE DOES THE HALO TAKE?")
print(f"{'='*70}")

print(f"""
   When a uniform background collapses gravitationally, the density
   profile of the resulting halo depends on whether collapse is:

   ISOTHERMAL (constant temperature analogy):
   → ρ(r) ∝ 1/r²  — gives perfectly FLAT rotation curves ✅
   → This is the "singular isothermal sphere" — standard in cosmology

   ADIABATIC (no heat exchange):
   → ρ(r) ∝ 1/(r(1+r)²)  — NFW profile
   → Gives rotation curves that PEAK then decline ⚠️

   For the eddy background, the relevant question is:
   Can eddy perturbations exchange energy with the background flow?

   If YES (eddy is "thermally coupled" to the 4D flow):
   → Isothermal collapse → 1/r² profile → flat rotation curves

   If NO (eddy perturbations are adiabatic):
   → NFW collapse → peaked rotation curves (like what we observe
     in simulations but need NFW correction in real galaxies)

   The eddy IS the background flow — perturbations in the eddy
   are embedded in what they're collapsing into. This suggests
   isothermal is the more natural regime.

   PREDICTION: If the eddy collapses isothermally,
   flat rotation curves emerge without any additional physics.
   The rotation velocity at large radius:

   v_flat = √(4πG × ρ_eddy_0 × r_scale²)

   where r_scale = Rs_galaxy / x₀  (the natural cutoff from the tilt)
""")

# Compute v_flat for Milky Way
M_MW     = 1e12 * const.M_sun
Rs_MW    = 2 * const.G * M_MW / const.c**2
r_scale  = Rs_MW / X_0
rho_0_eddy = RHO_EDDY
v_flat   = np.sqrt(4 * np.pi * const.G * rho_0_eddy * r_scale**2)

print(f"   Milky Way example:")
print(f"   M_MW = 10¹² M☉,  Rs = {Rs_MW:.3e} m")
print(f"   r_scale = Rs/x₀ = {r_scale/KPC_TO_M:.2e} kpc")
print(f"   v_flat = {v_flat/1e3:.1f} km/s")
print(f"   Observed MW flat velocity: ~220 km/s")
print(f"   Ratio: {v_flat/1e3/220:.3f}  "
      f"({'✅ same order' if 0.1 < v_flat/1e3/220 < 10 else '❌ wrong order of magnitude'})")

print(f"""
   NOTE on r_scale:
   The tilt formula naturally provides a characteristic halo radius
   r_halo = Rs/x₀ for each galaxy. This is where x_local(r) = x₀,
   i.e. where local curvature equals the cosmic background curvature.
   Inside this radius, the local tilt dominates.
   Outside, the background eddy dominates.
""")

# ============================================================================
# Section 5: The Sahara Analogy — Why Fluid Dynamics Works
# ============================================================================

print(f"{'='*70}")
print("SECTION 5: THE FLUID DYNAMICS ARGUMENT")
print(f"{'='*70}")

print(f"""
   SAND IN THE SAHARA:
   Individual grains: no attraction between them
   Collective: flows, forms dunes, creates patterns
   Mechanism: pressure from weight above + wind perturbation

   EDDY BACKGROUND IN SPACE:
   Individual eddy quanta: no direct particle interaction
   Collective: gravitational attraction between density perturbations
   Mechanism: ρ_eddy × G creates attraction + Hubble flow as perturbation

   JEANS INSTABILITY (the "pressure vs gravity" balance):

   At scale λ < λ_Jeans:  pressure wins → perturbations oscillate
   At scale λ > λ_Jeans:  gravity wins  → perturbations collapse

   This is identical to how Sahara sand forms dunes at the "right" scale:
   too small → grains bounce off each other
   too large → the whole desert doesn't move as one
   just right → dunes form at the scale where wind pressure equals
                the gravitational weight of the pile

   For the eddy background:
   ρ_eddy    = {RHO_EDDY:.4e} kg/m³
   t_ff      = {t_ff_gyr:.3f} Gyr
""")

# ============================================================================
# Section 6: Overall verdict
# ============================================================================

print(f"{'='*70}")
print("SECTION 6: VERDICT — DOES OPTION B WORK?")
print(f"{'='*70}")

# Find the best hypothesis
galactic_hits = [r for r in jeans_results if 1 < r[4] < 1000]

print(f"\n   Hypotheses producing galactic-scale Jeans lengths (1–1000 kpc):")
if galactic_hits:
    for name, cs, ratio, lj, lkpc, lmpc, desc in galactic_hits:
        M_j = (4*np.pi/3) * RHO_EDDY * (lj/2)**3 / const.M_sun
        print(f"   ✅ {name}: λ = {lkpc:.1f} kpc,  M = {M_j:.2e} M☉")
else:
    print(f"   ⚠️  No hypothesis produces galactic-scale collapse directly.")
    print(f"   The simplest sound speed candidates give too large a scale.")

print(f"""
   Free-fall time: {t_ff_gyr:.2f} Gyr
   {'✅ Within galaxy formation era' if t_ff_gyr < 5 else '⚠️  Longer than expected'}

   CONCLUSION:
""")

if galactic_hits and t_ff_gyr < 5:
    print("   ✅ Option B is physically viable.")
    print("   The eddy background has the right density and timescale.")
    print("   The sound speed hypothesis needs to be derived theoretically.")
    print("   Next step: derive c_s from the 4D tilt geometry.")
elif t_ff_gyr < 5:
    print("   ⚠️  Timescale is right but Jeans length needs refinement.")
    print("   The sound speed in the eddy is not yet derived from first principles.")
    print("   If c_s < c × x₀, the Jeans length shrinks to galactic scales.")
    print("   Next step: constrain c_s from the tilt geometry.")
else:
    print("   ⚠️  Both scale and timescale need refinement.")
    print("   The eddy density is correct — the sound speed is the unknown.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(f'Eddy Dark Matter: Jeans Analysis  (ρ_eddy = x₀ × ρ_crit, t_ff = {t_ff_gyr:.2f} Gyr)',
             fontsize=13, fontweight='bold')

# Plot 1: Jeans length vs sound speed
cs_range = np.logspace(-6, 0, 200) * const.c
lambda_range = cs_range * np.sqrt(np.pi / (const.G * RHO_EDDY))

ax = axes[0]
ax.loglog(cs_range/const.c, lambda_range/KPC_TO_M, 'black', linewidth=2.5)
ax.axhspan(1, 100, alpha=0.15, color='green', label='Galactic halo (1–100 kpc)')
ax.axhspan(100, 1000, alpha=0.1, color='yellow', label='Sub-cluster (0.1–1 Mpc)')
ax.axhspan(1000, 10000, alpha=0.1, color='orange', label='Cluster (1–10 Mpc)')

# Mark each hypothesis
colors_h = ['red', 'blue', 'green', 'purple', 'brown']
for i, (name, cs, ratio, lj, lkpc, lmpc, desc) in enumerate(jeans_results):
    ax.scatter([ratio], [lkpc], color=colors_h[i], s=150, zorder=5,
               label=f'{name}: {lkpc:.0f} kpc')

ax.set_xlabel('c_s / c  (sound speed fraction)', fontsize=11)
ax.set_ylabel('Jeans Length [kpc]', fontsize=11)
ax.set_title('Jeans Length vs Sound Speed', fontsize=11, fontweight='bold')
ax.legend(fontsize=7, loc='upper left'); ax.grid(alpha=0.3, which='both')

# Plot 2: Jeans mass vs sound speed
M_range = (4*np.pi/3) * RHO_EDDY * (lambda_range/2)**3 / const.M_sun

ax = axes[1]
ax.loglog(cs_range/const.c, M_range, 'black', linewidth=2.5)
ax.axhspan(1e8, 1e10, alpha=0.15, color='green', label='Dwarf–MW galaxy')
ax.axhspan(1e10, 1e13, alpha=0.1, color='yellow', label='Large galaxy')
ax.axhspan(1e13, 1e15, alpha=0.1, color='orange', label='Galaxy cluster')

for i, (name, cs, ratio, lj, lkpc, lmpc, desc) in enumerate(jeans_results):
    M_j = (4*np.pi/3) * RHO_EDDY * (lj/2)**3 / const.M_sun
    ax.scatter([ratio], [M_j], color=colors_h[i], s=150, zorder=5,
               label=f'{name}: {np.log10(max(M_j,1)):.1f}')

ax.set_xlabel('c_s / c', fontsize=11)
ax.set_ylabel('Jeans Mass [M☉]', fontsize=11)
ax.set_title('Jeans Mass vs Sound Speed', fontsize=11, fontweight='bold')
ax.legend(fontsize=7, loc='upper left'); ax.grid(alpha=0.3, which='both')

# Plot 3: What sound speed do we NEED for galactic halos?
target_scales_kpc = [1, 10, 50, 100, 200, 500, 1000]
target_scales_m   = [s * KPC_TO_M for s in target_scales_kpc]
required_cs       = [lj / np.sqrt(np.pi / (const.G * RHO_EDDY))
                     for lj in target_scales_m]

ax = axes[2]
ax.semilogy(target_scales_kpc, [c/const.c for c in required_cs],
            'black', linewidth=2.5, marker='o', markersize=6)
ax.axhspan(X_0**2, X_0, alpha=0.15, color='green',
           label=f'x₀² to x₀ ({X_0**2:.4f}–{X_0:.4f})')
ax.axhline(X_0,     color='red',    linewidth=2, linestyle='--', label=f'x₀ = {X_0:.4f}')
ax.axhline(X_0**2,  color='blue',   linewidth=2, linestyle='--', label=f'x₀² = {X_0**2:.4f}')
ax.axhline(OBS_NOW, color='purple', linewidth=2, linestyle=':', label=f'obs_now = {OBS_NOW:.4f}')
ax.axvspan(1, 100, alpha=0.1, color='green', label='Target galactic range')
ax.set_xlabel('Target Jeans length [kpc]', fontsize=11)
ax.set_ylabel('Required c_s / c', fontsize=11)
ax.set_title('Sound Speed Required\nfor Each Halo Scale', fontsize=11, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('jeans_length_eddy.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: jeans_length_eddy.png")
print("=" * 70)
