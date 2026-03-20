"""
test_nmax_investigation.py

Investigates whether N_MAX = 33.265 and B = 15.429 have
a natural physical interpretation.

QUESTION: Why does the tilt exponent in flat space equal 33.265?
Is this a meaningless fitting artifact, or does it relate to
known physical constants or geometric ratios?

APPROACH:
1. Check ratios of known constants and geometric quantities
2. Check if N_MAX and B are related to each other simply
3. Check if N_MAX connects to the Hubble radius / Schwarzschild geometry
4. Look for integer or simple fractional relationships
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Calibrated values
# ============================================================================

N_MAX = 33.2650
B     = 15.4291

# ============================================================================
# Known physical and geometric quantities
# ============================================================================

R_H         = const.c / const.H_0
r_universe  = 4.4e26
M_m87       = 6.5e9 * const.M_sun
Rs_m87      = 2 * const.G * M_m87 / const.c**2
r_photon    = 1.5 * Rs_m87

CURV_LOCAL  = Rs_m87 / r_photon       # 0.6667
CURV_COSM   = R_H / r_universe        # 0.3107

# ============================================================================
# Investigation 1: Simple relationships
# ============================================================================

print("=" * 70)
print("N_MAX INVESTIGATION: WHAT IS 33.265?")
print("=" * 70)

print(f"\nCalibrated values:")
print(f"   N_MAX = {N_MAX:.6f}")
print(f"   B     = {B:.6f}")

print(f"\n{'='*70}")
print("SECTION 1: SIMPLE RATIOS AND RELATIONSHIPS")
print(f"{'='*70}")

print(f"\n   N_MAX / 2π        = {N_MAX / (2*np.pi):.4f}  (= {N_MAX/(2*np.pi):.2f})")
print(f"   N_MAX / π         = {N_MAX / np.pi:.4f}  (= {N_MAX/np.pi:.2f})")
print(f"   N_MAX / e         = {N_MAX / np.e:.4f}  (= {N_MAX/np.e:.2f})")
print(f"   √N_MAX            = {np.sqrt(N_MAX):.4f}  (= {np.sqrt(N_MAX):.2f})")
print(f"   N_MAX²            = {N_MAX**2:.2f}")
print(f"   ln(N_MAX)         = {np.log(N_MAX):.4f}")
print(f"   log10(N_MAX)      = {np.log10(N_MAX):.4f}")
print(f"   1/N_MAX           = {1/N_MAX:.6f}")

print(f"\n   B / 2π            = {B / (2*np.pi):.4f}")
print(f"   B / π             = {B / np.pi:.4f}")
print(f"   B / e             = {B / np.e:.4f}")
print(f"   N_MAX / B         = {N_MAX/B:.4f}")
print(f"   B / N_MAX         = {B/N_MAX:.4f}")
print(f"   B / (2 × N_MAX)   = {B/(2*N_MAX):.4f}")

# ============================================================================
# Investigation 2: Connection to curvature scales
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: CONNECTION TO CURVATURE SCALES")
print(f"{'='*70}")

# At the activation point where n = 1, x = ln(N_MAX)/B
x_activation = np.log(N_MAX) / B
print(f"\n   n(x) = N_MAX × exp(-B × x) = 1  when:")
print(f"   x_activation = ln(N_MAX)/B = {x_activation:.4f}")
print(f"\n   Compare to:")
print(f"   CURVATURE_COSM  = {CURV_COSM:.4f}  (cosmological scale)")
print(f"   CURVATURE_LOCAL = {CURV_LOCAL:.4f}  (M87* photon sphere)")
print(f"   x_activation / CURV_COSM = {x_activation/CURV_COSM:.4f}")

# The formula n=1 is where the correction "turns on" meaningfully
# Check if this relates to the geometric mean of the two scales
geo_mean = np.sqrt(CURV_LOCAL * CURV_COSM)
print(f"\n   Geometric mean of scales: √(0.667 × 0.311) = {geo_mean:.4f}")
print(f"   x_activation / geo_mean = {x_activation/geo_mean:.4f}")

# Check N_MAX in terms of curvature ratio
print(f"\n   N_MAX as function of curvature ratio:")
print(f"   1/CURV_COSM     = {1/CURV_COSM:.4f}")
print(f"   1/CURV_LOCAL    = {1/CURV_LOCAL:.4f}")
print(f"   N_MAX × CURV_COSM  = {N_MAX * CURV_COSM:.4f}")
print(f"   N_MAX × CURV_LOCAL = {N_MAX * CURV_LOCAL:.4f}")

# ============================================================================
# Investigation 3: Geometric origin — could N_MAX = ln(r_universe/Rs) ?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: GEOMETRIC SCALE RATIOS")
print(f"{'='*70}")

# The most natural geometric quantity in the theory
# is the ratio of the largest scale to the smallest
# rs of sun, neutron star, stellar BH, M87*
Rs_sun = 2 * const.G * const.M_sun / const.c**2

candidates = {
    "ln(r_universe / Rs_m87)":  np.log(r_universe / Rs_m87),
    "ln(r_universe / R_H)":     np.log(r_universe / R_H),
    "ln(R_H / Rs_m87)":         np.log(R_H / Rs_m87),
    "ln(R_H / r_photon)":       np.log(R_H / r_photon),
    "log10(r_universe/Rs_m87)": np.log10(r_universe / Rs_m87),
    "log10(R_H / Rs_m87)":      np.log10(R_H / Rs_m87),
    "½ × ln(r_universe/Rs_m87)":0.5*np.log(r_universe / Rs_m87),
    "½ × ln(R_H / r_photon)":   0.5*np.log(R_H / r_photon),
    "ln(1/CURV_COSM)":          np.log(1/CURV_COSM),
    "ln(1/CURV_LOCAL)":         np.log(1/CURV_LOCAL),
    "ln(CURV_LOCAL/CURV_COSM)": np.log(CURV_LOCAL/CURV_COSM),
    "1/CURV_COSM - 1":          1/CURV_COSM - 1,
    "1/(CURV_COSM × CURV_LOCAL)":1/(CURV_COSM * CURV_LOCAL),
    "CURV_LOCAL/CURV_COSM²":    CURV_LOCAL/CURV_COSM**2,
}

print(f"\n   Searching for quantity ≈ N_MAX = {N_MAX:.4f}")
print(f"\n   {'Quantity':<35} {'Value':<12} {'Ratio to N_MAX'}")
print("   " + "-"*65)

matches = []
for name, val in candidates.items():
    ratio = val / N_MAX
    flag  = "  ← CLOSE!" if 0.95 < ratio < 1.05 else ""
    if 0.95 < ratio < 1.05:
        matches.append((name, val, ratio))
    print(f"   {name:<35} {val:<12.4f} {ratio:.4f}{flag}")

# ============================================================================
# Investigation 4: B in terms of same quantities
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: INVESTIGATING B = 15.429")
print(f"{'='*70}")

b_candidates = {
    "ln(r_universe / Rs_m87)":  np.log(r_universe / Rs_m87),
    "½ × ln(r_universe/Rs_m87)":0.5*np.log(r_universe / Rs_m87),
    "ln(R_H / r_photon)":       np.log(R_H / r_photon),
    "½ × ln(R_H / r_photon)":   0.5*np.log(R_H / r_photon),
    "log10(r_universe/R_H)":    np.log10(r_universe/R_H),
    "ln(r_universe/R_H)":       np.log(r_universe/R_H),
    "½ × ln(r_universe/R_H)":   0.5*np.log(r_universe/R_H),
    "N_MAX / 2":                N_MAX / 2,
    "N_MAX × CURV_COSM × π":   N_MAX * CURV_COSM * np.pi,
    "2π / CURV_COSM":           2*np.pi / CURV_COSM,
    "1/(CURV_COSM²)":           1/CURV_COSM**2,
    "ln(N_MAX) / CURV_COSM":    np.log(N_MAX) / CURV_COSM,
    "ln(N_MAX) / CURV_LOCAL":   np.log(N_MAX) / CURV_LOCAL,
}

print(f"\n   Searching for quantity ≈ B = {B:.4f}")
print(f"\n   {'Quantity':<35} {'Value':<12} {'Ratio to B'}")
print("   " + "-"*65)

b_matches = []
for name, val in b_candidates.items():
    ratio = val / B
    flag  = "  ← CLOSE!" if 0.95 < ratio < 1.05 else ""
    if 0.95 < ratio < 1.05:
        b_matches.append((name, val, ratio))
    print(f"   {name:<35} {val:<12.4f} {ratio:.4f}{flag}")

# ============================================================================
# Investigation 5: The ratio N_MAX/B
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 5: THE RATIO N_MAX/B AND ITS MEANING")
print(f"{'='*70}")

ratio_nb = N_MAX / B
print(f"\n   N_MAX / B = {ratio_nb:.4f}")
print(f"\n   This ratio equals x_activation (where n=1)")
print(f"   x_activation = ln(N_MAX)/B = {np.log(N_MAX)/B:.4f}")
print(f"\n   Hmm — is N_MAX/B ≈ ln(N_MAX)/B? No: ln(N_MAX)/B = {np.log(N_MAX)/B:.4f}")

# The characteristic curvature where n crosses 1
print(f"\n   Physical meaning of x_activation = {x_activation:.4f}:")
print(f"   Below this curvature: n > 1 (formula in 'strong' regime)")
print(f"   Above this curvature: n < 1 (formula in 'suppressed' regime)")
print(f"\n   Note: CURV_COSM = {CURV_COSM:.4f} < x_activation = {x_activation:.4f}")
print(f"         So cosmic scale is BELOW the activation threshold")
print(f"         meaning n > 1 cosmologically → very gentle tilt")

# ============================================================================
# Investigation 6: Does N_MAX relate to the number of e-folds?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 6: COSMOLOGICAL E-FOLDS CONNECTION")
print(f"{'='*70}")

# Number of e-folds of inflation ~ 60
# Number of e-folds of expansion since Big Bang
z_cmb  = 1100
efolds_expansion = np.log(1 + z_cmb)
efolds_inflation  = 60  # typical

print(f"\n   e-folds since CMB (ln(1+z_cmb)):  {efolds_expansion:.4f}")
print(f"   e-folds of inflation (typical):    {efolds_inflation}")
print(f"   N_MAX / e-folds CMB:               {N_MAX/efolds_expansion:.4f}")
print(f"   N_MAX / e-folds inflation:         {N_MAX/efolds_inflation:.4f}")
print(f"   N_MAX / (e-folds × ½):             {N_MAX/(efolds_inflation*0.5):.4f}")

# Check ln(N_MAX)
print(f"\n   ln(N_MAX) = {np.log(N_MAX):.4f}")
print(f"   This is the 'e-folding scale' of the formula")
print(f"   N_MAX × exp(-B × x) crosses 1 at x = ln(N_MAX)/B = {np.log(N_MAX)/B:.4f}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"\n   N_MAX = {N_MAX:.4f}")
print(f"   B     = {B:.4f}")
print(f"   N_MAX/B = {N_MAX/B:.4f}  (activation curvature)")
print(f"   ln(N_MAX) = {np.log(N_MAX):.4f}")

if matches:
    print(f"\n   ✅ N_MAX matches found:")
    for name, val, ratio in matches:
        print(f"   {name} = {val:.4f}  (ratio={ratio:.4f})")
else:
    print(f"\n   ⚠️  No exact matches found for N_MAX = {N_MAX:.4f}")
    print(f"   Closest candidates:")
    # Find top 3 closest
    all_cands = [(abs(v/N_MAX - 1), name, v)
                 for name, v in candidates.items() if v > 0]
    all_cands.sort()
    for diff, name, val in all_cands[:3]:
        print(f"   {name} = {val:.4f}  ({diff*100:.1f}% off)")

if b_matches:
    print(f"\n   ✅ B matches found:")
    for name, val, ratio in b_matches:
        print(f"   {name} = {val:.4f}  (ratio={ratio:.4f})")
else:
    print(f"\n   ⚠️  No exact matches found for B = {B:.4f}")
    all_b = [(abs(v/B - 1), name, v)
             for name, v in b_candidates.items() if v > 0]
    all_b.sort()
    print(f"   Closest candidates:")
    for diff, name, val in all_b[:3]:
        print(f"   {name} = {val:.4f}  ({diff*100:.1f}% off)")

print(f"\n   HONEST ASSESSMENT:")
print(f"   N_MAX = 33.265 is a calibration result, not yet derived")
print(f"   from first principles. Whether it relates to a known")
print(f"   geometric or physical quantity is still open.")
print(f"\n   The formula works. The physical meaning of N_MAX")
print(f"   is the next theoretical question to answer.")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('N_MAX Investigation: Physical Meaning of 33.265',
             fontsize=13, fontweight='bold')

x = np.linspace(0, 1, 500)

# Plot n(x) with annotations
ax = axes[0]
n_vals = N_MAX * np.exp(-B * x)
ax.plot(x, n_vals, 'purple', linewidth=3, label=f'n(x)={N_MAX:.3f}×exp(-{B:.3f}×x)')
ax.axhline(N_MAX, color='gray', linewidth=1, linestyle=':', alpha=0.7,
           label=f'N_MAX = {N_MAX:.3f} (flat space)')
ax.axhline(1.0,   color='orange', linewidth=1.5, linestyle='--',
           label='n = 1 (activation threshold)')
ax.axhline(0.5,   color='blue', linewidth=1, linestyle='--', alpha=0.5,
           label='n = 0.5')
ax.axvline(x_activation, color='orange', linewidth=1.5, linestyle=':',
           label=f'x_act = {x_activation:.3f}')
ax.axvline(CURV_COSM,    color='green',  linewidth=2, linestyle=':',
           label=f'Cosmic ({CURV_COSM:.3f})')
ax.axvline(CURV_LOCAL,   color='red',    linewidth=2, linestyle=':',
           label=f'M87* ({CURV_LOCAL:.3f})')
ax.set_xlabel('Curvature Ratio x', fontsize=12)
ax.set_ylabel('Dynamic n', fontsize=12)
ax.set_title('n(x): Shape and Key Points', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_xlim(0, 1); ax.set_ylim(0, min(N_MAX*1.1, 5))

# Plot scale ratios vs N_MAX
ax = axes[1]
scale_names = [
    'ln(R_H/r_photon)',
    'ln(r_univ/R_H)',
    'ln(r_univ/Rs_m87)',
    '½ln(r_univ/Rs_m87)',
    '½ln(R_H/r_photon)',
    'log₁₀(r_univ/Rs)',
    '1/x_cosm - 1',
]
scale_vals = [
    np.log(R_H/r_photon),
    np.log(r_universe/R_H),
    np.log(r_universe/Rs_m87),
    0.5*np.log(r_universe/Rs_m87),
    0.5*np.log(R_H/r_photon),
    np.log10(r_universe/Rs_m87),
    1/CURV_COSM - 1,
]
colors_bar = ['green' if 0.95 < v/N_MAX < 1.05 else 'steelblue'
              for v in scale_vals]

bars = ax.barh(scale_names, scale_vals, color=colors_bar, alpha=0.8,
               edgecolor='black')
ax.axvline(N_MAX, color='red', linewidth=2, linestyle='--',
           label=f'N_MAX = {N_MAX:.3f}')
ax.axvspan(N_MAX*0.95, N_MAX*1.05, alpha=0.15, color='red',
           label='±5% band')
ax.set_xlabel('Value', fontsize=12)
ax.set_title('Scale Ratios vs N_MAX', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig('nmax_investigation.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: nmax_investigation.png")
print("=" * 70)
