"""
test_planck_ruler.py

Tests whether N_MAX = 33.265 and B = 15.429 can be expressed as
combinations of Planck-scale ratios — using the Planck length as
a rigid, autonomous measuring stick.

THE INSIGHT:
Previous attempts used r_universe and Rs_m87 as reference scales.
Both are measured with rulers affected by the very phenomena ESTIF
is trying to describe — a bent ruler problem.

The Planck length l_P = √(ħG/c³) is built from fundamental constants
only. It cannot be bent by gravity or stretched by expansion.
It is the same everywhere in the universe at all times.

THE TEST:
Can N_MAX = 33.265 and B = 15.429 be expressed as:
    ln(some_scale / l_P)  or fractions thereof
where some_scale is a physically meaningful length?
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

N_MAX = 33.265
B     = 15.429

# ============================================================================
# Fundamental constants — the rigid rulers
# ============================================================================

hbar  = 1.054571817e-34   # Reduced Planck constant (J·s)
G     = const.G           # Gravitational constant
c     = const.c           # Speed of light

# Planck units — built from fundamentals only
l_P   = np.sqrt(hbar * G / c**3)    # Planck length  ~1.616e-35 m
m_P   = np.sqrt(hbar * c / G)       # Planck mass    ~2.176e-8 kg
t_P   = np.sqrt(hbar * G / c**5)    # Planck time    ~5.391e-44 s
E_P   = m_P * c**2                  # Planck energy  ~1.956e9 J

# ============================================================================
# Physical scales to test as numerators
# ============================================================================

R_H         = const.c / const.H_0           # Hubble radius
r_universe  = 4.4e26                         # Observable universe
M_m87       = 6.5e9 * const.M_sun
Rs_m87      = 2 * G * M_m87 / c**2          # M87* Schwarzschild radius
r_photon    = 1.5 * Rs_m87                  # M87* photon sphere

# Proton and electron for comparison
r_proton    = 8.41e-16                       # Proton charge radius (m)
r_electron  = 2.82e-15                       # Classical electron radius (m)

# ============================================================================
# Print Planck units
# ============================================================================

print("=" * 70)
print("PLANCK RULER TEST: AUTONOMOUS MEASURING STICKS")
print("=" * 70)

print(f"\nFundamental constants:")
print(f"   ħ = {hbar:.4e} J·s")
print(f"   G = {G:.4e} m³/(kg·s²)")
print(f"   c = {c:.4e} m/s")

print(f"\nPlanck units (rigid rulers):")
print(f"   l_P = √(ħG/c³) = {l_P:.4e} m")
print(f"   m_P = √(ħc/G)  = {m_P:.4e} kg")
print(f"   t_P = √(ħG/c⁵) = {t_P:.4e} s")

print(f"\nCalibrated targets:")
print(f"   N_MAX = {N_MAX:.6f}")
print(f"   B     = {B:.6f}")
print(f"   N_MAX/B = {N_MAX/B:.6f}")
print(f"   ln(N_MAX/0.5)/B = {np.log(N_MAX/0.5)/B:.6f} (GR crossover x)")

# ============================================================================
# Section 1: ln(scale / l_P) candidates for N_MAX
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 1: ln(scale / l_P) CANDIDATES FOR N_MAX = 33.265")
print(f"{'='*70}")

candidates_nmax = {
    "ln(r_universe / l_P)":        np.log(r_universe / l_P),
    "ln(R_H / l_P)":               np.log(R_H / l_P),
    "ln(Rs_m87 / l_P)":            np.log(Rs_m87 / l_P),
    "ln(r_photon / l_P)":          np.log(r_photon / l_P),
    "½ × ln(r_universe / l_P)":    0.5 * np.log(r_universe / l_P),
    "½ × ln(R_H / l_P)":           0.5 * np.log(R_H / l_P),
    "½ × ln(Rs_m87 / l_P)":        0.5 * np.log(Rs_m87 / l_P),
    "⅓ × ln(r_universe / l_P)":    (1/3) * np.log(r_universe / l_P),
    "⅓ × ln(R_H / l_P)":           (1/3) * np.log(R_H / l_P),
    "⅓ × ln(Rs_m87 / l_P)":        (1/3) * np.log(Rs_m87 / l_P),
    "¼ × ln(r_universe / l_P)":    0.25 * np.log(r_universe / l_P),
    "¼ × ln(R_H / l_P)":           0.25 * np.log(R_H / l_P),
    "ln(r_proton / l_P)":          np.log(r_proton / l_P),
    "ln(r_electron / l_P)":        np.log(r_electron / l_P),
    "½ × ln(r_proton / l_P)":      0.5 * np.log(r_proton / l_P),
    "ln(r_universe / r_proton)":   np.log(r_universe / r_proton),
    "½ × ln(r_universe/r_proton)": 0.5 * np.log(r_universe / r_proton),
    "ln(R_H / r_proton)":          np.log(R_H / r_proton),
    "½ × ln(R_H / r_proton)":      0.5 * np.log(R_H / r_proton),
    "ln(Rs_m87 / r_proton)":       np.log(Rs_m87 / r_proton),
}

print(f"\n   {'Quantity':<38} {'Value':<12} {'Ratio to N_MAX':<16} {'% off'}")
print("   " + "-"*75)

nmax_matches = []
for name, val in candidates_nmax.items():
    ratio = val / N_MAX
    pct   = abs(ratio - 1) * 100
    flag  = "  ✅ MATCH!" if pct < 1.0 else ("  ← close" if pct < 5.0 else "")
    if pct < 1.0:
        nmax_matches.append((name, val, ratio, pct))
    print(f"   {name:<38} {val:<12.4f} {ratio:<16.4f} {pct:.2f}%{flag}")

# ============================================================================
# Section 2: ln(scale / l_P) candidates for B
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 2: ln(scale / l_P) CANDIDATES FOR B = 15.429")
print(f"{'='*70}")

candidates_b = {
    "ln(r_universe / l_P)":        np.log(r_universe / l_P),
    "ln(R_H / l_P)":               np.log(R_H / l_P),
    "ln(Rs_m87 / l_P)":            np.log(Rs_m87 / l_P),
    "½ × ln(r_universe / l_P)":    0.5 * np.log(r_universe / l_P),
    "½ × ln(R_H / l_P)":           0.5 * np.log(R_H / l_P),
    "½ × ln(Rs_m87 / l_P)":        0.5 * np.log(Rs_m87 / l_P),
    "¼ × ln(r_universe / l_P)":    0.25 * np.log(r_universe / l_P),
    "¼ × ln(R_H / l_P)":           0.25 * np.log(R_H / l_P),
    "⅓ × ln(R_H / l_P)":           (1/3) * np.log(R_H / l_P),
    "½ × ln(r_proton / l_P)":      0.5 * np.log(r_proton / l_P),
    "⅓ × ln(r_proton / l_P)":      (1/3) * np.log(r_proton / l_P),
    "¼ × ln(r_proton / l_P)":      0.25 * np.log(r_proton / l_P),
    "ln(R_H / r_proton)":          np.log(R_H / r_proton),
    "½ × ln(R_H / r_proton)":      0.5 * np.log(R_H / r_proton),
    "⅓ × ln(R_H / r_proton)":      (1/3) * np.log(R_H / r_proton),
    "½ × ln(Rs_m87 / r_proton)":   0.5 * np.log(Rs_m87 / r_proton),
    "ln(r_universe / r_proton)":   np.log(r_universe / r_proton),
    "⅓ × ln(r_universe/r_proton)": (1/3) * np.log(r_universe / r_proton),
    "¼ × ln(r_universe/r_proton)": 0.25 * np.log(r_universe / r_proton),
    "N_MAX / 2":                   N_MAX / 2,
}

print(f"\n   {'Quantity':<38} {'Value':<12} {'Ratio to B':<16} {'% off'}")
print("   " + "-"*75)

b_matches = []
for name, val in candidates_b.items():
    ratio = val / B
    pct   = abs(ratio - 1) * 100
    flag  = "  ✅ MATCH!" if pct < 1.0 else ("  ← close" if pct < 5.0 else "")
    if pct < 1.0:
        b_matches.append((name, val, ratio, pct))
    print(f"   {name:<38} {val:<12.4f} {ratio:<16.4f} {pct:.2f}%{flag}")

# ============================================================================
# Section 3: Is N_MAX/B a clean Planck ratio?
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 3: THE RATIO N_MAX/B = 2.156")
print(f"{'='*70}")

ratio_nb = N_MAX / B
print(f"\n   N_MAX / B = {ratio_nb:.6f}")
print(f"\n   This is the activation curvature x_act = ln(N_MAX/0.5)/B")
print(f"   x_act = {np.log(N_MAX/0.5)/B:.6f}")
print(f"\n   Candidates for N_MAX/B = {ratio_nb:.4f}:")

ratio_candidates = {
    "ln(Rs_m87 / r_proton)":   np.log(Rs_m87 / r_proton),
    "ln(Rs_m87 / r_electron)": np.log(Rs_m87 / r_electron),
    "ln(l_P × c / ħ)":         np.log(l_P * c / hbar) if l_P*c/hbar > 0 else 0,
    "√(ln(R_H/l_P)/N_MAX)":    np.sqrt(np.log(R_H/l_P)/N_MAX),
    "2 × ln(2)":                2 * np.log(2),
    "ln(e²)":                   2.0,
    "√e":                       np.sqrt(np.e),
    "ln(π²)":                   np.log(np.pi**2),
    "π/ln(2)":                  np.pi/np.log(2),
}

for name, val in ratio_candidates.items():
    pct  = abs(val/ratio_nb - 1) * 100
    flag = "  ✅" if pct < 1 else ("  ←" if pct < 5 else "")
    print(f"   {name:<35} {val:.4f}   ({pct:.2f}% off){flag}")

# ============================================================================
# Section 4: The GR crossover x = 0.272 in Planck terms
# ============================================================================

print(f"\n{'='*70}")
print("SECTION 4: GR CROSSOVER x = 0.272 IN PLANCK TERMS")
print(f"{'='*70}")

x_cross = np.log(N_MAX / 0.5) / B
print(f"\n   x_crossover = {x_cross:.6f}")
print(f"\n   Physical meaning: curvature where ESTIF = GR time dilation")
print(f"   In terms of Rs/r: this is at r = Rs / {x_cross:.4f} = {1/x_cross:.4f} × Rs")
print(f"\n   For M87*:  r = {1/x_cross:.2f} × Rs_m87 = {1/x_cross * Rs_m87:.3e} m")
print(f"   For Sun:   r = {1/x_cross:.2f} × Rs_sun  = "
      f"{1/x_cross * 2*G*const.M_sun/c**2:.3e} m")
print(f"\n   Planck-based candidates for x_cross = {x_cross:.4f}:")

x_candidates = {
    "l_P / Rs_m87":             l_P / Rs_m87,
    "√(l_P / Rs_m87)":          np.sqrt(l_P / Rs_m87),
    "(l_P / Rs_m87)^(1/4)":     (l_P / Rs_m87)**(0.25),
    "1 / ln(Rs_m87/l_P)":       1 / np.log(Rs_m87/l_P),
    "1 / ln(R_H/l_P)":          1 / np.log(R_H/l_P),
    "ln(2) / ln(R_H/l_P)":      np.log(2) / np.log(R_H/l_P),
    "ln(2) / N_MAX":             np.log(2) / N_MAX,
    "1 / N_MAX":                 1 / N_MAX,
    "ln(N_MAX) / N_MAX":         np.log(N_MAX) / N_MAX,
    "√(1/N_MAX)":                np.sqrt(1/N_MAX),
}

x_matches = []
for name, val in x_candidates.items():
    pct  = abs(val/x_cross - 1) * 100
    flag = "  ✅ MATCH!" if pct < 1 else ("  ← close" if pct < 5 else "")
    if pct < 1:
        x_matches.append((name, val, pct))
    print(f"   {name:<35} {val:.6f}   ({pct:.2f}% off){flag}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"\n   Targets:  N_MAX={N_MAX:.4f}  B={B:.4f}  x_cross={x_cross:.4f}")

if nmax_matches:
    print(f"\n   ✅ N_MAX matches (within 1%):")
    for name, val, ratio, pct in nmax_matches:
        print(f"   {name} = {val:.4f}  ({pct:.3f}% off)")
else:
    print(f"\n   ⚠️  No exact N_MAX matches within 1%")
    all_n = sorted([(abs(v/N_MAX-1)*100, name, v)
                    for name, v in candidates_nmax.items() if v > 0])
    print(f"   Closest 3:")
    for pct, name, val in all_n[:3]:
        print(f"   {name} = {val:.4f}  ({pct:.2f}% off)")

if b_matches:
    print(f"\n   ✅ B matches (within 1%):")
    for name, val, ratio, pct in b_matches:
        print(f"   {name} = {val:.4f}  ({pct:.3f}% off)")
else:
    print(f"\n   ⚠️  No exact B matches within 1%")
    all_b = sorted([(abs(v/B-1)*100, name, v)
                    for name, v in candidates_b.items() if v > 0])
    print(f"   Closest 3:")
    for pct, name, val in all_b[:3]:
        print(f"   {name} = {val:.4f}  ({pct:.2f}% off)")

if x_matches:
    print(f"\n   ✅ x_crossover matches (within 1%):")
    for name, val, pct in x_matches:
        print(f"   {name} = {val:.6f}  ({pct:.3f}% off)")
else:
    print(f"\n   ⚠️  No exact x_crossover matches within 1%")

print(f"""
   HONEST ASSESSMENT:
   The Planck length provides a rigid, autonomous measuring stick that
   is immune to the dynamic ruler problem. Whether N_MAX and B can be
   expressed cleanly in Planck units determines whether the formula
   is truly fundamental or still contains bent-ruler artifacts.

   If matches are found: the formula has a natural Planck-scale origin
   and the GR crossover emerges from fundamental geometry.

   If no matches: either the formula needs a different reference scale,
   or the 7.5% gap is irreducible and points to missing physics.
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Planck Ruler Test: Can N_MAX and B be expressed in Planck units?',
             fontsize=13, fontweight='bold')

# Plot 1: All N_MAX candidates vs target
ax = axes[0]
names_n  = list(candidates_nmax.keys())
vals_n   = list(candidates_nmax.values())
colors_n = ['green' if abs(v/N_MAX-1)<0.01 else
            'orange' if abs(v/N_MAX-1)<0.05 else
            'steelblue' for v in vals_n]
y_pos = range(len(names_n))
ax.barh(y_pos, vals_n, color=colors_n, alpha=0.8, edgecolor='black')
ax.axvline(N_MAX, color='red', linewidth=2.5, linestyle='--',
           label=f'N_MAX = {N_MAX}')
ax.axvspan(N_MAX*0.99, N_MAX*1.01, alpha=0.2, color='red', label='±1%')
ax.set_yticks(y_pos)
ax.set_yticklabels(names_n, fontsize=7)
ax.set_xlabel('Value', fontsize=11)
ax.set_title('N_MAX Candidates', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(axis='x', alpha=0.3)

# Plot 2: All B candidates vs target
ax = axes[1]
names_b  = list(candidates_b.keys())
vals_b   = list(candidates_b.values())
colors_b = ['green' if abs(v/B-1)<0.01 else
            'orange' if abs(v/B-1)<0.05 else
            'steelblue' for v in vals_b]
y_pos_b = range(len(names_b))
ax.barh(y_pos_b, vals_b, color=colors_b, alpha=0.8, edgecolor='black')
ax.axvline(B, color='red', linewidth=2.5, linestyle='--',
           label=f'B = {B}')
ax.axvspan(B*0.99, B*1.01, alpha=0.2, color='red', label='±1%')
ax.set_yticks(y_pos_b)
ax.set_yticklabels(names_b, fontsize=7)
ax.set_xlabel('Value', fontsize=11)
ax.set_title('B Candidates', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig('planck_ruler_test.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Plot saved: planck_ruler_test.png")
print("=" * 70)
