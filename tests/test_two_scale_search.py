"""
test_two_scale_search.py

Searches combinations of TWO scales anchored to the Planck length
to find expressions for N_MAX = 33.265 and B = 15.429.

MOTIVATION FROM PREVIOUS RESULTS:
  N_MAX closest: ¼ × ln(R_H / l_P)     = 35.07  (5.4% off)
  B closest:     ⅓ × ln(r_proton / l_P) = 15.13  (1.9% off)

This suggests the formula might involve THREE scales simultaneously:
  Planck length (bottom)
  Proton (intermediate)
  Hubble radius (top)

A ratio or combination involving all three might close the gap.

SEARCH STRATEGY:
For each pair of scales (A, B), test:
  ln(A/l_P) × f + ln(B/l_P) × g  for small rational f, g
  ln(A/B) × f
  ln(A/l_P) / ln(B/l_P)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import itertools
import estif_ec_gr_constants as const

# ============================================================================
# Constants
# ============================================================================

hbar = 1.054571817e-34
G    = const.G
c    = const.c
l_P  = np.sqrt(hbar * G / c**3)
m_P  = np.sqrt(hbar * c / G)

N_MAX   = 33.265
B       = 15.429
x_cross = np.log(N_MAX / 0.5) / B   # 0.2721

# Physical scales
scales = {
    "r_universe":  4.4e26,
    "R_H":         c / const.H_0,
    "Rs_m87":      2*G*(6.5e9*const.M_sun)/c**2,
    "r_photon":    1.5 * 2*G*(6.5e9*const.M_sun)/c**2,
    "r_proton":    8.41e-16,
    "r_electron":  2.82e-15,
    "Rs_sun":      2*G*const.M_sun/c**2,
    "AU":          const.AU,
}

# Log ratios to Planck length
log_ratios = {name: np.log(val / l_P) for name, val in scales.items()}
log_ratios_pairs = {
    f"ln({a}/{b})": np.log(scales[a] / scales[b])
    for a in scales for b in scales if a != b and scales[a] > scales[b]
}

# Fractions to try
fractions = [1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8,
             2/3, 3/4, 2/5, 3/5, 4/5, 1/9, 1/10,
             2/7, 3/7, 4/7, 5/7]

# ============================================================================
# Search 1: Single scale × fraction
# ============================================================================

print("=" * 70)
print("TWO-SCALE COMBINATION SEARCH")
print("=" * 70)
print(f"\n   Targets: N_MAX={N_MAX}  B={B}  x_cross={x_cross:.4f}")

print(f"\n{'='*70}")
print("SEARCH 1: Single ln(scale/l_P) × fraction")
print(f"{'='*70}")

single_results_nmax = []
single_results_b    = []

for name, logr in log_ratios.items():
    for f in fractions:
        val  = logr * f
        if val <= 0:
            continue
        pct_n = abs(val/N_MAX - 1) * 100
        pct_b = abs(val/B    - 1) * 100
        if pct_n < 3.0:
            single_results_nmax.append((pct_n, f"{f:.4g}×ln({name}/l_P)", val))
        if pct_b < 2.0:
            single_results_b.append((pct_b, f"{f:.4g}×ln({name}/l_P)", val))

single_results_nmax.sort()
single_results_b.sort()

print(f"\n   Best single-scale candidates for N_MAX:")
for pct, name, val in single_results_nmax[:5]:
    print(f"   {name:<40} = {val:.4f}  ({pct:.3f}% off)")

print(f"\n   Best single-scale candidates for B:")
for pct, name, val in single_results_b[:5]:
    print(f"   {name:<40} = {val:.4f}  ({pct:.3f}% off)")

# ============================================================================
# Search 2: f×ln(A/l_P) + g×ln(B/l_P)
# ============================================================================

print(f"\n{'='*70}")
print("SEARCH 2: f×ln(A/l_P) + g×ln(B/l_P)")
print(f"{'='*70}")

combo_nmax = []
combo_b    = []

scale_names = list(log_ratios.keys())
for i, name_a in enumerate(scale_names):
    for j, name_b in enumerate(scale_names):
        if i >= j:
            continue
        la = log_ratios[name_a]
        lb = log_ratios[name_b]
        for f in fractions:
            for g in fractions:
                val = f*la + g*lb
                if val <= 0:
                    continue
                pct_n = abs(val/N_MAX - 1) * 100
                pct_b = abs(val/B    - 1) * 100
                label = f"{f:.3g}×ln({name_a}/l_P)+{g:.3g}×ln({name_b}/l_P)"
                if pct_n < 1.0:
                    combo_nmax.append((pct_n, label, val))
                if pct_b < 1.0:
                    combo_b.append((pct_b, label, val))

combo_nmax.sort()
combo_b.sort()

if combo_nmax:
    print(f"\n   ✅ MATCHES for N_MAX (within 1%):")
    for pct, name, val in combo_nmax[:5]:
        print(f"   {name}")
        print(f"   = {val:.6f}  ({pct:.4f}% off)")
else:
    print(f"\n   ⚠️  No combination matches for N_MAX within 1%")
    # Show top 3
    results = []
    for i, name_a in enumerate(scale_names):
        for j, name_b in enumerate(scale_names):
            if i >= j: continue
            la, lb = log_ratios[name_a], log_ratios[name_b]
            for f in fractions:
                for g in fractions:
                    val = f*la + g*lb
                    if val > 0:
                        results.append((abs(val/N_MAX-1)*100,
                                       f"{f:.3g}×{name_a}+{g:.3g}×{name_b}", val))
    results.sort()
    print(f"   Closest 3:")
    for pct, name, val in results[:3]:
        print(f"   {name} = {val:.4f}  ({pct:.3f}% off)")

if combo_b:
    print(f"\n   ✅ MATCHES for B (within 1%):")
    for pct, name, val in combo_b[:5]:
        print(f"   {name}")
        print(f"   = {val:.6f}  ({pct:.4f}% off)")
else:
    print(f"\n   ⚠️  No combination matches for B within 1%")
    results_b = []
    for i, name_a in enumerate(scale_names):
        for j, name_b in enumerate(scale_names):
            if i >= j: continue
            la, lb = log_ratios[name_a], log_ratios[name_b]
            for f in fractions:
                for g in fractions:
                    val = f*la + g*lb
                    if val > 0:
                        results_b.append((abs(val/B-1)*100,
                                         f"{f:.3g}×{name_a}+{g:.3g}×{name_b}", val))
    results_b.sort()
    print(f"   Closest 3:")
    for pct, name, val in results_b[:3]:
        print(f"   {name} = {val:.4f}  ({pct:.3f}% off)")

# ============================================================================
# Search 3: ln(A/B) without l_P — pure scale ratios
# ============================================================================

print(f"\n{'='*70}")
print("SEARCH 3: Pure scale ratios ln(A/B) × fraction")
print(f"{'='*70}")

ratio_nmax = []
ratio_b    = []

for pair_name, logr in log_ratios_pairs.items():
    for f in fractions:
        val   = logr * f
        if val <= 0:
            continue
        pct_n = abs(val/N_MAX - 1) * 100
        pct_b = abs(val/B    - 1) * 100
        label = f"{f:.4g}×{pair_name}"
        if pct_n < 1.5:
            ratio_nmax.append((pct_n, label, val))
        if pct_b < 1.5:
            ratio_b.append((pct_b, label, val))

ratio_nmax.sort()
ratio_b.sort()

print(f"\n   Best pure-ratio candidates for N_MAX (within 1.5%):")
for pct, name, val in ratio_nmax[:5]:
    print(f"   {name:<45} = {val:.4f}  ({pct:.3f}% off)")

print(f"\n   Best pure-ratio candidates for B (within 1.5%):")
for pct, name, val in ratio_b[:5]:
    print(f"   {name:<45} = {val:.4f}  ({pct:.3f}% off)")

# ============================================================================
# Search 4: N_MAX and B from SAME expression
# ============================================================================

print(f"\n{'='*70}")
print("SEARCH 4: N_MAX AND B FROM THE SAME BASE EXPRESSION")
print(f"{'='*70}")

print(f"\n   If N_MAX = f × X and B = g × X for same X,")
print(f"   then N_MAX/B = f/g = {N_MAX/B:.4f}")
print(f"\n   Looking for rational approximation to {N_MAX/B:.4f}:")
# Find rational approx
for num in range(1, 15):
    for den in range(1, 15):
        if abs(num/den - N_MAX/B) < 0.01:
            print(f"   {num}/{den} = {num/den:.4f}  ({abs(num/den-N_MAX/B)*100:.3f}% off)")

print(f"\n   If N_MAX/B ≈ 15/7 = {15/7:.4f}  ({abs(15/7-N_MAX/B)*100:.3f}% off)")
print(f"   Then: X = N_MAX/(15/7) = {N_MAX/(15/7):.4f}")
print(f"         X = B/(7/7) = {B:.4f}")
print(f"\n   Checking X = {N_MAX*7/15:.4f} against all scales:")

X_target = N_MAX * 7 / 15
for name, logr in log_ratios.items():
    for f in fractions:
        val = logr * f
        if val > 0 and abs(val/X_target - 1) < 0.02:
            print(f"   ✅ {f:.4g}×ln({name}/l_P) = {val:.4f}  "
                  f"({abs(val/X_target-1)*100:.3f}% off from X)")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

all_matches = combo_nmax + combo_b + \
              [r for r in ratio_nmax if r[0] < 0.5] + \
              [r for r in ratio_b    if r[0] < 0.5]

if all_matches:
    print(f"\n   ✅ Matches found within tight tolerance:")
    for pct, name, val in sorted(all_matches)[:10]:
        print(f"   {pct:.4f}%  {name} = {val:.4f}")
else:
    print(f"\n   ⚠️  No matches within 1% found in two-scale search.")
    print(f"\n   INTERPRETATION:")
    print(f"   The 7.5% gap in N_MAX and 1.9% gap in B suggest either:")
    print(f"   (a) The formula parameters are not expressible as simple")
    print(f"       logarithmic ratios of known scales — they may require")
    print(f"       a deeper geometric derivation from first principles.")
    print(f"   (b) The measurement uncertainties in r_universe, Rs_m87,")
    print(f"       and r_proton accumulate to produce the observed gaps.")
    print(f"   (c) A different class of rigid ruler is needed — perhaps")
    print(f"       based on dimensionless numbers like the fine structure")
    print(f"       constant α ≈ 1/137 or the proton/electron mass ratio.")

# ============================================================================
# Search 5: Dimensionless constants
# ============================================================================

print(f"\n{'='*70}")
print("SEARCH 5: DIMENSIONLESS CONSTANTS")
print(f"{'='*70}")

alpha       = 1/137.036          # Fine structure constant
mp_me       = 1836.15            # Proton/electron mass ratio
mp_mP       = const.M_sun*0/1 + 1.6726e-27 / m_P  # Proton/Planck mass ratio
proton_mass = 1.6726e-27         # kg
mp_mP_ratio = proton_mass / m_P

dimless = {
    "1/α (fine structure)":         1/alpha,
    "ln(1/α)":                      np.log(1/alpha),
    "mp/me (mass ratio)":           mp_me,
    "ln(mp/me)":                    np.log(mp_me),
    "ln(mp/mP)":                    np.log(proton_mass/m_P),
    "½×ln(mp/mP)":                  0.5*np.log(proton_mass/m_P),
    "⅓×ln(mp/mP)":                  (1/3)*np.log(proton_mass/m_P),
    "¼×ln(mp/mP)":                  0.25*np.log(proton_mass/m_P),
    "ln(mP/mp)":                    np.log(m_P/proton_mass),
    "½×ln(mP/mp)":                  0.5*np.log(m_P/proton_mass),
    "⅓×ln(mP/mp)":                  (1/3)*np.log(m_P/proton_mass),
    "¼×ln(mP/mp)":                  0.25*np.log(m_P/proton_mass),
    "ln(mP/M_sun)":                 np.log(m_P/const.M_sun),
    "ln(M_sun/mp)":                 np.log(const.M_sun/proton_mass),
    "½×ln(M_sun/mp)":               0.5*np.log(const.M_sun/proton_mass),
}

print(f"\n   {'Quantity':<35} {'Value':<12} {'vs N_MAX':<12} {'vs B'}")
print("   " + "-"*65)

dim_matches = []
for name, val in dimless.items():
    if val <= 0:
        continue
    pct_n = abs(val/N_MAX - 1) * 100
    pct_b = abs(val/B    - 1) * 100
    flag  = ""
    if pct_n < 2.0:
        flag = f"  ✅ N_MAX ({pct_n:.2f}%)"
        dim_matches.append((pct_n, f"N_MAX≈{name}", val))
    if pct_b < 2.0:
        flag = f"  ✅ B ({pct_b:.2f}%)"
        dim_matches.append((pct_b, f"B≈{name}", val))
    print(f"   {name:<35} {val:<12.4f} {pct_n:<12.2f}% {pct_b:.2f}%{flag}")

if dim_matches:
    print(f"\n   ✅ Dimensionless matches:")
    for pct, name, val in sorted(dim_matches):
        print(f"   {name} = {val:.4f}  ({pct:.3f}% off)")

print(f"\n✓ Search complete")
print("=" * 70)
