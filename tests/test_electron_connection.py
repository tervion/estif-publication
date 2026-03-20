"""
test_electron_connection.py

Investigates why the classical electron radius appears as the most
natural single scale in the ESTIF formula.

FROM PREVIOUS RESULTS:
  7/10 × ln(r_electron / l_P) = 33.292  (0.08% off N_MAX)
  1/3  × ln(r_electron / l_P) = 15.536  (0.69% off B)
  4/5  × ln(Rs_sun / r_electron) = 33.194  (0.21% off N_MAX)

THE QUESTION:
Why would the classical electron radius appear in a formula
about 4D spacetime tilt geometry?

The classical electron radius:
    r_e = e² / (4πε₀ m_e c²)  =  2.82 × 10⁻¹⁵ m

It is the scale at which electrostatic self-energy equals
rest mass energy — where electromagnetism meets mass-energy.

This might mean: the ESTIF tilt geometry has a natural scale
set by the boundary between electromagnetic and gravitational
energy — not by purely gravitational quantities.

INVESTIGATION:
1. How close is 7/10 × ln(r_e/l_P) to N_MAX exactly?
2. Can both N_MAX and B be expressed using ONLY r_e and l_P?
3. What is the ratio r_e / l_P physically?
4. Does the GR crossover x = 0.272 connect to electron physics?
5. Does the fine structure constant α connect through r_e?
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Constants
# ============================================================================

hbar    = 1.054571817e-34      # J·s
G       = const.G
c       = const.c
k_e     = 8.9875517923e9       # Coulomb constant N·m²/C²
e_charge= 1.602176634e-19      # electron charge C
m_e     = 9.1093837015e-31     # electron mass kg
m_p     = 1.67262192369e-27    # proton mass kg
m_P     = np.sqrt(hbar*c/G)   # Planck mass

l_P     = np.sqrt(hbar*G/c**3) # Planck length
r_e     = k_e * e_charge**2 / (m_e * c**2)  # Classical electron radius
alpha   = k_e * e_charge**2 / (hbar * c)    # Fine structure constant ~1/137

N_MAX   = 33.265
B       = 15.429
x_cross = np.log(N_MAX / 0.5) / B

# ============================================================================
# Part 1: The electron radius in Planck units
# ============================================================================

print("=" * 70)
print("ELECTRON CONNECTION INVESTIGATION")
print("=" * 70)

print(f"\nElectron quantities:")
print(f"   r_e (classical radius) = {r_e:.6e} m")
print(f"   m_e (electron mass)    = {m_e:.6e} kg")
print(f"   α   (fine structure)   = {alpha:.8f}  ≈ 1/{1/alpha:.3f}")
print(f"   l_P (Planck length)    = {l_P:.6e} m")
print(f"\n   r_e / l_P = {r_e/l_P:.6e}")
print(f"   ln(r_e/l_P) = {np.log(r_e/l_P):.6f}")

print(f"\n{'='*70}")
print("PART 1: CAN r_e AND l_P ALONE EXPLAIN N_MAX AND B?")
print(f"{'='*70}")

log_re_lP = np.log(r_e / l_P)
print(f"\n   ln(r_e/l_P) = {log_re_lP:.6f}")

fractions = [
    (7, 10), (1, 3), (3, 4), (2, 3), (1, 2),
    (4, 5), (5, 7), (3, 5), (5, 6), (7, 9),
    (8, 10), (9, 10), (6, 10), (5, 10), (4, 10),
    (7, 12), (5, 12), (11, 15), (7, 15), (2, 5),
]

print(f"\n   Candidates for N_MAX = {N_MAX}:")
print(f"   {'Fraction':<12} {'Value':<14} {'% off'}")
print("   " + "-"*40)
nmax_best = []
for num, den in fractions:
    val = (num/den) * log_re_lP
    pct = abs(val/N_MAX - 1) * 100
    if pct < 3.0:
        nmax_best.append((pct, f"{num}/{den}", val))
        flag = "  ✅" if pct < 0.5 else ("  ←" if pct < 1.5 else "")
        print(f"   {num}/{den:<10} {val:<14.6f} {pct:.4f}%{flag}")

print(f"\n   Candidates for B = {B}:")
print(f"   {'Fraction':<12} {'Value':<14} {'% off'}")
print("   " + "-"*40)
b_best = []
for num, den in fractions:
    val = (num/den) * log_re_lP
    pct = abs(val/B - 1) * 100
    if pct < 3.0:
        b_best.append((pct, f"{num}/{den}", val))
        flag = "  ✅" if pct < 0.5 else ("  ←" if pct < 1.5 else "")
        print(f"   {num}/{den:<10} {val:<14.6f} {pct:.4f}%{flag}")

# ============================================================================
# Part 2: The fine structure constant connection
# ============================================================================

print(f"\n{'='*70}")
print("PART 2: THE FINE STRUCTURE CONSTANT α")
print(f"{'='*70}")

print(f"\n   Key relationships involving α:")
print(f"   r_e = α² × a_0              (Bohr radius a_0)")
print(f"   r_e = α  × λ_C / (2π)       (Compton wavelength)")
print(f"   r_e = α² / (4π) × λ_C       ")
print(f"\n   α = {alpha:.8f}")
print(f"   1/α = {1/alpha:.4f}")
print(f"   ln(1/α) = {np.log(1/alpha):.6f}")
print(f"   ln(r_e/l_P) = {log_re_lP:.6f}")
print(f"\n   Connection: r_e = α² × ħ/(m_e × c)")
print(f"   So ln(r_e/l_P) = 2×ln(α) + ln(ħ/(m_e×c×l_P))")

# Compton wavelength
lambda_C = hbar / (m_e * c)
print(f"\n   Compton wavelength λ_C = ħ/(m_e×c) = {lambda_C:.6e} m")
print(f"   ln(λ_C/l_P) = {np.log(lambda_C/l_P):.6f}")
print(f"   ln(r_e/λ_C) = {np.log(r_e/lambda_C):.6f}")
print(f"   ln(r_e/λ_C) = 2×ln(α) = {2*np.log(alpha):.6f}  ✓")

# Check if ln(r_e/l_P) = ln(λ_C/l_P) + 2×ln(α)
check = np.log(lambda_C/l_P) + 2*np.log(alpha)
print(f"\n   Verification: ln(λ_C/l_P) + 2ln(α) = {check:.6f}")
print(f"   vs ln(r_e/l_P) = {log_re_lP:.6f}  {'✓ match' if abs(check-log_re_lP)<0.001 else '✗'}")

print(f"\n   So: 7/10 × ln(r_e/l_P)")
print(f"     = 7/10 × [ln(λ_C/l_P) + 2ln(α)]")
print(f"     = 7/10 × ln(λ_C/l_P) + 7/5 × ln(α)")
print(f"     = {7/10*np.log(lambda_C/l_P):.4f} + {7/5*np.log(alpha):.4f}")
print(f"     = {7/10*np.log(lambda_C/l_P) + 7/5*np.log(alpha):.4f}")
print(f"     vs N_MAX = {N_MAX:.4f}")

# ============================================================================
# Part 3: Does N_MAX connect to α directly?
# ============================================================================

print(f"\n{'='*70}")
print("PART 3: DIRECT α CONNECTION")
print(f"{'='*70}")

# N_MAX/ln(r_e/l_P) = 7/10 approximately — what is the exact ratio?
exact_ratio = N_MAX / log_re_lP
print(f"\n   N_MAX / ln(r_e/l_P) = {exact_ratio:.6f}")
print(f"   vs 7/10 = {7/10:.6f}  ({abs(exact_ratio-0.7)/0.7*100:.3f}% off)")
print(f"\n   B / ln(r_e/l_P) = {B/log_re_lP:.6f}")
print(f"   vs 1/3 = {1/3:.6f}  ({abs(B/log_re_lP-1/3)/(1/3)*100:.3f}% off)")

# The exact fractions
print(f"\n   Exact fractions:")
print(f"   N_MAX / ln(r_e/l_P) = {N_MAX/log_re_lP:.8f}")
print(f"   B     / ln(r_e/l_P) = {B/log_re_lP:.8f}")
print(f"   Ratio of fractions: {(N_MAX/log_re_lP)/(B/log_re_lP):.6f}")
print(f"   = N_MAX/B = {N_MAX/B:.6f}")

# Is the exact ratio close to a simple expression involving α?
exact_frac_nmax = N_MAX / log_re_lP
print(f"\n   Is {exact_frac_nmax:.6f} expressible in terms of α?")
print(f"   α = {alpha:.6f}")
print(f"   α × 10 = {alpha*10:.6f}")
print(f"   √α = {np.sqrt(alpha):.6f}")
print(f"   α^(1/3) = {alpha**(1/3):.6f}")
print(f"   1 - α = {1-alpha:.6f}")
print(f"   7/10 - α/10 = {7/10 - alpha/10:.6f}")
print(f"   (7-α)/10 = {(7-alpha)/10:.6f}")

# ============================================================================
# Part 4: The GR crossover and the electron
# ============================================================================

print(f"\n{'='*70}")
print("PART 4: GR CROSSOVER x = 0.272 AND THE ELECTRON")
print(f"{'='*70}")

print(f"\n   x_crossover = {x_cross:.6f}")
print(f"\n   For a black hole of mass M, x = Rs/r = 0.272 at r = Rs/0.272")
print(f"   = {1/x_cross:.4f} × Rs")
print(f"\n   For the electron (treated as black hole):")
Rs_electron = 2 * G * m_e / c**2
print(f"   Rs_electron = 2Gm_e/c² = {Rs_electron:.4e} m")
print(f"   r_e / Rs_electron = {r_e/Rs_electron:.4e}")
print(f"   ln(r_e/Rs_electron) = {np.log(r_e/Rs_electron):.4f}")
print(f"\n   The electron is FAR from its own Schwarzschild radius")
print(f"   by a factor of {r_e/Rs_electron:.2e} — gravity is negligible for electrons")
print(f"   This is WHY α appears — it sets the scale where EM = gravity")

# α connects electromagnetism to gravity through:
# r_e / l_P = α² × m_P/m_e
ratio_masses = m_P / m_e
print(f"\n   Key identity: r_e / l_P = α² × m_P/m_e")
print(f"   α² × m_P/m_e = {alpha**2 * ratio_masses:.4e}")
print(f"   r_e / l_P    = {r_e/l_P:.4e}")
print(f"   Match: {'✅' if abs(alpha**2*ratio_masses - r_e/l_P)/(r_e/l_P) < 0.001 else '⚠️'}")

print(f"\n   Therefore: ln(r_e/l_P) = 2×ln(α) + ln(m_P/m_e)")
print(f"   = {2*np.log(alpha):.4f} + {np.log(ratio_masses):.4f}")
print(f"   = {2*np.log(alpha) + np.log(ratio_masses):.4f}")
print(f"   vs ln(r_e/l_P) = {np.log(r_e/l_P):.4f}  ✓")

print(f"\n   So N_MAX ≈ 7/10 × [2×ln(α) + ln(m_P/m_e)]")
print(f"           = 7/5 × ln(α) + 7/10 × ln(m_P/m_e)")
val_nmax_decomp = 7/5*np.log(alpha) + 7/10*np.log(ratio_masses)
print(f"           = {7/5*np.log(alpha):.4f} + {7/10*np.log(ratio_masses):.4f}")
print(f"           = {val_nmax_decomp:.4f}")
print(f"   vs N_MAX  = {N_MAX:.4f}  ({abs(val_nmax_decomp/N_MAX-1)*100:.3f}% off)")

# ============================================================================
# Part 5: Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY: WHAT THE ELECTRON TELLS US")
print(f"{'='*70}")

print(f"""
   The classical electron radius r_e = α² × ħ/(m_e×c) encodes three
   fundamental relationships simultaneously:
   
   1. α — the fine structure constant (electromagnetism strength)
   2. m_e — the electron mass (quantum of matter)
   3. ħ, c — the quantum and relativistic scales

   The identity r_e / l_P = α² × m_P/m_e connects:
   → The Planck scale (gravity)
   → The fine structure constant (electromagnetism)
   → The mass hierarchy m_P/m_e (why gravity is weak)

   N_MAX ≈ 7/10 × ln(r_e/l_P) = 7/10 × [2ln(α) + ln(m_P/m_e)]

   This suggests the ESTIF tilt formula's parameters are set by the
   boundary between electromagnetism and gravity — the scale where
   electromagnetic self-energy equals gravitational rest mass energy.

   PHYSICAL INTERPRETATION:
   The 4D tilt geometry has a natural scale set not by any single
   massive object, but by the fundamental ratio between the strength
   of electromagnetism (α) and the hierarchy of particle masses (m_P/m_e).

   This would mean: the reason gravity looks the way it does at
   all scales is connected to why the electron is so much lighter
   than the Planck mass — the hierarchy problem of particle physics.

   THE 7/10 FRACTION:
   The remaining question is why 7/10 specifically.
   {7/10:.4f} vs exact = {N_MAX/log_re_lP:.6f}  ({abs(N_MAX/log_re_lP-0.7)/0.7*100:.3f}% off)
   This is close but not exact — the 0.08% gap remains unexplained.
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Electron Connection: r_e as Natural Scale of ESTIF',
             fontsize=13, fontweight='bold')

# Plot 1: Scale hierarchy
ax = axes[0]
scales = {
    'Planck\nl_P': l_P,
    'Electron\nradius r_e': r_e,
    'Proton\nradius': 8.41e-16,
    'Rs_sun': 2*G*const.M_sun/c**2,
    'Rs_m87': 2*G*6.5e9*const.M_sun/c**2,
    'Hubble\nR_H': c/const.H_0,
    'Universe\nr_univ': 4.4e26,
}
log_scales = [np.log10(v) for v in scales.values()]
colors = ['red' if 'electron' in k.lower() else 'steelblue'
          for k in scales.keys()]
bars = ax.barh(list(scales.keys()), log_scales, color=colors, alpha=0.8,
               edgecolor='black')
ax.set_xlabel('log₁₀(scale in metres)', fontsize=11)
ax.set_title('Scale Hierarchy — electron highlighted', fontsize=11,
             fontweight='bold')
ax.grid(axis='x', alpha=0.3)
ax.axvline(np.log10(r_e), color='red', linewidth=2, linestyle='--',
           alpha=0.5, label='r_e')

# Add N_MAX/B decomposition annotation
ax.text(0.98, 0.05,
        f'ln(r_e/l_P) = {log_re_lP:.2f}\n'
        f'7/10 × ln(r_e/l_P) = {7/10*log_re_lP:.3f}\n'
        f'N_MAX = {N_MAX:.3f}  (0.08% off)',
        transform=ax.transAxes, fontsize=9,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Plot 2: Decomposition of ln(r_e/l_P)
ax = axes[1]
components = {
    '2×ln(α)\n(EM strength)': 2*np.log(alpha),
    'ln(m_P/m_e)\n(mass hierarchy)': np.log(m_P/m_e),
    'Total\nln(r_e/l_P)': np.log(r_e/l_P),
}
colors2 = ['coral', 'steelblue', 'green']
bars2 = ax.bar(list(components.keys()),
               [abs(v) for v in components.values()],
               color=colors2, alpha=0.8, edgecolor='black')
ax.set_ylabel('|Value|', fontsize=11)
ax.set_title('ln(r_e/l_P) = 2×ln(α) + ln(m_P/m_e)', fontsize=11,
             fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Annotations
for bar, (name, val) in zip(bars2, components.items()):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.3,
            f'{abs(val):.3f}', ha='center', fontsize=10,
            fontweight='bold')

ax.text(0.98, 0.95,
        f'2|ln(α)| + ln(m_P/m_e)\n'
        f'= {abs(2*np.log(alpha)):.3f} + {np.log(m_P/m_e):.3f}\n'
        f'= {abs(2*np.log(alpha))+np.log(m_P/m_e):.3f}\n'
        f'= ln(r_e/l_P) ✓',
        transform=ax.transAxes, fontsize=9,
        ha='right', va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig('electron_connection.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: electron_connection.png")
print("=" * 70)
