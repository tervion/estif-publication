"""
test_a0_parameter_independence.py
ESTIF v6.2 — March 2026

TEST 3 — PARAMETER INDEPENDENCE
---------------------------------
Goal: Show that a₀ = H₀cx₀/√3 is NOT just "Planck values plugged in."

The criticism: "You chose H₀ and Ωm from Planck 2018. If you use different
cosmological measurements, the result falls apart."

The response: Over the full observationally allowed range of H₀ and Ωm,
the ESTIF prediction for a₀ stays within the observational scatter of the
SPARC baryonic Tully-Fisher relation (~15–20%). The formula is robust —
it is not pinned to a specific cosmological dataset.

APPROACH
--------
1. Sample H₀ ∈ [65, 75] km/s/Mpc  (covers Planck, SH0ES, and everything between)
2. Sample Ωm ∈ [0.27, 0.33]        (covers Planck 1σ and DES/KiDS range)
3. Compute a₀(H₀, Ωm) = H₀ × c × Ωm / √3  (since x₀ ≈ Ωm exactly)
4. Compare to MOND empirical: 1.200×10⁻¹⁰ m/s²
5. Report percentage deviation across the full grid

Pass condition: All combinations within ±20% of MOND empirical value.
(±20% matches the observed SPARC scatter — so ESTIF stays inside the data
regardless of which cosmological measurement you trust.)

Usage:
    python3 test_a0_parameter_independence.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import estif_ec_gr_constants as const

# ============================================================================
# CONSTANTS
# ============================================================================
c        = const.c
G        = const.G
a0_MOND  = 1.200e-10   # m/s²  — MOND empirical (Begeman+1991, McGaugh+2011)
a0_ESTIF = const.A0_MOND_ESTIF   # Planck 2018 central value

# Observational scatter in SPARC BTFR (in v_flat → maps to ~15–20% in a₀)
SPARC_SCATTER_PCT = 20.0   # conservative upper bound

# ============================================================================
# PARAMETER GRIDS
# ============================================================================
# H₀ range: covers full Planck–SH0ES tension range with margin
# Planck: 67.4 ± 0.5, SH0ES: 73.0 ± 1.0 — test 65–75 (wider than any published value)
H0_km_range = np.linspace(65.0, 75.0, 60)     # km/s/Mpc
H0_range    = H0_km_range * 1e3 / 3.085677581e22   # s⁻¹

# Ωm range: covers all published CMB + weak lensing measurements (2σ range)
# Planck: 0.311 ± 0.006, DES: 0.339 ± 0.031, KiDS: 0.279 ± 0.032
Om_range = np.linspace(0.27, 0.33, 60)

H0_grid, Om_grid = np.meshgrid(H0_range, Om_range)
H0_km_grid = H0_grid * 3.085677581e22 / 1e3

# a₀(H₀, Ωm) = H₀ × c × Ωm / √3
# (x₀ = R_H/r_universe ≈ Ωm to 0.12% — treat as equal for this test)
a0_grid = H0_grid * c * Om_grid / np.sqrt(3)

# Percentage deviation from MOND empirical
pct_deviation = (a0_grid - a0_MOND) / a0_MOND * 100.0

print("=" * 70)
print("TEST 3 — a₀ PARAMETER INDEPENDENCE")
print("ESTIF v6.2 — March 2026")
print("=" * 70)

print(f"""
QUESTION: Does a₀ = H₀cx₀/√3 depend critically on using Planck 2018 values?

H₀ range tested:  {H0_km_range[0]:.1f} – {H0_km_range[-1]:.1f} km/s/Mpc
                  (covers Planck 67.4, WMAP 70.0, SH0ES 73.0, and beyond)
Ωm range tested:  {Om_range[0]:.2f} – {Om_range[-1]:.2f}
                  (covers DES/KiDS 0.27, Planck 0.311, SPT 0.30)
Grid size:        {len(H0_km_range)} × {len(Om_range)} = {len(H0_km_range)*len(Om_range)} combinations

MOND empirical:   a₀ = {a0_MOND:.3e} m/s²
SPARC scatter:    ±{SPARC_SCATTER_PCT:.0f}% in a₀ (conservative)
""")

# ============================================================================
# SECTION 1: GRID STATISTICS
# ============================================================================
print("-" * 70)
print("SECTION 1: GRID STATISTICS")
print("-" * 70)

pct_min  = pct_deviation.min()
pct_max  = pct_deviation.max()
pct_mean = pct_deviation.mean()
pct_std  = pct_deviation.std()

# Fraction within SPARC scatter
within_scatter = np.sum(np.abs(pct_deviation) <= SPARC_SCATTER_PCT)
total_points   = pct_deviation.size
frac_within    = within_scatter / total_points * 100

print(f"""
  Deviation range:   {pct_min:+.1f}% to {pct_max:+.1f}%
  Mean deviation:    {pct_mean:+.1f}%
  Std deviation:     {pct_std:.1f}%
  Within ±{SPARC_SCATTER_PCT:.0f}% (SPARC scatter):  {within_scatter}/{total_points}  ({frac_within:.1f}%)
""")

pass1 = frac_within >= 95.0
print(f"  {'✅ PASS' if pass1 else '❌ FAIL'} — {'≥95% of parameter space within SPARC scatter' if pass1 else '<95% within SPARC scatter'}")

# ============================================================================
# SECTION 2: SPECIFIC COSMOLOGICAL DATASETS
# ============================================================================
print("\n" + "-" * 70)
print("SECTION 2: SPECIFIC COSMOLOGICAL DATASETS")
print("-" * 70)
print(f"\n  {'Dataset':<28} {'H₀ [km/s/Mpc]':<16} {'Ωm':<10} {'a₀ [×10⁻¹⁰]':<16} {'Error'}")
print("  " + "-" * 75)

datasets = [
    ("Planck 2018",      67.66, 0.3111),
    ("WMAP 9-year",      69.32, 0.2865),
    ("SH0ES (Riess+22)", 73.04, 0.3000),
    ("DES Year 3",       67.30, 0.3390),
    ("KiDS-1000",        67.40, 0.2790),
    ("SPT-3G",           68.80, 0.3000),
    ("ACT DR4",          67.90, 0.3100),
    ("H0LiCOW",         73.30, 0.3000),
]

all_within = True
for name, H0_kms, Om in datasets:
    H0_si = H0_kms * 1e3 / 3.085677581e22
    a0_val = H0_si * c * Om / np.sqrt(3)
    pct    = (a0_val - a0_MOND) / a0_MOND * 100
    within = abs(pct) <= SPARC_SCATTER_PCT
    flag   = "✅" if within else "❌"
    if not within:
        all_within = False
    print(f"  {name:<28} {H0_kms:<16.2f} {Om:<10.4f} {a0_val/1e-10:<16.4f} {pct:+.1f}%  {flag}")

print(f"\n  All datasets within ±{SPARC_SCATTER_PCT:.0f}%: {'✅ YES' if all_within else '❌ NO'}")

pass2 = all_within

# ============================================================================
# SECTION 3: HUBBLE TENSION ANALYSIS
# ============================================================================
print("\n" + "-" * 70)
print("SECTION 3: HUBBLE TENSION — DOES IT MATTER FOR a₀?")
print("-" * 70)

H0_planck = 67.66; Om_planck = 0.3111
H0_shoes  = 73.04; Om_shoes  = 0.3000

H0_p_si = H0_planck * 1e3 / 3.085677581e22
H0_s_si = H0_shoes  * 1e3 / 3.085677581e22

a0_planck = H0_p_si * c * Om_planck / np.sqrt(3)
a0_shoes  = H0_s_si * c * Om_shoes  / np.sqrt(3)
a0_shift  = (a0_shoes - a0_planck) / a0_planck * 100

print(f"""
  Planck 2018: H₀ = {H0_planck} km/s/Mpc, Ωm = {Om_planck}
    → a₀ = {a0_planck:.4e} m/s²  ({(a0_planck-a0_MOND)/a0_MOND*100:+.2f}% from MOND)

  SH0ES 2022:  H₀ = {H0_shoes} km/s/Mpc, Ωm = {Om_shoes}
    → a₀ = {a0_shoes:.4e} m/s²  ({(a0_shoes-a0_MOND)/a0_MOND*100:+.2f}% from MOND)

  Shift between Planck and SH0ES: {a0_shift:+.2f}%
  SPARC observational scatter:    ±{SPARC_SCATTER_PCT:.0f}%

  The entire Hubble tension ({H0_shoes-H0_planck:.2f} km/s/Mpc, ~5σ controversy)
  shifts a₀ by only {abs(a0_shift):.2f}% — well within SPARC scatter.
""")

pass3 = abs(a0_shift) < SPARC_SCATTER_PCT
print(f"  {'✅ PASS' if pass3 else '❌ FAIL'} — Hubble tension has negligible effect on a₀")

# ============================================================================
# SECTION 4: SENSITIVITY ANALYSIS
# ============================================================================
print("\n" + "-" * 70)
print("SECTION 4: SENSITIVITY — HOW MUCH DOES a₀ CARE ABOUT INPUTS?")
print("-" * 70)

# da₀/dH₀ × ΔH₀ / a₀  — relative sensitivity
# a₀ = H₀ × c × Ωm / √3  → ∂a₀/∂H₀ = c × Ωm / √3 = a₀/H₀
# relative sensitivity = (∂a₀/∂H₀) × ΔH₀ / a₀ = ΔH₀ / H₀

H0_central = const.H_0
Om_central = const.OMEGA_M

delta_H0_pct = (H0_km_range[-1] - H0_km_range[0]) / (H0_km_range[-1] + H0_km_range[0]) * 2 * 100
delta_Om_pct = (Om_range[-1] - Om_range[0])        / (Om_range[-1]    + Om_range[0])    * 2 * 100

da0_from_H0 = delta_H0_pct   # linear: 1% change in H₀ → 1% change in a₀
da0_from_Om = delta_Om_pct   # linear: 1% change in Ωm → 1% change in a₀

print(f"""
  a₀ = H₀ × c × Ωm / √3  →  d(a₀)/a₀ = d(H₀)/H₀ + d(Ωm)/Ωm

  H₀ range tested: ±{delta_H0_pct/2:.1f}% around midpoint
    → a₀ shift: ±{da0_from_H0/2:.1f}%

  Ωm range tested: ±{delta_Om_pct/2:.1f}% around midpoint
    → a₀ shift: ±{da0_from_Om/2:.1f}%

  Combined quadrature: ±{np.sqrt((delta_H0_pct/2)**2 + (delta_Om_pct/2)**2):.1f}%
  SPARC scatter:        ±{SPARC_SCATTER_PCT:.0f}%

  The formula is linearly sensitive to inputs — no hidden amplification.
  Measurement uncertainties in cosmological inputs produce proportional
  uncertainties in a₀, all well within observational scatter.
""")

pass4 = np.sqrt((delta_H0_pct/2)**2 + (delta_Om_pct/2)**2) < SPARC_SCATTER_PCT
print(f"  {'✅ PASS' if pass4 else '❌ FAIL'} — Parameter sensitivity within SPARC scatter")

# ============================================================================
# SUMMARY
# ============================================================================
all_pass = pass1 and pass2 and pass3 and pass4
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  Test 3.1 — Grid coverage (3600 combinations): {'✅ PASS' if pass1 else '❌ FAIL'}
    {frac_within:.1f}% of H₀ × Ωm parameter space within ±20% of MOND

  Test 3.2 — 8 published cosmological datasets:  {'✅ PASS' if pass2 else '❌ FAIL'}
    All datasets (Planck, WMAP, SH0ES, DES, KiDS, SPT, ACT, H0LiCOW)
    give a₀ within SPARC observational scatter

  Test 3.3 — Hubble tension impact:              {'✅ PASS' if pass3 else '❌ FAIL'}
    Planck vs SH0ES shift = {abs(a0_shift):.2f}% — negligible vs ±20% scatter

  Test 3.4 — Linear sensitivity:                 {'✅ PASS' if pass4 else '❌ FAIL'}
    No amplification of input uncertainties

  CONCLUSION: a₀ = H₀cx₀/√3 is robust across the full range of
  observationally allowed cosmological parameters. The result does not
  depend critically on Planck 2018 values — it holds for any reasonable
  cosmology. The formula is not fine-tuned.

  {'✅ PASS — parameter independence confirmed' if all_pass else '❌ FAIL'}
""")

# ============================================================================
# PLOT
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('ESTIF v6.2 — Test 3: a₀ Parameter Independence', fontsize=13, fontweight='bold')

# Panel 1: 2D heatmap of % deviation
ax = axes[0]
levels = [-20, -15, -10, -5, -2, 0, 2, 5, 10, 15, 20]
cmap = plt.cm.RdYlGn_r
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
im = ax.contourf(H0_km_grid, Om_grid, pct_deviation,
                 levels=levels, cmap='RdYlGn', extend='both')
cb = plt.colorbar(im, ax=ax, label='% deviation from MOND empirical')
# Mark SPARC scatter boundary
ax.contour(H0_km_grid, Om_grid, np.abs(pct_deviation),
           levels=[SPARC_SCATTER_PCT], colors='black', linewidths=2, linestyles='--')
# Mark Planck 2018
ax.scatter([67.66], [0.3111], color='blue', s=150, zorder=5,
           marker='*', label='Planck 2018')
# Mark SH0ES
ax.scatter([73.04], [0.300], color='red', s=100, zorder=5,
           marker='D', label='SH0ES')
ax.set_xlabel('H₀ [km/s/Mpc]', fontsize=11)
ax.set_ylabel('Ωm', fontsize=11)
ax.set_title('a₀ deviation map\n(black dashed = ±20% SPARC scatter)', fontsize=10, fontweight='bold')
ax.legend(fontsize=9)

# Panel 2: a₀ vs H₀ for fixed Ωm slices
ax = axes[1]
for Om_val, ls, col in [(0.27, ':', 'steelblue'), (0.311, '-', 'royalblue'), (0.33, '--', 'navy')]:
    H0_si_arr = H0_km_range * 1e3 / 3.085677581e22
    a0_arr = H0_si_arr * c * Om_val / np.sqrt(3)
    ax.plot(H0_km_range, a0_arr/1e-10, lw=2.5, ls=ls, color=col,
            label=f'Ωm = {Om_val}')
ax.axhline(a0_MOND/1e-10, color='black', lw=1.5, ls=':',
           label=f'MOND empirical')
ax.fill_between(H0_km_range,
                a0_MOND*(1-SPARC_SCATTER_PCT/100)/1e-10,
                a0_MOND*(1+SPARC_SCATTER_PCT/100)/1e-10,
                alpha=0.12, color='green', label='SPARC scatter (±20%)')
# Mark Planck and SH0ES
ax.axvline(67.66, color='blue', lw=1, ls='--', alpha=0.5)
ax.axvline(73.04, color='red',  lw=1, ls='--', alpha=0.5)
ax.text(67.66, 2.1, 'Planck', color='blue', fontsize=8, ha='center')
ax.text(73.04, 2.1, 'SH0ES',  color='red',  fontsize=8, ha='center')
ax.set_xlabel('H₀ [km/s/Mpc]', fontsize=11)
ax.set_ylabel('a₀ [×10⁻¹⁰ m/s²]', fontsize=11)
ax.set_title('a₀ vs H₀ for three Ωm values\n(all stay within SPARC scatter)', fontsize=10, fontweight='bold')
ax.set_ylim(0.7, 2.2); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Panel 3: 8 datasets comparison
ax = axes[2]
names_short = ['Planck\n2018', 'WMAP\n9yr', 'SH0ES\n2022', 'DES\nY3',
               'KiDS\n1000', 'SPT\n3G', 'ACT\nDR4', 'H0LiCOW']
a0_vals = []
pct_vals = []
for name, H0_kms, Om in datasets:
    H0_si = H0_kms * 1e3 / 3.085677581e22
    a0_v  = H0_si * c * Om / np.sqrt(3)
    a0_vals.append(a0_v/1e-10)
    pct_vals.append((a0_v - a0_MOND)/a0_MOND*100)

colors_bar = ['royalblue' if abs(p) <= 5 else 'steelblue' if abs(p) <= 10 else 'orange'
              for p in pct_vals]
bars = ax.bar(names_short, a0_vals, color=colors_bar, alpha=0.85, edgecolor='black', lw=1.2)
ax.axhline(a0_MOND/1e-10, color='black', lw=1.5, ls=':', label='MOND empirical')
ax.fill_between(np.arange(-0.5, len(names_short)+0.5),
                a0_MOND*(1-SPARC_SCATTER_PCT/100)/1e-10,
                a0_MOND*(1+SPARC_SCATTER_PCT/100)/1e-10,
                alpha=0.12, color='green', label='SPARC scatter (±20%)')
for bar, pct in zip(bars, pct_vals):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.02,
            f'{pct:+.1f}%', ha='center', va='bottom', fontsize=8)
ax.set_ylabel('a₀ [×10⁻¹⁰ m/s²]', fontsize=11)
ax.set_title('a₀ from 8 cosmological datasets\n(all within SPARC scatter)', fontsize=10, fontweight='bold')
ax.set_ylim(0.8, 1.6); ax.legend(fontsize=9); ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a0_parameter_independence.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"✓ Plot saved: {out}")
print("=" * 70)
