"""
test_lisa_tilt_scan.py

Tests LISA gravitational wave predictions using the EHT-constrained
tilt exponent range n = 0.05 to 0.215.

QUESTION:
With the same geometric suppression model that passes EHT (n ≈ 0.05–0.215),
does ESTIF still predict a LISA-detectable gravitational wave delay?

If yes → two independent predictions from one geometric parameter → strong science.
If no  → the model is consistent with EHT but silent on GW → weaker but not dead.

SUPPRESSION MODEL:
  sin(θ) = (Rs/r)^n
  β(r)   = cos(θ) = √(1 - (Rs/r)^(2n))

GW DELAY FORMULA (modified from gw_damping_delay):
  τ = (Rs/c) × β(Rs) × (Rs/r)
  At merger r = Rs:
  τ = (Rs/c) × β(Rs)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Detector Thresholds
# ============================================================================
LIGO_PRECISION_S  = 1.0e-3   # 1 ms
LISA_PRECISION_S  = 1.0e-5   # 10 μs  (launch ~2037)

# ============================================================================
# GW150914 Event Parameters
# ============================================================================
M_total_msun = 65.0
M_total      = M_total_msun * const.M_sun
Rs_merger    = 2 * const.G * M_total / const.c**2

# ============================================================================
# EHT-Constrained n Range (from tilt scan results)
# ============================================================================
N_CONSISTENT_LO  = 0.050   # Lower bound (1σ EHT consistent)
N_CONSISTENT_HI  = 0.215   # Upper bound (1σ EHT consistent)
N_BEST_FIT       = 0.050   # Best fit n from EHT scan
N_MODEL_A        = 0.500
N_MODEL_B        = 1.000

# ============================================================================
# Geometric Suppression Functions
# ============================================================================

def beta_tilt(r, Rs, n):
    """β(r) = cos(θ) = √(1 - (Rs/r)^(2n))"""
    ratio = (Rs / r) ** (2 * n)
    if ratio >= 1.0:
        return 0.0
    return np.sqrt(1.0 - ratio)


def gw_delay_geometric(M, n):
    """
    GW merger delay with geometric tilt suppression.

    GWs are generated during the inspiral which ends at the ISCO
    (Innermost Stable Circular Orbit) at r = 3Rs for a Schwarzschild BH.
    Evaluating β at r = Rs gives zero because the hypersurface is
    maximally tilted at the horizon — physically correct but the wrong
    radius for GW production.

    The friction correction the GW accumulates escaping from ISCO:

        τ = (Rs/c) × β(r_ISCO, n)

    where r_ISCO = 3 × Rs.
    """
    Rs     = 2 * const.G * M / const.c**2
    r_isco = 3 * Rs          # ISCO for Schwarzschild black hole
    tau    = (Rs / const.c) * beta_tilt(r_isco, Rs, n)
    return tau


def gw_delay_original(M):
    """
    Original gw_damping_delay formula (BETA_DRAG = 1.0):
        τ = (Rs/c) × BETA_DRAG × (Rs/Rs) = Rs/c × 1.0
    Kept for comparison.
    """
    Rs = 2 * const.G * M / const.c**2
    return (Rs / const.c) * const.BETA_DRAG * (Rs / Rs)


# ============================================================================
# Scan n and compute GW delay
# ============================================================================

n_values = np.linspace(0.01, 2.0, 2000)
delays   = np.array([gw_delay_geometric(M_total, n) for n in n_values])
snr_ligo = delays / LIGO_PRECISION_S
snr_lisa = delays / LISA_PRECISION_S

# Values at key n points
delay_best   = gw_delay_geometric(M_total, N_BEST_FIT)
delay_hi     = gw_delay_geometric(M_total, N_CONSISTENT_HI)
delay_A      = gw_delay_geometric(M_total, N_MODEL_A)
delay_B      = gw_delay_geometric(M_total, N_MODEL_B)
delay_orig   = gw_delay_original(M_total)

# ============================================================================
# Print Results
# ============================================================================

print("=" * 70)
print("LISA GW PREDICTIONS: EHT-CONSTRAINED TILT MODEL")
print("=" * 70)

print(f"\nEvent: GW150914-like binary merger")
print(f"   Total mass:        {M_total_msun:.0f} M_sun")
print(f"   Schwarzschild Rs:  {Rs_merger:.3e} m")
print(f"   Light-crossing:    {Rs_merger/const.c*1000:.4f} ms")

print(f"\nDetector Thresholds:")
print(f"   LIGO precision:    {LIGO_PRECISION_S*1000:.1f} ms")
print(f"   LISA precision:    {LISA_PRECISION_S*1e6:.0f} μs")

print(f"\nEHT-Constrained n Range: {N_CONSISTENT_LO:.3f} – {N_CONSISTENT_HI:.3f}")

print(f"\n{'='*70}")
print("GW DELAY PREDICTIONS")
print(f"{'='*70}")

rows = [
    ("EHT best fit   (n=0.050)", N_BEST_FIT,       delay_best),
    ("EHT upper bound(n=0.215)", N_CONSISTENT_HI,  delay_hi),
    ("Model A        (n=0.500)", N_MODEL_A,         delay_A),
    ("Model B        (n=1.000)", N_MODEL_B,         delay_B),
    ("No suppression (β=1.0)  ", None,              delay_orig),
]

print(f"\n{'Label':<32} {'n':<8} {'β at Rs':<10} {'Delay (s)':<14} {'LIGO S/N':<12} {'LISA S/N':<12} {'Status'}")
print("-" * 100)

for label, n, delay in rows:
    if n is not None:
        b = beta_tilt(Rs_merger, Rs_merger, n)
    else:
        b = const.BETA_DRAG
    snr_l  = delay / LIGO_PRECISION_S
    snr_li = delay / LISA_PRECISION_S

    if snr_li >= 3.0:
        status = "✅ LISA detectable"
    elif snr_li >= 1.0:
        status = "⚠️  LISA marginal"
    else:
        status = "❌ Below threshold"

    print(f"{label:<32} {str(n) if n else 'N/A':<8} {b:<10.4f} "
          f"{delay:<14.3e} {snr_l:<12.3e} {snr_li:<12.3e} {status}")

print(f"\n{'='*70}")
print("KEY FINDING")
print(f"{'='*70}")

snr_lisa_best = delay_best / LISA_PRECISION_S
snr_lisa_hi   = delay_hi   / LISA_PRECISION_S

print(f"\nWithin the EHT-consistent range (n = {N_CONSISTENT_LO} – {N_CONSISTENT_HI}):")
print(f"   LISA S/N at n={N_BEST_FIT}:  {snr_lisa_best:.2f}σ")
print(f"   LISA S/N at n={N_CONSISTENT_HI}: {snr_lisa_hi:.2f}σ")

if snr_lisa_best >= 3.0 and snr_lisa_hi >= 3.0:
    print(f"\n   ✅ STRONG RESULT:")
    print(f"   The entire EHT-consistent range predicts LISA-detectable GW delays.")
    print(f"   One geometric parameter (n) simultaneously satisfies:")
    print(f"   → EHT M87* shadow constraint")
    print(f"   → LISA gravitational wave detectability")
elif snr_lisa_best >= 1.0:
    print(f"\n   ⚠️  PARTIAL RESULT:")
    print(f"   Part of the EHT-consistent range predicts LISA-marginal signals.")
else:
    print(f"\n   ❌ NEGATIVE RESULT:")
    print(f"   EHT-consistent n values predict GW delays below LISA threshold.")

# ============================================================================
# Mass Dependence
# ============================================================================

print(f"\n{'='*70}")
print(f"MASS DEPENDENCE (at best-fit n = {N_BEST_FIT})")
print(f"{'='*70}")
print(f"\n{'Mass (M_sun)':<16} {'Rs (m)':<14} {'Delay (s)':<14} {'LISA S/N':<12} {'Status'}")
print("-" * 65)
for m_msun in [10, 30, 65, 100, 200, 500, 1000]:
    M      = m_msun * const.M_sun
    Rs_m   = 2 * const.G * M / const.c**2
    delay  = gw_delay_geometric(M, N_BEST_FIT)
    snr_li = delay / LISA_PRECISION_S
    status = "✅" if snr_li >= 3 else ("⚠️" if snr_li >= 1 else "❌")
    print(f"{m_msun:<16} {Rs_m:<14.3e} {delay:<14.3e} {snr_li:<12.2f} {status}")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(17, 6))
fig.suptitle('LISA GW Predictions: EHT-Constrained Tilt Model',
             fontsize=14, fontweight='bold')

# Shade EHT-consistent region on all plots
eht_band_kw = dict(alpha=0.2, color='green', label='EHT-consistent range')

# -- Plot 1: GW delay vs n --
axes[0].semilogy(n_values, delays * 1000, 'purple', linewidth=2.5)
axes[0].axhline(LIGO_PRECISION_S * 1000, color='blue',  linewidth=2,
                linestyle='--', label=f'LIGO: {LIGO_PRECISION_S*1000:.0f} ms')
axes[0].axhline(LISA_PRECISION_S * 1000, color='green', linewidth=2,
                linestyle='--', label=f'LISA: {LISA_PRECISION_S*1e6:.0f} μs')
axes[0].axvspan(N_CONSISTENT_LO, N_CONSISTENT_HI, **eht_band_kw)
axes[0].axvline(N_BEST_FIT, color='gold',  linewidth=2, linestyle='-',
                label=f'Best fit n={N_BEST_FIT}')
axes[0].axvline(N_MODEL_A,  color='blue',  linewidth=1.5, linestyle=':',
                alpha=0.7, label='Model A (n=0.5)')
axes[0].axvline(N_MODEL_B,  color='red',   linewidth=1.5, linestyle=':',
                alpha=0.7, label='Model B (n=1.0)')
axes[0].set_xlabel('Tilt Exponent n', fontsize=12)
axes[0].set_ylabel('GW Delay (ms)', fontsize=12)
axes[0].set_title('GW Delay vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)
axes[0].set_xlim(0.01, 2.0)

# -- Plot 2: LISA S/N vs n --
axes[1].semilogy(n_values, snr_lisa, 'purple', linewidth=2.5)
axes[1].axhline(3.0, color='green',  linewidth=2, linestyle='--',
                label='3σ detection threshold')
axes[1].axhline(1.0, color='orange', linewidth=2, linestyle='--',
                label='1σ marginal threshold')
axes[1].axvspan(N_CONSISTENT_LO, N_CONSISTENT_HI, **eht_band_kw)
axes[1].axvline(N_BEST_FIT, color='gold', linewidth=2, linestyle='-',
                label=f'EHT best fit n={N_BEST_FIT}')
axes[1].fill_between(n_values, 3.0, snr_lisa.max()*2,
                     where=(snr_lisa >= 3.0), alpha=0.15, color='green',
                     label='LISA-detectable zone')
axes[1].set_xlabel('Tilt Exponent n', fontsize=12)
axes[1].set_ylabel('LISA Signal-to-Noise Ratio', fontsize=12)
axes[1].set_title('LISA Detectability vs Tilt Exponent', fontsize=12, fontweight='bold')
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.3)
axes[1].set_xlim(0.01, 2.0)

# -- Plot 3: Delay vs merger mass at best-fit n --
masses_msun = np.logspace(0.8, 3.5, 200)
delays_mass = np.array([gw_delay_geometric(m * const.M_sun, N_BEST_FIT)
                        for m in masses_msun])

axes[2].loglog(masses_msun, delays_mass * 1000, 'purple', linewidth=2.5,
               label=f'ESTIF (n={N_BEST_FIT})')
axes[2].axhline(LIGO_PRECISION_S * 1000, color='blue',  linewidth=2,
                linestyle='--', label='LIGO: 1 ms')
axes[2].axhline(LISA_PRECISION_S * 1000, color='green', linewidth=2,
                linestyle='--', label='LISA: 10 μs')
axes[2].fill_between(masses_msun,
                     LISA_PRECISION_S * 1000, LIGO_PRECISION_S * 1000,
                     alpha=0.15, color='green', label='LISA-only range')
axes[2].fill_between(masses_msun,
                     LIGO_PRECISION_S * 1000, delays_mass.max() * 1000,
                     where=(delays_mass * 1000 >= LIGO_PRECISION_S * 1000),
                     alpha=0.15, color='blue', label='LIGO range')
axes[2].axvline(M_total_msun, color='red', linewidth=1.5, linestyle=':',
                label=f'GW150914 ({M_total_msun:.0f} M_sun)')
axes[2].set_xlabel('Binary Total Mass (M_sun)', fontsize=12)
axes[2].set_ylabel('GW Delay (ms)', fontsize=12)
axes[2].set_title(f'Delay vs Mass (n={N_BEST_FIT})', fontsize=12, fontweight='bold')
axes[2].legend(fontsize=9)
axes[2].grid(alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('lisa_tilt_scan_results.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n✓ Plot saved: lisa_tilt_scan_results.png")
print("=" * 70)
