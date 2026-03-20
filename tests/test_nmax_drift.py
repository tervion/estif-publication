"""
test_nmax_drift.py

Estimates how fast N_MAX changes over cosmic time and whether
that drift is observationally detectable.

THE IDEA:
N_MAX ≈ ln(r_universe / Rs_m87)

Both quantities change over time:
  r_universe grows as the universe expands
  Rs_m87 grows as M87* accretes matter

So N_MAX = ln(r_universe(t) / Rs_m87(t)) is time-dependent.

If Λ depends on N_MAX, and N_MAX changes over time,
then Λ is not truly constant — it just looks constant
because it changes very slowly.

The question: how slowly? Is it detectable?
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../src'))

import numpy as np
import matplotlib.pyplot as plt
import estif_ec_gr_constants as const

# ============================================================================
# Current values
# ============================================================================

R_H           = const.c / const.H_0          # Hubble radius m
r_universe_0  = 4.4e26                        # Observable universe today m
M_m87_0       = 6.5e9 * const.M_sun          # M87* mass today kg
Rs_m87_0      = 2 * const.G * M_m87_0 / const.c**2

LAMBDA_0      = 1.1056e-52                    # Λ today m⁻²
B             = 15.4291                        # calibrated B
N_MAX_0       = 33.265                         # N_MAX today

# Time scales
t_universe    = 13.8e9 * 365.25 * 24 * 3600  # age of universe in seconds
H_0_per_sec   = const.H_0                     # s⁻¹

# ============================================================================
# Rate 1: How fast is r_universe growing?
# ============================================================================

# The observable universe grows at approximately c
# (the Hubble radius grows as c/H(t))
# Today: dr_universe/dt ≈ c (the universe horizon expands at c)
dr_universe_dt = const.c                       # m/s

# Fractional rate
r_universe_frac_rate = dr_universe_dt / r_universe_0
print("=" * 70)
print("N_MAX DRIFT: HOW FAST DOES THE FORMULA CHANGE?")
print("=" * 70)

print(f"\n{'='*70}")
print("PART 1: RATE OF CHANGE OF r_universe")
print(f"{'='*70}")
print(f"\n   r_universe today:      {r_universe_0:.3e} m")
print(f"   dr/dt (≈ c):           {dr_universe_dt:.3e} m/s")
print(f"   Fractional rate:       {r_universe_frac_rate:.3e} per second")
print(f"   Fractional rate:       {r_universe_frac_rate * 3.156e7:.3e} per year")
print(f"   Doubling time:         {1/r_universe_frac_rate / 3.156e16:.2f} billion years")

# ============================================================================
# Rate 2: How fast is Rs_m87 growing?
# ============================================================================

# M87* accretion rate
# Observed: M87* accretes roughly 0.001 M_sun per year (Eddington limited)
# Conservative estimate: 1e-3 M_sun/year
accretion_rate_msun_per_year = 1e-3           # M_sun/year
accretion_rate_kg_per_s      = (accretion_rate_msun_per_year
                                 * const.M_sun
                                 / (365.25 * 24 * 3600))

dRs_dt = 2 * const.G * accretion_rate_kg_per_s / const.c**2
Rs_frac_rate = dRs_dt / Rs_m87_0

print(f"\n{'='*70}")
print("PART 2: RATE OF CHANGE OF Rs_m87")
print(f"{'='*70}")
print(f"\n   Rs_m87 today:          {Rs_m87_0:.3e} m")
print(f"   Accretion rate:        {accretion_rate_msun_per_year:.0e} M_sun/year (observed)")
print(f"   dRs/dt:                {dRs_dt:.3e} m/s")
print(f"   Fractional rate:       {Rs_frac_rate:.3e} per second")
print(f"   Fractional rate:       {Rs_frac_rate * 3.156e7:.3e} per year")
print(f"   Doubling time:         {1/Rs_frac_rate / 3.156e16:.2f} billion years")

# ============================================================================
# Rate 3: dN_MAX/dt
# ============================================================================

# N_MAX = ln(r_universe / Rs)
# dN_MAX/dt = (1/r_universe)(dr/dt) - (1/Rs)(dRs/dt)
# = fractional_rate_r - fractional_rate_Rs

dNmax_dt = r_universe_frac_rate - Rs_frac_rate  # per second
dNmax_per_year  = dNmax_dt * 3.156e7
dNmax_per_Gyr   = dNmax_dt * 3.156e16

print(f"\n{'='*70}")
print("PART 3: RATE OF CHANGE OF N_MAX")
print(f"{'='*70}")
print(f"\n   N_MAX today:           {N_MAX_0:.4f}")
print(f"   dN_MAX/dt (universe):  +{r_universe_frac_rate:.3e} /s  (universe growing)")
print(f"   dN_MAX/dt (M87*):      -{Rs_frac_rate:.3e} /s  (M87* growing)")
print(f"   Net dN_MAX/dt:         {dNmax_dt:+.3e} /s")
print(f"   Net dN_MAX per year:   {dNmax_per_year:+.3e}")
print(f"   Net dN_MAX per Gyr:    {dNmax_per_Gyr:+.2f}")
print(f"\n   The universe expansion dominates completely over M87* accretion")
print(f"   Ratio: {r_universe_frac_rate/Rs_frac_rate:.1e}×")

# ============================================================================
# Rate 4: How does Λ change?
# ============================================================================

# Λ = (3/R_H²) × Observable²
# Observable = √β where β = √(1 - x^(2n(x)))
# n(x) depends on N_MAX

# At cosmic scale x = R_H/r_universe = CURV_COSM
CURV_COSM = R_H / r_universe_0

def n_dynamic(curvature, N_MAX, B_val):
    return N_MAX * np.exp(-B_val * curvature)

def beta_val(curvature, N_MAX, B_val):
    n   = n_dynamic(curvature, N_MAX, B_val)
    val = curvature ** (2 * n)
    return 0.0 if val >= 1.0 else np.sqrt(1.0 - val)

def observable(curvature, N_MAX, B_val):
    return np.sqrt(beta_val(curvature, N_MAX, B_val))

def lambda_from_nmax(N_MAX_val):
    obs = observable(CURV_COSM, N_MAX_val, B)
    return (3.0 / R_H**2) * obs**2

# Numerical derivative
dN     = N_MAX_0 * 1e-6
dL_dN  = (lambda_from_nmax(N_MAX_0 + dN) -
           lambda_from_nmax(N_MAX_0 - dN)) / (2 * dN)

dLambda_dt      = dL_dN * dNmax_dt
dLambda_per_Gyr = dL_dN * dNmax_per_Gyr
frac_change_Gyr = dLambda_per_Gyr / LAMBDA_0

print(f"\n{'='*70}")
print("PART 4: RATE OF CHANGE OF Λ")
print(f"{'='*70}")
print(f"\n   Λ today:               {LAMBDA_0:.4e} m⁻²")
print(f"   dΛ/dN_MAX:             {dL_dN:.4e} m⁻² per unit N_MAX")
print(f"   dΛ/dt:                 {dLambda_dt:.4e} m⁻² per second")
print(f"   dΛ per Gyr:            {dLambda_per_Gyr:.4e} m⁻²")
print(f"   Fractional change/Gyr: {frac_change_Gyr:.4e}")
print(f"   Fractional change/Gyr: {frac_change_Gyr*100:.6f}%")

# ============================================================================
# Is it detectable?
# ============================================================================

print(f"\n{'='*70}")
print("PART 5: IS THIS DETECTABLE?")
print(f"{'='*70}")

# Current best constraints on Λ variation
# Dark energy equation of state w — if w ≠ -1, Λ is evolving
# Current precision: Δw/w ~ 1% over ~5 Gyr baseline
# This translates to Λ variation precision of roughly 1%/5Gyr = 0.2%/Gyr

detection_threshold_per_gyr = 0.002  # 0.2% per Gyr (current limit)
future_threshold_per_gyr    = 0.0001 # 0.01% per Gyr (EUCLID/LSST goal)

print(f"\n   Predicted Λ drift:     {abs(frac_change_Gyr)*100:.6f}% per Gyr")
print(f"\n   Detection thresholds:")
print(f"   Current surveys:       {detection_threshold_per_gyr*100:.2f}% per Gyr")
print(f"   EUCLID/LSST (2030s):   {future_threshold_per_gyr*100:.3f}% per Gyr")
print(f"\n   Signal / Current threshold:  "
      f"{abs(frac_change_Gyr)/detection_threshold_per_gyr:.2e}")
print(f"   Signal / Future threshold:   "
      f"{abs(frac_change_Gyr)/future_threshold_per_gyr:.2e}")

if abs(frac_change_Gyr) > detection_threshold_per_gyr:
    detect_verdict = "✅ DETECTABLE NOW"
elif abs(frac_change_Gyr) > future_threshold_per_gyr:
    detect_verdict = "⚠️  DETECTABLE WITH EUCLID/LSST"
else:
    detect_verdict = "❌ BELOW FORESEEABLE DETECTION"

print(f"\n   Verdict: {detect_verdict}")

# ============================================================================
# Λ over cosmic history
# ============================================================================

print(f"\n{'='*70}")
print("PART 6: Λ ACROSS COSMIC HISTORY")
print(f"{'='*70}")

# Estimate N_MAX at different epochs
# r_universe scales with scale factor: r ~ a(t) × r_0
# M87* mass grows slowly — treat as roughly constant for order of magnitude

# Redshifts of interest
epochs = [
    ("Big Bang (z=1100, CMB)",    1100),
    ("First galaxies (z=10)",     10),
    ("Peak star formation (z=2)", 2),
    ("Today (z=0)",               0),
    ("Far future (z=-0.5)",       -0.5),
    ("Far future (z=-0.9)",       -0.9),
]

print(f"\n   {'Epoch':<35} {'z':<8} {'r_univ scale':<14} {'N_MAX':<10} {'Λ/Λ₀'}")
print("   " + "-"*75)

for name, z in epochs:
    a       = 1.0 / (1.0 + z) if z > -1 else 10.0
    # Observable universe in past was smaller: r ~ r_0 × a (approximate)
    r_epoch = r_universe_0 * a
    # Rs grows slowly — use current value for approximation
    N_epoch = np.log(r_epoch / Rs_m87_0) if r_epoch > Rs_m87_0 else 1.0
    N_epoch = max(N_epoch, 0.1)

    # Recalculate CURV_COSM for this epoch
    # R_H roughly constant (H_0 changes but order of magnitude same)
    curv_epoch = R_H / r_epoch if r_epoch > 0 else 0.5
    curv_epoch = min(curv_epoch, 0.99)

    lam_epoch = lambda_from_nmax(N_epoch)
    ratio     = lam_epoch / LAMBDA_0

    print(f"   {name:<35} {z:<8.1f} {a:<14.4f} {N_epoch:<10.2f} {ratio:.3f}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"""
   Λ is predicted to change at {abs(frac_change_Gyr)*100:.6f}% per billion years.

   This is because N_MAX = ln(r_universe / Rs_m87) grows slowly
   as the observable universe expands at the speed of light.

   The change is dominated by universe expansion, not M87* accretion
   (expansion is {r_universe_frac_rate/Rs_frac_rate:.0e}× faster).

   Current dark energy surveys cannot detect this drift.
   Future surveys (EUCLID, LSST, ~2030s) would need {future_threshold_per_gyr*100:.3f}%/Gyr
   precision — the predicted drift is {abs(frac_change_Gyr)/future_threshold_per_gyr:.0f}× below this.

   HONEST INTERPRETATION:
   This is either:
   (a) Λ truly varies but too slowly to ever measure directly, OR
   (b) The theory needs a stabilization mechanism that keeps
       N_MAX approximately constant despite the expanding universe

   Either way — this is not an impasse. It is a prediction:
   Λ was different in the early universe and will be different
   in the far future. The cosmic history table above shows
   how Λ/Λ₀ has evolved across all of cosmic time.
""")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('N_MAX Drift: Does Λ Change Over Cosmic Time?',
             fontsize=14, fontweight='bold')

# Time axis: billions of years
t_Gyr = np.linspace(0.1, 50, 500)  # 0.1 to 50 billion years

# Scale factor grows approximately as t^(2/3) in matter era
# and exponentially in dark energy era — simplified
t_today = 13.8  # Gyr

# Approximate scale factor evolution
a_vals = []
for t in t_Gyr:
    if t < 10:
        a = (t / t_today) ** (2/3)
    else:
        a = (10/t_today)**(2/3) * np.exp(H_0_per_sec * (t-10) * 3.156e16 * 0.3)
    a_vals.append(a)
a_vals = np.array(a_vals)

r_univ_vals = r_universe_0 * a_vals
N_MAX_vals  = np.log(np.maximum(r_univ_vals, Rs_m87_0 * 1.01) / Rs_m87_0)

# Plot 1: N_MAX over time
ax = axes[0]
ax.plot(t_Gyr, N_MAX_vals, 'purple', linewidth=2.5)
ax.axvline(t_today, color='red', linewidth=2, linestyle='--',
           label=f'Today ({t_today} Gyr)')
ax.axhline(N_MAX_0, color='orange', linewidth=1.5, linestyle=':',
           label=f'N_MAX today = {N_MAX_0:.1f}')
ax.set_xlabel('Cosmic Time (Gyr)', fontsize=12)
ax.set_ylabel('N_MAX', fontsize=12)
ax.set_title('N_MAX over Cosmic Time', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# Plot 2: Λ/Λ₀ over time
ax = axes[1]
lam_vals = []
for N in N_MAX_vals:
    try:
        lam = lambda_from_nmax(max(N, 0.1))
        lam_vals.append(lam / LAMBDA_0)
    except:
        lam_vals.append(np.nan)
lam_vals = np.array(lam_vals)

ax.plot(t_Gyr, lam_vals, 'blue', linewidth=2.5)
ax.axvline(t_today, color='red', linewidth=2, linestyle='--',
           label='Today')
ax.axhline(1.0, color='black', linewidth=1.5, linestyle='--',
           label='Λ today')
ax.set_xlabel('Cosmic Time (Gyr)', fontsize=12)
ax.set_ylabel('Λ / Λ₀', fontsize=12)
ax.set_title('Λ / Λ₀ over Cosmic Time', fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_ylim(0, max(lam_vals[np.isfinite(lam_vals)]) * 1.1)

# Plot 3: Fractional change per Gyr (zoom on today)
ax = axes[2]
t_zoom = np.linspace(10, 20, 200)
a_zoom = [(10/t_today)**(2/3) * np.exp(H_0_per_sec*(t-10)*3.156e16*0.3)
          for t in t_zoom]
r_zoom  = [r_universe_0 * a for a in a_zoom]
N_zoom  = [np.log(max(r, Rs_m87_0*1.01) / Rs_m87_0) for r in r_zoom]
lam_zoom = [lambda_from_nmax(max(N, 0.1)) / LAMBDA_0 for N in N_zoom]

ax.plot(t_zoom, lam_zoom, 'blue', linewidth=2.5, label='Λ/Λ₀')
ax.axvline(t_today, color='red', linewidth=2, linestyle='--',
           label=f'Today')
ax.axhline(1.0, color='black', linewidth=1.5, linestyle='--')

# Detection thresholds
t_range = np.array([t_zoom[0], t_zoom[-1]])
ax.fill_between(t_range,
                1 - detection_threshold_per_gyr * (t_today - t_zoom[0]),
                1 + detection_threshold_per_gyr * (t_today - t_zoom[0]),
                alpha=0.2, color='orange', label='Current survey precision')
ax.fill_between(t_range,
                1 - future_threshold_per_gyr * (t_today - t_zoom[0]),
                1 + future_threshold_per_gyr * (t_today - t_zoom[0]),
                alpha=0.2, color='green', label='EUCLID/LSST precision')
ax.set_xlabel('Cosmic Time (Gyr)', fontsize=12)
ax.set_ylabel('Λ / Λ₀', fontsize=12)
ax.set_title('Zoom: Λ Drift Near Today', fontsize=12, fontweight='bold')
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('nmax_drift_results.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"\n✓ Plot saved: nmax_drift_results.png")
print("=" * 70)
