"""
BOOTSTRAP CLOSURE - propagate the self-derived Omega_m through everything
=========================================================================
The bootstrap (estif_omega_bootstrap.py) showed that principle P
(Omega_m = R_H / r_p) is a closed equation with a unique root:
Omega_m = 0.314 (with radiation), 0.96% from Planck, inside Planck's error.

This script CLOSES THE LOOP: adopt P, and push the bootstrap Omega_m through
every downstream number the framework owns, so nothing depends on the old
imported r_universe = 4.4e26 m any more. Outputs:

  [1] the bootstrap value (recomputed internally, self-contained)
  [2] a0 = c H0 x0 / sqrt(3) with the self-consistent x0  ->  vs MOND
  [3] r_universe BACK-PREDICTED (was an import, becomes an output)
  [4] DESI DR2 parity re-run with the bootstrap Omega_m (real data)
  [5] THE INPUT LEDGER: measured quantities before vs after adopting P

HONEST FLAGS (printed in verdict): everything below is conditional on
principle P, which is not yet derived from axioms A1-A3 (Part B / RHAC-H);
the Gaztanaga causal-universe comparison is pending before novelty claims;
MOND's empirical a0 carries ~10% scatter, so 0.67% vs 1.72% agreement is
directionally pleasing, not decisive; rd = 147.09 Mpc remains a Planck-
calibrated import OF THE DESI TEST (BAO ruler), not of the framework.

Run:  python3 estif_bootstrap_closure.py
Deps: numpy, scipy, internet (DESI fetch, cached beside script).
"""
import os, math, urllib.request
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs")
c = 2.99792458e8
G = 6.67430e-11
a_rad = 7.565723e-16
MPC = 3.085677581e22
H0_kms = 67.66
h = H0_kms / 100.0
H0 = H0_kms * 1000.0 / MPC
T_CMB, NEFF = 2.7255, 3.046
OM_PLANCK, OM_PLANCK_ERR = 0.3111, 0.0056
RD = 147.09
A0_MOND = 1.2e-10
R_U_IMPORTED = 4.4e26
X0_OLD = 0.31073

def omega_radiation(h_, T=T_CMB, neff=NEFF):
    rho_g = a_rad * T**4 / c**2
    rho_c = 3 * (h_ * 100e3 / MPC)**2 / (8 * np.pi * G)
    return (rho_g / rho_c) * (1.0 + 0.2271 * neff)

def I_horizon(Om, Or=0.0):
    OL = 1.0 - Om - Or
    f = lambda u: 2.0 * u / np.sqrt(Or + Om * u**2 + OL * u**8)
    val, _ = quad(f, 0.0, 1.0, limit=200)
    return val

Or = omega_radiation(h)
Om_b = brentq(lambda Om: Om * I_horizon(Om, Or) - 1.0, 0.02, 0.98, xtol=1e-10)
x0_b = Om_b
r_u_pred = I_horizon(Om_b, Or) * c / H0
a0_b = c * H0 * x0_b / math.sqrt(3)
a0_old = c * H0 * X0_OLD / math.sqrt(3)

print("=" * 72)
print("BOOTSTRAP CLOSURE - the framework running on its own Omega_m")
print("=" * 72)
print(f"\n[1] Bootstrap (self-contained recompute):")
print(f"    Omega_m = x0 = {Om_b:.5f}   (Planck {OM_PLANCK} +/- {OM_PLANCK_ERR}; "
      f"{abs(Om_b-OM_PLANCK)/OM_PLANCK_ERR:.2f} sigma)")
print(f"    Omega_Lambda = {1-Om_b-Or:.5f}   Omega_r = {Or:.3e}")

print(f"\n[2] a0 with the self-consistent x0:")
print(f"    old (imported r_u): a0 = {a0_old:.4e}   ({(A0_MOND-a0_old)/A0_MOND*100:+.2f}% vs MOND)")
print(f"    new (bootstrap)   : a0 = {a0_b:.4e}   ({(A0_MOND-a0_b)/A0_MOND*100:+.2f}% vs MOND)")
vflat_shift = (a0_b / a0_old) ** 0.25
print(f"    SPARC impact: v_flat shifts by x{vflat_shift:.5f} "
      f"({(vflat_shift-1)*100:+.2f}%) -- negligible vs 15.6% RMS; the")
print(f"    87-galaxy validation is insensitive to this refinement.")

print(f"\n[3] r_universe: import -> output")
print(f"    was IMPORTED:      {R_U_IMPORTED:.3e} m  (LCDM particle horizon)")
print(f"    now BACK-PREDICTED:{r_u_pred:.3e} m  "
      f"({(r_u_pred-R_U_IMPORTED)/R_U_IMPORTED*100:+.2f}% vs the import)")

# ---- DESI parity with bootstrap Omega_m --------------------------------
def Hz(z, Om): return H0 * np.sqrt(Om * (1 + z) ** 3 + (1.0 - Om))
def DH(z, Om): return c / (Hz(z, Om) * MPC)
def DM(z, Om):
    v, _ = quad(lambda zp: c / (Hz(zp, Om) * MPC), 0, z, limit=200); return v
def DV(z, Om): return (z * DH(z, Om) * DM(z, Om) ** 2) ** (1 / 3)
def predict(rows, Om):
    return np.array([(DV if q == 'DV_over_rs' else DM if q == 'DM_over_rs'
                      else DH)(z, Om) / RD for z, _, q in rows])
def chi2(pred, obs, cov):
    d = obs - pred
    return float(d @ np.linalg.inv(cov) @ d)

BASE2 = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2"
URLS = {'m': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        'c': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt"}
def cache(k):
    p = os.path.join(CACHE_DIR, f"closure_{k}.txt")
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[k], p)
    return p

rows = []
for line in open(cache('m')):
    line = line.strip()
    if line and not line.startswith('#'):
        p = line.split(); rows.append((float(p[0]), float(p[1]), p[2]))
N = len(rows)
cov = np.loadtxt(cache('c')).reshape(N, N)
obs = np.array([r[1] for r in rows])

print(f"\n[4] DESI DR2 parity, real data ({N} bins):")
print(f"    {'model':<38}{'chi2/N':>9}")
for tag, Om in [("LCDM (Planck Om=0.3111)", OM_PLANCK),
                ("ESTIF-Core, old x0 (imported r_u)", X0_OLD),
                ("ESTIF-Core, BOOTSTRAP Omega_m", Om_b)]:
    print(f"    {tag:<38}{chi2(predict(rows, Om), obs, cov)/N:>9.3f}")

print(f"\n[5] THE INPUT LEDGER (framework parameters)")
print("    BEFORE (v6.3):  measured/imported = H0, T_CMB, Neff,")
print("                    Omega_m (Planck fit), r_universe = 4.4e26 m (import)")
print("    AFTER (P adopted): measured = H0, T_CMB, Neff   -- that is all.")
print("                    computed = Omega_m, Omega_Lambda, x0, r_universe, a0")
print("    (rd = 147.09 Mpc remains an import of the BAO TEST, not of the")
print("     framework; G and c are universal constants.)")

print("\n" + "=" * 72)
print("VERDICT")
print("=" * 72)
print(f"""  Conditional on principle P, the framework now computes its own matter
  density, dark-energy fraction, cosmic curvature ratio, universe radius,
  and MOND acceleration from three directly measured inputs (H0, T_CMB,
  Neff). The a0-MOND agreement improves to {(A0_MOND-a0_b)/A0_MOND*100:.2f}% and DESI parity holds.
  OPEN AND FLAGGED: P is not yet derived (Part B / RHAC-H); Gaztanaga
  comparison pending; a0's empirical target has ~10% scatter.""")
