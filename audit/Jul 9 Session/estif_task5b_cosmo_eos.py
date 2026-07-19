"""
TASK 5b - Cosmological flow: what expansion history does DESI demand,
and can the ESTIF eddy reach it?
=====================================================================
Task 5 removed the circular ruler; the self-consistent ESTIF still sat at
chi^2/N = 3.35 vs LCDM's 1.92, with the residual concentrated at low z.
That localizes the problem to the SHAPE of the dark-energy term, not the
ruler. This script attacks that shape honestly.

What is DERIVED vs TESTED here
  DERIVED (Task 4 + the flow Friedmann analog):
    - flat expanding flow reproduces standard Friedmann H^2 = (8 pi G/3) rho
    - a CONSTANT eddy energy density => de Sitter => w = -1 exactly
      (the cosmological constant is the frozen-eddy limit of the flow)
    - so any w != -1 REQUIRES the eddy energy density to EVOLVE with z
  TESTED (the open physical question the framework does not yet fix):
    - HOW does the cosmic eddy energy dilute as the universe expands?
      Parameterize it by its equation of state w (constant), or w0,wa (CPL):
          rho_DE(z) = rho_DE0 * (1+z)^{3(1+w)}          [constant w]
      w is not a fudge factor: it is the physical dilution law of the eddy.
      w = -1 frozen (Lambda);  w = 0 matter-like;  w = -1/3 curvature-like.

The script then asks three questions against REAL DESI DR2 data:
  Q1  What single w best fits DESI in a flat flow cosmology? Is it a good fit?
  Q2  What w0,wa (evolving eddy) does DESI prefer, and how good is that fit?
  Q3  Does the ESTIF self-consistent tilt form (Task 5) sit near the DESI-
      preferred w(z) curve, or off it? i.e. is the tilt geometry pointing
      in the direction the data wants?

HONEST FRAMING (before the run): finding a w that fits is NOT itself a
derivation -- it is a CONSTRAINT the flow eddy must satisfy. The result is
useful either way: it tells you what the eddy must physically do, and
whether the tilt geometry already leans that way. It does not manufacture
a first-principles dark energy.

Run:  python3 estif_task5b_cosmo_eos.py
Deps: numpy, scipy, internet (fetches DESI data once, caches).
"""

import os
import urllib.request
import numpy as np
CACHE_DIR = os.path.dirname(os.path.abspath(__file__))  # portable: write beside this script
from scipy.integrate import quad

# ---- constants (Planck 2018, matching the project) -------------------
MPC = 3.085677581e22
c = 2.99792458e8
H0 = 67.66 * 1000.0 / MPC
OMEGA_M = 0.3111
OMEGA_DE = 1.0 - OMEGA_M
R_UNIV = 4.4e26
x0 = (c / H0) / R_UNIV
RD = 147.09
N_MAX, B = 33.265, 15.429


def observable(x):
    x = np.asarray(x, dtype=float)
    n = N_MAX * np.exp(-B * x)
    val = np.where(x > 0, x ** (2.0 * n), 0.0)
    beta = np.where(val >= 1.0, 0.0, np.sqrt(np.maximum(0.0, 1.0 - val)))
    return np.sqrt(beta)


OBS_NOW = float(observable(x0))


# ---- model H(z) families ---------------------------------------------
def H_lcdm(z):
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3 + OMEGA_DE)


def H_wconst(z, w):
    """Flat flow cosmology with a dark-energy component of constant EoS w.
    rho_DE(z) = rho_DE0 (1+z)^{3(1+w)}.  w=-1 recovers LCDM."""
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3
                        + OMEGA_DE * (1 + z) ** (3.0 * (1.0 + w)))


def H_cpl(z, w0, wa):
    """CPL evolving EoS w(z) = w0 + wa z/(1+z).
    rho_DE(z)/rho_DE0 = (1+z)^{3(1+w0+wa)} exp(-3 wa z/(1+z))."""
    fde = (1 + z) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * z / (1 + z))
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3 + OMEGA_DE * fde)


def H_estif_sc(z):
    """ESTIF self-consistent tilt (Task 5), no LCDM ruler."""
    matter = OMEGA_M * (1 + z) ** 3
    h = np.sqrt(matter + OMEGA_DE)
    for _ in range(500):
        x = x0 * (1 + z) / h
        obs_z = float(observable(x))
        om_tilt = OMEGA_DE * (OBS_NOW / obs_z) ** 2 if obs_z > 0 else OMEGA_DE
        h_new = np.sqrt(matter + om_tilt)
        if abs(h_new - h) < 1e-13:
            h = h_new
            break
        h = 0.5 * h + 0.5 * h_new
    return H0 * h


def w_estif_sc(z, dz=0.02):
    """Effective EoS of the ESTIF self-consistent dark-energy term."""
    def rho_de(zz):
        return (H_estif_sc(zz) / H0) ** 2 - OMEGA_M * (1 + zz) ** 3
    hi, lo = rho_de(z + dz), rho_de(max(z - dz, 1e-4))
    dln = (np.log(abs(hi) + 1e-30) - np.log(abs(lo) + 1e-30)) / (2 * dz)
    return -1.0 + (1.0 + z) / 3.0 * dln


# ---- BAO observables --------------------------------------------------
def DH(z, Hf):
    return c / (Hf(z) * MPC)


def DM(z, Hf):
    if z <= 0:
        return 0.0
    val, _ = quad(lambda zp: c / (Hf(zp) * MPC), 0, z, limit=200)
    return val


def DV(z, Hf):
    return (z * DH(z, Hf) * DM(z, Hf) ** 2) ** (1 / 3)


def predict(rows, Hf):
    out = []
    for z, _, qty in rows:
        if qty == 'DV_over_rs':
            out.append(DV(z, Hf) / RD)
        elif qty == 'DM_over_rs':
            out.append(DM(z, Hf) / RD)
        elif qty == 'DH_over_rs':
            out.append(DH(z, Hf) / RD)
    return np.array(out)


def chi2(pred, obs, cov):
    d = obs - pred
    try:
        return float(d @ np.linalg.inv(cov) @ d)
    except np.linalg.LinAlgError:
        return float(np.sum((d / np.sqrt(np.diag(cov))) ** 2))


# ---- data -------------------------------------------------------------
BASE2 = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2"
URLS = {'m': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        'c': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt"}


def cache(k):
    p = os.path.join(CACHE_DIR, f"dr2b_{k}.txt")
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[k], p)
    return p


def main():
    rows = []
    for line in open(cache('m')):
        line = line.strip()
        if line and not line.startswith('#'):
            a = line.split()
            rows.append((float(a[0]), float(a[1]), a[2]))
    n = len(rows)
    cov = np.loadtxt(cache('c')).reshape(n, n)
    obs = np.array([r[1] for r in rows])

    print("=" * 74)
    print("TASK 5b - COSMOLOGICAL EQUATION OF STATE vs REAL DESI DR2")
    print("=" * 74)
    print(f"  Real DESI DR2: {n} bins.  LCDM and ESTIF-SC shown as anchors.")
    print()

    c_lcdm = chi2(predict(rows, H_lcdm), obs, cov) / n
    c_estif = chi2(predict(rows, H_estif_sc), obs, cov) / n
    print(f"  ANCHORS:   LCDM chi2/N = {c_lcdm:.3f}    "
          f"ESTIF self-consistent tilt chi2/N = {c_estif:.3f}")
    print()

    # Q1: best constant w
    print("=" * 74)
    print("Q1  BEST CONSTANT-w FLOW COSMOLOGY")
    print("    (w = physical dilution law of the cosmic eddy; w=-1 is frozen)")
    print("=" * 74)
    ws = np.linspace(-1.6, -0.4, 241)
    cw = [chi2(predict(rows, lambda z, w=w: H_wconst(z, w)), obs, cov) / n
          for w in ws]
    iw = int(np.argmin(cw))
    w_best, c_wbest = ws[iw], cw[iw]
    print(f"  best constant w = {w_best:+.3f}   chi2/N = {c_wbest:.3f}")
    print(f"  (LCDM is w=-1 exactly: chi2/N = {c_lcdm:.3f})")
    for wv in [-1.0, -0.95, -0.9, w_best]:
        cc = chi2(predict(rows, lambda z, w=wv: H_wconst(z, w)), obs, cov) / n
        tag = '  <- best' if abs(wv - w_best) < 1e-9 else ''
        print(f"    w={wv:+.3f}:  chi2/N = {cc:.3f}{tag}")
    print()

    # Q2: best CPL (evolving eddy)
    print("=" * 74)
    print("Q2  BEST EVOLVING-w (CPL) FLOW COSMOLOGY")
    print("    w(z) = w0 + wa z/(1+z).  This is what DESI itself prefers.")
    print("=" * 74)
    w0s = np.linspace(-1.1, -0.3, 81)
    was = np.linspace(-2.5, 1.0, 71)
    best = (1e9, None, None)
    for w0 in w0s:
        for wa in was:
            cc = chi2(predict(rows, lambda z, a=w0, b=wa: H_cpl(z, a, b)), obs, cov)
            if cc < best[0]:
                best = (cc, w0, wa)
    c_cpl, w0b, wab = best[0] / n, best[1], best[2]
    print(f"  best CPL:  w0 = {w0b:+.3f}   wa = {wab:+.3f}   chi2/N = {c_cpl:.3f}")
    print(f"  DESI DR2 published (DESI+CMB+Union3): w0=-0.73, wa=-0.66")
    c_desipub = chi2(predict(rows, lambda z: H_cpl(z, -0.73, -0.66)), obs, cov) / n
    print(f"  at DESI's published w0,wa:  chi2/N = {c_desipub:.3f}")
    print()

    # Q3: where does ESTIF's derived w(z) sit vs the DESI-preferred curve?
    print("=" * 74)
    print("Q3  DOES THE ESTIF TILT LEAN THE WAY DESI WANTS?")
    print("=" * 74)
    print(f"  {'z':<6}{'w_ESTIF(z)':>12}{'w_CPLbest(z)':>14}{'w_DESIpub(z)':>14}")
    print("  " + "-" * 44)
    for z in [0.05, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        we = w_estif_sc(z)
        wc = w0b + wab * z / (1 + z)
        wd = -0.73 - 0.66 * z / (1 + z)
        print(f"  {z:<6.2f}{we:>12.3f}{wc:>14.3f}{wd:>14.3f}")
    print()
    w_e0, w_e1 = w_estif_sc(0.05), w_estif_sc(1.0)
    trend_estif = "increasing" if w_e1 > w_e0 else "decreasing"
    trend_desi = "increasing" if (-0.73) > (-0.73 - 0.66 * 0.5) else "decreasing"
    print(f"  ESTIF w(z) trend from z=0 to z=1: {trend_estif} "
          f"(w: {w_e0:+.3f} -> {w_e1:+.3f})")
    print(f"  DESI-preferred w(z) trend:        thawing/increasing toward z=0")
    print()

    print("=" * 74)
    print("VERDICT (runtime facts)")
    print("=" * 74)
    print(f"  LCDM (w=-1):                       chi2/N = {c_lcdm:.2f}")
    print(f"  ESTIF self-consistent tilt:        chi2/N = {c_estif:.2f}")
    print(f"  best constant-w flow (w={w_best:+.2f}):     chi2/N = {c_wbest:.2f}")
    print(f"  best evolving-w flow (CPL):        chi2/N = {c_cpl:.2f}")
    print(f"    at CPL w0={w0b:+.2f}, wa={wab:+.2f}")
    print()
    if c_cpl <= c_lcdm + 0.3:
        print("  KEY RESULT: an EVOLVING eddy EoS reaches (or beats) LCDM. The")
        print("  flow framework CAN fit DESI -- provided the cosmic eddy energy")
        print("  evolves with the w0,wa above. The remaining task is to DERIVE")
        print("  that evolution from the eddy physics, not to fit it.")
    elif c_wbest < c_estif:
        print("  KEY RESULT: even a simple constant-w eddy beats the current")
        print("  tilt form. The tilt functional (obs ratio)^2 is not the best")
        print("  the flow can do; a cleaner eddy dilution law fits better.")
    else:
        print("  KEY RESULT: neither constant nor CPL eddy clearly beats the")
        print("  tilt form here; the tension may be in the data covariances or")
        print("  need a genuinely different expansion shape.")
    print()
    print("  Interpretation of the DESI-preferred w(z): it is THAWING (w more")
    print("  negative in the past, rising toward now). Physically that means")
    print("  the cosmic eddy energy was MORE dominant / diluting SLOWER early")
    print("  on. Whether ESTIF's tilt geometry predicts exactly this thawing")
    print("  shape is the precise, now-isolated question for the derivation.")
    print()
    print("  HONEST BOUND: a good-fitting w0,wa is a CONSTRAINT the eddy must")
    print("  meet, not a derived result. Task 5b tells you the target and")
    print("  whether the flow can reach it; deriving the eddy's w(z) from the")
    print("  4D kinetic energy projection is the remaining cosmology work.")
    print("=" * 74)


if __name__ == "__main__":
    main()
