"""
TASK 5 - De-circularize x(z) and re-run the DESI test
======================================================
The old ESTIF cosmology defined its tilt term through
    x(z) = x0 (1+z) H0 / H_LCDM(z)      <-- uses LCDM as its own ruler
which is circular: ESTIF's expansion history is anchored to LCDM's. This
scored chi^2/N = 10.8 against DESI DR2 (LCDM: ~1.9) and its w_eff ~ -1.36
was falsified at 3.5 sigma.

This script replaces the LCDM ruler with ESTIF's OWN H(z), making x(z)
self-consistent. At each redshift the Friedmann relation becomes a
fixed-point equation solved for H(z) with no LCDM anywhere:

    x(z)          = x0 (1+z) H0 / H(z)                 [self-referential]
    omega_tilt(z) = Omega_L * ( obs(x0) / obs(x(z)) )^2
    H(z)^2        = H0^2 [ Omega_m (1+z)^3 + omega_tilt(z) ]

Only the RULER changed. The tilt geometry (obs = sqrt(beta), the calibrated
N_MAX,B) is untouched, so this isolates the effect of removing circularity.

It then downloads the REAL DESI DR1+DR2 BAO data (same CobayaSampler source
the project uses) and computes chi^2 for three models on identical data:
    - LCDM
    - ESTIF old (circular ruler)
    - ESTIF new (self-consistent ruler)   [with and without the z<2 cutoff]

HONEST FRAMING (stated before the run): de-circularizing is a CORRECTNESS
fix, not a guarantee of a good fit. Three outcomes are possible:
  (A) chi^2/N drops toward ~1-2  -> the circular ruler was the problem.
  (B) chi^2/N stays ~10          -> the tilt SHAPE is wrong, not the ruler.
  (C) somewhere between          -> partial; ruler mattered but shape needs work.
The script reports whichever occurs. It does not decide in advance.

Run:  python3 estif_task5_desi_selfconsistent.py
Deps: numpy, scipy, matplotlib, internet (fetches DESI data once, caches).
"""

import os
import urllib.request
import numpy as np
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs")  # caches live in tests/docs
from scipy.integrate import quad
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# Constants (Planck 2018, matching estif_ec_gr_constants / the model)
# ----------------------------------------------------------------------
MPC = 3.085677581e22            # metres per Mpc
c = 2.99792458e8               # m/s
H0 = 67.66 * 1000.0 / MPC      # s^-1  (67.66 km/s/Mpc, Planck 2018)
OMEGA_M = 0.3111
OMEGA_LAMBDA = 0.6889
R_UNIV = 4.4e26                # m  (observable-universe radius used by project)
x0 = (c / H0) / R_UNIV         # ~0.3107
RD = 147.09                    # Mpc, Planck 2018 sound horizon
Z_EFF_MAX = 2.0

N_MAX = 33.265                 # calibrated tilt parameters (unchanged)
B = 15.429


# ----------------------------------------------------------------------
# Tilt geometry (identical to estif_ec_gr_model.observable_combined)
# ----------------------------------------------------------------------
def observable(x):
    x = np.asarray(x, dtype=float)
    n = N_MAX * np.exp(-B * x)
    val = np.where(x > 0, x ** (2.0 * n), 0.0)
    beta = np.where(val >= 1.0, 0.0, np.sqrt(np.maximum(0.0, 1.0 - val)))
    return np.sqrt(beta)


OBS_NOW = float(observable(x0))


# ----------------------------------------------------------------------
# Model H(z): LCDM, ESTIF-old (circular), ESTIF-new (self-consistent)
# ----------------------------------------------------------------------
def H_lcdm(z):
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3 + OMEGA_LAMBDA)


def omega_tilt_circular(z):
    z_eff = min(z, Z_EFF_MAX)
    Hl = H0 * np.sqrt(OMEGA_M * (1 + z_eff) ** 3 + OMEGA_LAMBDA)
    x_z = x0 * (1 + z_eff) * H0 / Hl            # <-- LCDM ruler (circular)
    obs_z = float(observable(x_z))
    if obs_z <= 0:
        return OMEGA_LAMBDA
    return OMEGA_LAMBDA * (OBS_NOW / obs_z) ** 2


def H_estif_circular(z):
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3 + omega_tilt_circular(z))


def _h_selfconsistent(z, use_cutoff):
    """Return h = H/H0 solving the fixed point with ESTIF's own ruler.
    z_tilt is optionally clamped at 2 (parity with the old cutoff)."""
    z_tilt = min(z, Z_EFF_MAX) if use_cutoff else z
    matter = OMEGA_M * (1 + z) ** 3

    def g(h):
        # x uses the SAME h (ESTIF ruler), not LCDM
        x = x0 * (1 + z_tilt) * 1.0 / h        # = x0 (1+z) H0 / H
        obs_z = float(observable(x))
        omega_tilt = OMEGA_LAMBDA * (OBS_NOW / obs_z) ** 2 if obs_z > 0 else OMEGA_LAMBDA
        return np.sqrt(matter + omega_tilt)

    # fixed-point iteration with damping; robust and fast here
    h = np.sqrt(matter + OMEGA_LAMBDA)         # seed (value, not a ruler)
    for _ in range(500):
        h_new = g(h)
        if abs(h_new - h) < 1e-13:
            h = h_new
            break
        h = 0.5 * h + 0.5 * h_new              # damped
    return h


def H_estif_sc(z, use_cutoff=False):
    return H0 * _h_selfconsistent(z, use_cutoff)


# ----------------------------------------------------------------------
# BAO distance observables
# ----------------------------------------------------------------------
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


# ----------------------------------------------------------------------
# Data loading (real DESI DR1 + DR2, cached)
# ----------------------------------------------------------------------
BASE1 = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master"
BASE2 = f"{BASE1}/desi_bao_dr2"
URLS = {
    'dr1_mean': f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
    'dr1_cov':  f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
    'dr2_mean': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    'dr2_cov':  f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt",
}


def _cache(key):
    p = os.path.join(CACHE_DIR, f"{key}.txt")
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[key], p)
    return p


def load_mean(key):
    rows = []
    for line in open(_cache(key)):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[1]), parts[2]))
    return rows


def load_cov(key, n):
    return np.loadtxt(_cache(key)).reshape(n, n)


def main():
    print("=" * 74)
    print("TASK 5 - DE-CIRCULARIZED ESTIF vs DESI (real data)")
    print("=" * 74)
    print(f"  x0 = {x0:.5f}   obs(x0) = {OBS_NOW:.5f}   H0 = 67.66 km/s/Mpc")
    print(f"  Models: LCDM | ESTIF-old (circular) | ESTIF-new (self-consistent)")
    print()

    dr1, dr2 = load_mean('dr1_mean'), load_mean('dr2_mean')
    cov1, cov2 = load_cov('dr1_cov', len(dr1)), load_cov('dr2_cov', len(dr2))
    obs1 = np.array([r[1] for r in dr1])
    obs2 = np.array([r[1] for r in dr2])
    n1, n2 = len(dr1), len(dr2)
    print(f"  Loaded real data: DESI DR1 = {n1} bins, DR2 = {n2} bins")
    print()

    # sanity: self-consistent x(z) at a few z, confirm no LCDM used
    print("  Self-consistent solve check (no LCDM anywhere):")
    for z in [0.0, 0.5, 1.0, 2.0, 2.33]:
        h = _h_selfconsistent(z, use_cutoff=False)
        x = x0 * (1 + min(z, 1e9)) / h
        print(f"    z={z:<5} H/H0={h:.4f}  x(z)={x:.4f}  "
              f"omega_tilt={OMEGA_LAMBDA*(OBS_NOW/float(observable(x)))**2:.4f}")
    print()

    models = [
        ("LCDM", H_lcdm),
        ("ESTIF old (circular)", H_estif_circular),
        ("ESTIF new (self-consistent)", lambda z: H_estif_sc(z, use_cutoff=False)),
        ("ESTIF new (self-consist,+cutoff)", lambda z: H_estif_sc(z, use_cutoff=True)),
    ]

    print("  Computing predictions on real DESI data ...", flush=True)
    results = {}
    for name, Hf in models:
        p1 = predict(dr1, Hf)
        p2 = predict(dr2, Hf)
        results[name] = (chi2(p1, obs1, cov1), chi2(p2, obs2, cov2), p2)
    print("  Done.")
    print()

    print("=" * 74)
    print("CHI-SQUARED / N  (lower is better; ~1 is a good fit)")
    print("=" * 74)
    print(f"  {'Model':<36}{'DR1 chi2/N':>12}{'DR2 chi2/N':>12}")
    print("  " + "-" * 60)
    for name, _ in models:
        c1, c2, _ = results[name]
        print(f"  {name:<36}{c1/n1:>12.3f}{c2/n2:>12.3f}")
    print()

    # per-bin pulls for the new model vs DR2
    _, _, p2_new = results["ESTIF new (self-consistent)"]
    _, _, p2_old = results["ESTIF old (circular)"]
    _, _, p2_lcdm = results["LCDM"]
    err2 = np.sqrt(np.diag(cov2))
    print("=" * 74)
    print("PER-BIN PULLS vs DESI DR2  (model - data)/sigma")
    print("=" * 74)
    print(f"  {'z':<7}{'qty':<13}{'data':>9}{'LCDM':>8}{'old':>8}{'new':>8}"
          f"{'pull_new':>10}")
    print("  " + "-" * 62)
    n1s_new = n2s_new = 0
    for i, (z, ov, qty) in enumerate(dr2):
        pull_new = (p2_new[i] - ov) / err2[i]
        if abs(pull_new) < 1:
            n1s_new += 1
        if abs(pull_new) < 2:
            n2s_new += 1
        print(f"  {z:<7.3f}{qty:<13}{ov:>9.3f}{p2_lcdm[i]:>8.3f}"
              f"{p2_old[i]:>8.3f}{p2_new[i]:>8.3f}{pull_new:>+10.2f}")
    print(f"\n  New model: {n1s_new}/{n2} bins within 1sigma, "
          f"{n2s_new}/{n2} within 2sigma")
    print()

    # w_eff(z~0) for old vs new (numerical)
    def w_eff(Hf, z=0.05, dz=0.02):
        # w from d ln(rho_de)/d ln a via H: rho_de = 3H^2/8piG - matter
        def rho_de(zz):
            Hh = Hf(zz)
            return (Hh / H0) ** 2 - OMEGA_M * (1 + zz) ** 3
        hi, lo = rho_de(z + dz), rho_de(max(z - dz, 1e-4))
        dln = (np.log(abs(hi) + 1e-30) - np.log(abs(lo) + 1e-30)) / (2 * dz)
        return -1.0 + (1.0 + z) / 3.0 * dln

    w_old = w_eff(H_estif_circular)
    w_new = w_eff(lambda z: H_estif_sc(z, use_cutoff=False))
    print("=" * 74)
    print("EFFECTIVE w(z~0)   (DESI DR2 w0 = -0.73 +/- 0.10)")
    print("=" * 74)
    print(f"  ESTIF old (circular):        w_eff = {w_old:+.3f}")
    print(f"  ESTIF new (self-consistent): w_eff = {w_new:+.3f}")
    print(f"  LCDM value:                  w_eff = -1.000")
    print()

    # verdict
    c2_old = results["ESTIF old (circular)"][1] / n2
    c2_new = results["ESTIF new (self-consistent)"][1] / n2
    c2_lcdm = results["LCDM"][1] / n2
    print("=" * 74)
    print("VERDICT (runtime facts, not pre-judged)")
    print("=" * 74)
    print(f"  LCDM              chi2/N = {c2_lcdm:.2f}")
    print(f"  ESTIF old         chi2/N = {c2_old:.2f}")
    print(f"  ESTIF new         chi2/N = {c2_new:.2f}")
    improve = c2_old - c2_new
    print(f"  Change from de-circularizing: {improve:+.2f} in chi2/N "
          f"({'better' if improve > 0 else 'worse'})")
    print()
    if c2_new < 2.0:
        print("  OUTCOME A: self-consistent ruler brings ESTIF to an acceptable")
        print("  fit (chi2/N < 2). The circularity was the dominant problem.")
    elif c2_new < 0.6 * c2_old:
        print("  OUTCOME C: large improvement but not yet acceptable. The ruler")
        print("  mattered; the tilt SHAPE still needs work at high z.")
    else:
        print("  OUTCOME B: de-circularizing did NOT rescue the fit. The tilt")
        print("  SHAPE (obs = sqrt(beta) with these N_MAX,B), not the ruler, is")
        print("  the primary source of the DESI tension. This localizes the")
        print("  remaining problem to the geometry, not the self-reference.")
    print()
    print("  Either way: the circular dependency is REMOVED. x(z) now derives")
    print("  from ESTIF's own H(z). This is the correctness fix Task 5 targeted;")
    print("  the chi2 outcome above tells you where the next work must go.")
    print("=" * 74)


if __name__ == "__main__":
    main()
