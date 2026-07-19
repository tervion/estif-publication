"""
PATH ONE - Canonical cosmology test (ESTIF-Core)
================================================
THE CLAIM: the cosmic eddy energy density is CONSTANT. By the derived field
equation (Task 4), a constant effective energy density is exact de Sitter -- a
cosmological constant, w = -1. So ESTIF-Core's expansion history is

    H^2(z) = H0^2 [ Omega_m (1+z)^3 + Omega_Lambda ]

which is mathematically IDENTICAL to LCDM. Two honest consequences:
 (1) On expansion-history data ESTIF-Core is indistinguishable from LCDM, by
     construction -- a feature: as good as the standard model, reinterpreting
     dark energy as the frozen cosmic eddy, with no tilt formula, no fitted
     evolution, no circular ruler.
 (2) ESTIF adds a consistency relation: Omega_m = x0 = R_H/r_universe (0.12% from
     Planck). This is NOT an independent prediction -- r_universe is here the LCDM
     particle horizon, an integral that itself contains Omega_m (see CORRECTIONS
     v6.3.1, C2). An independent prediction requires deriving r_universe from the
     flow framework (RHAC Scenario H; partially addressed by the Omega_m bootstrap).

Run:  python3 estif_pathone_cosmology.py
Deps: numpy, scipy, internet (DESI fetch, cached beside script).
"""
import os, urllib.request
import numpy as np
from scipy.integrate import quad

CACHE_DIR = os.path.dirname(os.path.abspath(__file__))
MPC = 3.085677581e22
c = 2.99792458e8
H0 = 67.66 * 1000.0 / MPC
OMEGA_M_PLANCK = 0.3111
RD = 147.09
R_UNIV = 4.4e26
x0 = (c / H0) / R_UNIV
OMEGA_M_ESTIF = x0

def H_frozen(z, Om):
    return H0 * np.sqrt(Om * (1 + z) ** 3 + (1.0 - Om))
def H_lcdm(z): return H_frozen(z, OMEGA_M_PLANCK)
def H_estif_core(z): return H_frozen(z, OMEGA_M_ESTIF)

def DH(z, Hf): return c / (Hf(z) * MPC)
def DM(z, Hf):
    if z <= 0: return 0.0
    v, _ = quad(lambda zp: c / (Hf(zp) * MPC), 0, z, limit=200); return v
def DV(z, Hf): return (z * DH(z, Hf) * DM(z, Hf) ** 2) ** (1 / 3)
def predict(rows, Hf):
    out = []
    for z, _, q in rows:
        out.append((DV if q=='DV_over_rs' else DM if q=='DM_over_rs' else DH)(z, Hf) / RD)
    return np.array(out)
def chi2(pred, obs, cov):
    d = obs - pred
    try: return float(d @ np.linalg.inv(cov) @ d)
    except np.linalg.LinAlgError:
        return float(np.sum((d / np.sqrt(np.diag(cov))) ** 2))

BASE1 = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master"
BASE2 = f"{BASE1}/desi_bao_dr2"
URLS = {'dr1_m': f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
        'dr1_c': f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
        'dr2_m': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        'dr2_c': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt"}
def cache(k):
    p = os.path.join(CACHE_DIR, f"pathone_{k}.txt")
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[k], p)
    return p
def load_mean(k):
    rows = []
    for line in open(cache(k)):
        line = line.strip()
        if line and not line.startswith('#'):
            p = line.split(); rows.append((float(p[0]), float(p[1]), p[2]))
    return rows
def load_cov(k, n): return np.loadtxt(cache(k)).reshape(n, n)

def main():
    print("=" * 74)
    print("PATH ONE - CANONICAL COSMOLOGY TEST (ESTIF-Core)")
    print("Frozen cosmic eddy = cosmological constant")
    print("=" * 74)
    print(f"  Omega_m (Planck)     = {OMEGA_M_PLANCK:.4f}")
    print(f"  Omega_m (ESTIF = x0) = {OMEGA_M_ESTIF:.4f}   [consistency relation, see C2]")
    print(f"  |x0 - Om_Planck|     = {abs(x0-OMEGA_M_PLANCK)/OMEGA_M_PLANCK*100:.2f}%")
    print(f"  Omega_Lambda (ESTIF) = {1-OMEGA_M_ESTIF:.4f}   [= 1 - x0, predicted]")
    print()
    print("-" * 74); print("[1] frozen-eddy H(z) IS the LCDM form"); print("-" * 74)
    print(f"  {'z':<6}{'H_LCDM/H0':>12}{'H_ESTIF/H0':>13}")
    for z in [0.0, 0.5, 1.0, 2.0]:
        print(f"  {z:<6}{H_lcdm(z)/H0:>12.5f}{H_estif_core(z)/H0:>13.5f}")
    print("  -> same form; only Omega_m differs (x0 vs Planck, 0.12%).")
    print()
    dr1, dr2 = load_mean('dr1_m'), load_mean('dr2_m')
    cov1, cov2 = load_cov('dr1_c', len(dr1)), load_cov('dr2_c', len(dr2))
    o1, o2 = np.array([r[1] for r in dr1]), np.array([r[1] for r in dr2])
    n1, n2 = len(dr1), len(dr2)
    print("-" * 74); print("[2] Real DESI test (DR1 + DR2)"); print("-" * 74)
    print(f"  {'model':<40}{'DR1 chi2/N':>12}{'DR2 chi2/N':>12}")
    print("  " + "-" * 64)
    res = {}
    for name, Hf in [("LCDM (Om = Planck 0.3111)", H_lcdm),
                     ("ESTIF-Core frozen eddy (Om = x0)", H_estif_core)]:
        c1 = chi2(predict(dr1, Hf), o1, cov1) / n1
        c2 = chi2(predict(dr2, Hf), o2, cov2) / n2
        res[name] = (c1, c2)
        print(f"  {name:<40}{c1:>12.3f}{c2:>12.3f}")
    print()
    err2 = np.sqrt(np.diag(cov2)); p_core = predict(dr2, H_estif_core)
    n1s = n2s = 0
    print("-" * 74); print("[3] Per-bin pulls, ESTIF-Core vs DESI DR2"); print("-" * 74)
    print(f"  {'z':<7}{'qty':<13}{'data':>9}{'ESTIF':>9}{'pull':>9}")
    for i, (z, ov, q) in enumerate(dr2):
        pull = (p_core[i] - ov) / err2[i]; n1s += abs(pull) < 1; n2s += abs(pull) < 2
        print(f"  {z:<7.3f}{q:<13}{ov:>9.3f}{p_core[i]:>9.3f}{pull:>+9.2f}")
    print(f"\n  ESTIF-Core: {n1s}/{n2} within 1sigma, {n2s}/{n2} within 2sigma")
    print()
    c2_lcdm = res["LCDM (Om = Planck 0.3111)"][1]
    c2_core = res["ESTIF-Core frozen eddy (Om = x0)"][1]
    print("=" * 74); print("VERDICT"); print("=" * 74)
    print(f"  LCDM (Planck Om):          chi2/N = {c2_lcdm:.3f}")
    print(f"  ESTIF-Core (Om = x0):      chi2/N = {c2_core:.3f}")
    print(f"  Difference:                {abs(c2_core-c2_lcdm):.3f}")
    print()
    if abs(c2_core - c2_lcdm) < 0.1:
        print("  RESULT: ESTIF-Core ties LCDM on DESI (same expansion history).")
        print("  Path One cosmology CONFIRMED:")
        print("    - dark energy = frozen cosmic eddy = cosmological constant (w=-1);")
        print("    - expansion history indistinguishable from LCDM (by construction);")
        print("    - Omega_m = x0 is a CONSISTENCY RELATION (r_universe is the LCDM horizon, which itself depends on Omega_m); NOT an Omega_m-independent prediction;")
        print("    - NO tilt formula, NO fitted evolution, NO circular ruler.")
    else:
        print("  RESULT: unexpected difference -- investigate.")
    print()
    print("  HONEST SCOPE: cannot beat LCDM on expansion data (same H(z)). ESTIF's")
    print("  distinguishing content is the DERIVED gravity sector and, if derivable,")
    print("  the Path Two thawing correction. Path One cosmology provides parity +")
    print("  a geometric Omega_m CONSISTENCY relation + a physical dark-energy")
    print("  interpretation, with no fitted DARK-ENERGY parameters (w=-1 fixed).")
    print("=" * 74)

if __name__ == "__main__":
    main()
