#!/usr/bin/env python3
"""
PATH TWO - TARGET LOCK (pre-registration)
=========================================
Fits the CPL thawing shape w(a) = w0 + wa(1-a) on the EXACT Path One
DESI DR2 pipeline (same bins, same fixed ruler RD=147.09, same H0,
Omega_m FIXED at x0 -- the derivation gets no knob there either).

Purpose: freeze the bullseye BEFORE the vorticity derivation exists.
  - verifies the chi2/N ~ 0.66 target from runtime, not memory
  - records the (w0, wa) box the derived w(a) must land in
  - COMMIT THIS + ITS OUTPUT TO GIT BEFORE PHASE 2 BEGINS.

Convention: chi2/N with N = 13 bins (matches Path One reporting);
chi2/dof with dof = 11 also printed for honesty (2 fitted params here --
the DERIVATION will have zero, so its target is the chi2 value itself).
CAVEAT (C4): fixed rd and H0; marginalizing both could shift orderings.
"""
import os, numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
MPC = 3.085677581e22
c   = 2.99792458e8
H0  = 67.66 * 1000.0 / MPC
RD  = 147.09
R_UNIV = 4.4e26
x0  = (c / H0) / R_UNIV          # Om fixed = 0.3107, no refit

def E_cpl(z, w0, wa):
    a = 1.0 / (1.0 + z)
    de = (1.0 - x0) * a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
    return np.sqrt(x0 * (1.0 + z) ** 3 + de)

def Hf(z, w0, wa): return H0 * E_cpl(z, w0, wa)
def DHf(z, p): return c / (Hf(z, *p) * MPC)
def DMf(z, p):
    v, _ = quad(lambda zp: c / (Hf(zp, *p) * MPC), 0, z, limit=200); return v
def DVf(z, p): return (z * DHf(z, p) * DMf(z, p) ** 2) ** (1 / 3)

def load_mean(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if line and not line.startswith('#'):
            q = line.split(); rows.append((float(q[0]), float(q[1]), q[2]))
    return rows

rows = load_mean(os.path.join(HERE, "..", "docs", "pathone_dr2_m.txt"))
obs  = np.array([r[1] for r in rows])
cov  = np.loadtxt(os.path.join(HERE, "..", "docs", "pathone_dr2_c.txt")).reshape(len(rows), len(rows))
icov = np.linalg.inv(cov)

def predict(p):
    out = []
    for z, _, q in rows:
        out.append((DVf if q == 'DV_over_rs' else DMf if q == 'DM_over_rs' else DHf)(z, p) / RD)
    return np.array(out)

def chi2(p):
    d = obs - predict(p)
    return float(d @ icov @ d)

N = len(rows)
print("=" * 74)
print("PATH TWO - TARGET LOCK: CPL (w0, wa) on Path One DESI DR2 pipeline")
print(f"  bins N = {N} | Om = x0 = {x0:.4f} FIXED | RD = {RD} FIXED | H0 = 67.66 FIXED")
print("=" * 74)

c2_lcdm = chi2((-1.0, 0.0))
print(f"\n[0] Pipeline check, (w0,wa)=(-1,0) must equal Path One ESTIF-Core:")
print(f"    chi2 = {c2_lcdm:.3f}   chi2/N = {c2_lcdm/N:.3f}   (Path One printed 1.965)")

# coarse-to-fine grid: robust, no optimizer pathologies
w0g = np.linspace(-2.0, 0.0, 201)
wag = np.linspace(-4.0, 3.0, 141)
C = np.array([[chi2((w0, wa)) for wa in wag] for w0 in w0g])
i, j = np.unravel_index(np.argmin(C), C.shape)
# fine polish around minimum
w0f = np.linspace(max(-2, w0g[i]-0.06), min(0, w0g[i]+0.06), 61)
waf = np.linspace(wag[j]-0.30, wag[j]+0.30, 61)
Cf = np.array([[chi2((w0, wa)) for wa in waf] for w0 in w0f])
ii, jj = np.unravel_index(np.argmin(Cf), Cf.shape)
w0b, wab, c2b = w0f[ii], waf[jj], Cf[ii, jj]

d = C - c2b
in1 = d <= 2.30   # 1-sigma, 2 dof
print(f"\n[1] Best-fit thawing shape (2 knobs, for TARGET definition only):")
print(f"    w0 = {w0b:+.3f}   wa = {wab:+.3f}")
print(f"    chi2 = {c2b:.3f}   chi2/N = {c2b/N:.3f}   chi2/dof(11) = {c2b/(N-2):.3f}")
print(f"    Delta chi2 vs w=-1:  {c2_lcdm - c2b:.2f}   (AIC penalty for 2 params = 4)")
print(f"    w(z=0) = {w0b:+.3f}   w(z=1) = {w0b + wab*0.5:+.3f}")

print(f"\n[2] 1-sigma box (Delta chi2 <= 2.30) on the coarse grid:")
print(f"    w0 in [{w0g[in1.any(axis=1)].min():+.2f}, {w0g[in1.any(axis=1)].max():+.2f}]")
print(f"    wa in [{wag[in1.any(axis=0)].min():+.2f}, {wag[in1.any(axis=0)].max():+.2f}]")

print(f"""
[3] PRE-REGISTERED PASS/FAIL for the Phase 2 derivation (zero knobs):
    PASS  : derived w(a) lands inside the 1-sigma box above AND
            chi2 <= {c2b:.1f} + 2.3 on these exact bins
    SIGN  : w(z=0) > -1 with w rising toward z=0 (thawing), else FAIL
    NULL  : derivation returns exactly w = -1 -> Path Two closes,
            Path One stands as final (still publishable)
    CAVEAT: fixed rd, H0 (C4) -- ordering may shift under marginalization
""")
