import numpy as np

G      = 6.67430e-11
MSUN   = 1.98892e30
A0_EST = 1.192e-10          # ESTIF self-derived
A0_MND = 1.20e-10           # canonical MOND

# Mistele+2024 (ApJL 969 L3) Table 2: (log10 Mb [Msun], Vflat [km/s], sig_stat [km/s])
lens = {
 "All  R<300 ": [(10.10,137.3,11.5),(10.66,182.1,10.8),(10.96,204.7,8.2),(11.29,260.2,6.2)],
 "All  R<1000": [(10.11,135.8, 9.6),(10.66,175.3, 9.3),(10.96,211.0,6.6),(11.29,259.4,5.1)],
 "LTG  R<1000": [(10.08,131.8,10.4),(10.66,169.4,12.1),(10.95,197.5,10.1),(11.23,238.9,12.8)],
 "ETG  R<1000": [(10.68,199.4,14.5),(10.97,227.2, 8.7),(11.30,264.7,5.6)],
}
kin = [(8.69,54.2,1.0),(9.25,77.9,1.1),(9.67,89.1,1.3),(10.00,121.7,1.8),
       (10.46,160.7,1.7),(10.72,186.6,3.3),(11.03,219.3,1.8),(11.34,271.0,6.4)]

SYS_DEX = 0.10              # paper's stated M*/L systematic

def vpred(logM, a0):
    return (G * 10**logM * MSUN * a0)**0.25 / 1000.0

def report(name, pts, a0):
    print(f"\n{name}   [a0 = {a0:.4g} m/s^2]")
    print(f"{'logMb':>7} {'V_obs':>8} {'V_pred':>8} {'diff%':>8} {'n_sig':>7}")
    chi2 = chi2s = 0.0
    for logM, v, s in pts:
        vp = vpred(logM, a0)
        # systematic on Mb -> V uncertainty: dV/V = (1/4)*ln10*dlogM
        s_sys = v * 0.25 * np.log(10) * SYS_DEX
        s_tot = np.hypot(s, s_sys)
        chi2  += ((v-vp)/s)**2
        chi2s += ((v-vp)/s_tot)**2
        print(f"{logM:7.2f} {v:8.1f} {vp:8.1f} {100*(vp-v)/v:8.2f} {(v-vp)/s_tot:7.2f}")
    N = len(pts)
    print(f"  chi2/N (stat only)     = {chi2/N:.3f}")
    print(f"  chi2/N (stat + 0.1dex) = {chi2s/N:.3f}")
    return chi2s/N

def best_a0(pts):
    # deep-MOND: V^4 = G Mb a0  ->  a0_i = V^4/(G Mb); inverse-variance weight in log
    la, lw = [], []
    for logM, v, s in pts:
        a = (v*1000)**4 / (G * 10**logM * MSUN)
        rel = 4*s/v                      # frac err on a0 from stat
        rel = np.hypot(rel, np.log(10)*SYS_DEX)
        la.append(np.log(a)); lw.append(1/rel**2)
    la, lw = np.array(la), np.array(lw)
    m = np.sum(lw*la)/np.sum(lw)
    e = 1/np.sqrt(np.sum(lw))
    return np.exp(m), np.exp(m)*e

print("="*64)
print("WEAK-LENSING BTFR vs ESTIF-DERIVED a0   (zero free parameters)")
print("Data: Mistele, McGaugh, Lelli, Schombert, Li 2024, ApJL 969 L3")
print("="*64)

for k, v in lens.items():
    report(k, v, A0_EST)

print("\n" + "-"*64)
print("BASELINE: same test on kinematic BTFR (Lelli+2019 binned)")
report("Kinematic  ", kin, A0_EST)

print("\n" + "-"*64)
print("INVERTED: a0 implied by each dataset (deep-MOND BTFR)")
for k, v in list(lens.items()) + [("Kinematic  ", kin)]:
    a, e = best_a0(v)
    ns = (a - A0_EST)/e
    print(f"  {k}:  a0 = {a:.3e} +/- {e:.2e}   ({ns:+.2f} sigma from ESTIF)")

print("\n" + "-"*64)
print("Cross-check: ESTIF a0 vs canonical MOND a0 on primary sample")
report("All  R<1000 [MOND]", lens["All  R<1000"], A0_MND)
