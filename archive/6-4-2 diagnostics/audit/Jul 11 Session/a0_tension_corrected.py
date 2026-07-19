import numpy as np
G, MSUN, A0, SYS = 6.67430e-11, 1.98892e30, 1.192e-10, 0.10  # SYS: M*/L dex, CORRELATED

sets = {
 "Kinematic (Lelli+19)": [(8.69,54.2,1.0),(9.25,77.9,1.1),(9.67,89.1,1.3),(10.00,121.7,1.8),
                          (10.46,160.7,1.7),(10.72,186.6,3.3),(11.03,219.3,1.8),(11.34,271.0,6.4)],
 "Lensing All R<1000  ": [(10.11,135.8,9.6),(10.66,175.3,9.3),(10.96,211.0,6.6),(11.29,259.4,5.1)],
 "Lensing LTG R<1000  ": [(10.08,131.8,10.4),(10.66,169.4,12.1),(10.95,197.5,10.1),(11.23,238.9,12.8)],
 "Lensing ETG R<1000  ": [(10.68,199.4,14.5),(10.97,227.2,8.7),(11.30,264.7,5.6)],
}

print("a0 inversion, M*/L systematic FULLY CORRELATED across bins")
print(f"{'sample':<22} {'a0_implied':>11} {'stat_dex':>9} {'sys_dex':>8} {'tension':>9}")
for name, pts in sets.items():
    l = np.array([np.log10((v*1e3)**4/(G*10**m*MSUN)) for m,v,e in pts])
    w = np.array([(4*(e/v)/np.log(10))**-2 for m,v,e in pts])
    mean = np.sum(w*l)/np.sum(w)
    stat = 1/np.sqrt(np.sum(w))
    sig  = (mean - np.log10(A0))/np.hypot(stat, SYS)
    print(f"{name:<22} {10**mean:11.3e} {stat:9.4f} {SYS:8.2f} {sig:+8.2f}s")
print("\n(previous run: same offsets scored +2.67s / +2.12s -- error was")
print(" dividing the correlated systematic by sqrt(N) via per-bin weighting)")
