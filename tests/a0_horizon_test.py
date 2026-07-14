#!/usr/bin/env python3
"""
a0_horizon_test.py
------------------
Test whether the MOND acceleration scale a0 can be built from the
de Sitter / Lambda horizon (sqrt(Lambda) alone, matter genuinely gone)
vs. from the Hubble horizon (H0, which still hides matter via Friedmann).

Key relations
-------------
Friedmann (flat):        H0^2 = (8 pi G / 3) rho_total = (8 pi G/3) rho_m + Lambda c^2/3
de Sitter limit:         H_Lambda = c sqrt(Lambda/3)  =  H0 sqrt(Omega_Lambda)
                         => the Lambda-horizon rate is LOWER than H0 by sqrt(Omega_Lambda)
Omega_Lambda def:        Omega_Lambda = Lambda c^2 / (3 H0^2)

Observed:                a0 ~ 1.2e-10 m/s^2  (McGaugh/Lelli RAR/BTFR; lit. ~1.2-1.24e-10)

Author: (for Peter / ESTIF)  --  reproducible receipt for the P-2 sqrt(Lambda) test
"""

import math

# --- fundamental constants ---
c   = 2.99792458e8              # m/s
Mpc = 3.0856775814913673e22     # m

# --- observed MOND acceleration scale ---
a0_obs = 1.2e-10               # m/s^2

# --- cosmological parameter cases (H0 in km/s/Mpc, Omega_Lambda) ---
cases = {
    "Planck (H0=67.4, OmL=0.685)": (67.4, 0.685),
    "Local  (H0=73.0, OmL=0.685)": (73.0, 0.685),
}

print(f"Observed a0 = {a0_obs:.3e} m/s^2\n")
print(f"{'case':<32}{'cH0/2pi':>12}{'cH_L/2pi':>12}{'sqrtOmL':>9}")
print("-" * 65)

for name, (H0_kms, OmL) in cases.items():
    H0  = H0_kms * 1000 / Mpc          # s^-1
    H_L = H0 * math.sqrt(OmL)          # de Sitter (Lambda) horizon rate = c sqrt(Lambda/3)
    a0_H = c * H0  / (2 * math.pi)      # Hubble-horizon version
    a0_L = c * H_L / (2 * math.pi)      # Lambda-horizon version (matter genuinely gone)
    print(f"{name:<32}{a0_H:>12.3e}{a0_L:>12.3e}{math.sqrt(OmL):>9.3f}")

print()

# --- prefactor family, using Planck values ---
H0  = 67.4 * 1000 / Mpc
OmL = 0.685
H_L = H0 * math.sqrt(OmL)
Lambda = 3 * (H_L / c) ** 2             # m^-2   (standard value ~1.1e-52)
print(f"Lambda = {Lambda:.3e} m^-2   (standard value ~1.1e-52)\n")

forms = {
    "c^2 sqrt(Lambda/3) / 2pi"      : c**2 * math.sqrt(Lambda / 3) / (2 * math.pi),
    "c^2 sqrt(Lambda)   / 2pi"      : c**2 * math.sqrt(Lambda)     / (2 * math.pi),
    "c^2 sqrt(Lambda) * sqrt3 / 2pi": c**2 * math.sqrt(Lambda) * math.sqrt(3) / (2 * math.pi),
}
print("Lambda-only forms vs observed (ratio to 1.2e-10):")
for f, v in forms.items():
    print(f"  {f:<34} = {v:.3e}   ({v/a0_obs:.2f}x)")

# --- what single prefactor k in a0 = k * c^2 sqrt(Lambda) hits the data exactly? ---
k = a0_obs / (c**2 * math.sqrt(Lambda))
print(f"\nExact-match prefactor: a0 = {k:.4f} * c^2 sqrt(Lambda)")
print(f"   compare 1/(2pi)        = {1/(2*math.pi):.4f}")
print(f"   compare sqrt3/(2pi)    = {math.sqrt(3)/(2*math.pi):.4f}")
print(f"   ratio k / (1/2pi)      = {k*2*math.pi:.4f}   (~sqrt(OmL)={math.sqrt(OmL):.4f})")
