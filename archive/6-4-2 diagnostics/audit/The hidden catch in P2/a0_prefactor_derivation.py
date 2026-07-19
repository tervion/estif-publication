#!/usr/bin/env python3
"""
a0_prefactor_derivation.py
--------------------------
Step 1 of the P-2 test: try to DERIVE the a0 prefactor from black-hole-edge
geometry, and see whether it lands in the target band (~0.125) rather than the
bare horizon-thermodynamic value 1/(2pi)=0.159.

Contents
--------
  Part A: how big is the target, really?      -> it's a BAND, not the point "0.128"
  Part B: reference prefactors                -> 1/(2pi) and 1/6 both OVERSHOOT
  Part C: numerology floor                    -> ~1.3% of random O(1) forms hit by luck
  Part D: de Sitter Painleve-Gullstrand river -> honest first-principles attempt

Conclusion: the x_c=0.272 edge overshoots the needed 0.221*c*H_L by ~20%, and the
~20% is the ambiguity between acceleration DEFINITIONS. Decisive test needs the
ESTIF flow law + whether the a0-edge equals the mode-crossover edge.

Depends only on the Python standard library.  Run: python3 a0_prefactor_derivation.py
"""

import math, itertools

c   = 2.99792458e8
Mpc = 3.0856775814913673e22

# ------------------------------------------------------------------ Part A
def k_target(a0, H0_kms, OmL):
    """prefactor k such that a0 = k * c^2 sqrt(Lambda)."""
    H0  = H0_kms * 1000 / Mpc
    Lam = 3 * (H0 * math.sqrt(OmL) / c) ** 2
    return a0 / (c**2 * math.sqrt(Lam))

grid = [k_target(a0, H0, OmL)
        for a0  in (1.20e-10, 1.24e-10)
        for H0  in (67.4, 73.0)
        for OmL in (0.68, 0.69)]
lo, hi = min(grid), max(grid)
print("PART A  target band for k in a0 = k*c^2 sqrt(Lambda)")
print(f"   k_min={lo:.4f}  k_max={hi:.4f}  center~{sum(grid)/len(grid):.4f}")
print(f"   ('0.128' was just the corner a0=1.20,H0=67.4,OmL=0.685)\n")

# ------------------------------------------------------------------ Part B
refs = {
    "1/(2pi)  [Gibbons-Hawking/Unruh default]": 1/(2*math.pi),
    "1/6      [common a0=cH0/6 fit]"          : 1/6,
    "sqrt(OmL)/(2pi)"                          : math.sqrt(0.685)/(2*math.pi),
}
print("PART B  reference prefactors vs band")
for n, v in refs.items():
    tag = "IN" if lo <= v <= hi else ("ABOVE" if v > hi else "below")
    print(f"   {n:<42} = {v:.4f}  [{tag}]")
print()

# ------------------------------------------------------------------ Part C
xc = 0.272
blocks = {'1':1., '2':2., '1/2':.5, 'pi':math.pi, '2pi':2*math.pi,
          'sqrt3':math.sqrt(3), 'xc':xc, 'sqrtOmL':math.sqrt(0.685)}
names, hits, total = list(blocks), 0, 0
for nn in range(0, 3):
    for dd in range(1, 3):
        for num in itertools.combinations_with_replacement(names, nn):
            for den in itertools.combinations_with_replacement(names, dd):
                total += 1
                v = 1.0
                for x in num: v *= blocks[x]
                for x in den: v /= blocks[x]
                if lo <= v <= hi: hits += 1
print("PART C  numerology floor")
print(f"   {hits}/{total} = {100*hits/total:.1f}% of simple O(1) expressions hit the band by luck")
print("   -> an in-band match, by itself, is NOT evidence.\n")

# ------------------------------------------------------------------ Part D
H0 = 67.4*1000/Mpc; OmL = 0.685
HL = H0*math.sqrt(OmL); cHL = c*HL
a0_obs = 1.20e-10
need = a0_obs/cHL
print("PART D  de Sitter Painleve-Gullstrand river: accel at the v=x_c*c edge")
print(f"   c*H_L = {cHL:.3e},  needed fraction a0/cH_L = {need:.4f}")
edge = {
    "flow speed / flow gradient : x_c"          : xc,
    "proper accel (static)      : x_c/sqrt(1-x_c^2)": xc/math.sqrt(1-xc**2),
    "surface gravity (redshift) : x_c*sqrt(1-x_c^2)": xc*math.sqrt(1-xc**2),
    "kinetic 1/2 v^2            : x_c/2"          : xc/2,
}
for n, f in edge.items():
    print(f"   {n:<44} frac={f:.4f}  a0={f*cHL:.3e}  ({f*cHL/a0_obs:.2f}x)")
print(f"\n   inverse: law a0=c*H_L*x_c needs x_c={need:.4f}; you use 0.272 -> 23% high")
print("   => the a0-edge (0.221) and the mode-crossover edge (0.272) may be DIFFERENT surfaces.")
