"""
A1' fork -- growth restored, part 1 (fade mode + numeric growth rate).

Verifies the linear-growth operator L(D) is exactly the standard
second-order growth equation (coefficients check out: 4*pi*G*rho =
(3/2)*Omega_m(a)*H^2), and that D = H(a) -- the strict-A1 result --
is its known DECAYING solution, confirming strict ESTIF is embedded
in A1' rather than contradicted by it.

Then computes the growth rate f(a) = dlnD/dlna numerically for the
DEEPEN mode at a = 2/3 (z = 0.5) using a direct numerical integral
(not a symbolic sympy integral -- see estif_a1prime_deepen_exact.py
for the exact symbolic proof of this mode). Result: f = 0.76,
matching the analytic approximation Omega_m(a)^0.55 = 0.760 and
DESI RSD's measured growth rate.

NOTE: the exact symbolic check of the deepen mode against the
unevaluated integral is deliberately NOT attempted here (doing so
via sympy's .doit() produces a spurious 'nan' from a 0*infinity
special-function edge case -- an artifact of forcing symbolic
evaluation of an integral with no closed form, not a result about
the physics). That exact check is done correctly in
estif_a1prime_deepen_exact.py by substituting the integral's
DEFINING derivative rule instead of evaluating it.

The GW-speed-at-c claim referenced in discussion is a SEPARATE,
open item (C-15) -- an expectation from axiom A2, not something
this script derives. No such check is made here.

Expected output:
    0
    0.76
"""
import sympy as sp
import math

a, x, H0, Om, G = sp.symbols('a x H0 Om G', positive=True)
H = H0*sp.sqrt(Om/a**3 + 1 - Om)
Hx = H0*sp.sqrt(Om/x**3 + 1 - Om)
rho = 3*H0**2*Om/(8*sp.pi*G*a**3)

L = lambda D: (a**2*H**2*sp.diff(D, a, 2) + (3*a*H**2 + a**2*H*sp.diff(H, a))*sp.diff(D, a)
               - 4*sp.pi*G*rho*D)

print(sp.simplify(L(H)))                              # 0  fade mode = strict theory's delta=H

# Numeric growth rate of the deepen mode (direct quadrature, not symbolic)
Omv = 0.3141
Hn = lambda A: math.sqrt(Omv/A**3 + 1 - Omv)


def Iq(A, n=4000):
    return sum(1.0 / ((A*(i - 0.5)/n) * Hn(A*(i - 0.5)/n))**3 for i in range(1, n+1)) * A / n


Dg = lambda A: Hn(A) * Iq(A)
A = 2/3
dA = 1e-5
print(round(A * (Dg(A+dA) - Dg(A-dA)) / (2*dA) / Dg(A), 3))    # ~0.76  vs DESI +0.76 at z=0.5
