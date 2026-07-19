"""
Strict-A1 growth no-go theorem -- solution check.
Verifies: under the strict flat-slice law, a*H'/H = -(3/2)*Omega_m(a)
is an exact identity for a dust+Lambda background, which means
delta(a) = H(a)/H0 solves the (first-order, strict-A1) growth law
exactly -- the DECAYING mode only.

Also evaluates the growth rate f(a) = dlnD/dlna at a = 2/3 (z = 0.5)
with Omega_m = 0.3141: f = -0.911, opposite sign from DESI RSD's
measured f(z=0.5) = +0.76.

Status: verifies the solution and the sign flip; the underlying law
itself is derived in estif_growth_nogo_audit.py plus one pencil step
(continuity), stated as such there -- not re-derived here.

Expected output:
    0
    -0.911
"""
import sympy as sp

a, H0, Om, G = sp.symbols('a H0 Om G', positive=True)
H = H0 * sp.sqrt(Om*a**-3 + 1 - Om)            # Friedmann, dust + Lambda
rho = 3*H0**2*Om*a**-3 / (8*sp.pi*G)           # matter density

# derived flat-slice law:  dlnD/dlna = -(3/2)*Omega_m(a)
# claim: D(a) = H(a)/H0 solves it exactly
lhs = sp.simplify(a * sp.diff(H, a) / H)
rhs = sp.simplify(-sp.Rational(3, 2) * (8*sp.pi*G*rho) / (3*H**2))
print(sp.simplify(lhs - rhs))                  # -> 0   (D_ESTIF = H/H0, exact)

f = lhs                                        # growth rate f(a)
print(sp.nsimplify(0) + f.subs([(Om, 0.3141), (a, sp.Rational(2, 3))]).evalf(3))
                                               # -> about -0.91 at z = 0.5
                                               #    DESI RSD measures about +0.76
