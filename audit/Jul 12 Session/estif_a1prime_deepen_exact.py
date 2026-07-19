"""
A1' fork -- growth restored, part 2 (exact symbolic proof of the deepen mode).

This is the authoritative receipt for D+ = H(a) * Integral[da/(a*H)^3],
the growing-mode solution restored by the A1' amendment. Rather than
letting sympy evaluate the (closed-form-free) integral -- which
produces a spurious 'nan' -- this script substitutes the integral's
DEFINING derivative rule (I'(a) = 1/(a*H)^3) directly into the growth
operator L(D) from estif_a1prime_growth_restored.py. Every term
involving the integral's *value* then cancels by construction; only
terms involving its *derivative* remain, and those cancel exactly
against the background Friedmann identity.

This exact proof independently corroborates the numeric f = 0.76
result from estif_a1prime_growth_restored.py.

Expected output:
    0
"""
import sympy as sp

a, H0, Om, G = sp.symbols('a H0 Om G', positive=True)
H = H0*sp.sqrt(Om/a**3 + 1 - Om)
rho = 3*H0**2*Om/(8*sp.pi*G*a**3)
I = sp.Function('I')(a)                       # the integral, unevaluated
D = H*I
L = (a**2*H**2*sp.diff(D, a, 2)
     + (3*a*H**2 + a**2*H*sp.diff(H, a))*sp.diff(D, a)
     - 4*sp.pi*G*rho*D)
L = L.subs(sp.diff(I, a, 2), sp.diff(1/(a*H)**3, a)).subs(sp.diff(I, a), 1/(a*H)**3)
print(sp.simplify(L))                           # -> 0  deepen mode, exact
