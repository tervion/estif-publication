import sympy as sp

a,H0,Om,G = sp.symbols('a H0 Om G', positive=True)

H   = H0*sp.sqrt(Om/a**3 + 1 - Om)

rho = 3*H0**2*Om/(8*sp.pi*G*a**3)

I   = sp.Function('I')(a)                       # the integral, unevaluated

D   = H*I

L   = (a**2*H**2*sp.diff(D,a,2)

     + (3*a*H**2 + a**2*H*sp.diff(H,a))*sp.diff(D,a)

     - 4*sp.pi*G*rho*D)

L   = L.subs(sp.diff(I,a,2), sp.diff(1/(a*H)**3, a)).subs(sp.diff(I,a), 1/(a*H)**3)

print(sp.simplify(L))                           # -> 0  deepen mode, exact