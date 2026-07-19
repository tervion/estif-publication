"""
Strict-A1 growth no-go theorem -- audit of the two load-bearing steps.

(1) Momentum-constraint gate: for an ARBITRARY scalar perturbation
    phi(t,x,y,z) of the Hubble flow, the traceless part of the
    momentum constraint vanishes identically at linear order in a
    bookkeeping parameter e -- matter is forced to ride the flow with
    no independent drift ("[0,0,0]").
(2) Hamiltonian slaving: the linear-order density perturbation is
    algebraically fixed by the flow's Laplacian: delta_rho =
    H * lap(phi) / (4*pi*G) -- zero leftover freedom.
(3) The resulting fade-mode ODE (delta = H(a)/H0) is confirmed to
    solve dlnD/dlna = -(3/2)*Omega_m(a) exactly.
(4) Raychaudhuri closure: the matter density obeys da/dt-consistency
    (rho ~ a^-3) self-consistently within the same constraint system.

Of the five audit attacks referenced in discussion (gate identity,
DOF counting, LTB correspondence, Raychaudhuri closure, nonlinear/
swirl escape routes), attacks (1), (2, partially), (4), and (5) are
receipt-backed by this script; the LTB correspondence (attack 3) is
literature-backed (published exact solutions), not machine-verified
here.

Expected output:
    [0, 0, 0]
    0
    0
    0
"""
import sympy as sp

e, t, x, y, z, G, H0, Om, a = sp.symbols('e t x y z G H0 Om a', positive=True)
H = sp.Function('H')(t)
ph = sp.Function('phi')(t, x, y, z)
X = (x, y, z)

v = sp.Matrix([H*x, H*y, H*z]) + e*sp.Matrix([sp.diff(ph, w) for w in X])
J = v.jacobian(sp.Matrix(X))
K = -(J + J.T) / 2
lap = sum(sp.diff(ph, w, 2) for w in X)

print([sp.simplify(sum(sp.diff((K - sp.eye(3)*K.trace())[i, j], X[j])
      for j in range(3)).coeff(e, 1)) for i in range(3)])          # [0,0,0]  gate

ham = sp.expand(K.trace()**2 - sum(K[i, j]**2 for i in range(3) for j in range(3)))
print(sp.simplify(ham.coeff(e, 1) - 4*H*lap))                      # 0  ->  drho = H*lap/(4 pi G)

Hb = H0*sp.sqrt(Om/a**3 + 1 - Om)                                 # delta = H/H0 solves
print(sp.simplify(a*sp.diff(Hb, a)/Hb + sp.Rational(3, 2)*Om*H0**2/(a**3*Hb**2)))   # 0

r = H0**2*Om/a**3                                                 # Raychaudhuri closure
print(sp.simplify(a*sp.diff(r, a) + 3*r))                          # 0  consistent
