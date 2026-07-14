"""
Phase 2 -- Door 1 (passage-rate wobble) receipt.
Verdict: FORBIDDEN by axiom A2 (universal speed c); the ADM lapse
("rate dial") carries no stress on its own -- standard ADM fact,
asserted here, not derived by this script. What this script DOES
verify: the two banked Path One results (Schwarzschild clock rate,
Friedmann constraint) both live at rate == 1, so closing Door 1
costs zero banked physics.

Status: consistency check, not a derivation.
Expected output:
    static clock rate : sqrt(-2*G*M/r + 1)
    constraint density: 3*H**2/(8*pi*G)
    flow divergence   : 3*H
"""
import sympy as sp

G, M, H, r = sp.symbols('G M H r', positive=True)

# (a) PG-Schwarzschild: rate dial = 1, flow v^2 = 2GM/r (July engine result)
v = sp.sqrt(2*G*M/r)
print("static clock rate :", sp.sqrt(1 - v**2))      # -> sqrt(1 - 2*G*M/r)

# (b) Hubble flow v^i = H x^i: rate dial = 1, flat slice
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
X = sp.Matrix([x1, x2, x3])
V = H * X
J = V.jacobian(X)
K = -(J + J.T) / 2                       # extrinsic curvature at rate = 1
trK = K.trace()
KK = sum(K[i, j]**2 for i in range(3) for j in range(3))
print("constraint density:", sp.simplify((trK**2 - KK) / (16*sp.pi*G)))  # -> 3*H**2/(8*pi*G)
print("flow divergence   :", J.trace())              # -> 3*H
