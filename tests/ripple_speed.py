#!/usr/bin/env python3
"""
ripple_speed.py — does your flow's wobble travel at c?

Given an equation of motion (EoM) for a field it: (1) splits the field into
background + tiny ripple  phi -> phi0 + eps*delta;  (2) keeps only the O(eps)
linear part;  (3) plugs in a plane wave delta ~ exp(i(k x - omega t));
(4) solves the dispersion relation omega(k);  (5) reports phase speed omega/k,
group speed d(omega)/dk, and a verdict:
    - non-dispersive (all ripples one speed)?   <- what a real wave medium does
    - is that speed exactly c ?                 <- grounds "stiffness = c^2"

You edit ONE function, estif_eom, marked below. The rest is machinery, validated
on the ordinary wave equation and the massive (Klein-Gordon) equation.

1+1 D for clarity. In 3+1 D read d^2/dx^2 as the Laplacian and k as |k|; identical.
"""
import sympy as sp

t, x  = sp.symbols('t x', real=True)
w, k  = sp.symbols('omega k', positive=True)
c, cs = sp.symbols('c c_s', positive=True)
v0, m = sp.symbols('v_0 m', positive=True)
eps   = sp.symbols('epsilon', positive=True)

def d(u, var, n=1): return sp.diff(u, var, n)

# advection OPERATOR (d_t + v0 d_x) applied once — compose it for higher powers.
# (writing (d(u,t)+v0*d(u,x))**2 would be the ALGEBRAIC square = a nonlinear term,
#  whose linear part is zero. Operators must be COMPOSED, not squared.)
def adv(u):        return d(u,t) + v0*d(u,x)
def adv2(u):       return adv(adv(u))

# =====================================================================
#  MACHINERY (validated below; you should not need to touch this)
# =====================================================================
def analyze(eom, phi0, label):
    delta = sp.exp(sp.I*(k*x - w*t))
    lin   = sp.diff(eom(phi0 + eps*delta), eps).subs(eps, 0)   # O(eps) part
    char  = sp.simplify(sp.expand(lin/delta))                  # dispersion polynomial
    sols  = sp.solve(sp.Eq(char, 0), w)
    print("="*66); print(label); print("="*66)
    print(f"  dispersion:  {sp.nsimplify(char)} = 0")
    for s in sols:
        s   = sp.simplify(s)
        vph = sp.simplify(s/k); vgr = sp.simplify(sp.diff(s, k))
        nondisp = (sp.simplify(vph - vgr) == 0) and (k not in vph.free_symbols)
        print(f"  root omega = {s}")
        print(f"       phase speed  = {vph}")
        print(f"       group speed  = {vgr}")
        print(f"       {'NON-dispersive (one speed for all ripples)' if nondisp else 'DISPERSIVE (speed depends on k)'}")
    print()
    return sols

# ---- VALIDATION 1: wave equation  u_tt = c^2 u_xx  -> omega=c k, speed c ----
analyze(lambda u: d(u,t,2) - c**2*d(u,x,2), sp.Integer(0),
        "VALIDATION 1  wave equation  u_tt = c^2 u_xx")

# ---- VALIDATION 2: Klein-Gordon  u_tt = c^2 u_xx - m^2 u  -> dispersive ----
analyze(lambda u: d(u,t,2) - c**2*d(u,x,2) + m**2*u, sp.Integer(0),
        "VALIDATION 2  Klein-Gordon  u_tt = c^2 u_xx - m^2 u")

# ---- FLOW TEMPLATE: ripple on a moving medium (acoustic / Unruh) ----
#      (d_t + v0 d_x)^2 u = c_s^2 u_xx   ->   omega = v0 k +/- c_s k
#      i.e. ripple travels at c_s RELATIVE to the flow. ESTIF question: is c_s = c?
analyze(lambda u: adv2(u) - cs**2*d(u,x,2), sp.Integer(0),
        "FLOW TEMPLATE  (d_t + v0 d_x)^2 u = c_s^2 u_xx   [acoustic/Unruh]")

# =====================================================================
#  >>> ESTIF: PLUG YOUR FLOW'S EQUATION OF MOTION IN HERE <<<
#  Replace the body of estif_eom(u) with your actual EoM = 0, as a function of u
#  and its derivatives d(u,t), d(u,x), d(u,x,2), and/or the operators adv/adv2.
#  Set phi0 to your background (0, a constant, or a flow profile like v0*x).
#
#  READING THE OUTPUT:
#    non-dispersive AND phase speed = c  -> ripples travel at c: stiffness = c^2 is
#      grounded FROM your geometry, and the c*H acceleration floor follows. (These
#      ripples ARE gravitational waves — the same object as your LISA prediction.)
#    dispersive, or speed != c           -> mechanism not grounded; learned cheaply,
#      before any galaxy N-body.
#  Even a clean speed = c does NOT fix the O(1) prefactor (cH vs cH/2pi vs ~0.128).
#  That stays the separate, still-open problem.
# =====================================================================
def estif_eom(u):
    # ---- EDIT THIS (placeholder = wave eq so the file runs as-is) ----
    # A likely-relevant starting structure is the FLOW TEMPLATE above:
    #     return adv2(u) - cs**2*d(u,x,2)
    return d(u,t,2) - c**2*d(u,x,2)

analyze(estif_eom, sp.Integer(0), "ESTIF  (your EoM — edit estif_eom above)")
