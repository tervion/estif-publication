"""
A1' RE-CENSUS receipt -- the Phase-2 census re-executed under axiom A1'.
Written 19 Jul 2026 (v6.4.2) to close the receipt gap identified in
AUDIT_NOTES_JUL12 C-3: RHAC-006 and PHASE2_DECLARATION claimed a
"re-audit under A1'" citing a file that contains no such content.
This file IS that receipt (decision D4, option (a)).

Wording note: the census has FOUR doors (shape, rate, slosh, swirl) per
PHASE2_DECLARATION section 3; the earlier "five sectors" phrase is
amended in the same v6.4.2 round.

What A1' changed: strict A1 held the slice exactly even (no shape
deviations at all). A1' allows matter-SOURCED dents, with unsourced
evenness returning. The re-census question: does the newly opened
shape sector hide a Lambda-printer (a vacuum source of constant
positive energy density, w = -1)?

What this script verifies, in order:
  [1] Door 1 (rate): the three banked rate-dial facts unchanged
      (Schwarzschild clock rate, Friedmann constraint density, flow
      divergence). The lapse-carries-no-stress step remains an
      asserted ADM fact, exactly as in estif_p2_door1_rate_dial.py.
  [2] Door 2 (slosh): the pure-divergence theorem re-executed. It
      quantifies over an ARBITRARY scalar field, so the A1' dent
      field (scalar) is inside its quantification: cosmic average of
      the ledger entry = 0 for the whole class, momentum constraint
      passes identically. Verdict unchanged by A1'.
  [3] Door 3 (swirl): re-executed. Free vacuum tangle still forbidden
      (residual +w/2); sourced-tangle ledger still negative-definite,
      -(A^2+B^2+C^2)/(32*pi*G) < 0. A Lambda-printer needs a POSITIVE
      constant; a debt cannot print Lambda. Verdict unchanged by A1'.
  [4] Shape door, vacuum limit (the genuinely new A1' content):
      in the A1' growth operator L(D), the source term 4*pi*G*rho*D
      vanishes identically as Omega_m -> 0, and the vacuum ODE
      a^2 H0^2 D'' + 3 a H0^2 D' = 0 has general solution
      D = C1 + C2 * a^(-2): a constant relabeling plus a decaying
      branch, NO growing branch, no energy source. "Unsourced
      evenness returns" -- A1's own clause, machine-verified. The
      deepen kernel itself degenerates in vacuum: with Omega_m = 0,
      I'(a) = 1/(a H0)^3 integrates to D = H0*I proportional to
      const - a^(-2)/2 -- a mixture of the same two vacuum branches.
      Growth is strictly matter-sourced.
  [5] Sourced sanity: with Omega_m > 0 the deepen mode D+ = H*I is
      re-verified exact by the defining-derivative substitution
      (method of estif_a1prime_deepen_exact.py). The amendment opened
      exactly the matter-sourced dent, and only that.

VERDICT: four doors under A1' -- rate FORBIDDEN (A2, unchanged),
slosh ZERO (theorem covers A1' dents), swirl vacuum-forbidden /
sourced-negative (unchanged), shape vacuum-inert with matter-sourced
deepen open. NO LAMBDA-PRINTER. The RHAC-004 null stands under A1'.

Not receipted here (unchanged status): the ~1e11 suppression figure
remains an order-of-magnitude estimate; GW speed is C-15 (separate,
estif_C15_gw_sector.py); door-1 ADM lapse fact asserted, not derived.

Expected output ends with:  RE-CENSUS PASS 5/5 -- NULL STANDS UNDER A1'
Exit code 0 iff all five checks pass.
"""
import sys
import sympy as sp

results = []

# ---------------------------------------------------------------------
# [1] Door 1 re-execution (rate dial) -- from estif_p2_door1_rate_dial.py
# ---------------------------------------------------------------------
G, M, H, r = sp.symbols('G M H r', positive=True)
v = sp.sqrt(2*G*M/r)
clock = sp.sqrt(1 - v**2)
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
X = sp.Matrix([x1, x2, x3])
V = H * X
J = V.jacobian(X)
K = -(J + J.T) / 2
trK = K.trace()
KK = sum(K[i, j]**2 for i in range(3) for j in range(3))
c1 = (sp.simplify(clock - sp.sqrt(1 - 2*G*M/r)) == 0)
c2 = (sp.simplify((trK**2 - KK)/(16*sp.pi*G) - 3*H**2/(8*sp.pi*G)) == 0)
c3 = (sp.simplify(J.trace() - 3*H) == 0)
ok1 = c1 and c2 and c3
results.append(("door 1 (rate): banked facts unchanged", ok1))

# ---------------------------------------------------------------------
# [2] Door 2 re-execution (slosh) -- from estif_p2_door2_slosh_divergence.py
#     phi is ARBITRARY: the A1' dent field is inside this quantification.
# ---------------------------------------------------------------------
y1, y2, y3 = sp.symbols('y1 y2 y3', real=True)
Y = (y1, y2, y3)
phi = sp.Function('phi')(*Y)
S = sp.Matrix(3, 3, lambda i, j: sp.diff(phi, Y[i], Y[j]))
lap = S.trace()
resid = (3*H + lap)**2 - (3*H**2 + 2*H*lap +
                          sum(S[i, j]**2 for i in range(3) for j in range(3))) - 6*H**2
F = [4*H*sp.diff(phi, Y[i]) + sp.diff(phi, Y[i])*lap
     - sum(sp.diff(phi, Y[j])*S[i, j] for j in range(3)) for i in range(3)]
door2_zero = sp.simplify(resid - sum(sp.diff(F[i], Y[i]) for i in range(3)))
door2_mom = [sp.simplify(sum(sp.diff(S[i, j] - sp.eye(3)[i, j]*lap, Y[j])
                             for j in range(3))) for i in range(3)]
ok2 = (door2_zero == 0) and all(m == 0 for m in door2_mom)
results.append(("door 2 (slosh): pure-divergence theorem, arbitrary scalar", ok2))

# ---------------------------------------------------------------------
# [3] Door 3 re-execution (swirl) -- from estif_p2_door3_swirl_ledger.py
# ---------------------------------------------------------------------
x, y, z, A, B, C = sp.symbols('x y z A B C', positive=True)
XY = (x, y, z)
w = sp.Matrix([A*sp.sin(z) + C*sp.cos(y),
               B*sp.sin(x) + A*sp.cos(z),
               C*sp.sin(y) + B*sp.cos(x)])
div_w = sp.simplify(sum(sp.diff(w[i], XY[i]) for i in range(3)))
mom = sp.Matrix([-sp.Rational(1, 2)*sum(sp.diff(w[i], vv, 2) for vv in XY)
                 for i in range(3)])
free_tangle = sp.simplify(mom - w/2)
box = lambda e: sp.integrate(e, (x, 0, 2*sp.pi), (y, 0, 2*sp.pi),
                             (z, 0, 2*sp.pi)) / (2*sp.pi)**3
Ssw = sp.Matrix(3, 3, lambda i, j: -(sp.diff(w[i], XY[j]) + sp.diff(w[j], XY[i]))/2)
Kt = -H*sp.eye(3) + Ssw
led = (Kt.trace()**2 - sum(Kt[i, j]**2 for i in range(3) for j in range(3))
       - 6*H**2) / (16*sp.pi*G)
led_avg = sp.simplify(box(led))
ok3 = (div_w == 0
       and all(sp.simplify(free_tangle[i]) == 0 for i in range(3))
       and sp.simplify(led_avg + (A**2 + B**2 + C**2)/(32*sp.pi*G)) == 0)
results.append(("door 3 (swirl): free tangle forbidden; ledger negative", ok3))

# ---------------------------------------------------------------------
# [4] Shape door under A1', vacuum limit (the new content)
# ---------------------------------------------------------------------
a, H0, Om = sp.symbols('a H0 Om', positive=True)
Hbg = H0*sp.sqrt(Om/a**3 + 1 - Om)
rho = 3*H0**2*Om/(8*sp.pi*G*a**3)
D = sp.Function('D')(a)
Lop = (a**2*Hbg**2*sp.diff(D, a, 2)
       + (3*a*Hbg**2 + a**2*Hbg*sp.diff(Hbg, a))*sp.diff(D, a)
       - 4*sp.pi*G*rho*D)
source_vac = sp.simplify((4*sp.pi*G*rho).subs(Om, 0))
Lvac = sp.simplify(Lop.subs(Om, 0))
sol = sp.dsolve(sp.Eq(Lvac, 0), D)
basis_const = sp.simplify(Lvac.subs(sp.Derivative(D, (a, 2)), 0)
                          .subs(sp.Derivative(D, a), 0).subs(D, 1))
Dtest = a**-2
basis_decay = sp.simplify(
    Lvac.subs(sp.Derivative(D, (a, 2)), sp.diff(Dtest, a, 2))
        .subs(sp.Derivative(D, a), sp.diff(Dtest, a))
        .subs(D, Dtest))
I_vac = sp.integrate(1/(a*H0)**3, a)          # deepen kernel at Om = 0
D_vac_deepen = sp.simplify(H0*I_vac)          # -> -1/(2*a**2*H0**2): decaying branch
deepen_degenerate = sp.simplify(D_vac_deepen + 1/(2*a**2*H0**2)) == 0
sol_rhs = sol.rhs
structure_ok = sp.simplify(sp.expand(sol_rhs)
                           - (sol_rhs.coeff(a, 0) + sol_rhs.coeff(a, -2)*a**-2)) == 0
ok4 = (source_vac == 0 and basis_const == 0 and basis_decay == 0
       and deepen_degenerate and structure_ok)
results.append(("shape door, vacuum: source == 0; D = C1 + C2/a^2; no growth", ok4))

# ---------------------------------------------------------------------
# [5] Sourced sanity: deepen mode exact (method of estif_a1prime_deepen_exact.py)
# ---------------------------------------------------------------------
I = sp.Function('I')(a)
Dp = Hbg*I
L5 = (a**2*Hbg**2*sp.diff(Dp, a, 2)
      + (3*a*Hbg**2 + a**2*Hbg*sp.diff(Hbg, a))*sp.diff(Dp, a)
      - 4*sp.pi*G*rho*Dp)
L5 = L5.subs(sp.diff(I, a, 2), sp.diff(1/(a*Hbg)**3, a)).subs(sp.diff(I, a), 1/(a*Hbg)**3)
ok5 = (sp.simplify(L5) == 0)
results.append(("sourced deepen mode D+ = H*I exact (Om > 0)", ok5))

# ---------------------------------------------------------------------
print("A1' RE-CENSUS -- Phase-2 doors re-executed under the amended axiom")
print("-" * 68)
npass = 0
for name, ok in results:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    npass += int(ok)
print("-" * 68)
if npass == len(results):
    print(f"RE-CENSUS PASS {npass}/{len(results)} -- NULL STANDS UNDER A1'")
    sys.exit(0)
print(f"RE-CENSUS FAIL {npass}/{len(results)}")
sys.exit(1)
