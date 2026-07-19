"""
ESTIF T_munu derivation engine: Gauss-Codazzi / ADM symbolic calculator
========================================================================
Purpose
-------
Mechanizes the "pencil-and-paper" hypersurface calculation:
  1. Computes Christoffel symbols, Ricci tensor, Ricci scalar, and the
     Einstein tensor for an arbitrary metric (any dimension, symbolic).
  2. Computes ADM quantities (extrinsic curvature K_ij, its trace, the
     Hamiltonian constraint) for a t-foliation with zero shift.
  3. SELF-VALIDATES against two exact GR results before touching ESTIF:
       [V1] Flat FRW  -> Friedmann equation + pressure equation
       [V2] Global de Sitter (hyperboloid slicing) -> G_uv = -Lambda g_uv
       [V3] ADM Hamiltonian constraint == 2 * G_tt / N^2 (both cases)
  4. Runs an ESTIF-SHAPED ILLUSTRATIVE ANSATZ (clearly labeled, authored
     by Claude, NOT an ESTIF v6.2 result): a 3D hypersurface carried at
     bulk speed dW/dt through a 4D bulk whose spatial scale varies along
     the flow direction w via a profile b(w). The induced 4D metric is

         ds^2 = eps * (dW/dt)^2 dt^2 + b(W(t))^2 [dchi^2 + chi^2 dOmega_2^2]

     For eps = -1 (Lorentzian bulk direction) the engine extracts the
     effective energy density rho_eff and pressure p_eff that the
     Hamiltonian constraint and Einstein tensor FORCE, with no further
     postulate. For eps = +1 (Euclidean bulk) it reports the signature
     obstruction instead of faking an energy density.

What this script CANNOT do
--------------------------
It cannot choose ESTIF's physics. Two objects must come from the paper:
  (a) the signature mechanism (why the flow direction contributes -dt^2),
  (b) the flow law: W(t) and the bulk profile b(w), i.e. the precise
      mathematical statement of "expansion = 4D inward fall".
Insert them in the USER INPUT block below. Everything printed under
section [3] uses placeholder forms and is illustrative only. Peter's
publication rule applies: no illustrative output goes in the paper.

Run:  python3 estif_tmunu_gauss_codazzi.py
Deps: sympy (tested on 1.14.0)
"""

import sympy as sp

# ----------------------------------------------------------------------
# USER INPUT BLOCK  (replace with ESTIF v6.2 definitions before drawing
# any physical conclusion; defaults reproduce the illustrative run)
# ----------------------------------------------------------------------
EPSILON_CHOICE = -1          # -1 Lorentzian flow direction, +1 Euclidean
FLOW_AT_SPEED_C = True       # True: dW/dt = c (uniform inward fall)
                             # False: keep W(t) fully general/symbolic
# The bulk spatial profile b(w) stays a general symbolic function here.
# When ESTIF specifies b(w) explicitly, define it in estif_b_of_w below.
def estif_b_of_w(w_symbol):
    return None              # None means: keep b general (symbolic)
# ----------------------------------------------------------------------

t, chi, th, ph = sp.symbols('t chi theta phi', positive=True)
xx, yy, zz = sp.symbols('x y z')
L, c, G, w = sp.symbols('L c G w', positive=True)
eps = sp.Symbol('epsilon', real=True)


def christoffel(g, x):
    n = len(x)
    gi = g.inv()
    Gm = [[[sp.S(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for cc in range(b, n):
                expr = sum(
                    gi[a, d] * (sp.diff(g[d, b], x[cc])
                                + sp.diff(g[d, cc], x[b])
                                - sp.diff(g[b, cc], x[d]))
                    for d in range(n)) / 2
                expr = sp.simplify(expr)
                Gm[a][b][cc] = expr
                Gm[a][cc][b] = expr
    return Gm


def ricci_tensor(g, x):
    n = len(x)
    Gm = christoffel(g, x)
    Ric = sp.zeros(n, n)
    for b in range(n):
        for cc in range(b, n):
            e = sp.S(0)
            for a in range(n):
                e += sp.diff(Gm[a][b][cc], x[a]) - sp.diff(Gm[a][b][a], x[cc])
                for d in range(n):
                    e += Gm[a][a][d] * Gm[d][b][cc] - Gm[a][cc][d] * Gm[d][b][a]
            e = sp.simplify(e)
            Ric[b, cc] = e
            Ric[cc, b] = e
    return Ric


def einstein_tensor(g, x):
    Ric = ricci_tensor(g, x)
    gi = g.inv()
    n = len(x)
    Rs = sp.simplify(sum(gi[i, j] * Ric[i, j] for i in range(n) for j in range(n)))
    Ein = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            e = sp.simplify(Ric[i, j] - Rs * g[i, j] / 2)
            Ein[i, j] = e
            Ein[j, i] = e
    return Ein, Rs, Ric


def adm_hamiltonian(h3, x3, N):
    """Hamiltonian constraint scalar R3 + K^2 - K_ij K^ij for a
    t-foliation with lapse N, zero shift, 3-metric h3(t, x3)."""
    hi = h3.inv()
    K = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            K[i, j] = sp.diff(h3[i, j], t) / (2 * N)
    Ktr = sp.simplify(sum(hi[i, j] * K[j, i] for i in range(3) for j in range(3)))
    Kup = hi * K * hi
    KK = sp.simplify(sum(Kup[i, j] * K[i, j] for i in range(3) for j in range(3)))
    _, R3, _ = einstein_tensor(h3, x3)
    return sp.simplify(R3 + Ktr**2 - KK), R3, Ktr


def ssimp(e):
    return sp.simplify(sp.trigsimp(sp.simplify(sp.expand(e))))


def check(label, expr_should_be_zero):
    ok = ssimp(expr_should_be_zero) == 0
    print(('PASS  ' if ok else 'FAIL  ') + label)
    if not ok:
        print('       residual:', ssimp(expr_should_be_zero))
    return ok


print('=' * 72)
print('[1] ENGINE SELF-TEST: flat FRW -> Friedmann')
print('=' * 72)
a = sp.Function('a', positive=True)(t)
gF = sp.diag(-1, a**2, a**2, a**2)
xF = [t, xx, yy, zz]
EinF, RsF, _ = einstein_tensor(gF, xF)
rho_F = sp.simplify(EinF[0, 0] / (8 * sp.pi * G))
p_F = sp.simplify(EinF[1, 1] / (8 * sp.pi * G * a**2))
H = sp.diff(a, t) / a
ok1 = check('Friedmann:  rho = 3H^2 / 8piG',
            rho_F - 3 * H**2 / (8 * sp.pi * G))
ok2 = check('Pressure:   p = -(2 addot/a + H^2) / 8piG',
            p_F + (2 * sp.diff(a, t, 2) / a + H**2) / (8 * sp.pi * G))
h3F = sp.diag(a**2, a**2, a**2)
HamF, R3F, KtrF = adm_hamiltonian(h3F, [xx, yy, zz], sp.S(1))
ok3 = check('ADM identity:  R3 + K^2 - K_ij K^ij = 2 G_tt / N^2',
            HamF - 2 * EinF[0, 0])
print('    computed rho_eff  =', rho_F)
print('    computed p_eff    =', p_F)

print()
print('=' * 72)
print('[2] ENGINE SELF-TEST: global de Sitter (hyperboloid slicing)')
print('    embedding-derived metric; curved 3-sphere slices')
print('=' * 72)
adS = L * sp.cosh(t / L)
gdS = sp.diag(-1,
              adS**2,
              adS**2 * sp.sin(chi)**2,
              adS**2 * sp.sin(chi)**2 * sp.sin(th)**2)
xdS = [t, chi, th, ph]
EindS, RsdS, _ = einstein_tensor(gdS, xdS)
Lam = 3 / L**2
ok4 = all(
    ssimp((EindS[i, j] + Lam * gdS[i, j]).rewrite(sp.exp)) == 0
    for i in range(4) for j in range(4))
print(('PASS  ' if ok4 else 'FAIL  ') +
      'G_uv = -Lambda g_uv with Lambda = 3/L^2  (pure geometry -> vacuum energy)')
rho_dS = sp.simplify(EindS[0, 0] / (8 * sp.pi * G))
print('    computed rho_eff  =', ssimp(rho_dS.rewrite(sp.exp)),
      '   expected 3/(8 pi G L^2)')
h3dS = sp.diag(adS**2,
               adS**2 * sp.sin(chi)**2,
               adS**2 * sp.sin(chi)**2 * sp.sin(th)**2)
HamdS, R3dS, KtrdS = adm_hamiltonian(h3dS, [chi, th, ph], sp.S(1))
ok5 = check('ADM identity on curved slices:  Ham = 2 G_tt',
            (HamdS - 2 * EindS[0, 0]).rewrite(sp.exp))
print('    intrinsic 3-curvature of slice  R3 =', ssimp(R3dS.rewrite(sp.exp)))

print()
print('=' * 72)
print('[3] ESTIF-SHAPED ANSATZ  --  ILLUSTRATIVE ONLY (Claude ansatz,')
print('    NOT an ESTIF v6.2 result; replace W, b, eps with paper defs)')
print('=' * 72)
b = sp.Function('b', positive=True)

if FLOW_AT_SPEED_C:
    Wt = c * t
    Wdot = c
    print('flow law   : W(t) = c t          (uniform inward fall at c)')
else:
    Wfun = sp.Function('W')(t)
    Wt = Wfun
    Wdot = sp.diff(Wfun, t)
    print('flow law   : W(t) general symbolic')

user_b = estif_b_of_w(w)
if user_b is not None:
    bW = user_b.subs(w, Wt)
    print('bulk scale : b(w) =', user_b, '  (user supplied)')
else:
    bW = b(Wt)
    print('bulk scale : b(w) general symbolic  (ESTIF must supply)')

print('signature  : eps =', EPSILON_CHOICE)
print()

if EPSILON_CHOICE == +1:
    print('SIGNATURE OBSTRUCTION (structural finding, not an error):')
    print('  g_tt = + (dW/dt)^2 > 0  ->  induced metric is (+,+,+,+).')
    print('  Riemannian worldvolume: no timelike direction, no light cones,')
    print('  no causal structure; calling G_tt/8piG an energy density is')
    print('  unjustified. ESTIF must specify the mechanism that produces')
    print('  Lorentzian signature from the flow (this is input (a)).')
else:
    gE = sp.diag(-Wdot**2,
                 bW**2,
                 bW**2 * chi**2,
                 bW**2 * chi**2 * sp.sin(th)**2)
    xE = [t, chi, th, ph]
    EinE, RsE, _ = einstein_tensor(gE, xE)

    N2 = Wdot**2
    rho_eff = sp.simplify(EinE[0, 0] / (8 * sp.pi * G * N2))
    p_eff = sp.simplify(EinE[1, 1] / (8 * sp.pi * G * bW**2))

    offdiag_ok = all(
        ssimp(EinE[i, j]) == 0 for i in range(4) for j in range(4) if i != j)
    print(('PASS  ' if offdiag_ok else 'FAIL  ') +
          'momentum sector: all off-diagonal G_uv = 0 (homogeneous flow)')

    h3E = sp.diag(bW**2, bW**2 * chi**2, bW**2 * chi**2 * sp.sin(th)**2)
    HamE, R3E, KtrE = adm_hamiltonian(h3E, [chi, th, ph], sp.sqrt(N2))
    okE = check('ADM identity:  Ham = 2 G_tt / N^2',
                HamE - 2 * EinE[0, 0] / N2)

    print()
    print('FORCED effective sources (normal-observer frame), no postulates:')
    print('  rho_eff =', rho_eff)
    print('  p_eff   =', p_eff)
    print()
    print('Friedmann-analog read-off:')
    print('  proper time      d tau = |dW/dt| dt / c')
    aa = bW
    Hphys = sp.simplify(sp.diff(aa, t) / (aa * sp.sqrt(N2)) * c)
    print('  H (proper time)  =', Hphys)
    print('  i.e.  H = (flow speed) x (d ln b / dw) : expansion rate equals')
    print('        the bulk log-gradient traversed per unit proper time.')
    print()
    print('Consistency: 8 pi G rho_eff - 3 (H/c)^2  =',
          ssimp(8 * sp.pi * G * rho_eff - 3 * (Hphys / c)**2))

print()
print('=' * 72)
print('[4] WHAT MUST COME FROM ESTIF v6.2 (the physics the engine cannot')
print('    supply; insert in the USER INPUT block and rerun)')
print('=' * 72)
print(' (a) signature mechanism: justification that the flow direction')
print('     contributes with eps = -1 (Lorentzian). If v6.2 does not state')
print('     one explicitly, that gap is itself a finding of this exercise.')
print(' (b) flow law: is dW/dt = c exact (inward fall at c), or W(t) with')
print('     its own dynamics? One line of the paper should fix this.')
print(' (c) bulk profile b(w): the precise geometric statement of the')
print('     inward-fall law. This is where x enters. Given b(w), the')
print('     Hamiltonian constraint above yields H(w), hence H(z) and x(z)')
print('     as the solution of an ODE  --  dynamically determined, not')
print('     self-referential. That is the de-circularization path.')
print(' (d) eddies (gravity sector): promote W(t) -> W(t, r) and/or')
print('     b(w) -> b(w, r); the same engine computes the perturbed')
print('     constraints unchanged. Background first, eddies second.')
print()
print('All numeric/symbolic values printed above are actual runtime output')
print('of the validated engine. Section [3] uses illustrative ansatz forms')
print('and must not be quoted in the publication until (a)-(c) are replaced')
print('by ESTIF v6.2 definitions.')
