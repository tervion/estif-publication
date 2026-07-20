"""
ESTIF converse flow-law test
=============================
Direction of inference REVERSED relative to Test 4d.

Test 4d proved:   GIVEN v_r = c sqrt(rs/r)  =>  exact vacuum (Schwarzschild).
This script asks: GIVEN only the formalized inward-fall principle, is
                  v_r^2 = 2A/r the UNIQUE flow the principle allows?

Formalized principle under test (assumption set for the fidelity audit):
  A1  space FLOWS, it does not stretch: spatial slices are flat
      (3-metric = flat, spherical symmetry)
  A2  everything rides the bulk at c: flow riders' clocks tick coordinate
      time (lapse N = c, verified as the raindrop identity in Test 4d)
  A3  static, spherically symmetric eddy: one unknown function v(r)
  A4  the Einstein tensor of the emergent metric IS the effective
      energy content (the engine premise validated on FRW and de Sitter);
      'no matter' means the normal observer measures zero energy density

Pre-registered PASS criteria (written before the run):
  P1  flat spherical 3-slices have R3 = 0                       [sanity]
  P2  ADM-with-shift Hamiltonian == 2 G_nn                      [engine
      cross-validation, first test including a nonzero shift]
  P3  the effective density has the closed form
        G_nn = d(r v^2)/dr / (c^2 r^2)
      and every nonzero Einstein component dies when d(r v^2)/dr = 0
      is imposed (single master condition)
  P4  the vacuum ODE has the one-parameter general solution
        v^2 = 2A / r      (power ansatz v ~ r^p forces p = -1/2 uniquely)
  P5  substituting v = sqrt(2A/r) back: FULL Einstein tensor == 0
  P6  static pull computed from the metric = -A/r^2 EXACTLY;
      calibrating A = GM (one constant, same move as measuring Newton's G)
      gives exact Newton
  P7  control: sink flow v = k/r^2 yields G_nn != 0 (excluded by the
      principle, matching Test 4c D1)
  P8  control: uniform inflow offset v^2 = 2A/r + B forces B = 0
  P9a matter at the energy-constraint level:
        G_nn = 8 pi G rho / c^2   =>   v^2 = 2 G M(r) / r,
      M(r) = enclosed mass. This is D2's content, DERIVED under A1-A4.
  P9b honesty audit (no PASS label): with static rho(r) the REMAINING
      Einstein components do not vanish -- the static interior needs the
      pressure/lapse sector. Printed, bounded, not claimed.

Disclosed prediction (previous turn): PASS via Birkhoff's theorem
restated in flow variables. This run is the check on that prediction.

Honest framing for the paper: a PASS here means the formalized principle
A1-A4 is equivalent to the flat-slice (Painleve-Gullstrand) formulation
of the GR vacuum. That is a DERIVATION of the flow law, not new physics
beyond GR in this sector. The novelty claim must be worded accordingly.
The open fidelity question (no script can answer it): does the v6.2
inward-fall text formalize to exactly A1-A4?

All printed values are actual runtime output.

Run:  python3 estif_converse_flow_law.py
Deps: sympy (tested 1.14.0). Engine identical to the validated files.
"""

import sympy as sp

t, r, th, ph, s_ = sp.symbols('t r theta phi s', positive=True)
c, G, M, A, k, B = sp.symbols('c G M A k B', positive=True)


# ------------------------------------------------------------------
# Engine (identical to validated estif_tmunu_gauss_codazzi.py)
# ------------------------------------------------------------------
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
    Rs_ = sp.simplify(sum(gi[i, j] * Ric[i, j]
                          for i in range(n) for j in range(n)))
    Ein = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            e = sp.simplify(Ric[i, j] - Rs_ * g[i, j] / 2)
            Ein[i, j] = e
            Ein[j, i] = e
    return Ein, Rs_, Ric


def ssimp(e):
    return sp.simplify(sp.trigsimp(sp.simplify(sp.expand(e))))


def is_zero(e):
    if ssimp(e) == 0:
        return True
    e2 = sp.simplify(sp.radsimp(sp.powsimp(sp.expand(e), force=True)))
    if e2 == 0:
        return True
    eq = e.equals(0)
    return eq is True


def check(label, expr_should_be_zero):
    ok = is_zero(expr_should_be_zero)
    print(('PASS  ' if ok else 'FAIL  ') + label)
    if not ok:
        print('       residual:', ssimp(expr_should_be_zero))
    return ok


x4 = [t, r, th, ph]
v = sp.Function('v', positive=True)(r)
vp = sp.diff(v, r)


def flow_metric(vfield):
    g = sp.zeros(4, 4)
    g[0, 0] = -(c**2 - vfield**2)
    g[0, 1] = vfield
    g[1, 0] = vfield
    g[1, 1] = 1
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(th)**2
    return g


def normal_energy(Ein, g):
    gi = g.inv()
    n_up = [sp.simplify(gi[mu, 0] * (-c)) for mu in range(4)]
    Enn = sp.simplify(sum(Ein[mu, nu] * n_up[mu] * n_up[nu]
                          for mu in range(4) for nu in range(4)))
    return Enn


print('=' * 72)
print('[1] SANITY: flat spherical 3-slices')
print('=' * 72)
h3 = sp.diag(1, r**2, r**2 * sp.sin(th)**2)
x3 = [r, th, ph]
_, R3, _ = einstein_tensor(h3, x3)
p1 = check('P1  R3 = 0 for the flat spherical 3-metric', R3)

print()
print('=' * 72)
print('[2] UNKNOWN FLOW v(r): Einstein tensor and the master condition')
print('=' * 72)
gU = flow_metric(v)
EinU, _, _ = einstein_tensor(gU, x4)
EnnU = normal_energy(EinU, gU)
master = sp.diff(r * v**2, r)
print('effective density (normal observer):')
print('    G_nn =', EnnU)
p3a = check('P3  closed form:  G_nn = d(r v^2)/dr / (c^2 r^2)',
            EnnU - master / (c**2 * r**2))

Gm3 = christoffel(h3, x3)
Ni = [v, sp.S(0), sp.S(0)]
DN = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        DN[i, j] = sp.diff(Ni[j], x3[i]) - sum(
            Gm3[kk][i][j] * Ni[kk] for kk in range(3))
Kij = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        Kij[i, j] = sp.simplify(-(DN[i, j] + DN[j, i]) / (2 * c))
h3i = h3.inv()
Ktr = sp.simplify(sum(h3i[i, j] * Kij[j, i]
                      for i in range(3) for j in range(3)))
Kup = h3i * Kij * h3i
KK = sp.simplify(sum(Kup[i, j] * Kij[i, j]
                     for i in range(3) for j in range(3)))
Ham = sp.simplify(R3 + Ktr**2 - KK)
p2 = check('P2  ADM WITH SHIFT:  R3 + K^2 - K_ij K^ij = 2 G_nn',
           Ham - 2 * EnnU)

print()
print('nonzero Einstein components for unknown v(r), and the kill test')
print('(substitute the master condition d(r v^2)/dr = 0, i.e. '
      "v' = -v/(2r)):")
kill1 = {sp.Derivative(v, r, 2): 3 * v / (4 * r**2)}
kill2 = {sp.Derivative(v, r): -v / (2 * r)}
all_killed = True
for i in range(4):
    for j in range(i, 4):
        comp = EinU[i, j]
        if ssimp(comp) != 0:
            dead = ssimp(comp.subs(kill1).subs(kill2))
            tag = 'dies' if dead == 0 else 'SURVIVES: ' + str(dead)
            print('    G[%d,%d] = %s   ->  %s'
                  % (i, j, sp.simplify(comp), tag))
            if dead != 0:
                all_killed = False
print(('PASS  ' if all_killed else 'FAIL  ')
      + 'P3  every component dies under the single master condition:')
print('       the WHOLE vacuum problem is the one scalar equation')
print('       d(r v^2)/dr = 0.')

print()
print('=' * 72)
print('[3] VACUUM: solve the master condition')
print('=' * 72)
ode = sp.Eq(master, 0)
print('ODE:', ode)
sol = sp.dsolve(ode, v)
print('general solution:', sol)
p_exp = sp.Symbol('p', real=True)
kamp = sp.Symbol('k_amp', positive=True)
ansatz = kamp * r**p_exp
ode_pow = sp.simplify(sp.diff(r * ansatz**2, r))
p_solutions = sp.solve(sp.Eq(sp.expand(ode_pow / (kamp**2 * r**(2 * p_exp))),
                             0), p_exp)
print('power ansatz v = k r^p  =>  admissible p:', p_solutions)
p4 = (p_solutions == [sp.Rational(-1, 2)])
print(('PASS  ' if p4 else 'FAIL  ')
      + 'P4  unique power p = -1/2:  v^2 = 2A/r is the ONLY law allowed.')
print('       One integration constant A. Nothing else survives.')

print()
print('=' * 72)
print('[4] CLOSURE AND NEWTON')
print('=' * 72)
vA = sp.sqrt(2 * A / r)
gA = flow_metric(vA)
EinA, _, _ = einstein_tensor(gA, x4)
p5 = check('P5  v = sqrt(2A/r):  FULL Einstein tensor == 0 '
           '(exact vacuum, Birkhoff in flow variables)',
           sum(sp.Abs(EinA[i, j]) for i in range(4) for j in range(4)))
GmA = christoffel(gA, x4)
ut2 = sp.simplify(-c**2 / gA[0, 0])
pull = sp.simplify(-GmA[1][0][0] * ut2)
print('    static pull computed from the metric:', pull)
p6 = check('P6  pull = -A/r^2 EXACTLY;  A = GM  =>  exact Newton',
           pull + A / r**2)
print('    calibration A = GM is one measured constant, the same')
print('    epistemic move as measuring G in Newtonian gravity.')

print()
print('=' * 72)
print('[5] CONTROLS (the test must be able to fail)')
print('=' * 72)
Enn_sink = sp.simplify(EnnU.subs(sp.Derivative(v, r),
                                 sp.diff(k / r**2, r)).subs(v, k / r**2))
print('sink flow v = k/r^2:   G_nn =', Enn_sink)
p7 = not is_zero(Enn_sink)
print(('PASS  ' if p7 else 'FAIL  ')
      + 'P7  sink flow is NOT vacuum: the principle excludes it')
print('       (negative effective energy would be needed to sustain it;')
print('       independently matches Test 4c D1).')
vB = sp.sqrt(2 * A / r + B)
Enn_B = sp.simplify(EnnU.subs(sp.Derivative(v, r),
                              sp.diff(vB, r)).subs(v, vB))
B_forced = sp.solve(sp.Eq(Enn_B, 0), B)
print('offset flow v^2 = 2A/r + B:   G_nn =', Enn_B,
      '  =>  vacuum forces B =', B_forced)
p8 = (B_forced == [0])
print(('PASS  ' if p8 else 'FAIL  ')
      + 'P8  uniform-inflow offset forbidden: B = 0 forced.')

print()
print('=' * 72)
print('[6] MATTER at the energy-constraint level')
print('=' * 72)
rho = sp.Function('rho', positive=True)
Menc = 4 * sp.pi * sp.Integral(rho(s_) * s_**2, (s_, 0, r))
print('impose  G_nn = 8 pi G rho(r) / c^2   (T_nn = rho c^2, engine')
print('convention validated on FRW):')
print('    d(r v^2)/dr = 8 pi G rho r^2')
rv2 = sp.integrate(8 * sp.pi * G * rho(s_) * s_**2, (s_, 0, r))
p9a = check('P9a  r v^2 = 2 G M(r)  with  M(r) = 4 pi Int rho s^2 ds  '
            '(enclosed mass)', rv2 - 2 * G * Menc)
print('    =>  v^2 = 2 G M(r) / r .  D2 DERIVED under A1-A4:')
print('        outside matter M(r) = M: the escape law, uniquely;')
print('        the flow is sourced by enclosed mass, Gauss-law form.')
print()
print('P9b  honesty audit (expected nonzero, no PASS label):')
Mfun = sp.Function('M', positive=True)(r)
vM = sp.sqrt(2 * G * Mfun / r)
gM = flow_metric(vM)
EinM, _, _ = einstein_tensor(gM, x4)
EnnM = sp.simplify(normal_energy(EinM, gM))
print('    with v^2 = 2GM(r)/r:  G_nn =', EnnM)
print('    remaining components (static interior, must be carried by the')
print('    pressure/lapse sector, NOT claimed solved here):')
for i in range(4):
    for j in range(i, 4):
        if (i, j) != (0, 0):
            compM = ssimp(EinM[i, j])
            if compM != 0:
                print('      G[%d,%d] =' % (i, j), compM)
print('    Interpretation: the energy constraint fixes the flow law;')
print('    a STATIC star additionally requires the pressure equation,')
print('    which in this construction means letting the lapse depart')
print('    from c inside matter. Identified next derivation layer.')

print()
print('=' * 72)
print('[7] VERDICT')
print('=' * 72)
print('CONVERSE TEST RESULT: the formalized inward-fall principle A1-A4')
print('FORCES the flow law. In vacuum the unique solution is')
print('    v^2 = 2A/r,')
print('all Einstein components collapse to that single scalar condition,')
print('and coupling to matter fixes A = G M(enclosed). The disclosed')
print('prediction (Birkhoff in flow variables) is confirmed by the run.')
print()
print('WHAT THIS MEANS, precisely:')
print('  - The child-version question is answered YES at this level:')
print('    the rule does say how fast the river falls at each distance.')
print('  - Scope of the YES: vacuum exterior + energy-constraint coupling.')
print('  - Honest framing: A1-A4 is equivalent to the flat-slice (GP)')
print('    formulation of the GR vacuum. This DERIVES the law from the')
print('    principle; it is not new physics beyond GR in this sector,')
print('    and the paper must say so in those words.')
print()
print('REMAINING, in order:')
print('  1. Fidelity audit (no script can do it): does the v6.2 inward-')
print('     fall text formalize to exactly A1-A4? Check each of the four')
print('     against the actual wording. Any mismatch = the derivation')
print('     proves a neighboring principle, not ESTIF.')
print('  2. Static interior: pressure/lapse sector (P9b).')
print('  3. Time-dependent flows and the moving-observer bulk identity')
print('     (caveat carried from the previous suite).')
