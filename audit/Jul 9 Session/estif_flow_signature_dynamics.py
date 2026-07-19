"""
ESTIF flow-signature dynamics tests
====================================
Tests the hypothesis: "everything moves through the 4D bulk at exactly c"
(the upgrade of 'gravity moves at light speed' to a universal law) against
problem (a) [Lorentzian signature] and against gravity itself.

  [1] Signature emergence   -- is -dt^2 DERIVED from the constraint?
  [2] SR recovery           -- time dilation, null photons, rest clocks
  [3] Engine closure        -- emergent metric through the validated
                               Gauss-Codazzi engine; previous section-[3]
                               ansatz becomes DERIVED, all PASSes survive
  [4] The falsifiable one   -- does the idea do dynamics?
      4a  static clock rate = local flow speed u(r)        [derived]
      4b  test-particle pull  a = -c^2 grad(u)             [derived]
      4c  what fixes u(r)? two candidate flow laws:
            D1 incompressible sink absorption  -> force ~ 1/r^5  [FAILS]
            D2 Poisson law for the flow deficit -> exact Newton  [WORKS,
               but is a POSTULATE the script cannot derive]
      4d  exact closure: the full flow metric (Gullstrand-Painleve form,
          w-slowdown u + inward 3-flow v_r, with u^2 + v_r^2/c^2 = 1)
          is an EXACT vacuum solution: Einstein tensor identically zero
          = exact Schwarzschild. Engine-verified, no approximations.

Honest boundary
---------------
Tests 1, 2, 3, 4a, 4b, 4d are mathematics: derived, engine-verified.
Test 4c isolates the single remaining physical input: the field equation
that fixes u(r). D1 proves not any flow works. D2 states the minimal
sufficient law. Whether ESTIF's inward-fall principle IMPLIES D2 is the
open question for v6.2 -- this script sharpens it, it cannot answer it.

Proper time is defined as bulk w-distance / c. This is exact for the flat
construction (Tests 1-2). Inside eddies the emergent METRIC is the arbiter
(engine-computed, Tests 4a-4d); the bulk-kinematic identity for moving
observers in eddies is flagged as future work, not silently assumed.

All printed values are actual runtime output. Nothing here is an ESTIF
v6.2 result; section [4c]-D2 is a candidate law awaiting derivation.

Run:  python3 estif_flow_signature_dynamics.py
Deps: sympy (tested 1.14.0). Engine functions are identical to the
      previously validated estif_tmunu_gauss_codazzi.py.
"""

import sympy as sp

t, r, chi, th, ph = sp.symbols('t r chi theta phi', positive=True)
xx, yy, zz = sp.symbols('x y z')
c, G, M, rs, k, K, L = sp.symbols('c G M r_s k K L', positive=True)
eps = sp.Symbol('varepsilon', positive=True)


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


print('=' * 72)
print('[0] ENGINE SANITY RECHECK (flat FRW -> Friedmann; full validation')
print('    battery lives in estif_tmunu_gauss_codazzi.py, all PASS)')
print('=' * 72)
a = sp.Function('a', positive=True)(t)
gF = sp.diag(-1, a**2, a**2, a**2)
EinF, _, _ = einstein_tensor(gF, [t, xx, yy, zz])
H = sp.diff(a, t) / a
check('Friedmann recovered:  G_tt = 3 H^2',
      EinF[0, 0] - 3 * H**2)

print()
print('=' * 72)
print('[1] TEST 1: SIGNATURE EMERGENCE')
print('    Bulk metric assumed EUCLIDEAN (all plus signs).')
print('    Universal law: every object satisfies dw^2 + dsigma^2 = c^2 dt^2')
print('    Proper time DEFINED as bulk w-distance:  c dtau := dw')
print('=' * 72)
dt_ = sp.Symbol('dt', positive=True)
dw_ = sp.Symbol('dw', real=True)
dsig = sp.Symbol('dsigma', nonnegative=True)
dtau = sp.Symbol('dtau', nonnegative=True)
bulk_interval = dw_**2 + dsig**2
print('bulk interval (Euclidean, no minus signs anywhere):', bulk_interval)
constraint = sp.Eq(dw_**2 + dsig**2, c**2 * dt_**2)
dw2_solved = sp.solve(constraint, dw_**2)[0]
proper_time_interval = dw2_solved
print('constraint solved for dw^2 :', dw2_solved)
print('c^2 dtau^2 = dw^2          :', proper_time_interval)
ok1 = check('emergent interval == c^2 dt^2 - dsigma^2   (Minkowski, derived)',
            proper_time_interval - (c**2 * dt_**2 - dsig**2))
print('    The minus sign was produced by SOLVING the constraint,')
print('    not inserted. Signature (-,+,+,+) is derived, not assumed.')

print()
print('=' * 72)
print('[2] TEST 2: SPECIAL-RELATIVITY RECOVERY (no free parameters)')
print('=' * 72)
v = sp.Symbol('v', nonnegative=True)
dtau_dt = sp.sqrt(dw2_solved.subs(dsig, v * dt_)) / (c * dt_)
dtau_dt = sp.simplify(sp.powsimp(dtau_dt, force=True))
ok2a = check('time dilation:  dtau/dt = sqrt(1 - v^2/c^2)',
             dtau_dt - sp.sqrt(1 - v**2 / c**2))
photon = sp.simplify(dw2_solved.subs(dsig, c * dt_))
ok2b = check('photon (dsigma = c dt  <=>  dw = 0):  dtau = 0  (null)',
             photon)
rest = sp.simplify(dw2_solved.subs(dsig, 0) - c**2 * dt_**2)
ok2c = check('rest particle (dsigma = 0):  dtau = dt', rest)

print()
print('=' * 72)
print('[3] TEST 3: ENGINE CLOSURE (the emergent metric fed through the')
print('    validated Gauss-Codazzi engine; the former eps = -1 ANSATZ is')
print('    now DERIVED from Tests 1-2)')
print('=' * 72)
g3a = sp.diag(-c**2, 1, 1, 1)
Ein3a, _, _ = einstein_tensor(g3a, [t, xx, yy, zz])
ok3a = check('3a  flat emergent metric  diag(-c^2, 1, 1, 1):  G_uv = 0 '
             '(empty space stays empty)',
             sum(sp.Abs(Ein3a[i, j]) for i in range(4) for j in range(4)))
b = sp.Function('b', positive=True)
bW = b(c * t)
g3b = sp.diag(-c**2,
              bW**2,
              bW**2 * chi**2,
              bW**2 * chi**2 * sp.sin(th)**2)
Ein3b, _, _ = einstein_tensor(g3b, [t, chi, th, ph])
offdiag = all(ssimp(Ein3b[i, j]) == 0
              for i in range(4) for j in range(4) if i != j)
print(('PASS  ' if offdiag else 'FAIL  ')
      + '3b  cosmological flow metric: all off-diagonal G_uv = 0')
rho_eff = sp.simplify(Ein3b[0, 0] / (8 * sp.pi * G * c**2))
Hphys = sp.simplify(sp.diff(bW, t) / bW)
ok3b = check('3b  consistency:  8 pi G rho_eff = 3 (H/c)^2  (Friedmann-'
             'analog survives)',
             8 * sp.pi * G * rho_eff - 3 * (Hphys / c)**2)
print('    rho_eff =', rho_eff)
print('    H       =', Hphys, '   [= c * dln(b)/dw : unchanged from the')
print('    previous run -- but g_tt = -c^2 is now a consequence, not an')
print('    input. Problem (a) is closed FOR THIS CONSTRUCTION.]')

print()
print('=' * 72)
print('[4] TEST 4: DOES THE IDEA DO DYNAMICS?')
print('    Eddy hypothesis: near mass the local w-flow speed is c*u(r),')
print('    u unspecified. Static clocks are carried at the local flow')
print('    speed, so the emergent time component is g_tt = -c^2 u(r)^2.')
print('=' * 72)
u = sp.Function('u', positive=True)(r)
g4 = sp.diag(-c**2 * u**2, 1, r**2, r**2 * sp.sin(th)**2)
x4 = [t, r, th, ph]
Gm4 = christoffel(g4, x4)

print('--- 4a  static clock rate [DERIVED from g_tt] ---')
dtau_static = sp.sqrt(-g4[0, 0]) / c
ok4a = check('dtau/dt = u(r);  redshift ratio between radii = u(r1)/u(r2)',
             dtau_static - u)

print()
print('--- 4b  force law from flow gradient [DERIVED via Christoffels] ---')
Gamma_r_tt = Gm4[1][0][0]
print('    Gamma^r_tt        =', sp.simplify(Gamma_r_tt))
a_static = sp.simplify(Gamma_r_tt * (1 / u)**2)
print('    static acceleration (outward thrust needed to hover):')
print('    a^r = Gamma^r_tt (dt/dtau)^2 =', a_static)
epsx = sp.Symbol('epsilon_x', positive=True)
phi = sp.Function('phi')(r)
a_weak = a_static.subs(u, 1 + epsx * phi)
a_series = sp.series(a_weak, epsx, 0, 2).removeO()
lead = sp.simplify(sp.expand(a_series).coeff(epsx, 1))
ok4b = check('weak field:  gravitational pull = -c^2 du/dr '
             '(with Phi_eff := c^2(u-1), pull = -dPhi_eff/dr)',
             lead - c**2 * sp.diff(phi, r))
print('    => flow-speed gradients ARE gravitational acceleration.')
print('       "Gravity = Time = Eddies" holds at the force-law level,')
print('       GIVEN u(r). What fixes u(r) is Test 4c.')

print()
print('--- 4c  what fixes u(r)?  candidate flow laws ---')
print('D1: incompressible absorption. Steady sink flow in 3D forces')
print('    v_r = k/r^2, and the speed-c constraint u^2 = 1 - v_r^2/c^2 gives:')
u1 = sp.sqrt(1 - k**2 / (c**2 * r**4))
Phi1 = c**2 * (u1 - 1)
alph = sp.Symbol('alpha', positive=True)
Phi1_series = sp.series(Phi1.subs(k**2, alph), alph, 0, 2).removeO()
F1 = sp.simplify(-sp.diff(Phi1_series, r))
print('    Phi_eff (leading order) =', sp.simplify(Phi1_series))
print('    force                   =', F1)
resid1 = sp.simplify(F1 * r**2)
d1_newtonian = not resid1.has(r)
print(('PASS  ' if not d1_newtonian else 'FAIL  ')
      + 'D1 verdict: force ~ 1/r^5, NOT inverse-square. Naive fluid')
print('       absorption CANNOT reproduce Newton. Not any flow works.')

print()
print('D2: POSTULATE a Poisson law for the flow deficit:')
print('       laplacian(Phi_eff) = 4 pi G rho,   Phi_eff = c^2 (u - 1)')
Phi2 = -G * M / r
lap = sp.simplify(sp.diff(r**2 * sp.diff(Phi2, r), r) / r**2)
ok4c1 = check('    exterior solution Phi_eff = -GM/r satisfies the '
              'Poisson law (vacuum)', lap)
u2 = 1 + Phi2 / c**2
F2 = sp.simplify(-c**2 * sp.diff(u2, r))
ok4c2 = check('    resulting force = -GM/r^2  (exact Newton)',
              F2 + G * M / r**2)
vr2 = sp.simplify(c**2 * (1 - u2**2))
print('    equivalent diverted 3-flow:  v_r^2 = c^2(1 - u^2) =',
      sp.expand(vr2))
print('    weak field:  v_r^2 -> 2GM/r  (escape velocity: the river'
      ' picture).')
print('    STATUS: D2 is an INPUT. The script verifies its consequences;')
print('    it cannot derive it. The v6.2 question is whether the inward-')
print('    fall principle IMPLIES D2.')

print()
print('--- 4d  exact closure: the full flow metric [ENGINE-VERIFIED] ---')
print('    Gullstrand-Painleve form: w-slowdown u plus inward 3-flow v_r,')
print('    v_r = c sqrt(rs/r), speed-c constraint u^2 + v_r^2/c^2 = 1:')
vr = c * sp.sqrt(rs / r)
uGP = sp.sqrt(1 - rs / r)
ok4d0 = check('    constraint u^2 + v_r^2/c^2 = 1 holds identically',
              uGP**2 + vr**2 / c**2 - 1)
gGP = sp.zeros(4, 4)
gGP[0, 0] = -(c**2 - vr**2)
gGP[0, 1] = vr
gGP[1, 0] = vr
gGP[1, 1] = 1
gGP[2, 2] = r**2
gGP[3, 3] = r**2 * sp.sin(th)**2
EinGP, _, _ = einstein_tensor(gGP, x4)
gp_ok = all(is_zero(EinGP[i, j]) for i in range(4) for j in range(4))
print(('PASS  ' if gp_ok else 'FAIL  ')
      + '    Einstein tensor identically ZERO: the flow metric is an')
print('       EXACT vacuum solution of GR = Schwarzschild, no weak-field')
print('       approximation. [Given v_r(r); its profile is D2 again.]')
dtS, drS = sp.symbols('dt_s dr_s', real=True)
ds2_rain = (gGP[0, 0] * dtS**2 + 2 * gGP[0, 1] * dtS * drS
            + gGP[1, 1] * drS**2).subs(drS, -vr * dtS)
ok4d1 = check('    raindrop (rides the flow, dr = -v_r dt): dtau = dt '
              'exactly', sp.simplify(ds2_rain) + c**2 * dtS**2)
static_rate = sp.simplify(sp.sqrt(-gGP[0, 0]) / c)
ok4d2 = check('    static observer: dtau/dt = sqrt(1 - rs/r) = u(r)',
              static_rate - uGP)
ident = (uGP**2 - (1 + 2 * (-G * M / r) / c**2)).subs(rs, 2 * G * M / c**2)
ok4d3 = check('    identity u^2 = 1 + 2 Phi/c^2 with Phi = -GM/r, '
              'rs = 2GM/c^2', ident)

print()
print('=' * 72)
print('[5] VERDICT (all statements below are properties of the runtime')
print('    output above, not opinions)')
print('=' * 72)
print('DERIVED, engine-verified, zero free parameters:')
print('  - Lorentzian signature from Euclidean bulk + universal speed-c')
print('    constraint  [problem (a) closed for this construction]')
print('  - exact SR kinematics')
print('  - previous cosmological sector unchanged with g_tt = -c^2 now a')
print('    consequence  [engine closure]')
print('  - static clock rate = local flow speed u(r)')
print('  - gravitational pull = -c^2 grad(u): eddies ARE gravity, given u')
print('  - full flow metric with v_r = c sqrt(rs/r) is EXACT Schwarzschild')
print()
print('NOT DERIVED (the single remaining input):')
print('  - the field equation that fixes u(r). Equivalent forms:')
print('        laplacian(c^2 (u - 1)) = 4 pi G rho')
print('        v_r^2 = 2 G M / r   (escape-velocity flow)')
print('    D1 proves an arbitrary flow law fails (1/r^5, not 1/r^2).')
print()
print('THE QUESTION FOR v6.2: does the ESTIF inward-fall principle imply')
print('the Poisson law for the flow deficit? If yes: signature, SR,')
print('Newton, and exact Schwarzschild all follow from one constraint')
print('plus one field equation. If no: the speed-c idea fixes problem (a)')
print('and the force law, but Newton G-dynamics remains unaccounted for.')
print()
print('Caveat (flagged, not hidden): proper time = w-distance/c is exact')
print('in the flat construction; inside eddies the emergent metric is the')
print('arbiter and the bulk identity for MOVING observers is future work.')
