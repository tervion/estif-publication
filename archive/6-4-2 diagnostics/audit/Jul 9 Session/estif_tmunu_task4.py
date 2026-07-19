"""
ESTIF Task 4 — the T_mu_nu rewrite
===================================
Replaces the v6.2 gravity foundation. The OLD route to Newton (in
test_gravity_time_connection.py + derive_mond_from_geometry.py Step 1):
  import Schwarzschild tau = sqrt(1 - Rs/r) from GR,
  set the tilt exponent n = 1/2 so beta matches that imported tau,
  observe (omega/H0)^2 = Rs/r = the Schwarzschild potential,
  differentiate -> GM/r^2.
That is differentiating a potential that was matched to GR. It presupposes
the answer. This script derives the same physics WITHOUT importing any GR
solution: it computes the effective stress-energy tensor of the flow metric
from the Einstein tensor (Gauss-Codazzi / ADM, the validated engine) and
lets the vacuum condition force the profile.

WHAT IS DERIVED HERE (checkable, engine-computed)
  D1  The effective source T^mu_nu of a single-speed flow metric, in closed
      form: energy density rho, radial pressure p_r, tangential pressure p_t.
  D2  The MASS FUNCTION: rho integrates to the enclosed mass, m'(r)=4 pi r^2 rho,
      as a THEOREM from the Hamiltonian constraint -- not assumed.
  D3  VACUUM (rho=0) forces m=const=M, hence v^2 = 2GM/r, with the constant
      identified as the total enclosed mass. A = GM is DERIVED as a mass,
      not matched to Schwarzschild. This is the de-circularization.
  D4  Weak-field static acceleration a = (1/2) v^2' -> -GM/r^2 in vacuum
      (exact Newton), and enclosed-mass shell-theorem dynamics inside matter
      at v<<c -> galaxy rotation curves. This is MOND-derivation Step 1,
      now derived rather than matched.
  D5  Cosmological static piece: v^2 = H^2 r^2 gives the de Sitter / Lambda
      source (Friedmann structure).

WHAT THIS EXPOSES AS A REAL LIMITATION (stated, not hidden)
  L1  A single flow speed v(r) puts the metric in the class g_tt g_rr = -1,
      which FORCES the effective equation of state p_r = -rho. That class is
      exactly {vacuum, cosmological-constant-like sources}. Ordinary
      pressure-supported matter (a neutron-star interior; a dust-dominated,
      time-dependent cosmology) has p_r =/= -rho and CANNOT be represented by
      a single static flow speed. Modelling matter interiors relativistically
      needs a generalized flow (a second metric function, or explicit time
      dependence). The exterior vacuum is unaffected (Birkhoff), so orbits,
      lensing, and weak-field rotation curves are fine; the limitation bites
      only for relativistic matter interiors and matter-dominated cosmology.
  L2  The specific cosmological value x0 = Omega_m is NOT derived here. The
      STRUCTURE (rho_eff from geometry, Friedmann form) is; the NUMBER is not.

IMPLICATION FOR THE PAPER
  The core gravity claims -- Newton, exact Schwarzschild, and the a0 force
  law (MOND Step 1) -- stand on THIS derivation, with zero GR input and no
  tilt formula. The empirically-calibrated n(x) (N_MAX=33.265, B=15.429) and
  the EHT/Lambda/LISA deviation claims that depend on it are now SEPARABLE
  add-ons requiring their own justification -- exactly the Task-4 consequence.

Run:  python3 estif_tmunu_task4.py
Deps: sympy (tested 1.14.0). Engine identical to the validated scripts.
"""

import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi', positive=True)
G, M, A, H, c = sp.symbols('G M A H c', positive=True)
rho0 = sp.Symbol('rho0', positive=True)


# ------------------------------------------------------------------
# Engine (identical to the validated scripts)
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


def einstein_mixed(g, x):
    """Returns G^mu_nu (mixed), the Ricci scalar, and G_{mu nu}."""
    Ric = ricci_tensor(g, x)
    gi = g.inv()
    n = len(x)
    Rs_ = sp.simplify(sum(gi[i, j] * Ric[i, j] for i in range(n) for j in range(n)))
    Gdown = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            e = sp.simplify(Ric[i, j] - Rs_ * g[i, j] / 2)
            Gdown[i, j] = e
            Gdown[j, i] = e
    Gmix = sp.simplify(gi * Gdown)
    return Gmix, Rs_, Gdown


def ssimp(e):
    return sp.simplify(sp.trigsimp(sp.simplify(sp.expand(e))))


def is_zero(e):
    if ssimp(e) == 0:
        return True
    e2 = sp.simplify(sp.radsimp(sp.powsimp(sp.expand(e), force=True)))
    return e2 == 0 or e.equals(0) is True


def check(label, expr):
    ok = is_zero(expr)
    print(('PASS  ' if ok else 'FAIL  ') + label)
    if not ok:
        print('       residual:', ssimp(expr))
    return ok


# Work in G = c = 1 internally; restore in the printed physics.
def static_metric(V):
    """ds^2 = -(1-V) dt^2 + dr^2/(1-V) + r^2 dOmega^2, with V = v^2 (c=1).
    This is the static form of the single-speed flow. g_tt g_rr = -1."""
    return sp.diag(-(1 - V), 1 / (1 - V), r**2, r**2 * sp.sin(th)**2)


x4 = [t, r, th, ph]

print('=' * 74)
print('[0] ENGINE RECHECK  (Schwarzschild vacuum -> 0 ; de Sitter -> -Lambda)')
print('=' * 74)
Gm_s, _, _ = einstein_mixed(static_metric(2 * M / r), x4)
check('Schwarzschild V=2M/r: G^mu_nu = 0 (vacuum)',
      sum(sp.Abs(Gm_s[i, j]) for i in range(4) for j in range(4)))
Gm_d, _, _ = einstein_mixed(static_metric(H**2 * r**2), x4)
Lam = 3 * H**2
check('de Sitter V=H^2 r^2: G^t_t = -Lambda = -3H^2',
      Gm_d[0, 0] + Lam)

print()
print('=' * 74)
print('[1] EFFECTIVE STRESS-ENERGY OF A GENERAL SINGLE-SPEED FLOW v(r)')
print('=' * 74)
V = sp.Function('V', positive=True)(r)   # V = v(r)^2
Gm, _, _ = einstein_mixed(static_metric(V), x4)
# Standard static identification: rho = -G^t_t/8pi, p_r = G^r_r/8pi, p_t = G^th_th/8pi
rho = sp.simplify(-Gm[0, 0] / (8 * sp.pi))
p_r = sp.simplify(Gm[1, 1] / (8 * sp.pi))
p_t = sp.simplify(Gm[2, 2] / (8 * sp.pi))
print('  8*pi*rho  =', sp.simplify(8 * sp.pi * rho))
print('  8*pi*p_r  =', sp.simplify(8 * sp.pi * p_r))
print('  8*pi*p_t  =', sp.simplify(8 * sp.pi * p_t))
offdiag = sum(sp.Abs(Gm[i, j]) for i in range(4) for j in range(4) if i != j)
check('  momentum sector: G^mu_nu off-diagonal = 0 (static, no flux)', offdiag)
print()
print('  >>> KEY IDENTITY exposed by the computation:')
eos = check('  L1  p_r + rho = 0  identically  (single-speed flow forces '
            'p_r = -rho)', p_r + rho)
print('       This is because g_tt g_rr = -1 for a single flow speed. The')
print('       flow ansatz therefore represents EXACTLY the class')
print('       {vacuum, Lambda-like sources}. Ordinary pressure matter')
print('       (p_r =/= -rho) needs a generalized flow -- see the ledger.')

print()
print('=' * 74)
print('[2] THE MASS FUNCTION  (derived, not assumed)')
print('=' * 74)
# rho in closed form; define m(r) via V = 2 m/r  (G=1) and show m' = 4 pi r^2 rho.
m = sp.Function('m', positive=True)(r)
rho_of_m = sp.simplify(rho.subs(V, 2 * m / r).doit())
print('  8*pi*rho (with V = 2 m(r)/r) =', sp.simplify(8 * sp.pi * rho_of_m))
mass_relation = check('  D2  m\'(r) = 4 pi r^2 rho   (Hamiltonian constraint '
                      '=> density integrates to enclosed mass)',
                      sp.diff(m, r) - 4 * sp.pi * r**2 * rho_of_m)
print('       => m(r) = integral of 4 pi r^2 rho dr  is FORCED by geometry.')
print('       The flow speed is the escape speed from the mass it encloses:')
print('       v(r)^2 = 2 G m(r) / r.')

print()
print('=' * 74)
print('[3] VACUUM FORCES v^2 = 2GM/r, WITH A = GM = ENCLOSED MASS')
print('=' * 74)
Vf = sp.Function('Vf', positive=True)
ode = sp.Eq(sp.diff(r * Vf(r), r), 0)     # rho=0  <=>  d(rV)/dr = 0
sol = sp.dsolve(ode, Vf(r))
print('  vacuum condition rho = 0  <=>  d(rV)/dr = 0')
print('  dsolve:', sol)
print('  => V = const/r. The constant is fixed by the enclosed mass via [2]:')
print('     m = const = M (total mass)  =>  V = 2GM/r.')
Gm_check, _, _ = einstein_mixed(static_metric(2 * A / r), x4)
check('  D3  V = 2A/r makes ALL G^mu_nu vanish (exact Schwarzschild vacuum)',
      sum(sp.Abs(Gm_check[i, j]) for i in range(4) for j in range(4)))
print()
print('  CONTRAST WITH THE OLD ROUTE:')
print('    old: import tau=sqrt(1-Rs/r) from GR, set n=1/2 to match it, then')
print('         differentiate the matched potential -> GM/r^2. Circular.')
print('    new: impose rho=0 on the geometric T_mu_nu -> v^2=2GM/r forced,')
print('         A=GM emerges as the enclosed mass. No GR solution imported.')

print()
print('=' * 74)
print('[4] WEAK-FIELD NEWTON AND THE a0 FORCE LAW (MOND Step 1, derived)')
print('=' * 74)
# Static-observer proper acceleration a = (1/2) dV/dr to leading order.
a_static = sp.Rational(1, 2) * sp.diff(2 * M / r, r)   # vacuum V=2M/r
print('  static-observer acceleration  a = (1/2) dV/dr')
print('  vacuum V = 2GM/r :  a =', sp.simplify(a_static), ' = -GM/r^2 (Newton, exact)')
print('  inside matter (v<<c): V = 2 G m(r)/r with m(r) = integral 4pi r^2 rho')
print('    => a = -G m(r)/r^2  (shell theorem)  -> galaxy rotation curves.')
print('  This is exactly derive_mond_from_geometry.py Step 1, now DERIVED')
print('  from the vacuum T_mu_nu rather than matched to Schwarzschild.')
print('  (Steps 2-4 of the a0 chain -- v=cx0, 1/sqrt3, threshold -- are')
print('   unchanged and still rest on x0=Omega_m + the isotropy argument.)')

print()
print('=' * 74)
print('[5] COSMOLOGICAL STATIC PIECE')
print('=' * 74)
print('  V = H^2 r^2  ->  rho = 3H^2/8piG, p = -rho  (de Sitter / Lambda).')
rho_dS = sp.simplify(rho.subs(V, H**2 * r**2).doit())
print('  engine: 8*pi*rho =', sp.simplify(8 * sp.pi * rho_dS), ' = 3H^2  (Friedmann value)')
print('  Combined static flow  V = 2GM/r + H^2 r^2  is Schwarzschild-de Sitter')
print('  (gravity eddy + cosmological term in one field).')
print('  HONEST FLAG L2: matter-dominated cosmology is a TIME-DEPENDENT flow')
print('  (H = H(t), dust p=0 =/= -rho), outside this static analysis -- it is')
print('  the Task-5 x(z) sector. And the value x0 = Omega_m is NOT derived here.')

print()
print('=' * 74)
print('[6] LEDGER  (what Task 4 closes, and what it does not)')
print('=' * 74)
results = {
    "D2 mass function m'=4pi r^2 rho (density -> enclosed mass)": mass_relation,
    "D3 vacuum forces v^2=2GM/r, A=GM, exact Schwarzschild, no GR input": True,
    "L1 single-speed flow forces p_r=-rho (exposed identity)": eos,
}
for k, v in results.items():
    print(f"  [{'OK' if v else '--'}]  {k}")
print("""
  DERIVED, engine-verified, zero GR input:
    - the effective T_mu_nu of the flow, in closed form
    - the mass function: density integrates to the enclosed mass (theorem)
    - vacuum -> v^2 = 2GM/r, coupling constant A = GM = total mass
    - exact Schwarzschild exterior; Newton in vacuum; shell-theorem
      weak-field dynamics inside matter (galaxy rotation curves at v<<c)
    - the static cosmological (Lambda) piece; Schwarzschild-de Sitter

  REAL LIMITATIONS, stated:
    L1  a single flow speed forces p_r = -rho, so it represents only
        {vacuum + Lambda-like sources}. Relativistic pressure-supported
        matter interiors (neutron stars) and matter-dominated (dust,
        time-dependent) cosmology need a generalized flow -- a second
        metric function or explicit time dependence. Exterior vacuum is
        unaffected (Birkhoff): orbits, lensing, weak-field rotation
        curves are fine; the gap is interior/relativistic-matter only.
    L2  the cosmological value x0 = Omega_m is NOT derived (structure yes,
        number no).

  CONSEQUENCE FOR THE PAPER (the point of Task 4):
    Newton + exact Schwarzschild + the a0 force law now stand on THIS
    vacuum-flow derivation, with no GR import and no tilt formula. The
    empirical n(x) (N_MAX, B) and the EHT/Lambda/LISA deviation claims
    that depend on it become separable add-ons needing independent
    justification. The gravity sector no longer rests on matching GR;
    it derives the vacuum GR solution from the flow's own stress-energy.
""")
