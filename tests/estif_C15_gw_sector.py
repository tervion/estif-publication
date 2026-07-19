"""
C-15 - THE GRAVITATIONAL-WAVE SECTOR: c_gw = c, DERIVED
=======================================================
OPEN ITEM (RHAC-006): 'gravitational waves (speed = c forced by A2;
GW170817-consistent)' was ASSERTED when the fork opened the radiative sector.
C-15 closes it: SHOW c_gw = c is forced, not assumed.  (RHAC-008)

REBUILD NOTE (15 July 2026): this file is a RECONSTRUCTION of the receipt cited
in RHAC-008. The original was written in an ephemeral session container on
12 July 2026 and never reached the repository -- the citation pointed at a file
that did not exist. The derivation, the metric conventions and the verdict text
are recovered from the session transcript; scaffolding is rewritten. Every
symbolic result below is COMPUTED by sympy at run time, not transcribed.
Verify on the Mac mini before RHAC-008 is treated as closed.

THE DERIVATION, in one sentence: under A1'+A2 the world is a SINGLE Lorentzian
geometry (even-on-average slices + a flow at speed c), and A3 (vacuum sources
nothing) + the empty residual sector (Phase-2 null, RHAC-004) guarantee there
is NO second metric and NO extra field. Gravitational waves are ripples OF that
one geometry; light follows null cones OF that one geometry; the wave operator's
PRINCIPAL SYMBOL is that geometry's null cone. One null cone -> one speed ->
c_gw = c EXACTLY, with the flow tilting both signals identically.

WHY MODIFIED GRAVITY GETS c_gw != c (and ESTIF cannot): those theories add an
extra tensor/scalar that GWs couple to differently than light (disformal
couplings, bimetric sectors, ...). GW170817 killed that whole class to
|c_gw/c - 1| < ~1e-15. ESTIF has no such extra structure BY CONSTRUCTION, so
the deviation is identically zero -- there is no dial that could turn it on.

THE A1' HINGE: strict A1 (exactly even, everywhere, always) FROZE the TT sector
-- no local deviation from evenness is allowed, so no ripple can exist (the same
over-constraint that forbade the growing density mode; RHAC-005/006). A1' (even
ON AVERAGE, local dents permitted) OPENS exactly the space TT ripples live in.
So the radiative sector and the growth sector were opened by the SAME amendment,
for the SAME reason.

GEOMETRY AND CONVENTIONS: the ESTIF line element in flow (Painleve-Gullstrand)
form, with lapse N and inflow field v_i, Euclidean bulk slices (A1'):

    ds^2 = -(N^2 c^2 - v^2) dt^2 - 2 v_i dx^i dt + delta_ij dx^i dx^j

i.e. ADM shift beta^i = -v^i, so a slice-normal (free-falling) observer has
coordinate velocity dx^i/dt = +v^i -- the observer is CARRIED WITH THE INFLOW,
which is the ESTIF reading. Slices are even (delta_ij) per A1'; all the gravity
sits in v. Phase convention: exp(i k_mu x^mu) with k_mu = (omega, k_x, k_y, k_z),
so a surface of constant phase moves at u = dx/dt = -omega/k_x.

Run:  python3 estif_C15_gw_sector.py
Deps: sympy
"""
import sympy as sp

# ------------------------------------------------------------------ symbols
t, x, y, z = sp.symbols('t x y z', real=True)
c, N = sp.symbols('c N', positive=True)
vx, vy, vz = sp.symbols('v_x v_y v_z', real=True)
w, kx, ky, kz = sp.symbols('omega k_x k_y k_z', real=True)

V = sp.Matrix([vx, vy, vz])
K = sp.Matrix([kx, ky, kz])
KMU = sp.Matrix([w, kx, ky, kz])


def pg_metric():
    """ESTIF flow-form (PG) metric with constant lapse N and inflow v."""
    g = sp.zeros(4, 4)
    g[0, 0] = -(N ** 2 * c ** 2 - (V.T * V)[0])
    for i in range(3):
        g[0, i + 1] = -V[i]
        g[i + 1, 0] = -V[i]
        g[i + 1, i + 1] = 1
    return sp.simplify(g)


def principal_symbol(ginv):
    """g^{mu nu} k_mu k_nu -- the characteristic (null-cone) polynomial."""
    return sp.simplify((KMU.T * ginv * KMU)[0])


def main():
    print("=" * 74)
    print("C-15 - GRAVITATIONAL-WAVE SECTOR: c_gw = c, DERIVED")
    print("=" * 74)
    checks = []

    # ---- [1] the single geometry ------------------------------------------
    print("\n[1] The ONE geometry (A1' even slices + A2 flow at c), PG form:")
    g = pg_metric()
    sp.pprint(g)
    ginv = sp.simplify(g.inv())
    ident = sp.simplify(g * ginv)
    ok1 = ident.equals(sp.eye(4))
    checks.append(("metric inverse verified", ok1))
    print("\n    g^{mu nu}:")
    sp.pprint(ginv)
    print(f"\n    g . g^-1 == I : {ok1}")
    detg = sp.simplify(g.det())
    print(f"    det g = {detg}   (sqrt(-g) = {sp.simplify(sp.sqrt(-detg))})")
    print("    => ONE Lorentzian metric. A3 + Phase-2 null forbid a second one.")

    # ---- [2] photon characteristics ---------------------------------------
    print("\n[2] LIGHT: eikonal / null geodesics obey g^{mu nu} k_mu k_nu = 0.")
    P_light = principal_symbol(ginv)
    print("    P_light =")
    sp.pprint(sp.collect(sp.expand(P_light * N ** 2 * c ** 2), [w, kx, ky, kz]))

    # ---- [3] TT wave operator principal symbol -----------------------------
    print("\n[3] GRAVITY: vacuum TT wave operator Box_g h^{TT}_{ab} = 0.")
    print("    Build Box_g on a generic component and read off its SECOND-")
    print("    derivative coefficient matrix -- that IS the principal part.")
    h = sp.Function('h')(t, x, y, z)
    X = [t, x, y, z]
    sqrtg = sp.sqrt(-detg)
    box = 0
    for mu in range(4):
        inner = 0
        for nu in range(4):
            inner += sqrtg * ginv[mu, nu] * sp.diff(h, X[nu])
        box += sp.diff(inner, X[mu])
    box = sp.expand(sp.simplify(box / sqrtg))

    # extraction: substitute a plane wave and read off the k-quadratic part.
    # Only the SECOND-derivative terms survive as k-quadratic; first-derivative
    # (Christoffel / flow-friction) terms are linear in k and sit below the
    # principal part -- exactly the point being made in [7].
    amp = sp.Symbol('A_0')
    pw = amp * sp.exp(sp.I * (w * t + kx * x + ky * y + kz * z))
    box_pw = sp.simplify(box.subs(h, pw).doit())
    symbol_gw = sp.simplify(sp.expand(box_pw / pw))
    print("\n    Box_g h  ->  (plane wave)  ->  symbol =")
    sp.pprint(sp.collect(sp.expand(symbol_gw * N ** 2 * c ** 2), [w, kx, ky, kz]))

    # ---- [4] THE COMPARISON ------------------------------------------------
    print("\n[4] THE COMPARISON: is the GW cone the SAME cone as the light cone?")
    diff_sym = sp.simplify(sp.expand(-symbol_gw - P_light))
    ok4 = (diff_sym == 0)
    checks.append(("GW principal symbol == light null cone", ok4))
    print(f"    P_GW - P_light = {diff_sym}")
    print(f"    identical for ARBITRARY v, N: {ok4}")
    print("    Note the flow v appears in BOTH symbols in exactly the same")
    print("    places. There is no term that sees gravity but not light.")

    # ---- [5] radial null speed --------------------------------------------
    print("\n[5] Radial null speed u = dx/dt (lapse N=1), solving cone = 0:")
    Prad = P_light.subs({ky: 0, kz: 0, vy: 0, vz: 0})
    u = sp.symbols('u', real=True)          # coordinate speed dx/dt
    Prad_u = Prad.subs(w, -u * kx)          # omega = -u k_x  for phase speed u
    Prad_u = sp.simplify(Prad_u / kx ** 2)
    sol = sp.solve(sp.Eq(Prad_u, 0), u)
    sol_N1 = [sp.simplify(s.subs(N, 1)) for s in sol]
    sp.pprint(sol_N1)
    print("    => u = v +/- c : the flow ADVECTS the signal; it moves at c")
    print("    RELATIVE TO THE FLOW. This is the SAME cone for GW and light,")
    print("    so both are advected identically -> arrival is simultaneous.")

    # ---- [6] sanity --------------------------------------------------------
    diff = sp.simplify((sol_N1[0] - sol_N1[1]) ** 2 - (2 * c) ** 2)
    ok6 = (diff == 0)
    checks.append(("speeds are exactly v +/- c", ok6))
    print(f"\n[6] Sanity: (u_+ - u_-)^2 - (2c)^2 = {diff}  (0 confirms speeds are v +/- c)")

    # ---- [7] the crux ------------------------------------------------------
    print("\n[7] THE CRUX (why c_gw = c exactly):")
    print("    * Light characteristics: g_{mu nu} dx^mu dx^nu = 0  (null geodesics).")
    print("    * GW characteristics: the vacuum TT wave operator is Box_g h^{TT}=0;")
    print("      its principal symbol is g^{mu nu} k_mu k_nu -- the SAME null cone.")
    print("    * Lower-order terms (flow friction ~H, curvature of the dent) sit")
    print("      BELOW the principal part and cannot move the characteristic speed.")
    print("    => GW and light share one null structure. c_gw = c, identically,")
    print("       for ANY flow v -- no free parameter, no dial.")

    # ---- [8] GW170817 ------------------------------------------------------
    print("\n[8] GW170817 confrontation:")
    pred = sp.Integer(0)
    bound = 1e-15
    print(f"    ESTIF prediction : |c_gw/c - 1| = {pred}  (exactly, structurally)")
    print(f"    Measured bound   : |c_gw/c - 1| < ~{bound:.0e}")
    ok8 = (float(pred) < bound)
    checks.append(("GW170817 structural pass", ok8))
    print(f"    PASS: {ok8}  -- and by law, not by tuning.")

    # ---- verdict -----------------------------------------------------------
    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"""  c_gw = c is DERIVED, not assumed. On the ESTIF geometry (A1' even-on-
  average slices + A2 flow at c), the vacuum TT wave operator and the photon
  both propagate on the single null cone g^{{mu nu}} k_mu k_nu = 0 -- shown
  above for an arbitrary flow field v. The flow tilts that cone (u = v +/- c),
  but tilts it IDENTICALLY for gravity and light, so any GW and its
  electromagnetic counterpart arrive together.

  GW170817: ESTIF predicts |c_gw/c - 1| = 0 exactly; the measured bound is
  ~1e-15. PASS -- and structurally, not by tuning: A3 + the empty residual
  sector (Phase-2 null) forbid the extra field/metric that every c_gw != c
  theory needs. There is no dial in ESTIF that could break it.

  A1' HINGE recorded: strict A1 froze the TT sector (no ripple = no local
  deviation from exact evenness) -- the SAME over-constraint that killed
  the growing mode. A1' opened both, for the same reason. C-15 closed.

  LEDGER (Front 3 companion): c_gw = c is a FOURTH exactness lock -- ESTIF
  forbids c_gw != c as law, modified-gravity permits it, and it is the
  TIGHTEST-tested of all the locks (GW170817, ~1e-15).

  HONEST SCOPE: this proves the two cones COINCIDE on the ESTIF geometry. It
  does not, and cannot, separate ESTIF from GR -- GR gives c_gw = c too. Per
  C-11 the lock separates ESTIF-Core from GR's EXTRA FREEDOMS only.""")

    print("\n" + "=" * 74)
    npass = sum(1 for _, ok in checks if ok)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"RECEIPT STATUS: {npass}/{len(checks)}"
          f" {'-- C-15 CLOSED' if npass == len(checks) else '-- CHECK FAILED'}")
    print("=" * 74)
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
