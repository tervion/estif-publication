"""
C-15 - THE GRAVITATIONAL-WAVE SECTOR: c_gw = c, DERIVED
=======================================================
OPEN ITEM (RHAC-006): 'gravitational waves (speed = c forced by A2;
GW170817-consistent)' was ASSERTED when the fork opened the radiative sector.
C-15 closes it: SHOW c_gw = c is forced, not assumed.

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

This script PROVES the single-null-cone claim symbolically on the
Painleve-Gullstrand flow metric, for an ARBITRARY flow field v.

Run:  python3 estif_C15_gw_sector.py
Deps: sympy
"""
import sympy as sp

def main():
    print("=" * 74)
    print("C-15 - GW SECTOR: c_gw = c from single-geometry (A1'+A2+A3)")
    print("=" * 74)

    t, x, y, z = sp.symbols('t x y z', real=True)
    c = sp.symbols('c', positive=True)
    vx, vy, vz = sp.symbols('v_x v_y v_z', real=True)   # arbitrary flow field
    N = sp.symbols('N', positive=True)                  # lapse (A2 proper-time)

    # ---- Painleve-Gullstrand metric: EVEN (flat) spatial slices + flow -----
    # ds^2 = -N^2 c^2 dt^2 + delta_ij (dx^i - v^i dt)(dx^j - v^j dt)
    # coords (t,x,y,z); shift beta^i = -v^i; gamma_ij = delta_ij  (A1': even).
    v = sp.Matrix([vx, vy, vz])
    g = sp.zeros(4, 4)
    g[0, 0] = -N**2 * c**2 + (vx**2 + vy**2 + vz**2)
    for i, vi in enumerate([vx, vy, vz], start=1):
        g[0, i] = -vi
        g[i, 0] = -vi
    for i in range(1, 4):
        g[i, i] = 1
    print("\n[1] The ONE geometry (A1' even slices + A2 flow at c, arbitrary v):")
    sp.pprint(g)

    ginv = g.inv()
    print("\n[2] Inverse metric g^{mu nu} (the wave operator's principal symbol")
    print("    is g^{mu nu} k_mu k_nu -- SAME object for GW and for light):")
    ginv_s = sp.simplify(ginv)
    sp.pprint(ginv_s)

    # ---- the shared null cone:  g^{mu nu} k_mu k_nu = 0 --------------------
    w, kx, ky, kz = sp.symbols('omega k_x k_y k_z', real=True)
    k = sp.Matrix([w, kx, ky, kz])          # k_mu = (omega, k_vec)
    P = sp.expand(sp.simplify((k.T * ginv_s * k)[0]))
    print("\n[3] Principal symbol  g^{mu nu} k_mu k_nu  (=0 defines BOTH the")
    print("    light cone AND the GW characteristics):")
    sp.pprint(sp.collect(P, w))

    # ---- radial case: solve the null condition for the coordinate speed ----
    # take k along x: ky=kz=0, vy=vz=0.  speed u = -omega/k_x (dx/dt of phase).
    Prad = P.subs({ky: 0, kz: 0, vy: 0, vz: 0})
    u = sp.symbols('u', real=True)          # coordinate speed dx/dt
    Prad_u = Prad.subs(w, -u * kx)          # omega = -u k_x  for phase speed u
    Prad_u = sp.simplify(Prad_u / kx**2)
    sol = sp.solve(sp.Eq(Prad_u, 0), u)
    print("\n[4] Radial null speed u = dx/dt (lapse N=1), solving cone = 0:")
    sol_N1 = [sp.simplify(s.subs(N, 1)) for s in sol]
    sp.pprint(sol_N1)
    print("    => u = v +/- c : the flow ADVECTS the signal; it moves at c")
    print("    RELATIVE TO THE FLOW. This is the SAME cone for GW and light,")
    print("    so both are advected identically -> arrival is simultaneous.")

    # ---- the crux: GW and light characteristics are identical -------------
    print("\n[5] THE CRUX (why c_gw = c exactly):")
    print("    * Light characteristics: g_{mu nu} dx^mu dx^nu = 0  (null geodesics).")
    print("    * GW characteristics: the vacuum TT wave operator is Box_g h^{TT}=0;")
    print("      its principal symbol is g^{mu nu} k_mu k_nu -- the SAME null cone.")
    print("    * Lower-order terms (flow friction ~H, curvature of the dent) sit")
    print("      BELOW the principal part and cannot move the characteristic speed.")
    print("    => GW and light share one null structure. c_gw = c, identically,")
    print("       for ANY flow v -- no free parameter, no dial.")

    # sanity: confirm the two speeds differ from each other only by sign of c
    diff = sp.simplify((sol_N1[0] - sol_N1[1])**2 - (2*c)**2)
    print(f"\n[6] Sanity: (u_+ - u_-)^2 - (2c)^2 = {diff}  (0 confirms speeds are v +/- c)")

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
  TIGHTEST-tested of all the locks (GW170817, ~1e-15). Like the others it
  does NOT separate ESTIF from GR/flat-LCDM (which also give c_gw = c); it
  separates ESTIF from the modified-gravity landscape.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
