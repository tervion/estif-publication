"""
FRONT 2 - THE FIRST-BLACK-HOLE RECIPE (computable inside the A1' fork)
======================================================================
THE QUESTION (author's): when does the first hole form, and by which route?

WHY IT IS NOW COMPUTABLE: strict A1 forbade the growing density mode -- no
dent ever deepens, no patch ever folds back, no hole ever forms (a retro-kill
the fork dodged, RHAC-005/006). Under A1' the recipe is four steps:

    seed (zeta ~ 1e-5, imported)  ->  deepen along D+(z) [Front 1]
    ->  FOLD BACK when the linearly-carried contrast reaches delta_c
    ->  rare dense patches become the first collapsed objects -> holes.

FOLD-BACK RULE (step 3, the only new physics input, and it is standard):
    a patch of mass M sitting nu standard deviations above average folds
    back (detaches from the mean flow and collapses) at the z solving
        nu * sigma(M, z=0) * D(z)/D(0) = delta_c = 1.686
    delta_c is the spherical-collapse threshold (EdS value; exact to <0.5%
    at the matter-dominated epochs relevant here). Since the A1' deepen
    mode IS GR's growth (Front 1), the recipe's numbers are LCDM's --
    inherited, not invented. ESTIF's contribution is that the recipe now
    EXISTS inside the framework, in the dent sector.

CHANNELS (step 4; astrophysics imports, flagged, Barkana & Loeb 2001):
    (i)  Pop III remnant: smallest gas-cooling patch (H2 cooling,
         T_vir ~ 2200 K, mu = 1.22) -> massive first star, ~3 Myr life,
         leaves a ~1e2 Msun hole. The "star channel".
    (ii) Direct collapse (the author's NO-STAR channel): atomic-cooling
         patch (T_vir ~ 1e4 K) folds back and the gas runs straight to a
         1e4-1e6 Msun hole without a stellar detour.
    (iii) Primordial: CLOSED under the zeta = 1e-5 passport -- horizon-
         reentry holes need contrasts ~1e4-1e5 times larger; no such
         boost is derived or assumed. Honest ledger.

PASSPORTS: Om = 0.3141 (banked), H0 = 67.66 (Path One convention),
    amplitude/shape of sigma(M) = colossus 'planck18' (imported, same
    baseline frozen in JWST_TEST_SPEC.md). ESTIF growth offset vs that
    baseline is g = 0.998 (Front 1): shifts every z below by ~0.1 -- quoted
    once, ignored after. delta_c, T_vir formula, Salpeter time: textbook
    imports, not ESTIF-derived.

Run:  python3 estif_front2_first_hole_recipe.py
Deps: numpy, scipy, colossus (same as Front 1)
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

try:
    from colossus.cosmology import cosmology
    from colossus.lss import mass_function, peaks
except ImportError:
    raise SystemExit("needs colossus:  pip install colossus --break-system-packages")

if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz

cosmo = cosmology.setCosmology('planck18')
OM = 0.3141                      # banked (bootstrap); background + D+ run on this
H0_KMSMPC = 67.66
h = cosmo.H0 / 100.0
DELTA_C = 1.686
MPC_KM = 3.0857e19
H0_S = H0_KMSMPC / MPC_KM        # 1/s
GYR_S = 3.1557e16
C_KMS = 2.99792458e5

# ---------------- Front 1 machinery: A1' deepen mode (Heath form), Om banked
def _gint(ap, Om=OM):
    return ap ** 1.5 / (Om + (1.0 - Om) * ap ** 3) ** 1.5

def D_norm(z, Om=OM):
    def Draw(a):
        v, _ = quad(_gint, 0.0, a, args=(Om,), limit=200)
        return np.sqrt(Om / a ** 3 + 1 - Om) * v
    return Draw(1.0 / (1.0 + z)) / Draw(1.0)

# ---------------- background clock and horizon (frozen-Lambda form, banked Om)
def t_of_z(z):
    """cosmic time, Gyr; analytic matter+Lambda"""
    a = 1.0 / (1.0 + z)
    OL = 1.0 - OM
    t = (2.0 / (3.0 * H0_S * np.sqrt(OL))) * np.arcsinh(np.sqrt(OL / OM) * a ** 1.5)
    return t / GYR_S

def comoving_horizon_mpc():
    f = lambda a: 1.0 / (a ** 0.5 * np.sqrt(OM + (1.0 - OM) * a ** 3))
    v, _ = quad(f, 0.0, 1.0, limit=200)
    return (C_KMS / H0_KMSMPC) * v

# ---------------- sigma(M) and abundances on the frozen baseline
def sigma_M0(M_msun):
    R = peaks.lagrangianR(M_msun * h)            # Mpc/h, M in Msun/h
    return cosmo.sigma(R, 0.0)

def z_foldback(M_msun, nu):
    """z where nu * sigma(M,0) * D(z) = delta_c; None if patch never folds"""
    target = DELTA_C / (nu * sigma_M0(M_msun))
    if target >= 1.0:
        return None
    return brentq(lambda z: D_norm(z) - target, 0.0, 199.0, xtol=1e-3)

def n_cum(Mmin_msun, z):
    lnM = np.linspace(np.log(Mmin_msun * h), np.log(1e15), 400)
    M = np.exp(lnM)
    dndlnM = mass_function.massFunction(M, z, mdef='fof', model='sheth99',
                                        q_in='M', q_out='dndlnM')
    return np.trapezoid(dndlnM, lnM)             # (h/Mpc)^3 comoving

def z_first_in_volume(M_msun, V_h3):
    """earliest z with expected count N(>M) * V = 1 in the observable volume"""
    f = lambda z: np.log(max(n_cum(M_msun, z) * V_h3, 1e-300))
    return brentq(f, 2.0, 90.0, xtol=1e-2)

# ---------------- channel mass thresholds (Barkana & Loeb 2001, mu = 1.22)
def M_vir(T_K, z):
    """halo mass (Msun) at virial temperature T; neutral primordial gas"""
    mu = 1.22
    M8_h = (T_K / (1.98e4 * (mu / 0.6))) ** 1.5 * (10.0 / (1.0 + z)) ** 1.5
    return M8_h * 1e8 / h

def solve_channel(T_K, nu):
    """self-consistent fold-back z for the nu-sigma patch at threshold T_vir"""
    def g(z):
        zc = z_foldback(M_vir(T_K, z), nu)
        return (zc if zc is not None else -1.0) - z
    return brentq(g, 0.5, 80.0, xtol=1e-2)

def solve_channel_first(T_K, V_h3):
    """earliest-in-observable-volume z at the (z-dependent) threshold mass"""
    def g(z):
        return z_first_in_volume(M_vir(T_K, z), V_h3) - z
    return brentq(g, 2.0, 85.0, xtol=1e-2)

# ============================================================== the receipt
def main():
    print("=" * 74)
    print("FRONT 2 - FIRST-BLACK-HOLE RECIPE (A1' dent sector, D+ from Front 1)")
    print("=" * 74)
    Rh = comoving_horizon_mpc()
    V_h3 = (4.0 * np.pi / 3.0) * (Rh * h) ** 3
    print(f"  Om = {OM}   H0 = {H0_KMSMPC}   delta_c = {DELTA_C}")
    print(f"  comoving horizon = {Rh / 1e3:.2f} Gpc   "
          f"observable volume = {V_h3:.3e} (Mpc/h)^3")
    print(f"  age today = {t_of_z(0):.2f} Gyr   (matter+Lambda clock, banked Om)")
    print()

    print("-" * 74)
    print("[1] FOLD-BACK CHART: collapse z (and cosmic time) for a nu-sigma patch")
    print("-" * 74)
    masses = [1e5, 1e6, 1e7, 1e8, 1e9]
    nus = [3.0, 4.0, 5.0]
    hdr = "".join(f"{f'nu={n:.0f}':>16}" for n in nus)
    print(f"  {'M [Msun]':<11}{'sigma(M,0)':>11}{hdr}")
    for M in masses:
        s = sigma_M0(M)
        row = ""
        for n in nus:
            zc = z_foldback(M, n)
            row += f"{f'{zc:5.1f} ({t_of_z(zc) * 1e3:4.0f} Myr)':>16}"
        print(f"  {M:<11.0e}{s:>11.2f}{row}")
    print()

    print("-" * 74)
    print("[2] EARLIEST IN THE OBSERVABLE UNIVERSE (expected count = 1, ST)")
    print("    extreme-peak extrapolation of the registered machinery; indicative")
    print("-" * 74)
    for M in [1e5, 1e6, 1e7]:
        zf = z_first_in_volume(M, V_h3)
        print(f"  first patch of M >= {M:.0e} Msun:  z = {zf:5.1f}"
              f"   (t = {t_of_z(zf) * 1e3:4.0f} Myr)")
    print()

    print("-" * 74)
    print("[3] THE CHANNELS (thresholds: Barkana & Loeb 2001, mu = 1.22 -- import)")
    print("-" * 74)
    print("  (i) STAR CHANNEL - H2-cooling patch, T_vir = 2200 K:")
    z3 = solve_channel(2200.0, 3.0)
    zf = solve_channel_first(2200.0, V_h3)
    M3, Mf = M_vir(2200.0, z3), M_vir(2200.0, zf)
    print(f"      typical (3-sigma): z = {z3:.1f}, M = {M3:.1e} Msun,"
          f" t = {t_of_z(z3) * 1e3:.0f} Myr")
    print(f"      earliest in volume: z = {zf:.1f}, M = {Mf:.1e} Msun,"
          f" t = {t_of_z(zf) * 1e3:.0f} Myr")
    print(f"      + ~3 Myr massive-star lifetime -> FIRST ~1e2 Msun HOLE at")
    print(f"      t ~ {t_of_z(zf) * 1e3 + 3:.0f} Myr (earliest) /"
          f" {t_of_z(z3) * 1e3 + 3:.0f} Myr (typical)")
    print()
    print("  (ii) NO-STAR CHANNEL (direct collapse) - atomic patch, T_vir = 1e4 K:")
    z3a = solve_channel(1e4, 3.0)
    zfa = solve_channel_first(1e4, V_h3)
    print(f"      typical (3-sigma): z = {z3a:.1f}, M = {M_vir(1e4, z3a):.1e} Msun,"
          f" t = {t_of_z(z3a) * 1e3:.0f} Myr")
    print(f"      earliest in volume: z = {zfa:.1f}, M = {M_vir(1e4, zfa):.1e} Msun,"
          f" t = {t_of_z(zfa) * 1e3:.0f} Myr")
    print(f"      gas runs straight to a 1e4-1e6 Msun hole (seed mass:")
    print(f"      astrophysics import, order of magnitude)")
    print()
    print("  (iii) PRIMORDIAL - CLOSED. Horizon-reentry holes need contrasts")
    print("      O(0.1-1); the zeta = 1e-5 passport is 1e4-1e5 short. No boost")
    print("      is derived; none is assumed.")
    print()

    print("-" * 74)
    print("[4] THE GROWTH CLOCK - why the no-star channel matters (Salpeter import:")
    print("    45 Myr per e-fold at Eddington, eps = 0.1). Rare quasars grow from")
    print("    RARE seeds: start = 5-sigma fold-back at each channel threshold.")
    print("-" * 74)
    tS = 0.045
    z5s = solve_channel(2200.0, 5.0)
    z5a = solve_channel(1e4, 5.0)
    print(f"  rare-seed starts: star channel z = {z5s:.1f}"
          f" ({t_of_z(z5s) * 1e3:.0f} Myr), no-star z = {z5a:.1f}"
          f" ({t_of_z(z5a) * 1e3:.0f} Myr)")
    for zt in [10.0, 7.0]:
        avail_l = t_of_z(zt) - t_of_z(z5s) - 0.003
        avail_h = t_of_z(zt) - t_of_z(z5a)
        need_l = np.log(1e9 / 1e2) * tS
        need_h = np.log(1e9 / 1e5) * tS
        print(f"  reach 1e9 Msun by z = {zt:.0f} (t = {t_of_z(zt):.2f} Gyr):")
        print(f"    light seed 1e2 (star channel, from z={z5s:.0f}):  need"
              f" {need_l:.2f} Gyr, have {avail_l:.2f} Gyr"
              f"  -> {'OK' if avail_l > need_l else 'STRAINED'}")
        print(f"    heavy seed 1e5 (no-star,     from z={z5a:.0f}):  need"
              f" {need_h:.2f} Gyr, have {avail_h:.2f} Gyr"
              f"  -> {'OK' if avail_h > need_h else 'STRAINED'}")
    print()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"""  The recipe is COMPUTABLE and computed. Under A1', seeds deepen along
  GR's D+ and fold back on schedule:
    star channel:   typical 3-sigma z = {z3:.0f} -> first ~1e2 Msun hole at
                    t ~ {t_of_z(z3) * 1e3 + 3:.0f} Myr; rarest-in-volume z = {zf:.0f} -> t ~ {t_of_z(zf) * 1e3 + 3:.0f} Myr.
    no-star channel: typical 3-sigma z = {z3a:.0f} (t ~ {t_of_z(z3a) * 1e3:.0f} Myr); rare 5-sigma
                    z = {z5a:.0f} (t ~ {t_of_z(z5a) * 1e3:.0f} Myr); 1e4-1e6 Msun holes, no stellar detour.
    primordial:     closed under the zeta = 1e-5 passport.
  The clock says the no-star channel is the comfortable route to 1e9 Msun
  by z = 7 (heavy seed: OK; light seed: strained even from the rarest
  start) -- the author's channel lands where the strain is.
  Every number is LCDM's, inherited through the GR-equivalent deepen mode
  on the banked Om and one imported amplitude. ESTIF-specific content:
  (a) under strict A1 NO hole ever forms -- the fork rescued the question
  itself; (b) fold-back is geometric: a dent deepening until the patch
  detaches from the mean inward fall. Distinguishing numbers do not live
  here; they remain with the curvature kill-shot and the C-11 hunt.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
