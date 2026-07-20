"""
FRONT 1 - THE SHARPENED PENCIL: sigma8 / JWST growth numbers under A1'
======================================================================
INPUT (derived 11 July 2026, RHAC-006): the A1' deepen mode

    D+(a) = H(a) * INTEGRAL_0^a da' / (a' H(a'))^3        [Heath 1977 form]

verified symbolically in tests/test_UKN2.py to solve the amended wrinkle law
d2(delta)/dt2 + 2H d(delta)/dt = 4 pi G rho delta EXACTLY on the frozen-Lambda
background H^2 = H0^2 [Om/a^3 + (1-Om)]. This is GR's linear growth factor:
no ESTIF-specific enhancement exists at linear order (PATH_ONE_CHECKLIST B-6).

THIS SCRIPT (straight numerics, laptop):
  [1] numeric echo of the fork receipt: f(z=0.5) from D+ must return 0.76
  [2] the clumping history: D(z), f(z), sigma8(z), f*sigma8(z), z = 0..12
  [3] check against DESI measured growth data (verified, cited below)
  [4] JWST go/no-go per the FROZEN spec docs/report/JWST_TEST_SPEC.md v1.1:
      g(9.1) = D_ESTIF/D_LCDM fed to the Sheth-Tormen abundance machinery

PASSPORTS (honest ledger):
  Om_ESTIF  = 0.3141   banked (Principle P bootstrap; closed transcendental eq.)
  H0        = 67.66    Path One convention (= colossus 'planck18', Planck+BAO)
  sigma8_0  = colossus 'planck18' value (imported amplitude; same passport as
              Lambda and the 1e-5 seed zeta -- value not derived by ESTIF)
  LCDM ref  = colossus 'planck18' (Om = 0.3111): the SAME baseline cosmology
              frozen in JWST_TEST_SPEC.md / estif_jwst_growth_spec.py

DESI growth data used in [3] (only values verified against the sources):
  P1: f*sigma8(z_eff=0.07) = 0.450 +/- 0.055
      DESI DR1 Peculiar Velocity consensus (Qin et al. 2026, A&A 708, A219;
      Lai et al. 2026, JCAP 04, 026). Their quoted Planck+LCDM prediction at
      this z: 0.449 +/- 0.008.
  P2: f*sigma8(BGS, z_eff=0.295) = 0.38 +/- 0.09
      DESI DR1 Full-Shape RSD, BGS bin (DESI Collaboration 2025e), as quoted
      in the DESI PV correlation-function paper (arXiv:2512.03230).
  A1: sigma8 = 0.841 +/- 0.034
      DESI DR1 Full-Shape + BAO + BBN amplitude (DESI 2024 V, arXiv:2411.12021).
  NOTE: the full six-bin ShapeFit f*sigma8 vector (DESI 2024 V, Sec. 7.1 /
  App. A; bins z_eff = 0.295/0.510/0.706/0.919/1.317/1.491) can be added by
  transcribing the published table -- values are NOT embedded here from memory.
  P1 and P2 share BGS data; pulls below are diagonal and indicative, not a
  covariance-correct likelihood (lesson of the BTFR correlated-error fix).

Run:  python3 estif_front1_growth_sigma8_jwst.py
Deps: numpy, scipy, colossus (already installed for estif_jwst_growth_spec.py)
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

try:
    from colossus.cosmology import cosmology
    from colossus.lss import mass_function
except ImportError:
    raise SystemExit("needs colossus:  pip install colossus --break-system-packages")

if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz

# ------------------------------------------------------------------ inputs
cosmo = cosmology.setCosmology('planck18')      # frozen-spec baseline
OM_LCDM = cosmo.Om0                              # 0.3111 (Planck+BAO)
OM_ESTIF = 0.3141                                # banked bootstrap value
SIG8_0 = cosmo.sigma8                            # imported amplitude
h = cosmo.H0 / 100.0
fb = cosmo.Ob0 / cosmo.Om0

# --------------------------------------------- A1' deepen mode, Heath form
def _grow_integrand(ap, Om):
    # 1/(a E)^3 rewritten smooth at a=0:  a^1.5 / (Om + (1-Om) a^3)^1.5
    return ap ** 1.5 / (Om + (1.0 - Om) * ap ** 3) ** 1.5

def E_of_a(a, Om):
    return np.sqrt(Om / a ** 3 + (1.0 - Om))

def I_of_a(a, Om):
    v, _ = quad(_grow_integrand, 0.0, a, args=(Om,), limit=200)
    return v

def D_raw(a, Om):
    return E_of_a(a, Om) * I_of_a(a, Om)          # D+ = H * int da/(aH)^3, H0 cancels

def D_norm(z, Om):
    a = 1.0 / (1.0 + z)
    return D_raw(a, Om) / D_raw(1.0, Om)

def f_growth(z, Om):
    # f = dlnD/dlna = -(3/2) Om(a)  +  a * g(a) / [ (aE)^3 I(a) ]   (exact)
    a = 1.0 / (1.0 + z)
    E2 = Om / a ** 3 + (1.0 - Om)
    om_a = (Om / a ** 3) / E2
    term = a * _grow_integrand(a, Om) / I_of_a(a, Om)
    return -1.5 * om_a + term

def sigma8_z(z, Om):
    return SIG8_0 * D_norm(z, Om)

def fs8(z, Om):
    return f_growth(z, Om) * sigma8_z(z, Om)

# --------------------------------- Sheth-Tormen abundance (frozen-spec code)
def _mf(Mmin_Msun, z, want):
    Mmin_h = Mmin_Msun * h
    lnM = np.linspace(np.log(Mmin_h), np.log(1e15), 400)
    M = np.exp(lnM)
    dndlnM = mass_function.massFunction(M, z, mdef='fof', model='sheth99',
                                        q_in='M', q_out='dndlnM')
    if want == 'n':
        return np.trapezoid(dndlnM, lnM)
    return np.trapezoid(M * dndlnM, lnM)

def n_cum(Mmin, z):    return _mf(Mmin, z, 'n')
def rho_coll(Mmin, z): return _mf(Mmin, z, 'rho')

Dc = lambda z: cosmo.growthFactor(z)             # colossus D for the boost map

def z_eff_for_boost(g, z):
    # spec helper, extended: g may sit below 1 (root then lies above z)
    target = g * Dc(z)
    if target >= Dc(0.0):
        return 0.0
    return brentq(lambda zz: Dc(zz) - target, 0.0, 30.0, xtol=1e-4)

# ============================================================== the receipt
def main():
    print("=" * 74)
    print("FRONT 1 - A1' GROWTH HISTORY: sigma8(z), DESI growth check, JWST verdict")
    print("=" * 74)
    print(f"  Om (ESTIF, banked)  = {OM_ESTIF:.4f}     Om (LCDM ref) = {OM_LCDM:.4f}")
    print(f"  sigma8_0 (imported) = {SIG8_0:.4f}     h = {h:.4f}   f_b = {fb:.4f}")
    print()

    print("-" * 74)
    print("[1] NUMERIC ECHO of the fork receipt (test_UKN2.py, symbolic -> 0)")
    print("-" * 74)
    f05 = f_growth(0.5, OM_ESTIF)
    print(f"  f(z=0.5) from D+ = H*int da/(aH)^3 with Om = {OM_ESTIF}:  {f05:.4f}")
    print(f"  fork receipt / DESI RSD anchor:                          0.76")
    print(f"  match: {'PASS' if abs(f05 - 0.76) < 0.005 else 'FAIL'}")
    a05 = 1 / 1.5
    om05 = (OM_ESTIF / a05 ** 3) / (OM_ESTIF / a05 ** 3 + 1 - OM_ESTIF)
    print(f"  strict-A1 counterfactual (delta = H/H0): f(0.5) = {-1.5 * om05:.2f}"
          f"  [the healed sign flip]")
    print()

    print("-" * 74)
    print("[2] THE CLUMPING HISTORY (ESTIF = A1' deepen mode, zero free parameters")
    print("    beyond the three passports above)")
    print("-" * 74)
    print(f"  {'z':<7}{'D/D0':>9}{'f':>9}{'sigma8':>10}{'f*sigma8':>11}")
    print("  " + "-" * 46)
    zs = [0.0, 0.07, 0.295, 0.5, 0.706, 0.919, 1.317, 1.491,
          2.0, 3.0, 5.0, 8.0, 9.1, 10.0, 12.0]
    for z in zs:
        print(f"  {z:<7.3f}{D_norm(z, OM_ESTIF):>9.4f}{f_growth(z, OM_ESTIF):>9.4f}"
              f"{sigma8_z(z, OM_ESTIF):>10.4f}{fs8(z, OM_ESTIF):>11.4f}")
    print()

    print("-" * 74)
    print("[3] CHECK vs DESI MEASURED GROWTH (verified points; diagonal pulls)")
    print("-" * 74)
    data = [("PV consensus  z=0.070", 0.070, 0.450, 0.055),
            ("FS BGS (RSD)  z=0.295", 0.295, 0.380, 0.090)]
    chi2 = 0.0
    print(f"  {'point':<24}{'measured':>12}{'ESTIF':>9}{'pull':>8}")
    for name, z, val, err in data:
        pred = fs8(z, OM_ESTIF)
        pull = (pred - val) / err
        chi2 += pull ** 2
        print(f"  {name:<24}{val:>7.3f} +/- {err:<5.3f}{pred:>7.3f}{pull:>+8.2f}")
    print(f"\n  diagonal chi2/N (N=2, indicative only) = {chi2 / 2:.2f}")
    s8_pred, s8_meas, s8_err = sigma8_z(0.0, OM_ESTIF), 0.841, 0.034
    print(f"  amplitude: sigma8(0) carried = {s8_pred:.3f} (import)  vs  DESI DR1")
    print(f"  FS+BAO+BBN sigma8 = {s8_meas:.3f} +/- {s8_err:.3f}   "
          f"pull = {(s8_pred - s8_meas) / s8_err:+.2f}")
    print("  (an import-vs-data consistency check, inherited unchanged from LCDM)")
    print()

    print("-" * 74)
    print("[4] JWST GO/NO-GO per the FROZEN spec (JWST_TEST_SPEC.md v1.1, Sec. 8)")
    print("-" * 74)
    print(f"  {'z':<7}{'D_E/D_E0':>10}{'D_L/D_L0':>10}{'g = ratio':>11}")
    for z in [8.0, 9.1, 10.0, 12.0]:
        gE, gL = D_norm(z, OM_ESTIF), D_norm(z, OM_LCDM)
        print(f"  {z:<7.1f}{gE:>10.5f}{gL:>10.5f}{gE / gL:>11.5f}")
    z0 = 9.1
    g91 = D_norm(z0, OM_ESTIF) / D_norm(z0, OM_LCDM)
    ze = z_eff_for_boost(g91, z0)
    b11 = n_cum(1e11, ze) / n_cum(1e11, z0)
    b31 = n_cum(3e11, ze) / n_cum(3e11, z0)
    print(f"\n  g(9.1) = {g91:.5f}   ->  abundance factor n(>1e11): {b11:.3f}x,"
          f"  n(>3e11): {b31:.3f}x")
    print(f"  (the ~0.2% growth offset is the banked Om = 0.3141 vs the ref 0.3111;")
    print(f"   second-order Om effects on the mass function, <1%, are not")
    print(f"   propagated -- the frozen spec's g-knob is the registered currency)")
    print(f"  spec target for relief: g >= 1.13 (5x reservoir)."
          f"  Shortfall: {1.13 - g91:+.3f}")
    print(f"  eps_required(ESTIF)/eps_required(LCDM) = {1.0 / b11:.3f}"
          f"  (baryon-conversion strain essentially identical)")
    print(f"  Sec. 5a two-sided filter: no boost persists to z=0"
          f" (sigma8(0) = import) -> PASS (trivially: no enhancement exists)")
    E91 = E_of_a(1 / (1 + z0), OM_ESTIF)
    print(f"  strict-A1 counterfactual: delta = H/H0 DECAYS; contrast at z=9.1")
    print(f"  was {E91:.1f}x today's and shrinking -> outcome 5 (falsified);")
    print(f"  the A1' fork is what rescued this test from a kill.")
    print()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"""  Spec Sec. 8, outcome 4 -- HONEST NULL:
    g(9.1) = {g91:.4f} (~= 1.0). ESTIF's linear growth is GR's; it does NOT
    resolve the JWST too-big-too-early tension; the tension is INHERITED from
    LCDM, shared, not worsened. This matches the spec's own registered prior
    expectation ("if the generalized flow reduces to standard Poisson at
    linear order, outcome 4 is forced").
  Growth history vs DESI: f(0.5) = {f05:.2f} (= measured 0.76); f*sigma8 pulls
    {(fs8(0.070, OM_ESTIF) - 0.450) / 0.055:+.2f} / {(fs8(0.295, OM_ESTIF) - 0.380) / 0.090:+.2f} sigma at z = 0.07 / 0.295. The clumping history is
    LCDM's, carried on the banked Om = 0.3141 and one imported amplitude.
  Distinguishing content remains where registered: mean spatial curvature = 0
    exactly at all epochs (kill-shot), and the C-11 discriminator hunt.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
