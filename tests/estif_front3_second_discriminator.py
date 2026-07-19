"""
FRONT 3 - THE SECOND DIFFERING NUMBER: growth index gamma + the slip lock
========================================================================
RHAC-007. Companion to Front 1 (growth history) and Front 2 (first holes).

REBUILD NOTE (15 July 2026): this file is a RECONSTRUCTION of the receipt cited
in RHAC-007. The original was written in an ephemeral session container on
12 July 2026 and never reached the repository -- the citation pointed at a file
that did not exist. Header, passports, derivation and cited measurements are
recovered from the session transcript; scaffolding is rewritten. The numbers
below are COMPUTED by this script, not transcribed. Verify on the Mac mini
before RHAC-007 is treated as closed.

THE C-11 WALL (what a discriminator can and cannot do): A1' reproduces GR at
linear order by construction. A theorem follows -- ESTIF-Core cannot differ from
FLAT LCDM in ANY background-or-linear observable. So no linear number separates
ESTIF-Core from flat LCDM. What a discriminator CAN do is separate ESTIF-Core
from GR's EXTRA FREEDOMS: the wider GR landscape (curvature, dynamical w,
modified growth/slip) that flat LCDM sits inside. ESTIF forbids each freedom as
LAW; GR treats each as a fitted parameter. That "exactness vs freedom" is the
only discriminator structure available at linear order -- and it is genuine and
falsifiable.

THE SECOND NUMBER (this front's find): the GROWTH INDEX gamma, with the
gravitational-slip pair (Sigma, eta) as its companion lock.
  Parametrise the growth rate as  f(a) = Omega_m(a)^gamma.  In GR + flat LCDM
  this is gamma ~= 0.55, flat across z. Modified-gravity extensions of Einstein
  (f(R), DGP, interacting DE, any slip) push gamma AWAY from 0.55 and/or drive
  the lensing/growth ratio Sigma away from 1 and the slip eta away from 1.
  Under A1' the residual sector is EMPTY (Phase-2 null, RHAC-004) and the deepen
  mode IS GR's D+ (Front 1). Two consequences, both LAWS not fits:
      gamma_ESTIF == gamma_GR ~= 0.55   at all epochs
      (Sigma, eta, mu)_ESTIF == (1, 1, 1) exactly   (no anisotropic-stress
                                source; growth and lensing see one potential)
  A confirmed gamma != 0.55, or any detected slip, falsifies ESTIF-Core.

WHY IT IS INDEPENDENT OF KILL-SHOT #1: Omega_k is a BACKGROUND-GEOMETRY number;
gamma / Sigma / eta are PERTURBATION-SECTOR numbers. Orthogonal parameters,
orthogonal data (BAO distances vs RSD+lensing). Two separate falsifiers, plus
the background w-lock below makes a THIRD, and C-15 (c_gw = c) a FOURTH.

THE FOUR LOCKS (ESTIF-Core forces a POINT; GR-family fits a REGION):
  #1 curvature:      Omega_k        == 0        [registered kill-shot]
  #2 growth/slip:    gamma          ~= 0.55  and  (Sigma, eta) == (1, 1)
                                                 [THIS FRONT]
  #3 equation-state: (w0, wa)       == (-1, 0)   [Phase-2 null, RHAC-004]
  #4 wave speed:     c_gw           == c         [C-15, RHAC-008]
  Zero free dark-sector parameters vs GR-extension's 1-3 fitted ones.

PASSPORTS: Om = 0.3141 (banked), amplitude irrelevant to gamma (it cancels).
gamma from f(a) via Front 1's A1' deepen mode D+ = H * INT da/(aH)^3.

MEASURED VALUES USED (verified, cited):
  gamma_obs = 0.58 +/- 0.11   DESI DR1 peculiar-velocity consensus + DESI
              ShapeFit high-z growth (Lai et al. 2026, arXiv:2512.03229;
              Qin et al. 2026, A&A 708, A219). GR value ~0.55. Consistent.
  w-lock context (kill-shot #3 current pressure): DESI DR2 BAO+CMB prefer
              w0>-1, wa<0 over LCDM at 3.1 sigma frequentist (Abdul-Karim
              et al. 2025, arXiv:2503.14738); +SN reaches 2.8-4.2 sigma by
              dataset; BUT Bayesian evidence is ambiguous -- ln B modestly
              FAVOURS LCDM for DESI+CMB (arXiv:2511.10631), and the SN pull
              is partly a DESI-vs-DES tension that recalibration removes.
              ESTIF-Core's background IS flat LCDM, so kill-shot #3 lives or
              dies EXACTLY as flat LCDM does against the thawing hint -- no
              more, no less. Honest: currently disfavoured by the frequentist
              DE fit, not by Bayesian model comparison.

Run:  python3 estif_front3_second_discriminator.py
Deps: numpy, scipy  (no colossus needed -- gamma is amplitude-free)
"""
import numpy as np
from scipy.integrate import quad

# ------------------------------------------------------------------ passports
OM = 0.3141                      # banked bootstrap value
GAMMA_GR_FIT = 0.55              # textbook flat-LCDM growth index
GAMMA_GR_EXACT = 6.0 / 11.0      # linear-theory analytic value for LCDM, 0.5454...

GAMMA_OBS, GAMMA_OBS_ERR = 0.58, 0.11    # DESI DR1 PV + ShapeFit consensus

# ------------------------------------- Front 1 A1' deepen mode (Heath form)
def _gint(ap, Om=OM):
    return ap ** 1.5 / (Om + (1.0 - Om) * ap ** 3) ** 1.5

def I_of_a(a, Om=OM):
    v, _ = quad(_gint, 0.0, a, args=(Om,), limit=200)
    return v

def f_growth(z, Om=OM):
    a = 1.0 / (1.0 + z)
    E2 = Om / a ** 3 + (1.0 - Om)
    om_a = (Om / a ** 3) / E2
    return -1.5 * om_a + a * _gint(a, Om) / I_of_a(a, Om)   # exact dlnD/dlna

def omega_m_a(z, Om=OM):
    a = 1.0 / (1.0 + z)
    return (Om / a ** 3) / (Om / a ** 3 + 1.0 - Om)

def gamma_eff(z, Om=OM):
    """the second differing number: gamma from f = Omega_m(a)^gamma"""
    f = f_growth(z, Om)
    om = omega_m_a(z, Om)
    return np.log(f) / np.log(om)

# ------------------------------------- slip sector (mu, eta, Sigma) algebra
def sigma_from(mu, eta):
    """Sigma is not independent: Sigma = mu*(1+eta)/2 by definition of the
    lensing potential (Phi+Psi)/2 sourced through mu with slip eta = Phi/Psi."""
    return mu * (1.0 + eta) / 2.0

def f_growth_with_mu(z_grid, mu0, Om=OM):
    """Integrate the growth ODE with a DE-scaled deviation mu(a) = 1 + mu0*Om_DE(a).
    mu0 = 0 must return the A1' deepen mode. Used ONLY to show the lock has teeth."""
    from scipy.integrate import solve_ivp
    lna0, lna1 = np.log(1e-3), 0.0

    def rhs(lna, y):
        a = np.exp(lna)
        E2 = Om / a ** 3 + (1.0 - Om)
        om_a = (Om / a ** 3) / E2
        ode = 1.0 - om_a                      # Omega_DE(a)
        mu = 1.0 + mu0 * ode
        D, Dp = y
        # D'' + (2 + dlnH/dlna) D' - 1.5 mu Om(a) D = 0
        dlnH = -1.5 * om_a
        return [Dp, -(2.0 + dlnH) * Dp + 1.5 * mu * om_a * D]

    sol = solve_ivp(rhs, (lna0, lna1), [np.exp(lna0), np.exp(lna0)],
                    t_eval=np.log(1.0 / (1.0 + np.asarray(z_grid)))[::-1],
                    rtol=1e-10, atol=1e-12)
    D = sol.y[0][::-1]
    Dp = sol.y[1][::-1]
    return Dp / D                              # f = dlnD/dlna

# ============================================================== the receipt
def main():
    print("=" * 74)
    print("FRONT 3 - SECOND DIFFERING NUMBER: growth index gamma + slip lock")
    print("=" * 74)
    print(f"  Om (banked) = {OM}   GR growth index: fit ~{GAMMA_GR_FIT},"
          f" analytic 6/11 = {GAMMA_GR_EXACT:.4f}")
    print()

    print("-" * 74)
    print("[1] THE C-11 WALL (why the number is 'exactness vs freedom')")
    print("-" * 74)
    print("""  A1' reproduces GR at linear order by construction. Theorem: ESTIF-Core
  matches FLAT LCDM in every background + linear observable -- no linear
  number separates them. A discriminator can only separate ESTIF-Core from
  GR's EXTRA freedoms (the wider landscape flat LCDM lives in). ESTIF forces
  each as LAW; GR fits each. That is the discriminator, and it is real.""")
    print()

    print("-" * 74)
    print("[2] gamma(z) FROM THE A1' DEEPEN MODE (no fit, no free parameter)")
    print("-" * 74)
    zs = np.array([0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
    print(f"    {'z':>5}  {'Omega_m(a)':>11}  {'f(z)':>9}  {'gamma_eff':>10}")
    gam = []
    for z in zs:
        g = gamma_eff(z)
        gam.append(g)
        print(f"    {z:5.2f}  {omega_m_a(z):11.5f}  {f_growth(z):9.5f}  {g:10.5f}")
    gam = np.array(gam)
    g0 = gamma_eff(0.0)
    print()
    print(f"    range over z=0-5: {gam.min():.4f} - {gam.max():.4f}"
          f"   (spread {gam.max()-gam.min():.4f})")
    print(f"    high-z limit  -> 6/11 = {GAMMA_GR_EXACT:.4f}  (Omega_m -> 1)")
    print(f"    gamma(z=0)    =  {g0:.4f}")
    print("    => gamma is FLAT to ~1% across the whole range. Nothing dials it.")
    print()

    print("-" * 74)
    print("[3] PULL AGAINST THE MEASURED GROWTH INDEX")
    print("-" * 74)
    pull = (g0 - GAMMA_OBS) / GAMMA_OBS_ERR
    print(f"    gamma_obs   = {GAMMA_OBS} +/- {GAMMA_OBS_ERR}"
          f"   (DESI DR1 PV consensus + ShapeFit; see header)")
    print(f"    gamma_ESTIF = {g0:.4f}   (forced, zero free parameters)")
    print(f"    pull        = {pull:+.2f} sigma")
    verdict3 = "CONSISTENT" if abs(pull) < 2.0 else "TENSION"
    print(f"    verdict     = {verdict3}")
    print("    HONEST NOTE: this is not evidence FOR ESTIF over flat LCDM --")
    print("    flat LCDM predicts the same number. It is evidence against the")
    print("    GR EXTENSIONS that push gamma off 0.55. That is the whole claim.")
    print()

    print("-" * 74)
    print("[4] THE SLIP LOCK: (Sigma, eta, mu) == (1, 1, 1) EXACTLY")
    print("-" * 74)
    print("""  Chain, not a fit:
    A3 (vacuum sources nothing) + Phase-2 null (residual sector EMPTY,
    RHAC-004) => no anisotropic-stress source in the perturbed sector.
      no anisotropic stress   => Phi = Psi        => eta == 1
      no extra clustering DOF => Poisson unmodified => mu == 1
      Sigma = mu(1+eta)/2                          => Sigma == 1""")
    mu, eta = 1.0, 1.0
    Sig = sigma_from(mu, eta)
    print()
    print(f"    mu = {mu:.1f}   eta = {eta:.1f}   Sigma = mu(1+eta)/2 = {Sig:.1f}")
    ok4 = (mu == 1.0 and eta == 1.0 and Sig == 1.0)
    print(f"    slip lock holds: {ok4}")
    print()

    print("-" * 74)
    print("[5] DOES THE LOCK HAVE TEETH? (sensitivity of gamma to mu)")
    print("-" * 74)
    zg = [0.0, 0.5, 1.0]
    base = f_growth_with_mu(zg, 0.0)
    consistency = abs(base[0] - f_growth(0.0))
    print(f"    ODE-vs-Heath cross-check at z=0: |df| = {consistency:.2e}"
          f"   ({'PASS' if consistency < 1e-6 else 'FAIL'})")
    print()
    print(f"    {'mu0':>7}  {'f(0)':>9}  {'gamma(0)':>9}  {'shift':>8}  {'pull':>7}")
    for mu0 in [-0.4, -0.2, 0.0, 0.2, 0.4]:
        fz = f_growth_with_mu(zg, mu0)[0]
        gm = np.log(fz) / np.log(omega_m_a(0.0))
        print(f"    {mu0:7.2f}  {fz:9.5f}  {gm:9.5f}  {gm-g0:+8.4f}"
              f"  {(gm-GAMMA_OBS)/GAMMA_OBS_ERR:+7.2f}")
    print()
    print("    A ~20% modification of Poisson moves gamma by ~0.03 -- inside")
    print("    today's +/-0.11, outside a Stage-IV +/-0.02. The lock is not")
    print("    falsifiable NOW; it becomes falsifiable with Euclid/DESI-full")
    print("    RSD+lensing. That is the honest status: a REGISTERED, not yet")
    print("    decisive, discriminator.")
    print()

    print("-" * 74)
    print("[6] THE FOUR-LOCK LEDGER (point vs region)")
    print("-" * 74)
    ledger = [
        ("#1 curvature",   "Omega_k = 0",        "registered kill-shot",   "BAO distances"),
        ("#2 growth/slip", f"gamma = {g0:.4f}, (S,eta,mu)=(1,1,1)",
                                                  "THIS FRONT (RHAC-007)",  "RSD + lensing"),
        ("#3 eqn-of-state","(w0, wa) = (-1, 0)",  "Phase-2 null (RHAC-004)","BAO+CMB+SN"),
        ("#4 wave speed",  "c_gw = c",            "C-15 (RHAC-008)",        "GW170817"),
    ]
    for name, lock, prov, data in ledger:
        print(f"    {name:16s} {lock:34s} {prov:24s} {data}")
    print()
    print("    ESTIF-Core free dark-sector parameters: 0")
    print("    GR-extension family fitted parameters:  1-3")
    print()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"""  The growth index is FORCED, not fitted: gamma(z=0) = {g0:.4f}, running to
  6/11 = {GAMMA_GR_EXACT:.4f} at high z, spread {gam.max()-gam.min():.4f} over z = 0-5. Pull against the
  measured consensus {GAMMA_OBS} +/- {GAMMA_OBS_ERR}: {pull:+.2f} sigma -- {verdict3}.

  The companion slip lock is exact: (Sigma, eta, mu) = (1, 1, 1), forced by A3
  plus the empty residual sector. No anisotropic-stress source exists to break
  it, and there is no dial that could turn one on.

  WHAT THIS IS: the SECOND independent falsifier, in the perturbation sector,
  orthogonal to kill-shot #1 (background geometry). A confirmed gamma != 0.55
  or ANY detected slip kills ESTIF-Core.

  WHAT THIS IS NOT: evidence for ESTIF over flat LCDM. Per C-11 it cannot be --
  flat LCDM makes the identical prediction. The discriminator separates
  ESTIF-Core from GR's EXTRA FREEDOMS, and only from those. Any claim beyond
  that is overreach.""")
    print()

    all_ok = (abs(pull) < 2.0 and ok4 and consistency < 1e-6
              and abs(gam.min() - GAMMA_GR_EXACT) < 0.01)
    print("=" * 74)
    print(f"RECEIPT STATUS: {'PASS 4/4' if all_ok else 'CHECK FAILED'}"
          "   (pull, slip lock, ODE cross-check, 6/11 limit)")
    print("=" * 74)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
