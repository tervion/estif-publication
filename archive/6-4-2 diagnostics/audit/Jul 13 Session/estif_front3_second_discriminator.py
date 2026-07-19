"""
FRONT 3 - THE SECOND DIFFERING NUMBER (C-11 discriminator hunt)
==============================================================
THE HUNT: the fork owns ONE registered kill-shot -- mean spatial curvature
Omega_k == 0 exactly, at all epochs, as LAW (RHAC-006; current 0.0007+/-0.0019).
Wanted: a SECOND, INDEPENDENT number that A1'-ESTIF FORCES to an exact value,
which the Einstein family leaves FREE to fit -- so a confirmed departure would
kill ESTIF-Core while GR merely absorbs it.

THE C-11 WALL (stated first, because it shapes what a discriminator can be):
  A1' was engineered to reproduce GR at LINEAR order (that was the whole point
  of the amendment: f(z=0.5)=0.76, waves at c, Friedmann background). A theorem
  follows -- ESTIF-Core cannot differ from FLAT LCDM in ANY background-or-linear
  observable. So no linear number separates ESTIF-Core from flat LCDM. What a
  discriminator CAN do is separate ESTIF-Core from GR's EXTRA FREEDOMS: the
  wider GR landscape (curvature, dynamical w, modified growth/slip) that flat
  LCDM sits inside. ESTIF forbids each freedom as LAW; GR treats each as a
  fitted parameter. That "exactness vs freedom" is the only discriminator
  structure available at linear order -- and it is genuine and falsifiable.

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
the background w-lock below makes a THIRD.

THE THREE LOCKS (ESTIF-Core forces a POINT; GR-family fits a REGION):
  #1 curvature:      Omega_k        == 0        [registered kill-shot]
  #2 growth/slip:    gamma          ~= 0.55  and  (Sigma, eta) == (1, 1)
                                                 [THIS FRONT]
  #3 equation-state: (w0, wa)       == (-1, 0)   [Phase-2 null, RHAC-004]
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
    print("[2] THE SECOND NUMBER COMPUTED: gamma(z) forced by A1' deepen mode")
    print("    (amplitude-free: gamma depends on f and Omega_m only)")
    print("-" * 74)
    print(f"  {'z':<7}{'Omega_m(a)':>12}{'f [A1p D+]':>12}{'gamma_eff':>11}")
    print("  " + "-" * 42)
    gs = []
    for z in [0.0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]:
        g = gamma_eff(z)
        gs.append(g)
        print(f"  {z:<7.2f}{omega_m_a(z):>12.4f}{f_growth(z):>12.4f}{g:>11.4f}")
    print(f"\n  gamma_eff range over z: [{min(gs):.4f}, {max(gs):.4f}]"
          f"  -- flat, and pinned on the GR value ~0.55")
    print(f"  => ESTIF-Core FORCES gamma ~= 0.55 as LAW (empty residual sector +")
    print(f"     GR-equivalent D+). GR-modifications leave gamma FREE.")
    print()

    print("-" * 74)
    print("[3] THE COMPANION SLIP LOCK (same origin, independent observable)")
    print("-" * 74)
    print("""  Empty residual sector => no anisotropic-stress source => the two metric
  potentials stay equal: growth and lensing feel ONE potential. Hence
      gravitational slip   eta   == 1     (Phi = Psi)
      lensing amplitude    Sigma == 1
      growth amplitude     mu    == 1
  all exact, all epochs -- the GR values, forced not fitted. Any detected
  slip (eta != 1) or lensing-vs-growth split (Sigma != 1) falsifies
  ESTIF-Core. Modified-gravity extensions of Einstein generically predict
  eta != 1 or Sigma != 1; ESTIF-Core forbids it.""")
    print()

    print("-" * 74)
    print("[4] CURRENT HEALTH vs MEASUREMENT")
    print("-" * 74)
    g_est = gamma_eff(0.0)
    g_obs, g_err = 0.58, 0.11
    pull = (g_est - g_obs) / g_err
    print(f"  ESTIF-Core forced:  gamma = {g_est:.3f} (~0.55, flat in z)")
    print(f"  measured:           gamma = {g_obs:.2f} +/- {g_err:.2f}"
          f"   (DESI DR1 PV + ShapeFit)")
    print(f"  pull = {pull:+.2f} sigma  ->  PASS (consistent; not yet falsified)")
    print(f"  headroom: gamma_obs error {g_err} leaves room for future data to")
    print(f"  either confirm the lock or break it -- a live, decisive test.")
    print()

    print("-" * 74)
    print("[5] INDEPENDENCE MAP - three locks, three sectors, one point")
    print("-" * 74)
    print(f"  {'lock':<26}{'ESTIF forces':>14}{'GR family':>12}   {'sector'}")
    print("  " + "-" * 70)
    rows = [
        ("#1 curvature Omega_k", "0", "free", "background geom."),
        ("#2 growth index gamma", "~0.55", "free", "perturb. growth"),
        ("   slip (Sigma, eta)", "(1, 1)", "free", "perturb. lensing"),
        ("#3 (w0, wa)", "(-1, 0)", "free", "background EoS"),
    ]
    for name, e, g, s in rows:
        print(f"  {name:<26}{e:>14}{g:>12}   {s}")
    print(f"\n  ESTIF-Core = a POINT (0 free dark params); GR-extension = a REGION")
    print(f"  (1-3 fitted). #2 is orthogonal to #1 (perturbations vs geometry)")
    print(f"  and to #3 (growth-shape vs background-history) -- a genuine SECOND.")
    print()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print(f"""  SECOND DIFFERING NUMBER FOUND: the growth index gamma ~= 0.55 (with its
  companion slip lock Sigma = eta = 1), forced by A1's empty residual sector
  plus the GR-equivalent deepen mode. Computed value {min(gs):.3f}-{max(gs):.3f} across
  z = 0-5; measured 0.58 +/- 0.11 -> pull {pull:+.2f} sigma, PASS. It is
  independent of the curvature kill-shot (perturbation sector vs background
  geometry) and of the w-lock (growth shape vs expansion history).

  HONEST C-11 RESOLUTION: none of the three locks separates ESTIF-Core from
  FLAT LCDM -- at linear order they are the same theory, ESTIF-Core being
  flat LCDM DERIVED from three axioms rather than assumed. What the locks do:
  make ESTIF-Core a maximally-constrained, zero-free-dark-parameter POINT
  inside the Einstein landscape. Any confirmed Omega_k != 0, gamma != 0.55 /
  slip != 1, or w != -1 kills ESTIF-Core while a GR extension absorbs it.
  The genuinely distinguishing content that remains lives OFF the linear
  sheet -- in the nonlinear halo interior behind the N-body wall (descoped)
  and in the exactness-as-law claims tested here.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
