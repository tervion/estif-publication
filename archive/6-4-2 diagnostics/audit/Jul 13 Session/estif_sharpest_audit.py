"""
THE SHARPEST AUDIT - ESTIF 6.4.0 PARAMETER HONESTY LEDGER
========================================================
QUESTION (three, answered together):
  (a) How many hardcoded, externally imported parameters does ESTIF use?
  (b) Does ESTIF still use LCDM as a knob?
  (c) What does A1' produce if it refuses to borrow LCDM -- keeping ONLY
      the dimensional anchor(s) plus the three axioms + Principle P?

METHOD: strip the colossus 'planck18' baseline out entirely. Re-derive every
number ESTIF claims as its own from scratch, tag each line NATIVE (falls out
of axioms + P + H0) or IMPORT (no ESTIF mechanism produces it). Count the
imports. Reproduce the bootstrap from first code, in BOTH modes, so the
input-dependence of the crown-jewel a0 number is explicit, not buried.

AXIOMS (v6.4.0): A1' even-on-average slices (Omega_k == 0 by law) - A2
universal speed c - A3 vacuum sources nothing. Plus the load-bearing extra:
  PRINCIPLE P:  Omega_m = R_H / r_p(Omega_m)   [matter fraction = Hubble
                radius / particle horizon] -> closed transcendental bootstrap.
  *** P IS NOT DERIVED FROM A1-A3. *** It is an added postulate (repo's own
  flag: estif_omega_bootstrap.py, "Part B / RHAC-H"). So "Omega_m derived"
  means "derived FROM P", not "derived from the three axioms". This audit
  treats P as a fourth principle and says so in every relevant line.

Run:  python3 estif_sharpest_audit.py
Deps: numpy, scipy  (NO colossus -- that is the whole point)
"""
import numpy as np
import math
from scipy.integrate import quad
from scipy.optimize import brentq

# ------------------- ONLY these dimensional inputs are permitted ------------
c     = 2.99792458e8          # universal constant (not a parameter)
G     = 6.67430e-11           # universal constant (not a parameter)
MPC   = 3.085677581e22
H0_kms = 67.66                # THE dimensional anchor (import #1)
h     = H0_kms / 100.0
H0    = H0_kms * 1e3 / MPC
# radiation inputs -- measured constants of nature, needed ONLY for the
# sub-percent bootstrap mode (import #2, #3 if that mode is used):
T_CMB = 2.7255
NEFF  = 3.046
a_rad = 7.565723e-16
# external target for comparison only (NOT an ESTIF input):
A0_MOND = 1.2e-10             # empirical, ~10% intrinsic scatter

# ================= NATIVE BLOCK 1 : the Principle-P bootstrap ===============
def Or_rad(h_=h, T=T_CMB, neff=NEFF):
    rho_g = a_rad * T ** 4 / c ** 2
    rho_c = 3 * (h_ * 100e3 / MPC) ** 2 / (8 * np.pi * G)
    return (rho_g / rho_c) * (1.0 + 0.2271 * neff)

def I_horizon(Om, Or=0.0):
    OL = 1.0 - Om - Or
    f = lambda u: 2.0 * u / np.sqrt(Or + Om * u ** 2 + OL * u ** 8)
    v, _ = quad(f, 0.0, 1.0, limit=200)          # = INT_0^inf dz/E(z), H0-free
    return v

def bootstrap(Or=0.0):
    return brentq(lambda Om: Om * I_horizon(Om, Or) - 1.0, 0.02, 0.98, xtol=1e-12)

# ================= NATIVE BLOCK 2 : growth index (A1' deepen mode) ==========
def _gint(ap, Om):
    return ap ** 1.5 / (Om + (1.0 - Om) * ap ** 3) ** 1.5
def _I(a, Om):
    v, _ = quad(_gint, 0.0, a, args=(Om,), limit=200); return v
def f_growth(z, Om):
    a = 1.0 / (1.0 + z)
    om = (Om / a ** 3) / (Om / a ** 3 + 1 - Om)
    return -1.5 * om + a * _gint(a, Om) / _I(a, Om)
def gamma_eff(z, Om):
    a = 1.0 / (1.0 + z)
    om = (Om / a ** 3) / (Om / a ** 3 + 1 - Om)
    return math.log(f_growth(z, Om)) / math.log(om)

# ============================================================ the audit
def main():
    print("=" * 76)
    print("THE SHARPEST AUDIT - ESTIF 6.4.0 running on axioms + P + H0, no planck18")
    print("=" * 76)

    print("-" * 76)
    print("[BOOTSTRAP] Omega_m from Principle P -- reproduced from first code")
    print("-" * 76)
    Om1 = bootstrap(0.0)
    Or  = Or_rad()
    Om2 = bootstrap(Or)
    a0_1 = c * H0 * Om1 / math.sqrt(3)
    a0_2 = c * H0 * Om2 / math.sqrt(3)
    ru1  = I_horizon(Om1) * c / H0
    ru2  = I_horizon(Om2, Or) * c / H0
    print(f"  MODE A  matter+Lambda, ZERO inputs beyond H0:")
    print(f"    Omega_m = {Om1:.5f}   (Planck 0.3111+/-0.0056: {(Om1-0.3111)/0.0056:+.2f} sigma)")
    print(f"    a0      = {a0_1:.4e}   ({abs(a0_1/A0_MOND-1)*100:.2f}% from MOND)")
    print(f"    r_univ  = {ru1:.4e} m  (old hardcode 4.4e26)")
    print(f"  MODE B  + radiation, needs T_CMB + N_eff (measured constants):")
    print(f"    Omega_m = {Om2:.5f}   (Planck: {(Om2-0.3111)/0.0056:+.2f} sigma)  <- the 0.3141 headline")
    print(f"    a0      = {a0_2:.4e}   ({abs(a0_2/A0_MOND-1)*100:.2f}% from MOND)  <- the 0.66% headline")
    print(f"    r_univ  = {ru2:.4e} m  ({(ru2/4.4e26-1)*100:+.2f}% vs hardcode)")
    print(f"  HONEST: the sub-percent a0 (MODE B) is NOT zero-input -- it needs")
    print(f"  two measured radiation constants. The truly parameter-free number")
    print(f"  (MODE A) sits at {abs(a0_1/A0_MOND-1)*100:.1f}%. Both are NATIVE (not fits); neither is")
    print(f"  0 free parameters. And BOTH stand on Principle P, not on A1-A3.")
    print()

    Om = Om2                     # adopt MODE B for downstream (the repo's choice)
    print("-" * 76)
    print(f"[DERIVED DOWNSTREAM] everything below uses the bootstrap Om={Om:.4f} only")
    print("-" * 76)
    print(f"  Omega_Lambda = 1 - Om - Or = {1 - Om - Or:.5f}   (NATIVE: flatness law)")
    print(f"  Omega_k      = 0 exactly              (NATIVE: A1' law, kill-shot #1)")
    print(f"  Lambda value = OL * 3H0^2/8piG        (NATIVE given H0; NO dynamical")
    print(f"                                         mechanism -- cc problem intact)")
    gs = [gamma_eff(z, Om) for z in (0.0, 0.5, 1.0, 2.0)]
    print(f"  growth index gamma(z=0..2) = {min(gs):.4f}-{max(gs):.4f}  (NATIVE: A1' D+,")
    print(f"                                         pinned on GR 0.55, Front 3)")
    print(f"  H(z), D+(z), f(z)          = GR/flat-LCDM, independently DERIVED")
    print(f"  slip (Sigma, eta, mu)      = (1,1,1) exactly (NATIVE: empty residual")
    print(f"                                         sector, Phase-2 null)")
    print()

    print("-" * 76)
    print("[THE IMPORT LEDGER] what ESTIF CANNOT produce and must borrow")
    print("-" * 76)
    print("""  ESTIF has NO early-universe mechanism (no inflation, no fluctuation
  generation, no recombination physics of its own). So it must import the
  initial conditions and the linear transfer physics wholesale:
    IMPORT  amplitude  As / sigma8      -- primordial fluctuation power
    IMPORT  tilt       ns               -- primordial spectral slope
    IMPORT  baryons    Omega_b          -- sets the transfer function shape
    IMPORT  transfer function / recombination (the planck18 baseline itself)
    IMPORT  MOND shape mu(x)            -- the a0 SCALE is native, but the
                                           RAR CURVE is not (flagged limitation;
                                           n(x) collapses g_obs/g_N at galaxy
                                           accelerations)
  ANSWER TO 'does ESTIF still use LCDM as a knob?': YES -- at the perturbation
  level, entirely. The BACKGROUND no longer borrows LCDM (Om is bootstrapped),
  but the primordial power spectrum + transfer function ARE the planck18
  baseline. Every Front script ran on it.""")
    print()

    print("-" * 76)
    print("[THE COUNT]")
    print("-" * 76)
    print(f"""  GENUINE IMPORTED PARAMETERS (framework, late-time predictions):
      1. H0                 dimensional anchor (sets all scales)
      2. As / sigma8        primordial amplitude       }} no ESTIF
      3. ns                 primordial tilt            }} mechanism
      4. Omega_b            baryon fraction (transfer) }} exists
    = 4 core imports.
    + 2 measured radiation constants (T_CMB, N_eff) IF the 0.66% bootstrap
      mode is used (MODE A avoids them at 3.75%).
    + 1 imported FUNCTION: the MOND interpolation shape mu(x).

  NATIVE OUTPUTS (computed, given Principle P + H0 [+ radiation]):
      Omega_m, Omega_Lambda, Omega_k(=0), Lambda's value, r_universe,
      a0 (acceleration scale), H(z), D+(z), f(z), gamma(~0.55), (Sigma,eta,mu).

  NOT framework parameters (test-only / universal): r_d = 147.09 Mpc is a
      Planck-calibrated import OF THE BAO TEST, not of ESTIF; delta_c,
      Barkana-Loeb thresholds, Salpeter time are application constants; c, G.""")
    print()

    print("=" * 76)
    print("VERDICT - THE HONEST BOTTOM LINE")
    print("=" * 76)
    print(f"""  vs standard flat LCDM (6-param base: H0, omega_b, omega_c, As, ns, tau):
  ESTIF-Core REMOVES ONE fitted cosmological knob -- omega_c / Omega_m
  graduates from a Planck fit to a computed number ({Om1:.4f} zero-input /
  {Om2:.4f} with radiation) -- and DEMOTES curvature from 'fitted ~= 0' to
  '0 by law'. It does not compute the CMB, so tau does not enter its late-time
  calc. Net: one fewer fitted knob, plus a genuine extra prediction (a0) that
  LCDM has no analogue for.

  THE PRICE, stated plainly:
    * The reduction hinges ENTIRELY on Principle P, which is NOT derived from
      A1-A3 -- it is a fourth postulate. 'Zero-input Omega_m' is really
      'Omega_m from one unproven principle'. Deriving P is the real open work.
    * The crown-jewel 0.66% a0 needs radiation inputs; parameter-free it is
      3.75%. Honest either way, but not both 'sub-percent' AND 'zero-input'.
    * P's root ({Om1:.4f}) overlaps Gaztanaga's causal-universe scale (~0.3176);
      novelty is NOT established until that comparison is done (repo flag).
    * a0's empirical target carries ~10% scatter -- 0.66% is directionally
      pleasing, not decisive.
    * Structure (amplitude, shape, transfer) is still 100% imported LCDM;
      ESTIF earns the geometry and the a0 scale, not the initial conditions.

  So: ESTIF-Core imports ~4 parameters (+2 radiation constants, +1 function),
  derives what LCDM fits for Om and curvature, and its one irreducible native
  win over LCDM is the a0 acceleration scale. Everything else at linear order
  is flat LCDM re-derived -- provided Principle P survives derivation.""")
    print("=" * 76)

if __name__ == "__main__":
    main()
