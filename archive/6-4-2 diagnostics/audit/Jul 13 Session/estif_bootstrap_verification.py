"""
BOOTSTRAP VERIFICATION - the calculations behind the sharpest-audit numbers
==========================================================================
WHY THIS FILE EXISTS: in the parameter-audit turn I quoted the Principle-P
bootstrap numbers (the naive 0.304-0.322 cutoff spread; the two-mode cross-check
against the repo's exact functions: 0.3043 / 3.75% and 0.31408 / 0.66%) from
throwaway probes without surfacing the code. This receipt reproduces EVERY one of
those numbers so they are reviewable and runnable before any of it is written into
the docs. Nothing here is new physics; it is the transparency this turn owed.

WHAT PRINCIPLE P SAYS:  Omega_m = R_H / r_p   (matter fraction = Hubble radius /
particle horizon). Since r_p contains Omega_m, this closes into a self-consistency
equation  Omega_m * I(Omega_m) = 1,  I = r_p / R_H = INT_0^inf dz/E(z).  H0 cancels.

THE INVESTIGATION, IN THREE STEPS (as it actually happened):
  STEP 1  from-scratch, naive: integrate the particle horizon in z. Discover the
          root DEPENDS ON THE CUTOFF (z->inf vs z=1100) -> 0.304 to 0.322. This
          cutoff-sensitivity is what made me stop and go read the repo.
  STEP 2  repo-exact form: the substitution a=u^2 turns the horizon integral into
          INT_0^1 2u du / sqrt(Or + Om u^2 + OL u^8), and it ADMITS RADIATION (Or).
          Matter+Lambda (Or=0) = Case 1; +radiation = Case 2.
  STEP 3  cross-check: naive matter+Lambda (STEP 1, z->inf) MUST equal repo-exact
          matter+Lambda (STEP 2, Or=0). If they agree, the gap was purely the
          radiation inputs, and the 0.66% headline is the +radiation case.

Run:  python3 estif_bootstrap_verification.py
Deps: numpy, scipy   (no colossus)
"""
import numpy as np
import math
from scipy.integrate import quad
from scipy.optimize import brentq

# ---- only these dimensional inputs are permitted --------------------------
c       = 2.99792458e8
G       = 6.67430e-11
MPC     = 3.085677581e22
H0_kms  = 67.66
h       = H0_kms / 100.0
H0      = H0_kms * 1e3 / MPC
T_CMB   = 2.7255
NEFF    = 3.046
a_rad   = 7.565723e-16
A0_MOND = 1.2e-10          # empirical target (~10% intrinsic scatter), NOT an input

# =========================================================================
# STEP 1 - NAIVE FROM-SCRATCH: particle horizon integrated in z (matter+Lambda)
#          I(Om) = INT_0^zmax dz / sqrt(Om (1+z)^3 + (1-Om))
# =========================================================================
def I_naive(Om, zmax):
    f = lambda z: 1.0 / np.sqrt(Om * (1 + z) ** 3 + (1 - Om))
    v, _ = quad(f, 0.0, zmax, limit=400)
    return v

def root_naive(zmax):
    return brentq(lambda Om: Om * I_naive(Om, zmax) - 1.0, 0.05, 0.95, xtol=1e-10)

# =========================================================================
# STEP 2 - REPO-EXACT FORM: substitution a = u^2 gives
#          I(Om,Or) = INT_0^1  2u du / sqrt(Or + Om u^2 + OL u^8),  OL = 1-Om-Or
#          (this is the same particle horizon, reparametrised, admitting radiation)
# =========================================================================
def Or_rad(h_=h, T=T_CMB, neff=NEFF):
    rho_g = a_rad * T ** 4 / c ** 2
    rho_c = 3 * (h_ * 100e3 / MPC) ** 2 / (8 * np.pi * G)
    return (rho_g / rho_c) * (1.0 + 0.2271 * neff)

def I_repo(Om, Or=0.0):
    OL = 1.0 - Om - Or
    f = lambda u: 2.0 * u / np.sqrt(Or + Om * u ** 2 + OL * u ** 8)
    v, _ = quad(f, 0.0, 1.0, limit=200)
    return v

def root_repo(Or=0.0):
    return brentq(lambda Om: Om * I_repo(Om, Or) - 1.0, 0.02, 0.98, xtol=1e-12)

# ---- helpers --------------------------------------------------------------
def a0_of(Om):        return c * H0 * Om / math.sqrt(3)      # = c^2/(r_univ sqrt3)
def offset(Om):       return abs(a0_of(Om) / A0_MOND - 1) * 100
def sigma_planck(Om): return (Om - 0.3111) / 0.0056
def r_univ(I):        return I * c / H0

# ============================================================ report
def main():
    print("=" * 74)
    print("BOOTSTRAP VERIFICATION  (H0 = %.2f km/s/Mpc, c/H0 = %.4e m)" % (H0_kms, c / H0))
    print("=" * 74)

    print("\nSTEP 1 - naive from-scratch (matter+Lambda), root depends on cutoff:")
    r_inf  = root_naive(np.inf)
    r_1100 = root_naive(1100.0)
    print(f"  z -> inf :  Om = {r_inf:.6f}   a0 off {offset(r_inf):.2f}%   "
          f"r_univ = {r_univ(I_naive(r_inf, np.inf)):.4e} m")
    print(f"  z = 1100 :  Om = {r_1100:.6f}   a0 off {offset(r_1100):.2f}%")
    print(f"  => cutoff-sensitivity 0.304-0.322: THIS is why I stopped and read the repo.")

    print("\nSTEP 2 - repo-exact form (substitution a=u^2), admits radiation:")
    Or = Or_rad()
    Om1 = root_repo(0.0)          # Case 1: matter + Lambda, zero inputs
    Om2 = root_repo(Or)           # Case 2: + radiation (T_CMB, N_eff, h)
    print(f"  radiation density from T_CMB={T_CMB}, N_eff={NEFF}, h={h:.4f}:  Or = {Or:.4e}")
    print(f"  CASE 1  matter+Lambda (0 inputs):  Om = {Om1:.5f}   "
          f"a0 off {offset(Om1):.2f}%   {sigma_planck(Om1):+.2f} sigma vs Planck")
    print(f"  CASE 2  + radiation (3 inputs)  :  Om = {Om2:.5f}   "
          f"a0 off {offset(Om2):.2f}%   {sigma_planck(Om2):+.2f} sigma vs Planck")
    print(f"  CASE 2  r_universe back-predicted = {r_univ(I_repo(Om2, Or)):.4e} m "
          f"(hardcode 4.4e26, {(r_univ(I_repo(Om2,Or))/4.4e26-1)*100:+.2f}%)")

    print("\nSTEP 3 - cross-check: naive(z->inf, matter+Lambda) vs repo-exact(Or=0):")
    print(f"  naive  = {r_inf:.6f}")
    print(f"  repo   = {Om1:.6f}")
    print(f"  |diff| = {abs(r_inf - Om1):.2e}   -> {'AGREE' if abs(r_inf-Om1)<1e-4 else 'DISAGREE'}")
    print(f"  => the two parametrisations are the SAME horizon. The whole gap between")
    print(f"     my 3.75% and the repo's 0.66% was exactly the radiation inputs.")

    print("\n" + "=" * 74)
    print("WHAT THESE NUMBERS ESTABLISH (the honest summary the audit reported)")
    print("=" * 74)
    print(f"""  * The bootstrap MECHANISM is real: a unique, H0-free self-consistency root
    that back-predicts r_universe ({r_univ(I_repo(Om2,Or)):.3e} m vs the 4.4e26 hardcode).
  * The sub-percent a0 (0.66%, Om={Om2:.4f}) is the +RADIATION case -- it needs
    T_CMB, N_eff, h. The truly zero-input number (Om={Om1:.4f}) gives a0 off
    {offset(Om1):.1f}%. Sub-percent OR zero-input, not both.
  * Every quantity here is computed from {{H0, T_CMB, N_eff}} + the P closure.
    None of it derives P itself -- that remains the postulate the audit and the
    P-obstruction receipt address. These are the verification numbers only.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
