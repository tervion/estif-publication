"""
PRINCIPLE P - DERIVATION ATTEMPT, AND THE OBSTRUCTION
====================================================
GOAL: derive Principle P  ( Omega_m = R_H / r_p )  from axioms A1'-A2-A3,
so that Omega_m graduates from 'postulated via P' to 'theorem'.

RESULT (negative, and sharp): P CANNOT be an axiom-derived law, because it is
NOT an identity that holds at all epochs. It holds at ONE epoch -- ours -- and
fails everywhere else. A law derivable from the three axioms would have to hold
always; P does not; therefore no such derivation exists. P is a TODAY-CONDITION
(a coincidence-of-our-epoch), of the same character as the cosmic-coincidence
problem -- not a consequence of A1-A3.

--------------------------------------------------------------------------
EQUIVALENT FORMS OF P (all algebraically identical)
--------------------------------------------------------------------------
  (P1)  Omega_m           = R_H / r_p            [matter fraction = Hubble
                                                  radius / particle horizon]
  (P2)  Omega_m * I        = 1,   I = INT_0^inf dz/E(z)   [the bootstrap]
  (P3)  (4 pi G / 3) rho_m r_p  = c H0 / 2       [matter-gravity at the
                                                  horizon = half cH0]
  (P4)  rho_m * r_p        = rho_crit * R_H       [matter x particle-horizon
                                                  = critical x Hubble-radius]
  Note (P3)/(P4): the CRITICAL-density version is AUTOMATIC --
      (4 pi G/3) rho_crit R_H == c H0 / 2   is just the Friedmann equation,
  true at every epoch. P is the claim that (rho_matter, particle-horizon)
  REPRODUCES that automatic critical identity. That reproduction requires
  Omega_m(a) * [r_p/R_H] = 1, which -- as the table shows -- holds only at a=1.

--------------------------------------------------------------------------
Run:  python3 estif_P_derivation_attempt.py
Deps: numpy, scipy
"""
import numpy as np
from scipy.integrate import quad

Om0 = 0.31408                 # bootstrap value (matter+Lambda)
OL = 1.0 - Om0

def E(a):
    return np.sqrt(Om0 * a ** -3 + OL)

def chi_p(a):                 # comoving particle horizon to a, units c/H0
    f = lambda ap: 1.0 / (ap ** 2 * E(ap))
    v, _ = quad(f, 1e-9, a, limit=400)
    return v

def main():
    print("=" * 74)
    print("PRINCIPLE P - DERIVATION ATTEMPT: does P hold at all epochs?")
    print("=" * 74)
    print("""  A law from A1-A3 must hold at every epoch. Test whether P's core
  equality -- instantaneous matter fraction Omega_m(a) = R_H(a)/d_p(a) --
  is an identity (holds always) or a coincidence (holds only now).
""")
    print(f"  {'a':>6}{'z':>9}{'Om_inst(a)':>13}{'R_H/d_p(a)':>13}{'mismatch':>12}")
    print("  " + "-" * 60)
    rows = [3.0, 2.0, 1.5, 1.0, 0.7, 0.5, 0.3, 0.1, 0.03, 0.01]
    for a in rows:
        z = 1.0 / a - 1.0
        Om_inst = Om0 * a ** -3 / E(a) ** 2
        RH_dp = 1.0 / (E(a) * a * chi_p(a))
        rel = abs(Om_inst - RH_dp) / Om_inst * 100
        tag = "MATCH" if rel < 2 else f"{rel:.0f}%"
        print(f"  {a:>6.2f}{z:>9.2f}{Om_inst:>13.4f}{RH_dp:>13.4f}{tag:>12}")
    print(f"""
  READING: both curves fall monotonically through cosmic history --
  Omega_m(a): 1 (high z) -> 0.31 (now) -> 0 (future)
  R_H/d_p(a): 1/2 (high z) -> 0.31 (now) -> 0 (future)
  They CROSS once, at a ~= 1. P is exactly 'we live at the crossing'.
  High-z analytic limit confirms divergence: Om_inst -> 1 while R_H/d_p ->
  1/2 (the EdS Schwarzschild-horizon value). NOT an identity.
""")
    print("=" * 74)
    print("VERDICT - WHAT THIS MEANS FOR THE DERIVATION PROGRAM")
    print("=" * 74)
    print("""  P is a TODAY-CONDITION, not an all-epoch law. Consequences:

  1. DERIVE-AS-LAW is impossible. You cannot derive from three time-symmetric
     axioms an identity that holds at only one time. The obstruction is
     structural, not a gap in cleverness.

  2. THREE ROUTES REMAIN, none of which is a derivation from A1-A3:
       (a) ATTRACTOR: add flow DYNAMICS (a 4th axiom) that drive the universe
           to, and observers/structure to peak at, the crossing epoch. This is
           a new research program, not a consequence of the current axioms.
       (b) ANTHROPIC: argue our epoch is selected. Legitimate but external;
           and it is the SAME unsolved cosmic-coincidence problem LCDM faces
           (why Om ~ OL now) -- P does not solve it, it restates it.
       (c) ACCEPT AS POSTULATE: keep P as a PREDICTIVE postulate. This is
           honest and defensible: P converts 'Om is whatever we measure' into
           'Om MUST be 0.314', a parameter-free number matching Planck at
           0.53 sigma. That predictive content is real and more than LCDM
           offers for Om -- but it is a successful postulate, NOT a theorem.

  3. HONEST RELABEL for the paper and the audit: 'Om derived from three
     axioms' OVERSTATES it twice over -- Om follows from P (shown earlier),
     and P itself is not derivable as a law (shown here). Correct statement:
     'Om is fixed by a predictive coincidence-principle P that holds at our
     epoch; deriving or dynamically explaining P is open (arguably the
     framework's deepest open problem).'

  4. PARAMETER-COUNT IMPACT: the 'one fewer fitted knob' claim survives only
     if P is granted. Strictly, P trades a data-fitted Om for a postulate-
     fixed Om -- a lateral epistemic move, redeemed by P's falsifiable
     prediction. Not a free reduction; a more predictive assumption.""")
    print("=" * 74)

if __name__ == "__main__":
    main()
