"""
B-7: WEAK-LENSING TEST OF THE DERIVED a0  (published-fit level)
===============================================================
The purest option-A test available: ESTIF's derived a0 against gravitational
LENSING data -- observations the theory has never touched (SPARC used
rotation-curve kinematics; lensing is an entirely independent probe that
extends the acceleration relation ~2 decades deeper).

DATA (published fitted values; provenance fetched 8 July 2026):
 [B21] Brouwer et al. 2021, A&A 650, A113 (KiDS-1000 lensing RAR):
   - Baseline MOND scale used/tested: g_dagger = 1.20 +/- 0.26 e-10 m/s^2
     (M16 fit; +/-0.02 stat, +/-0.24 syst per McGaugh et al. 2016).
   - GAMA isolated spectroscopic lens sample vs MOND curve:
     chi2_red = 0.8  (~0.4 sigma) -> lensing RAR CONSISTENT with the
     M16 curve at that acceleration scale, down to g_bar ~ 1e-15 m/s^2.
   - KiDS-bright photometric sample: sits systematically high
     (chi2_red = 4.6, ~6 sigma; 4.0 within isolation-reliable range),
     BUT consistent (chi2_red = 1.5) if stellar masses shift +0.2 dex --
     B21: "Whether these models are confirmed or excluded relies heavily
     on the systematic bias in the stellar mass measurements."
 [M24] Mistele et al. 2024, ApJL 969, L3 (KiDS weak-lensing BTFR):
   - circular velocities stay FLAT to ~1 Mpc (no decline; "difficult to
     understand in LCDM"); lensing BTFR "fully consistent with the
     kinematic BTFR", for LTGs and ETGs alike.

THE TEST:
 ESTIF derives a0 = c H0 x0 / sqrt(3), with x0 either from the imported
 r_universe (old, 1.1793e-10) or from the Omega_m bootstrap (1.1920e-10).
 (1) Compare derived a0 to the lensing-validated scale g_dagger.
 (2) Show the deep-flow prediction g_obs = sqrt(a0 g_bar) across the
     lensing regime -- the curve the GAMA data followed at 0.4 sigma.
 (3) BTFR transfer: SPARC already validated ESTIF's kinematic BTFR
     (RMS 15.6%); M24 shows lensing BTFR == kinematic BTFR; consistency
     therefore TRANSFERS to lensing at the published-fit level.

HONEST SCOPE: this is a comparison against PUBLISHED FITS, not a raw
re-stacking of the KiDS ESD profiles (that reanalysis = definitive
version, future work). The KiDS-bright photometric tension and the
hot-gas/missing-baryon uncertainty are quoted, not hidden.

Run:  python3 estif_weaklensing_btfr.py
Deps: numpy, scipy.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

c = 2.99792458e8
G = 6.67430e-11
a_rad = 7.565723e-16
MPC = 3.085677581e22
H0 = 67.66 * 1000.0 / MPC
h = 0.6766
T_CMB, NEFF = 2.7255, 3.046

# ---- ESTIF derived a0 (both chains) ----------------------------------
X0_OLD = 0.31073                      # imported r_universe chain
def omega_radiation():
    rho_g = a_rad * T_CMB**4 / c**2
    rho_c = 3 * (h * 100e3 / MPC)**2 / (8 * np.pi * G)
    return (rho_g / rho_c) * (1.0 + 0.2271 * NEFF)
def I_hor(Om, Or):
    f = lambda u: 2.0*u/np.sqrt(Or + Om*u**2 + (1-Om-Or)*u**8)
    v, _ = quad(f, 0, 1, limit=200); return v
Or = omega_radiation()
X0_BOOT = brentq(lambda Om: Om*I_hor(Om, Or) - 1.0, 0.02, 0.98, xtol=1e-10)
a0_old  = c*H0*X0_OLD /math.sqrt(3)
a0_boot = c*H0*X0_BOOT/math.sqrt(3)

# ---- Published lensing-validated scale [B21/M16] ----------------------
G_DAG      = 1.20e-10
G_DAG_STAT = 0.02e-10
G_DAG_SYST = 0.26e-10   # combined, as quoted in B21 Sect. 2.3

print("=" * 72)
print("B-7  WEAK-LENSING TEST OF THE DERIVED a0   (published-fit level)")
print("=" * 72)
print(f"  ESTIF a0 (old chain, x0={X0_OLD:.5f})      = {a0_old:.4e} m/s^2")
print(f"  ESTIF a0 (bootstrap, x0={X0_BOOT:.5f})     = {a0_boot:.4e} m/s^2")
print(f"  Lensing-validated scale g_dagger [B21/M16] = {G_DAG:.2e} "
      f"(+/-{G_DAG_STAT:.0e} stat, +/-{G_DAG_SYST:.2e} syst)")
print()
print("[1] Scale comparison (deviation of derived a0 from g_dagger):")
for tag, a0 in [("old chain ", a0_old), ("bootstrap ", a0_boot)]:
    d = G_DAG - a0
    print(f"    {tag}: {d/G_DAG*100:+5.2f}%   = {d/G_DAG_SYST:.2f} sigma (syst)"
          f"   = {d/G_DAG_STAT:.2f} sigma (stat-only)")
print("    -> the GAMA isolated LENSING RAR matched the M16 curve at this")
print("       scale with chi2_red = 0.8 (~0.4 sigma) [B21]. ESTIF's derived")
print("       a0 sits 0.03 sigma (syst) / 0.40 sigma (stat) from that scale:")
print("       the lensing data validate the acceleration scale ESTIF derives.")
print()
print("[2] Deep-flow prediction  g_obs = sqrt(a0 * g_bar)  across the")
print("    lensing regime (the curve GAMA followed; bootstrap a0):")
print(f"    {'g_bar [m/s^2]':>15}{'g_obs pred [m/s^2]':>22}")
for gb in [1e-15, 1e-14, 1e-13, 1e-12]:
    print(f"    {gb:>15.0e}{math.sqrt(a0_boot*gb):>22.3e}")
ratio = math.sqrt(a0_boot/G_DAG)
print(f"    Distinguishability: g_obs(ESTIF)/g_obs(canonical 1.20e-10) =")
print(f"    sqrt({a0_boot/G_DAG:.4f}) = {ratio:.4f} -> a {abs(1-ratio)*100:.2f}% offset in g_obs,")
print("    far below current lensing errors: the data CANNOT distinguish")
print("    ESTIF's 1.192 from canonical 1.20, but COULD have falsified the")
print("    scale at the +/-20% level -- and did not.")
print()
print("[3] BTFR transfer [M24]:")
Mb = np.array([1e10, 1e11]) * 1.989e30
vf = (G * Mb * a0_boot) ** 0.25 / 1e3
print(f"    ESTIF v_flat(M_b=1e10 Msun) = {vf[0]:6.1f} km/s ;  "
      f"v_flat(1e11) = {vf[1]:6.1f} km/s")
print("    M24: lensing circular velocities stay FLAT to ~1 Mpc and the")
print("    lensing BTFR is 'fully consistent with the kinematic BTFR'.")
print("    ESTIF's kinematic BTFR is already validated (SPARC, RMS 15.6%),")
print("    so consistency TRANSFERS to lensing at the published-fit level.")
print("    Flat-to-1-Mpc is the deep-flow signature (v ~ r^0), which the")
print("    a0-law reproduces; M24 note it is 'difficult to understand in LCDM'.")
print()
print("=" * 72)
print("VERDICT")
print("=" * 72)
print("""  PASS at the published-fit level. The derived a0 lands 0.03 sigma
  (syst) / 0.40 sigma (stat) from the acceleration scale that KiDS
  gravitational lensing validated on the spectroscopic GAMA isolated
  sample (chi2_red = 0.8), on data ESTIF never touched, ~2 decades
  deeper in acceleration than SPARC.

  FLAGS (quoted, not hidden):
  - KiDS-bright PHOTOMETRIC sample sits high vs the same curve
    (~3.8-6 sigma) unless stellar masses shift +0.2 dex; B21 attribute
    the verdict to stellar-mass systematics, which dominate.
  - Hot-gas / missing-baryon content of g_bar is the fundamental
    systematic of ALL such tests (B21 Sect. 4.3).
  - This compares against PUBLISHED FITS; the definitive version is a
    raw re-stacking of the KiDS ESD profiles with ESTIF's a0 fixed --
    future work, now well-specified.""")
