"""
JWST EARLY-STRUCTURE-FORMATION TEST - baseline + target specification
=====================================================================
This computes the LCDM baseline for the JWST "too massive too early" tension
and quantifies exactly how much ESTIF must enhance linear structure growth to
relieve it. It is a READY SPECIFICATION: run it to see the number ESTIF must
beat and how the required growth enhancement maps to massive-halo abundance.

KEY PHYSICS (why this is laptop-tractable, not N-body-blocked):
  The ABUNDANCE of massive halos at high z (the JWST observable) is set by
  LINEAR growth + a collapse threshold via the Press-Schechter/Sheth-Tormen
  mass function -- all analytic. N-body is needed only for the NONLINEAR
  internal structure of halos (rotation curves), a SEPARATE question. So the
  population-level JWST test can be done on a laptop, given one ESTIF input.

THE ONE ESTIF INPUT REQUIRED:
  D_ESTIF(z), the linear growth factor under the ESTIF field equation. This is
  the natural extension of the derived field equation (Task 4) to LINEAR
  PERTURBATIONS on an FRW background -- a perturbed Gauss-Codazzi calculation,
  hard but analytic (NOT compute-blocked). Until it is derived, this script
  treats the growth enhancement g = D_ESTIF/D_LCDM as a knob and reports the
  abundance response, so the target is explicit.

CORRECTED PATH ASSIGNMENT (honest):
  At z~9 the universe is matter-dominated (dark energy ~0.2% of the density),
  so BOTH the frozen eddy (Path One, H=LCDM) and the thawing (Path Two,
  negligible at z~9) give essentially LCDM's H(z). The growth boost therefore
  CANNOT come from the expansion history. It must come from ESTIF's MODIFIED
  GRAVITY (a stronger effective source for perturbations) -- i.e. the derived
  field-equation sector, not the dark-energy sector.

Run:  python3 estif_jwst_growth_spec.py
Deps: numpy, scipy, colossus  (pip install colossus --break-system-packages)
"""
import sys
import numpy as np
from scipy.optimize import brentq

try:
    from colossus.cosmology import cosmology
    from colossus.lss import mass_function
except ImportError:
    print("This script needs colossus:  pip install colossus --break-system-packages")
    sys.exit(1)

if not hasattr(np, "trapezoid"):        # numpy<2 compatibility
    np.trapezoid = np.trapz

cosmo = cosmology.setCosmology('planck18')
h = cosmo.H0 / 100.0
Om = cosmo.Om0
fb = cosmo.Ob0 / Om
rho_crit0 = 2.77536627e11               # h^2 Msun/Mpc^3
rho_m_com = Om * rho_crit0              # h^2 Msun/Mpc^3
D = lambda z: cosmo.growthFactor(z)


def _mf(Mmin_Msun, z, want):
    Mmin_h = Mmin_Msun * h
    lnM = np.linspace(np.log(Mmin_h), np.log(1e15), 400)
    M = np.exp(lnM)
    dndlnM = mass_function.massFunction(M, z, mdef='fof', model='sheth99',
                                        q_in='M', q_out='dndlnM')  # (h/Mpc)^3
    if want == 'n':
        return np.trapezoid(dndlnM, lnM)             # (h/Mpc)^3
    return np.trapezoid(M * dndlnM, lnM)             # h^2 Msun/Mpc^3


def n_cum(Mmin, z):   return _mf(Mmin, z, 'n')
def rho_coll(Mmin, z): return _mf(Mmin, z, 'rho')


def z_eff_for_boost(g, z):
    """Redshift at which LCDM growth equals g x D(z): emulates a growth boost."""
    target = g * D(z)
    if target >= D(0.0):
        return 0.0
    return brentq(lambda zz: D(zz) - target, 0.0, z, xtol=1e-4)


print("=" * 74)
print("JWST EARLY-STRUCTURE-FORMATION TEST - LCDM baseline + ESTIF target")
print("=" * 74)
print(f"  Cosmology: Planck18   h={h:.4f}  Om={Om:.4f}  f_b={fb:.4f}")
print()

print("-" * 74)
print("[1] LCDM BASELINE - cumulative comoving number density of halos")
print("    (the JWST observable: how many massive halos exist at high z)")
print("-" * 74)
print(f"  {'z':<6}{'n(>1e10)':>12}{'n(>1e11)':>12}{'n(>3e11)':>12}{'n(>1e12)':>12}")
print("  " + "-" * 56 + "   [physical Mpc^-3]")
zlist = [8.0, 9.1, 10.0, 12.0]
for z in zlist:
    vals = [n_cum(M, z) * h**3 for M in (1e10, 1e11, 3e11, 1e12)]
    print(f"  {z:<6.1f}{vals[0]:>12.2e}{vals[1]:>12.2e}{vals[2]:>12.2e}{vals[3]:>12.2e}")
print()
print("  Maximum available stellar mass density (eps=1 ceiling, physical):")
print(f"  {'z':<6}{'rho*(>1e10)':>16}{'rho*(>1e11)':>16}   [Msun/Mpc^3]")
for z in zlist:
    r10 = fb * rho_coll(1e10, z) * h**2
    r11 = fb * rho_coll(1e11, z) * h**2
    print(f"  {z:<6.1f}{r10:>16.2e}{r11:>16.2e}")
print()
print("  Reference JWST observation (Labbe+2023 / Boylan-Kolchin 2023):")
print("  the two most massive z~7.5-9.1 candidates (M* ~ 10^10.5-10^11) imply")
print("  a stellar-mass density that pushes the required star-formation")
print("  efficiency toward or beyond plausible values in LCDM -- i.e. LCDM")
print("  must convert an implausibly large fraction of baryons to stars.")
print()

print("-" * 74)
print("[2] THE TARGET - how a linear-growth boost raises massive-halo abundance")
print("    (growth enhancement g = D_ESTIF/D_LCDM at z = 9.1)")
print("-" * 74)
z0 = 9.1
print(f"  D_LCDM(9.1) = {D(z0):.4f}")
print(f"  {'g':<6}{'z_eff':>8}{'boost n(>1e11)':>16}{'boost n(>3e11)':>16}")
print("  " + "-" * 46)
for g in [1.00, 1.05, 1.10, 1.13, 1.20, 1.30, 1.50]:
    ze = z_eff_for_boost(g, z0)
    b11 = n_cum(1e11, ze) / n_cum(1e11, z0)
    b31 = n_cum(3e11, ze) / n_cum(3e11, z0)
    print(f"  {g:<6.2f}{ze:>8.3f}{b11:>16.2f}{b31:>16.2f}")
print()
print("  Inverse - growth enhancement needed for a given 'reservoir boost'")
print("  (boost in collapsed mass above 1e11 Msun at z = 9.1):")
for target in [2, 5, 10]:
    g_needed = brentq(lambda g: rho_coll(1e11, z_eff_for_boost(g, z0)) / rho_coll(1e11, z0) - target,
                      1.0, 3.0, xtol=1e-3)
    print(f"    {target:>2}x reservoir  <-  growth boost g = {g_needed:.3f}  "
          f"({(g_needed-1)*100:.0f}% enhancement of D at z=9.1)")
print()

print("-" * 74)
print("[3] HARD FILTER - sigma8/S8 (the low-z structure constraint)")
print("    A growth boost that PERSISTS to z=0 overproduces present-day structure.")
print("-" * 74)
sigma8_planck = 0.8111; err8 = 0.0060; S8_lens = 0.766
print(f"  Planck sigma8 = {sigma8_planck:.4f} +/- {err8:.4f};  KiDS lensing S8 ~ {S8_lens}")
print("  If the growth enhancement g at z=9 PERSISTS unchanged to z=0 (scale-indep):")
print(f"  {'g':<6}{'sigma8(z=0)':>14}{'sigma above Planck':>20}")
for g in [1.05, 1.13, 1.20]:
    s8 = sigma8_planck * g; nsig = (s8 - sigma8_planck) / err8
    print(f"  {g:<6.2f}{s8:>14.3f}{nsig:>18.0f}")
print("  -> a persistent 13% boost gives sigma8 ~ 0.92 (~18 sigma above Planck) and")
print("     the WRONG sign for lensing (which mildly prefers SUPPRESSED low-z growth).")
print("  REQUIREMENT: D_ESTIF(z) enhancement must be TRANSIENT (high-z only, -> 1 by")
print("  z <~ 2) and/or SCALE-DEPENDENT (small-scale/high-k only, leaving the 8 Mpc/h")
print("  scale that sets sigma8 ~ unchanged). This is a TWO-SIDED constraint: enough")
print("  early growth for JWST, ~standard late growth for sigma8. It sharpens the")
print("  D_ESTIF target -- more specific, not merely larger.")
print()
print("=" * 74)
print("SPECIFICATION SUMMARY")
print("=" * 74)
print("""  OBSERVABLE TO PREDICT:
    cumulative comoving number density n(>M_halo, z) of massive halos at
    z ~ 8-12, equivalently the maximum stellar-mass density; compared to the
    JWST massive-galaxy candidates.

  LCDM NUMBER TO BEAT (z = 9.1, physical):
    n(>1e11 Msun) = {:.2e} Mpc^-3     n(>3e11 Msun) = {:.2e} Mpc^-3

  ESTIF INPUT REQUIRED (the one missing piece):
    D_ESTIF(z), the linear growth factor from the ESTIF field equation
    linearized on an FRW background (perturbed Gauss-Codazzi). Analytic,
    laptop-tractable, NOT N-body. This is the natural extension of the
    derived field equation (Task 4) to linear perturbations.

  TARGET (two-sided):
    D_ESTIF(9.1) must exceed D_LCDM(9.1) by ~13% to yield a 5x reservoir
    boost at high z (relieving JWST), WHILE the enhancement must be transient
    or scale-dependent so that sigma8(z=0) stays ~Planck (a persistent 13%
    boost gives sigma8 ~ 0.92, ~18 sigma excluded). So the target is a
    high-z-localised or small-scale growth enhancement, not a monotonic one.

  WHAT IS LAPTOP-TRACTABLE vs BLOCKED:
    Laptop:  D_ESTIF(z) derivation (perturbed field eq) + spherical-collapse
             threshold delta_c + Sheth-Tormen mass function -> abundance.
    Blocked: nonlinear INTERNAL halo structure (rotation curves, delta~1e5) -
             that is the SEPARATE N-body question, not this test.

  PATH ASSIGNMENT (corrected, honest):
    The growth boost comes from ESTIF's MODIFIED GRAVITY (perturbation source),
    NOT the expansion history (identical to LCDM at z~9 for both paths). So the
    JWST test is a GRAVITY-SECTOR prediction tied to the derived field equation
    - closer to Path One's core than to Path Two's dark-energy work.
""".format(n_cum(1e11, z0) * h**3, n_cum(3e11, z0) * h**3))
print("=" * 74)
