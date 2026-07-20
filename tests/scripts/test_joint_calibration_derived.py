"""
test_joint_calibration_derived.py

JOINT CALIBRATION WITH DERIVED PARAMETERS - the non-circular version.
=====================================================================
Created 15 July 2026. Companion to `test_joint_calibration.py`. Reads only;
writes no plots; modifies nothing. Run both and compare.

THE PROBLEM WITH THE FITTED VERSION
-----------------------------------
`test_joint_calibration.py` solves a 2-equation, 2-unknown system: it fits
N_MAX and B by fsolve so that
    equation 1:  EHT shadow tension  == 0 sigma exactly
    equation 2:  Lambda ratio        == 1.0 exactly
It then prints "EHT PASS (-0.000 sigma)", "Lambda PASS (ratio=1.0000)", and
"With no free parameters, simultaneously describes ...".

Those two PASSes are not results. They are the fit targets. A 2-parameter solve
hitting its own 2 targets exactly cannot fail, and reporting that as agreement
is circular. The phrase "no free parameters" is false as the script stands:
there are exactly two, and they are fitted to the two numbers being "predicted".
Only the LISA leg (TEST 3) is a genuine output of that script, because nothing
was fitted to it.

WHAT THIS SCRIPT DOES INSTEAD
-----------------------------
Take N_MAX and B from the INDEPENDENT derivation chain -- the electron
connection, `test_electron_connection.py` / `test_multiplier_derivation.py`:

    L      = ln(r_e / l_P)          (electron radius over Planck length)
    B      = L / 3                  (derived)
    N_MAX  = (5/7) * L              (derived, CONDITIONAL ON x_c -- see caveat)

Nothing in that chain touches the EHT shadow or Lambda. Feed those numbers into
the SAME formula functions (copied verbatim from test_joint_calibration.py) and
report what comes out. Now EHT and Lambda are real tests: they could miss.

CAVEAT ON N_MAX (carry this, do not drop it)
--------------------------------------------
TEST_INDEX.md records `test_multiplier_derivation.py` as: "B = L/3 derived;
N_MAX = 5/7*L CONDITIONAL ON x_c". So B is a clean derivation and N_MAX is not
fully unconditional. This script is therefore a ZERO-PARAMETER test of the EHT
and Lambda sectors only in so far as that conditionality holds. It is still
strictly stronger than fitting both to the answers. State it this way in any
write-up; do not upgrade it to "no free parameters" without closing x_c.

WHAT THIS SCRIPT IS NOT
-----------------------
It does not defend the LISA leg. `delay = (Rs/c) * sqrt(beta)` is the source's
Schwarzschild light-crossing time scaled by sqrt(beta) at ISCO. It carries NO
distance term, so it is a source-side quantity, not a propagation effect -- which
means C-15 / RHAC-008 (c_gw = c exactly, |c_gw/c - 1| = 0 for arbitrary flow)
does NOT forbid it. But nothing in either script derives why that quantity is
something LISA would clock as a delay. It is unjustified, not falsified. Open
item; do not cite the LISA sigma as a prediction until it has a derivation.

Run:  python3 test_joint_calibration_derived.py
Deps: numpy, scipy   (no matplotlib -- deliberately writes no figures)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))

import numpy as np
from scipy.optimize import fsolve
import estif_ec_gr_constants as const

# ============================================================================
# Constants -- mirrored EXACTLY from test_joint_calibration.py
# ============================================================================

M_m87           = 6.5e9 * const.M_sun
Rs_m87          = 2 * const.G * M_m87 / const.c**2
r_photon        = 1.5 * Rs_m87
CURVATURE_LOCAL = Rs_m87 / r_photon

R_H             = const.c / const.H_0
r_universe      = 4.4e26
CURVATURE_COSM  = R_H / r_universe

LAMBDA_MEASURED = 1.1056e-52
SHADOW_OBS      = 42.0
SHADOW_ERR      = 3.0
DISTANCE_M87    = 16.8 * 3.086e22

N_MAX_OLD = 1.5940
B_OLD     = 3.7311

# ============================================================================
# Electron connection -- the independent derivation of N_MAX and B
# ============================================================================

def _get(name, fallback, label):
    """Prefer the repo's own constant; fall back to CODATA with a loud note."""
    if hasattr(const, name):
        return getattr(const, name), f"const.{name}"
    return fallback, f"CODATA fallback ({label}) -- NOT from estif_ec_gr_constants"

r_e, r_e_src = _get('r_e', 2.8179403262e-15, 'classical electron radius')
l_P, l_P_src = _get('l_P', 1.616255e-35, 'Planck length')

L_ELECTRON = np.log(r_e / l_P)
B_DERIVED     = L_ELECTRON / 3.0
N_MAX_DERIVED = (5.0 / 7.0) * L_ELECTRON

# ============================================================================
# Formula Functions -- copied VERBATIM from test_joint_calibration.py
# ============================================================================

def n_dynamic(curvature, N_MAX, B):
    return N_MAX * np.exp(-B * curvature)

def beta_val(curvature, N_MAX, B):
    n   = n_dynamic(curvature, N_MAX, B)
    val = curvature ** (2 * n)
    if np.isscalar(val):
        return 0.0 if val >= 1.0 else np.sqrt(1.0 - val)
    return np.where(val >= 1.0, 0.0, np.sqrt(np.maximum(0.0, 1.0 - val)))

def observable(curvature, N_MAX, B):
    return np.sqrt(beta_val(curvature, N_MAX, B))

def shadow_sigma(N_MAX, B):
    """EHT tension in sigma for given H4 parameters."""
    R_shadow_gr  = np.sqrt(27) * Rs_m87
    theta_gr_uas = (R_shadow_gr / DISTANCE_M87) * 206265 * 1e6
    curv         = CURVATURE_LOCAL
    obs          = observable(curv, N_MAX, B)
    theta_gr_pt  = (4 * const.G * M_m87) / (const.c**2 * r_photon)
    theta_estif  = theta_gr_pt * (1 + obs * (Rs_m87 / (2 * r_photon)))
    shadow_pred  = theta_gr_uas * (theta_estif / theta_gr_pt)
    return (shadow_pred - SHADOW_OBS) / SHADOW_ERR

def shadow_uas(N_MAX, B):
    R_shadow_gr  = np.sqrt(27) * Rs_m87
    theta_gr_uas = (R_shadow_gr / DISTANCE_M87) * 206265 * 1e6
    obs          = observable(CURVATURE_LOCAL, N_MAX, B)
    theta_gr_pt  = (4 * const.G * M_m87) / (const.c**2 * r_photon)
    theta_estif  = theta_gr_pt * (1 + obs * (Rs_m87 / (2 * r_photon)))
    return theta_gr_uas * (theta_estif / theta_gr_pt)

def lambda_ratio(N_MAX, B):
    """Lambda ratio for given H4 parameters."""
    obs  = observable(CURVATURE_COSM, N_MAX, B)
    lam  = (3.0 / R_H**2) * obs**2
    return lam / LAMBDA_MEASURED

def lisa(N_MAX, B):
    """Source-side timing quantity. NOT a propagation delay. Undefended."""
    M_gw    = 65 * const.M_sun
    Rs_gw   = 2 * const.G * M_gw / const.c**2
    r_isco  = 3 * Rs_gw
    curv_gw = Rs_gw / r_isco
    d = (Rs_gw / const.c) * observable(curv_gw, N_MAX, B)
    return d, d / 1e-5

def equations(params):
    return [shadow_sigma(params[0], params[1]),
            lambda_ratio(params[0], params[1]) - 1.0]

# ============================================================================
# The receipt
# ============================================================================

def main():
    print("=" * 78)
    print("JOINT CALIBRATION -- DERIVED PARAMETERS (the non-circular version)")
    print("=" * 78)

    print("\n[0] PROVENANCE OF THE INPUTS")
    print("-" * 78)
    print(f"    r_e  = {r_e:.6e} m   <- {r_e_src}")
    print(f"    l_P  = {l_P:.6e} m   <- {l_P_src}")
    print(f"    L    = ln(r_e/l_P)   = {L_ELECTRON:.4f}")
    print(f"    B     = L/3          = {B_DERIVED:.4f}     [derived]")
    print(f"    N_MAX = (5/7)*L      = {N_MAX_DERIVED:.4f}    [derived, CONDITIONAL on x_c]")
    print("    Neither input touches the EHT shadow or Lambda. EHT and Lambda")
    print("    below are therefore genuine tests -- they are free to miss.")
    if 'fallback' in r_e_src or 'fallback' in l_P_src:
        print("\n    !! WARNING: one or more inputs came from a CODATA fallback, not")
        print("       from estif_ec_gr_constants. Add them to the constants module")
        print("       so this receipt reads the repo's own numbers.")

    print("\n[1] THE FITTED PARAMETERS (what test_joint_calibration.py solves for)")
    print("-" * 78)
    fit = fsolve(equations, [N_MAX_DERIVED, B_DERIVED], full_output=False)
    N_FIT, B_FIT = float(fit[0]), float(fit[1])
    print(f"    fsolve -> N_MAX = {N_FIT:.4f}   B = {B_FIT:.4f}")
    print(f"    These are fitted to force EHT = 0 sigma and Lambda ratio = 1.0.")
    print(f"    Any PASS they produce is a tautology, not a result.")

    print("\n[2] FITTED vs DERIVED -- do the two chains agree?")
    print("-" * 78)
    dN = abs(N_MAX_DERIVED - N_FIT) / N_FIT * 100
    dB = abs(B_DERIVED - B_FIT) / B_FIT * 100
    print(f"    {'':10} {'fitted':>10} {'derived':>10} {'apart':>9}")
    print(f"    {'N_MAX':10} {N_FIT:10.4f} {N_MAX_DERIVED:10.4f} {dN:8.2f}%")
    print(f"    {'B':10} {B_FIT:10.4f} {B_DERIVED:10.4f} {dB:8.2f}%")
    print()
    print("    This is the substantive claim: two numbers fitted to EHT+Lambda")
    print("    land close to two numbers derived from r_e/l_P. That agreement is")
    print("    the finding -- NOT the 0.000 sigma the fitted run prints.")

    print("\n[3] THE HONEST NUMBERS (derived parameters, nothing fitted)")
    print("-" * 78)
    rows = [("FITTED (circular)", N_FIT, B_FIT),
            ("DERIVED ln(r_e/l_P)", N_MAX_DERIVED, B_DERIVED),
            ("OLD (uncalibrated)", N_MAX_OLD, B_OLD)]
    print(f"    {'':21} {'N_MAX':>9} {'B':>9} {'EHT sig':>9} {'shadow':>9} {'Lam ratio':>10}")
    for lab, N, B in rows:
        print(f"    {lab:21} {N:9.4f} {B:9.4f} {shadow_sigma(N,B):9.4f}"
              f" {shadow_uas(N,B):9.3f} {lambda_ratio(N,B):10.4f}")
    print(f"\n    EHT observed: {SHADOW_OBS} +/- {SHADOW_ERR} uas")

    s_der = shadow_sigma(N_MAX_DERIVED, B_DERIVED)
    r_der = lambda_ratio(N_MAX_DERIVED, B_DERIVED)
    eht_ok = abs(s_der) < 1.0
    lam_ok = abs(r_der - 1.0) < 0.01

    print("\n[4] LISA LEG -- reported, NOT claimed")
    print("-" * 78)
    d_der, snr_der = lisa(N_MAX_DERIVED, B_DERIVED)
    d_fit, snr_fit = lisa(N_FIT, B_FIT)
    print(f"    derived: {d_der*1e6:7.1f} us  (S/N {snr_der:5.1f})")
    print(f"    fitted : {d_fit*1e6:7.1f} us  (S/N {snr_fit:5.1f})")
    print(f"    archived ESTIF-Gravity fork README claimed 32.03 us / 3.2 sigma")
    print(f"    -> superseded by recalibration, by a factor {d_der/32.03e-6:.1f}x.")
    print("    NOT counted in the verdict below: the quantity (Rs/c)*sqrt(beta)")
    print("    has no derivation. C-15 does not forbid it (no distance term, so")
    print("    not a propagation effect), but nothing justifies it either.")

    print("\n" + "=" * 78)
    print("VERDICT -- derived parameters, EHT and Lambda only")
    print("=" * 78)
    print(f"\n    EHT M87* shadow:  {'PASS' if eht_ok else 'FAIL'}"
          f"  ({s_der:+.4f} sigma, predicted {shadow_uas(N_MAX_DERIVED,B_DERIVED):.3f} uas)")
    print(f"    Cosmological Lam: {'PASS' if lam_ok else 'FAIL'}"
          f"  (ratio {r_der:.4f}, {abs(r_der-1)*100:.2f}% off; gate is 1%)")
    print(f"    LISA:             NOT ASSESSED (undefended quantity)")
    print(f"""
    Read this straight. With N_MAX and B taken from the electron connection --
    nothing fitted, nothing touching the EHT measurement -- the M87* shadow
    lands at {s_der:+.4f} sigma. That is a real result: it was free to be
    anything, and it wasn't. The fitted run's "0.000 sigma" carries none of
    that weight, because 0.000 was what it solved for.

    Lambda comes out {abs(r_der-1)*100:.2f}% off, which {'passes' if lam_ok else 'MISSES'} the script's own 1% gate.
    That is a real residual, not a failure of nerve. Report it.

    HONEST SUMMARY: {int(eht_ok) + int(lam_ok)}/2 on the derived parameters, with the LISA leg
    unassessed pending a derivation. This is a WEAKER headline than the fitted
    script's "3/3, no free parameters" and a STRONGER claim, because it is not
    circular. Do not restore the old wording.""")

    print("\n" + "=" * 78)
    print("This script wrote no files and modified nothing.")
    print("=" * 78)
    return 0 if eht_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
