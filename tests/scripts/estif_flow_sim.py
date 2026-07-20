#!/usr/bin/env python3
"""
estif_flow_sim.py  —  the river-model test, both combination rules, honest.

PICTURE
-------
Space is a fluid flowing inward toward masses (drains). A star is a boat carried
by the current; the pull it feels is the material acceleration a = (v.grad)v —
how fast the current it rides speeds up as it drifts (Peter's tides/gradient choice).

  galaxy inflow (one drain):  v_gal = -sqrt(2GM/r) r_hat     (escape-velocity river)
  background ("the moment"):  v_bg  = +H r      r_hat        (de Sitter / Hubble flow)

TWO RULES (the whole point of "test both"):
  RULE A  velocities add:     v = v_gal+v_bg, then a=(v.grad)v      [flows INTERACT]
  RULE B  accelerations add:  a=(v_gal.grad)v_gal+(v_bg.grad)v_bg   [flows DON'T interact]
  A - B  =  the cross term  (v_gal.grad)v_bg + (v_bg.grad)v_gal.

For a single drain the algebra closes exactly (inward-positive):
  RULE A:  a = GM/r^2  +  (H/2)sqrt(2GM/r)  -  H^2 r
              |Newton|    |__cross (INWARD)_|    |deSitter (OUTWARD)|
  RULE B:  a = GM/r^2                        -  H^2 r

TESTS (the real ones, not "hit 0.221"):
  (1) does a bend appear?   (2) is its scale MASS-INDEPENDENT?  (make-or-break)
  (3) does it TRACK H?      (4) does its MAGNITUDE match observed a0 = 1.2e-10?

Nothing about a0 or any switch is inserted. numpy + matplotlib only. ~seconds.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G, Msun = 6.674e-11, 1.989e30
kpc = 3.0856775814913673e19; Mpc = 1000*kpc
c   = 2.99792458e8
H0  = 67.4*1000/Mpc
a0_obs = 1.2e-10

# ---- single-drain closed forms (inward positive) --------------------
def a_newton(r, M):     return G*M/r**2
def a_cross (r, M, H):  return 0.5*H*np.sqrt(2*G*M/r)          # extra INWARD term of rule A
def a_deSit (r, H):     return H**2*r                          # OUTWARD (repulsive)
def a_ruleA (r, M, H):  return a_newton(r,M) + a_cross(r,M,H) - a_deSit(r,H)
def a_ruleB (r, M, H):  return a_newton(r,M)                   - a_deSit(r,H)

def cross_scale(M, H):
    """radius & accel where rule-A's EXTRA INWARD term equals Newton (its natural scale)."""
    r_c = (np.sqrt(2*G*M)/H)**(2/3)                            # 2^{1/3}(GM)^{1/3}H^{-2/3}
    return r_c, a_newton(r_c, M)

# =====================================================================
#  VALIDATION — 1-D radial finite differences vs closed form (fine grid)
# =====================================================================
def validate():
    print("="*70)
    print("VALIDATION: 1-D radial finite-difference vs closed form (single drain)")
    print("="*70)
    M, H = 1e11*Msun, H0
    r = np.logspace(np.log10(kpc), np.log10(3000*kpc), 40000)
    v = H*r - np.sqrt(2*G*M/r)                 # radial flow, outward +
    a_fd = np.abs(v*np.gradient(v, r))         # |material accel| by finite diff
    a_cf = np.abs(a_ruleA(r, M, H))            # closed form
    good = (r > 5*kpc) & (r < 2000*kpc)
    err  = np.median(np.abs(a_fd[good]-a_cf[good])/a_cf[good])
    print(f"  median |FD - closedform| / closedform = {err*100:.3f}%")
    print(f"  => the (v.grad)v formula is confirmed to <1%. Sweeps below use closed form.\n")

# =====================================================================
#  TEST 1 — single drain: rotation curve, acceleration profile, RAR
# =====================================================================
def test_single_drain():
    M, H = 1e11*Msun, H0
    r  = np.logspace(np.log10(0.3*kpc), np.log10(3000*kpc), 700)
    aN = a_newton(r, M); aC = a_cross(r, M, H); aD = a_deSit(r, H)
    aA = np.abs(a_ruleA(r, M, H)); aB = np.abs(a_ruleB(r, M, H))
    vN = np.sqrt(aN*r)/1e3; vA = np.sqrt(aA*r)/1e3

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.7))

    ax[0].loglog(r/kpc, vN, lw=2, label="Newton")
    ax[0].loglog(r/kpc, vA, lw=2, label="Rule A")
    ax[0].set(xlabel="radius [kpc]", ylabel="circular speed [km/s]",
              title="Rotation curve (M=1e11 Msun)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")

    ax[1].loglog(r/kpc, aN, lw=2, label="Newton GM/r^2")
    ax[1].loglog(r/kpc, aC, lw=1.5, label="cross term (extra INWARD)")
    ax[1].loglog(r/kpc, aD, lw=1.5, label="de Sitter (OUTWARD)")
    ax[1].loglog(r/kpc, aA, "k-",  lw=2.2, label="Rule A total")
    ax[1].loglog(r/kpc, aB, "k--", lw=1.5, label="Rule B total")
    ax[1].axhline(a0_obs, color="crimson", ls="-.", lw=1.2, label="observed a0")
    ax[1].set(xlabel="radius [kpc]", ylabel="acceleration [m/s^2]",
              title="Acceleration pieces\n(cross term is real & inward, but tiny)")
    ax[1].legend(fontsize=7.5); ax[1].grid(alpha=.3, which="both")
    ax[1].set_ylim(1e-16, 1e-9)

    ax[2].loglog(aN, aA, lw=2, label="Rule A")
    ax[2].loglog(aN, aB, "--", lw=2, label="Rule B")
    ln = np.array([aN.min(), aN.max()])
    ax[2].loglog(ln, ln, "k:", lw=1, label="Newton (slope 1)")
    ax[2].loglog(ln, np.sqrt(ln*a0_obs), "r:", lw=1, label="MOND (slope 1/2)")
    ax[2].axvline(a0_obs, color="grey", ls="-.", lw=.8)
    ax[2].set(xlabel="a_baryon = GM/r^2 [m/s^2]", ylabel="a_obs [m/s^2]",
              title="Radial Acceleration Relation")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3, which="both")

    fig.tight_layout(); fig.savefig("fig1_single_drain.png", dpi=110)
    plt.close(fig)

# =====================================================================
#  TEST 2 — mass independence   (a0 vs galaxy mass)
# =====================================================================
def test_mass_independence():
    H = H0
    M = np.logspace(8, 12.5, 40)*Msun
    a_c = np.array([cross_scale(m, H)[1] for m in M])
    slope = np.polyfit(np.log10(M/Msun), np.log10(a_c), 1)[0]
    fig, ax = plt.subplots(figsize=(6.8, 5))
    ax.loglog(M/Msun, a_c, "o-", label="Rule A extra-term scale a0(M)")
    ax.axhline(a0_obs, color="k", ls="-.", label="observed a0 (mass-INDEPENDENT)")
    ax.set(xlabel="galaxy mass [Msun]", ylabel="crossover accel a0 [m/s^2]",
           title=f"MASS-INDEPENDENCE TEST\nRule A: a0 proportional to M^{slope:+.2f}  (need 0.00  ->  FAILS)")
    ax.legend(fontsize=9); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig("fig2_mass_independence.png", dpi=110)
    plt.close(fig)
    return slope

# =====================================================================
#  TEST 3 — H tracking   (a0 vs expansion rate)
# =====================================================================
def test_H_tracking():
    M = 1e11*Msun
    H = H0*np.logspace(-0.6, 0.6, 30)
    a_c = np.array([cross_scale(M, h)[1] for h in H])
    slope = np.polyfit(np.log10(H), np.log10(a_c), 1)[0]
    fig, ax = plt.subplots(figsize=(6.8, 5))
    ax.loglog(H/H0, a_c, "s-", label=f"Rule A a0(H) prop. to H^{slope:.2f}")
    ax.loglog(H/H0, c*H, "r--", label="c*H  (what observed a0 tracks)")
    ax.set(xlabel="H / H0", ylabel="acceleration [m/s^2]",
           title="H-TRACKING TEST\na0 DOES move with the 'moment' (good) —\nbut scale is H*v_gal, not c*H (magnitude wrong)")
    ax.legend(fontsize=9); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig("fig3_H_tracking.png", dpi=110)
    plt.close(fig)
    return slope

# =====================================================================
#  TEST 4 — "right direction, wrong magnitude" + the c-vs-v_gal gap
# =====================================================================
def test_magnitude_gap():
    H = H0
    r = np.logspace(np.log10(kpc), np.log10(200*kpc), 400)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for M in (1e9, 1e10, 1e11, 1e12):
        frac = a_cross(r, M*Msun, H)/a_newton(r, M*Msun)*100
        ax[0].semilogx(r/kpc, frac, label=f"M={M:.0e} Msun")
    ax[0].axhline(100, color="k", ls=":", lw=1, label="100% (needed to flatten curves)")
    ax[0].set(xlabel="radius [kpc]", ylabel="extra inward pull / Newton  [%]",
              title="Extra inward pull is REAL but ~100x too weak\ninside galaxy radii")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3); ax[0].set_ylim(0, 30)

    M = np.logspace(8, 12.5, 30)*Msun
    a_sim = np.array([cross_scale(m, H)[1] for m in M])
    ax[1].loglog(M/Msun, a_sim, "o-", label="Rule A a0 (uses H*v_gal)")
    ax[1].axhline(a0_obs, color="k", ls="-.", label="observed a0")
    ax[1].axhline(c*H0/(2*np.pi), color="r", ls="--", label="c*H0/2pi (uses c)")
    ax[1].set(xlabel="galaxy mass [Msun]", ylabel="a0 [m/s^2]",
              title="Absolute scale: sim undershoots by ~1000x\nbecause a0 wants c, not v_gal")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig("fig4_magnitude_gap.png", dpi=110)
    plt.close(fig)

# =====================================================================
def main():
    validate()
    test_single_drain()
    sM = test_mass_independence()
    sH = test_H_tracking()
    test_magnitude_gap()

    M = 1e11*Msun
    r_c, a_c = cross_scale(M, H0)
    v_c = np.sqrt(2*G*M/r_c)
    frac30 = a_cross(30*kpc, M, H0)/a_newton(30*kpc, M)*100

    print("="*70)
    print("VERDICT   (both rules, de Sitter background, M=1e11 Msun reference)")
    print("="*70)
    print(f"(1) bend appears?         YES.")
    print(f"    - Rule A adds an extra INWARD term (right MOND-like direction).")
    print(f"    - Rule B (no cross term) adds only the repulsive de Sitter push: no boost.")
    print(f"(2) mass-independent?     NO.  a0 prop. to M^{sM:+.2f}  (need 0.00).  FAILS key test.")
    print(f"(3) tracks H?             YES. a0 prop. to H^{sH:+.2f}  -> a0 tied to the MOMENT,")
    print(f"                          not a fixed medium constant. Peter's framing survives.")
    print(f"(4) magnitude?            NO. extra pull is only {frac30:.1f}% of Newton at 30 kpc")
    print(f"                          (need ~100%); absolute a0={a_c:.1e} vs obs {a0_obs:.1e}")
    print(f"                          -> ~{a0_obs/a_c:.0f}x too small.")
    print("-"*70)
    print("THE ROOT CAUSE (and the real result):")
    print(f"  Rule A's extra pull is (H/2)*v_gal — H times the GALAXY inflow speed")
    print(f"  ({v_c/1e3:.0f} km/s here). Observed a0 ~ c*H uses the SPEED OF LIGHT.")
    print(f"  Gap = c/v_gal ~ {c/v_c:.0f}. To get real a0 the background flow must be ~c,")
    print( "  which happens only at the COSMIC HORIZON — not locally near a galaxy.")
    print("-"*70)
    print("SO: the horizon instinct in ESTIF is CORRECT, but it proves a0 is NON-LOCAL.")
    print("    a0 cannot be a purely local galaxy calculation. No simple combination")
    print("    rule + de Sitter background gives mass-independence, and all undershoot")
    print("    by ~1000x. The next honest move is a horizon-referenced background,")
    print("    not more N-body. Figures: fig1..fig4 .png")

if __name__ == "__main__":
    main()
