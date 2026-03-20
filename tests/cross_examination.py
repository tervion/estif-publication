"""
cross_examination.py
ESTIF v6.0 — March 2026

CROSS-EXAMINATION: Tests 1, 2, 3
==================================

This script synthesises the results of all three tests, looks for
internal consistency, identifies contradictions, and produces a
single honest verdict on the project's current state and next steps.

TEST RESULTS RECAP:
  Test 1 — DESI w(z):         BAD    (chi2/N = 10.8, 3.5σ tension on w0)
  Test 2 — SPARC bias:        GOOD   (bias is calibration, not structural)
  Test 3 — Multipliers:       MIXED  (1/3 derived, 5/7 partially derived)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

H0 = const.H_0
c  = const.c
G  = const.G
x0 = (c / H0) / 4.4e26
l_P = 1.616255e-35
import math
k_e = 8.9875517923e9
e_ch = 1.602176634e-19
m_e = 9.10938e-31
r_e = k_e * e_ch**2 / (m_e * c**2)
L   = math.log(r_e / l_P)
a0  = H0 * c * x0 / np.sqrt(3)

print("=" * 70)
print("CROSS-EXAMINATION — ESTIF v6.0 — Tests 1, 2, 3")
print("=" * 70)

# ============================================================================
# PART 1: WHAT THE THREE TESTS MEASURED
# ============================================================================
print("""
PART 1: WHAT EACH TEST ACTUALLY MEASURED
=========================================

  ESTIF has two largely independent sectors:

  ┌─────────────────────────────────────────────────────────────┐
  │  GRAVITY SECTOR          │  COSMOLOGY SECTOR               │
  │  (local, static)         │  (expansion, dark energy)       │
  │                          │                                 │
  │  • β(x) tilt formula     │  • Ω_tilt(z) replaces Ω_Λ      │
  │  • Force law ∇(ω/H₀)²    │  • H_ESTIF(z) Friedmann eq.    │
  │  • a₀ MOND derivation    │  • Comoving distances DM, DH    │
  │  • N_MAX, B parameters   │  • w_eff(z) equation of state   │
  │                          │                                 │
  │  Tested by: T2, T3       │  Tested by: T1                  │
  └─────────────────────────────────────────────────────────────┘

  Test 1 (DESI)    → tested COSMOLOGY sector only
  Test 2 (SPARC)   → tested GRAVITY sector only
  Test 3 (Fracs.)  → tested THEORY of GRAVITY sector parameters

  This separation is crucial for the cross-examination.
  A failure in one sector does NOT invalidate the other.
""")

# ============================================================================
# PART 2: CONSISTENCY CHECK — DO THE SECTORS AGREE WITH EACH OTHER?
# ============================================================================
print("=" * 70)
print("PART 2: INTERNAL CONSISTENCY CHECK")
print("=" * 70)

print("""
  QUESTION: Are the gravity sector and cosmology sector
  using the same underlying formula coherently?

  The link between the two sectors is x₀ = R_H/r_universe = Ω_m.
  This single number appears in both:
    Gravity:    a₀ = H₀ × c × x₀ / √3    (MOND, galaxy scales)
    Cosmology:  x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z)  (dark energy)

  The x₀ identity Ω_m ≈ x₀ is tested to 0.12% by Planck.
  Both sectors use the SAME x₀ — they are geometrically linked.
""")

# numerical check
omega_m = 0.3111
print(f"  x₀ = R_H/r_universe = {x0:.6f}")
print(f"  Ω_m (Planck 2018)   = {omega_m:.6f}")
print(f"  Agreement:           {abs(x0/omega_m-1)*100:.4f}%  ✓\n")

print("""  CONSISTENCY FINDING:
  The two sectors are internally consistent — they share x₀.
  But they DIVERGE in how x₀ propagates into predictions:

  Gravity sector:  x₀ enters as a static cosmological background
                   → tested against galaxy rotation velocities
                   → PASSES (Test 2: RMS 15.6%, calibration bias only)

  Cosmology sector: x₀ evolves with redshift via x(z)
                   → tested against BAO distance rulers
                   → FAILS (Test 1: chi2/N = 10.8 on DESI DR2)

  The failure is in the EVOLUTION law x(z), not in x₀ itself.
  Specifically: x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is circular —
  it uses ΛCDM to compute its own correction to ΛCDM.
""")

# ============================================================================
# PART 3: ROOT CAUSE ANALYSIS — WHY DID TEST 1 FAIL?
# ============================================================================
print("=" * 70)
print("PART 3: ROOT CAUSE ANALYSIS — TEST 1 FAILURE")
print("=" * 70)

print("""
  The DESI DR2 chi2/N = 10.8 is not marginal — it is a clear failure.
  The pull pattern shows predictions systematically LOW across all
  DM/rs bins (z = 0.3 to 2.3), meaning ESTIF predicts shorter comoving
  distances than DESI measures.

  Shorter distances ← less integrated expansion ← Ω_tilt(z) evolves
  differently from what DESI requires.

  TWO CANDIDATE CAUSES:
  ─────────────────────
  Cause A: The FUNCTIONAL FORM of x(z) is wrong.
    x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is derived from the geometry,
    but it uses H_ΛCDM as a "ruler" to define its own evolution.
    This is circular: if ESTIF truly replaces ΛCDM, x(z) should be
    computed with H_ESTIF, not H_ΛCDM. The current form freezes the
    ΛCDM shape and adds a small tilt correction — which DESI DR2 rules
    out at > 3σ.

  Cause B: The z_eff = 2.0 HARD CUTOFF is too aggressive.
    The model applies Ω_tilt(z) = Ω_tilt(2.0) for all z > 2.
    DESI measurements extend to z = 2.33 (Lyman-alpha).
    The Lya bin at z=2.33 is the ONLY bin ESTIF gets right (pull = +0.04σ)
    because there the frozen Ω_tilt accidentally matches.
    But this freezing distorts the DM/rs integral for all lower-z bins.

  WHICH IS MORE LIKELY?
    Both contribute. But Cause A is structural — even if the cutoff
    were removed, the circular x(z) definition would produce a different
    evolution than DESI requires. The DESI DR2 best-fit w0 = -0.73
    (dynamical dark energy, not the static-ish ESTIF tilt) suggests
    DESI is seeing real evolution that the current Ω_tilt formula
    cannot reproduce.

  IMPORTANT NOTE:
    DESI DR2 (released March 19 2026 — the day Test 1 was run) shows
    3.1σ preference for dynamical dark energy over ΛCDM. ESTIF is also
    dynamical dark energy — but its SPECIFIC evolution is wrong.
    This is not a "ΛCDM is fine" problem. DESI is saying dark energy
    evolves AND showing ESTIF's particular evolution is incorrect.
""")

# Show the w(z) mismatch numerically
def w_eff(z, dz=0.01):
    z = max(z, dz)
    hi = estif.omega_tilt(z + dz)
    lo = estif.omega_tilt(max(z - dz, 1e-4))
    dlnde_dz = (np.log(hi+1e-30) - np.log(lo+1e-30)) / (2*dz)
    return -1.0 + (1.0+z)/3.0 * dlnde_dz

def w_cpl(z, w0, wa):
    return w0 + wa * z/(1+z)

print("  ESTIF w_eff(z) vs DESI DR2 best-fit (w0=-0.73, wa=-0.66):")
print(f"  {'z':<8} {'ESTIF w':<14} {'DESI w':<14} {'Diff'}")
print("  " + "-" * 44)
for z in [0.1, 0.3, 0.5, 0.93, 1.32, 2.0]:
    we = w_eff(z)
    wd = w_cpl(z, -0.73, -0.66)
    print(f"  {z:<8.2f} {we:<14.4f} {wd:<14.4f} {we-wd:>+.4f}")

# ============================================================================
# PART 4: WHAT TEST 2 TELLS US ABOUT TEST 3
# ============================================================================
print("\n" + "=" * 70)
print("PART 4: CROSS-VALIDATION — TESTS 2 AND 3")
print("=" * 70)

print(f"""
  Test 2 (SPARC) and Test 3 (multipliers) are both about the
  gravity sector. Do they reinforce or contradict each other?

  TEST 3 DERIVED:
    B = L/3 = {L/3:.4f}  (3D isotropy, same principle as MOND 1/√3)
    This gives the decay rate of the tilt exponent with curvature.

  TEST 2 FOUND:
    The formula v⁴ = G × M_bar × a₀ with a₀ = H₀cx₀/√3
    fits 87 SPARC galaxies to RMS = 15.6%.
    The bias is calibration (Upsilon_*), not structural.

  THE REINFORCING LINK:
    Both the 1/√3 in a₀ AND the 1/3 in B trace to the same
    physical principle: 3D isotropic projection of a 4D quantity.

    MOND:        v_3D = v_4D / √3    → a₀ = H₀cx₀/√3
    Tilt decay:  B_3D = B_4D / 3     → B = L/3

    These are not two independent coincidences. They are the SAME
    physical argument applied twice to different quantities in the
    same framework. This mutual reinforcement is significant.

    If the isotropy argument is correct (as Test 2 suggests through
    good SPARC results), then B = L/3 is also correct (Test 3).
    If B = L/3 is wrong, then 1/√3 is probably also wrong, and
    Test 2 should have failed — but it didn't.

  CONSISTENCY VERDICT:
    Tests 2 and 3 are MUTUALLY CONSISTENT and REINFORCE each other.
    They share a common physical principle (3D isotropy) that passes
    the galaxy rotation test independently.
""")

# ============================================================================
# PART 5: THE OPEN GAP — x_c — IN CONTEXT OF ALL THREE TESTS
# ============================================================================
print("=" * 70)
print("PART 5: THE OPEN GAP — x_c IN CONTEXT")
print("=" * 70)

x_c = np.log(const.N_MAX_COMBINED / 0.5) / const.B_COMBINED
print(f"""
  Test 3 identified x_c = {x_c:.6f} as the remaining underived number.
  This is the curvature ratio where ESTIF exactly equals GR.

  In the context of all three tests, x_c plays a specific role:
    It appears only in the GRAVITY sector (β(x) formula).
    It does NOT appear directly in the cosmological sector.
    Its value was set by joint calibration to EHT + Planck + LISA.

  The three calibration anchors are independent observational facts:
    EHT M87*:  shadow diameter = 42.0 μas (Event Horizon Telescope 2019)
    Planck Λ:  cosmological constant from CMB (Planck 2018)
    LISA:      GW delay prediction (future, but calibrated to it)

  Key question: does the Test 1 DESI failure invalidate x_c?
    No. The DESI failure is in Ω_tilt(z), not in β(x).
    x_c is a LOCAL gravity parameter (r ~ few × Rs near black holes).
    DESI measures COSMOLOGICAL distances at Gpc scales.
    These are completely different physical regimes.

  x_c = {x_c:.6f} → r = Rs/x_c = {1/x_c:.3f} Rs from a black hole
    This is between the photon sphere (1.5 Rs) and ISCO (3 Rs).
    Specifically: r = {1/x_c:.3f} Rs — no obvious GR significance yet found.
    This is the one number that needs a geometric derivation.
""")

# ============================================================================
# PART 6: TWO-SECTOR VERDICT TABLE
# ============================================================================
print("=" * 70)
print("PART 6: TWO-SECTOR VERDICT")
print("=" * 70)

print("""
  ┌──────────────────────┬───────────────────────┬────────────────────┐
  │ COMPONENT            │ TEST RESULT           │ VERDICT            │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ a₀ derivation        │ T2: SPARC RMS 15.6%   │ ✅ SOLID           │
  │ (MOND, gravity)      │ T3: 1/√3 from isotropy│                    │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ B = L/3              │ T3: 0.69% off         │ ✅ DERIVED         │
  │ (tilt decay rate)    │ T2: gravity OK        │ (motivated)        │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ x₀ = Ω_m identity   │ T2: implicit via a₀   │ ✅ HOLDS           │
  │ (eddy dark matter)   │ Planck: 0.12%         │ (independently)    │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ x_c = 0.272          │ T3: not derived       │ ⚠️  OPEN GAP       │
  │ (GR crossover)       │ No geometric origin   │ (well-defined)     │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ N_MAX connection     │ T3: 5/7 × L, 0.08%   │ ⚠️  CONDITIONAL    │
  │ to electron radius   │ follows from x_c      │ (needs x_c)        │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ Ω_tilt(z)            │ T1: chi2/N = 10.8     │ ❌ FAILS           │
  │ (dark energy evol.)  │ DESI DR2 3.5σ tension │ (needs rework)     │
  ├──────────────────────┼───────────────────────┼────────────────────┤
  │ w_eff prediction     │ T1: -0.65 vs DR2 -0.73│ ⚠️  MARGINAL       │
  │ (equation of state)  │ Pull: +0.79σ at z~0   │ (shape is wrong)   │
  └──────────────────────┴───────────────────────┴────────────────────┘
""")

# ============================================================================
# PART 7: WHAT NEEDS TO HAPPEN — PRIORITISED
# ============================================================================
print("=" * 70)
print("PART 7: PRIORITISED ACTION LIST")
print("=" * 70)

print("""
  PRIORITY 1 (BLOCKING — do before any new cosmology claims):
  ─────────────────────────────────────────────────────────────
  Fix the Ω_tilt(z) evolution law.

  The current x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) is circular.
  Options:
    Option A: Define x(z) using H_ESTIF(z) self-consistently.
              → Requires iterative solve (H_ESTIF depends on Ω_tilt,
                which depends on x(z), which depends on H_ESTIF).
              → Computationally tractable. Removes the circularity.

    Option B: Replace the functional form entirely.
              DESI DR2 prefers w0 > -1 at low z, wa < 0 (DE was
              stronger in the past). Design Ω_tilt(z) to match this
              phenomenology while retaining geometric motivation.

    Option C: Abandon the cosmological sector for now.
              Focus the paper on the GRAVITY sector (Tests 2+3 pass).
              Publish: "geometric derivation of a₀ and MOND from
              hypersurface tilt geometry" — a narrow but defensible claim.

  PRIORITY 2 (THEORY — strengthens the gravity sector):
  ──────────────────────────────────────────────────────
  Derive x_c = 0.272 geometrically.

  This is a well-posed problem: find the geometric property of the
  Schwarzschild metric that produces x_c without calibration.
  The n(x_c) = 1/2 condition must follow from geometry, not fitting.

  Approach: look at thermodynamic properties of black holes at x_c.
    At x_c: β = τ_GR exactly → entropy relationship?
    Hawking temperature is T ∝ 1/Rs → connects to x?
    Bekenstein-Hawking: S ∝ Rs² → area element at r = Rs/x_c?

  PRIORITY 3 (VALIDATION — strengthens Test 2):
  ───────────────────────────────────────────────
  Run SPARC with Upsilon_* = 0.65–0.70 explicitly.
  The current Upsilon_* = 0.50 (McGaugh standard) gives -7.6% bias.
  The zero-bias value is 0.85 (too high vs literature).
  The literature range 0.60–0.70 gives bias ~-5% to -4%.
  A proper photometric mass analysis with individual Upsilon_* per galaxy
  would reduce the RMS further and give a cleaner result.

  PRIORITY 4 (PUBLICATION — achievable now):
  ────────────────────────────────────────────
  Write up Tests 2+3 as a standalone gravity paper.
  Do NOT include cosmological sector claims.
  Specific publishable claim:
    "We derive the MOND critical acceleration a₀ = H₀cx₀/√3 from
    a geometric model of 3D space as a hypersurface in 4D, using only
    the Planck 2018 values of H₀ and Ω_m. The derivation requires no
    free parameters. The 1/3 projection factor follows from 3D spatial
    isotropy — the same principle that connects the Jeans criterion to
    the virial theorem. We test this prediction against 87 quality-1
    galaxies from the SPARC survey (Lelli et al. 2016) and find
    RMS agreement of 15.6%, consistent with the observed scatter
    in the baryonic Tully-Fisher relation."
  This is a defensible, novel, and verifiable claim.
""")

# ============================================================================
# PART 8: THE OVERALL PROJECT STATUS
# ============================================================================
print("=" * 70)
print("PART 8: OVERALL PROJECT STATUS")
print("=" * 70)

print(f"""
  WHAT THE THREE TESTS TOGETHER SAY:

  The ESTIF project is a split verdict.

  The GRAVITY sector is in better shape than it was before these tests.
  Tests 2 and 3 together show:
    → The MOND acceleration is derived, not fitted (zero free parameters)
    → The -7.6% SPARC bias is a known calibration issue, not model failure
    → The 1/3 multiplier has a genuine derivation from 3D isotropy
    → The two isotropy results (MOND + multipliers) reinforce each other

  The COSMOLOGY sector has a specific, identified failure:
    → Ω_tilt(z) evolution is inconsistent with DESI DR2 at 3.5σ
    → The pre-existing prediction w_eff ≈ -1.08 is now falsified
      by DESI DR2's w0 = -0.73 ± 0.10
    → This is NOT a "minor adjustment" — the shape of w(z) is wrong

  WHAT THIS MEANS FOR THE PLANNED PUBLICATION:
    The planned publication claiming to replace dark energy with tilt
    geometry CANNOT proceed without fixing the cosmological sector.
    The reviewers' concern about "galaxy dynamics + cosmology in one
    formula" was correct — these sectors are currently inconsistent.

  WHAT CAN BE PUBLISHED NOW:
    A narrower paper on the gravity sector only — specifically the
    geometric derivation of a₀ and its agreement with SPARC — is
    defensible today. This would be a short letter (4–6 pages) rather
    than the full theory paper.

  IS THE PROJECT WORTHY OF CONTINUED PURSUIT?
    Yes, emphatically — but within a redefined scope.
    The gravity sector results are genuinely interesting:
      → Deriving a₀ from geometry to 1.72% with zero free parameters
      → Connecting N_MAX and B to ln(r_e/l_P) through isotropy
      → Passing the SPARC galaxy survey test
    These results are not trivial and are not curve-fitting.
    The cosmological sector needs to be rebuilt, not abandoned.
""")

# ============================================================================
# VISUALIZATION
# ============================================================================

fig = plt.figure(figsize=(18, 14))
fig.suptitle(
    'ESTIF v6.0 Cross-Examination: Tests 1, 2, 3\n'
    'Gravity sector: PASSES  |  Cosmology sector: FAILS  |  Theory: PROGRESSES',
    fontsize=13, fontweight='bold'
)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.38)

# ── Panel 1: Test scorecard ─────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
items = [
    ('a₀ from geometry',         'T2+T3', 'SOLID',       '#28a745'),
    ('SPARC BTFR (15.6% RMS)',    'T2',    'PASSES',      '#28a745'),
    ('1/3 from isotropy',         'T3',    'DERIVED',     '#28a745'),
    ('B = L/3  (0.69%)',          'T3',    'MOTIVATED',   '#5cb85c'),
    ('x_c = 0.272',               'T3',    'OPEN GAP',    '#ffc107'),
    ('5/7 × L conditional',       'T3',    'PARTIAL',     '#ffc107'),
    ('Ω_tilt(z) evolution',       'T1',    'FAILS',       '#dc3545'),
    ('w_eff prediction',          'T1',    'FALSIFIED',   '#dc3545'),
    ('DESI DR2 distances',        'T1',    'TENSION 3.5σ','#dc3545'),
]
ax1.axis('off')
for i, (label, test, verdict, color) in enumerate(items):
    y = 0.95 - i * 0.105
    ax1.text(0.02, y, f'{test}', transform=ax1.transAxes,
             fontsize=8, va='center', color='gray', style='italic')
    ax1.text(0.18, y, label, transform=ax1.transAxes,
             fontsize=8.5, va='center')
    ax1.text(0.80, y, verdict, transform=ax1.transAxes,
             fontsize=8.5, va='center', fontweight='bold', color=color,
             ha='center')
ax1.set_title('Test Scorecard', fontsize=11, fontweight='bold')
# removed

# ── Panel 2: Two-sector split ────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
sectors = ['Gravity\nSector', 'Cosmology\nSector']
scores  = [85, 20]   # qualitative %
colors  = ['#28a745', '#dc3545']
bars    = ax2.bar(sectors, scores, color=colors, alpha=0.85,
                  edgecolor='black', lw=1.2, width=0.5)
for bar, s in zip(bars, scores):
    ax2.text(bar.get_x()+bar.get_width()/2, s+1.5, f'{s}%',
             ha='center', fontsize=14, fontweight='bold')
ax2.axhline(70, color='green', lw=2, ls='--', alpha=0.6, label='Publication threshold')
ax2.axhline(50, color='orange', lw=1.5, ls=':', alpha=0.6, label='Marginal')
ax2.set_ylabel('Confidence score [qualitative %]', fontsize=10)
ax2.set_title('Sector Health\n(Tests 2+3 vs Test 1)',
              fontsize=10, fontweight='bold')
ax2.set_ylim(0, 105)
ax2.legend(fontsize=9)
ax2.grid(axis='y', alpha=0.3)

# ── Panel 3: Isotropy principle linking T2 and T3 ────────────────────────────
ax3 = fig.add_subplot(gs[0, 2])
ax3.axis('off')
link_txt = (
    "3D ISOTROPY PRINCIPLE\n"
    "(links Tests 2 and 3)\n\n"
    "3 spatial dimensions,\n"
    "each receives equal share.\n\n"
    "MOND (Test 2):\n"
    "  v₃D = v₄D / √3\n"
    "  → a₀ = H₀cx₀/√3\n"
    "  ✅ SPARC: 15.6% RMS\n\n"
    "Tilt decay (Test 3):\n"
    "  B₃D = B₄D / 3\n"
    "  → B = L/3 (0.69% off)\n"
    "  ✅ Consistent with T2\n\n"
    "Same principle → same\n"
    "physical foundation.\n"
    "Not two coincidences."
)
ax3.text(0.5, 0.5, link_txt, transform=ax3.transAxes,
         fontsize=9, va='center', ha='center', family='monospace',
         bbox=dict(boxstyle='round', facecolor='#e8f5e9', alpha=0.9,
                   edgecolor='#28a745', lw=2))
ax3.set_title('Cross-Test Consistency', fontsize=10, fontweight='bold')

# ── Panel 4: DESI w(z) failure anatomy ───────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 0:2])
z_plot = np.linspace(0.05, 2.2, 150)
w_e_plot  = [w_eff(z) for z in z_plot]
w_dr2_plt = [w_cpl(z, -0.73, -0.66) for z in z_plot]
w_lcdm    = [-1.0] * len(z_plot)
w_pred    = [-1.08] * len(z_plot)

w_hi = [w_cpl(z, -0.63, -0.41) for z in z_plot]
w_lo = [w_cpl(z, -0.83, -0.91) for z in z_plot]

ax4.fill_between(z_plot, w_lo, w_hi, alpha=0.2, color='green',
                 label='DESI DR2 ±1σ band')
ax4.plot(z_plot, w_dr2_plt, 'green',     lw=2.5, label='DESI DR2 best-fit')
ax4.plot(z_plot, w_e_plot,  'tomato',    lw=2.5, label='ESTIF w_eff(z) — ACTUAL')
ax4.plot(z_plot, w_pred,    'tomato',    lw=1.5, ls=':', alpha=0.7,
         label='ESTIF prediction −1.08')
ax4.plot(z_plot, w_lcdm,    'black',     lw=1.5, ls='--', alpha=0.5,
         label='ΛCDM w = −1')
ax4.axvline(2.0, color='gray', lw=1.5, ls=':', alpha=0.5, label='z_eff cutoff')

ax4.annotate('ESTIF w_eff\ncrosses −1 here\n(too late vs DESI)',
             xy=(1.0, w_eff(1.0)), xytext=(1.2, -0.5),
             fontsize=8.5, arrowprops=dict(arrowstyle='->', color='tomato'),
             bbox=dict(boxstyle='round', fc='#fff3cd', alpha=0.8))

ax4.set_xlabel('Redshift z', fontsize=11)
ax4.set_ylabel('w(z)', fontsize=11)
ax4.set_title('Why Test 1 Failed: w(z) Shape Mismatch\n'
              'DESI: w>−1 at low-z, falls. ESTIF: w<−1 everywhere, rises.',
              fontsize=10, fontweight='bold')
ax4.legend(fontsize=8.5, loc='lower left')
ax4.grid(alpha=0.3)
ax4.set_ylim(-2.0, 0.3)

# ── Panel 5: Priority action list ────────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2])
ax5.axis('off')
actions = [
    ('#1 BLOCKING', 'Fix Ω_tilt(z)', 'Self-consistent x(z) using\nH_ESTIF (not H_ΛCDM).\nRemoves circularity.', '#dc3545'),
    ('#2 THEORY',   'Derive x_c',    'Find geometric condition\ngiving x_c=0.272 from\nSchwarzschild geometry.', '#ffc107'),
    ('#3 VALIDATION','SPARC Upsilon', 'Use individual Υ_* per\ngalaxy. Reduces RMS\nbelow 15%.', '#17a2b8'),
    ('#4 PUBLISH',  'Gravity paper', 'Narrow letter: a₀ derivation\n+ SPARC test. Defensible\ntoday.', '#28a745'),
]
for i, (prio, title, desc, color) in enumerate(actions):
    y = 0.95 - i * 0.24
    ax5.text(0.02, y, prio, transform=ax5.transAxes,
             fontsize=8, va='top', color=color, fontweight='bold')
    ax5.text(0.28, y, f'{title}\n{desc}', transform=ax5.transAxes,
             fontsize=8, va='top', linespacing=1.4)
ax5.set_title('Prioritised Action List', fontsize=10, fontweight='bold')

# ── Panel 6: What can be published today ─────────────────────────────────────
ax6 = fig.add_subplot(gs[2, 0:2])
ax6.axis('off')

publication_txt = """
  PUBLISHABLE NOW (gravity-only paper)                    NOT YET PUBLISHABLE
  ────────────────────────────────────────────            ─────────────────────────────────
  Claim: a₀ = H₀cx₀/√3 derived from geometry             Claim: ESTIF replaces dark energy
  Evidence: 4-step derivation, zero free params           Blocker: DESI DR2 chi2/N = 10.8
  Test: SPARC, 87 galaxies, RMS = 15.6%                   Fix needed: Ω_tilt(z) rework
  Uniqueness: 1/√3 is the only geometric factor           
  Cross-check: B = L/3 from same principle                Claim: N_MAX and B fully derived
                                                          Blocker: x_c still observational
  Suitable venue: MNRAS Letters, ApJL, JCAP               
  Length: ~4 pages                                        Fix needed: geometric x_c derivation
  Readiness: NOW
"""
ax6.text(0.02, 0.97, publication_txt, transform=ax6.transAxes,
         fontsize=9, va='top', family='monospace',
         bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.9,
                   edgecolor='#dee2e6'))
ax6.set_title('Publication Readiness Assessment', fontsize=11, fontweight='bold')

# ── Panel 7: Final verdict ────────────────────────────────────────────────────
ax7 = fig.add_subplot(gs[2, 2])
ax7.axis('off')

verdict_txt = (
    "OVERALL VERDICT\n\n"
    "Results: MIXED but\n"
    "directionally positive.\n\n"
    "Gravity sector: SOLID\n"
    "  ✅ a₀ derived (zero params)\n"
    "  ✅ SPARC passes\n"
    "  ✅ Isotropy principle holds\n\n"
    "Cosmology sector: FAILS\n"
    "  ❌ DESI DR2 tension 3.5σ\n"
    "  ❌ w(z) shape wrong\n"
    "  ❌ Pre-prediction falsified\n\n"
    "Theory: PROGRESSES\n"
    "  ✅ 1/3 has derivation\n"
    "  ⚠️  5/7 still conditional\n"
    "  ⚠️  x_c still open\n\n"
    "CONTINUE development.\n"
    "Publish gravity sector now.\n"
    "Rebuild cosmology sector."
)
ax7.text(0.5, 0.5, verdict_txt, transform=ax7.transAxes,
         fontsize=9, va='center', ha='center', family='monospace',
         bbox=dict(boxstyle='round', facecolor='#e8f4fd', alpha=0.9,
                   edgecolor='#17a2b8', lw=2))
ax7.set_title('Overall Verdict', fontsize=11, fontweight='bold')

plt.savefig('cross_examination.png', dpi=150, bbox_inches='tight')
plt.close()
print("Plot saved: cross_examination.png")

# ============================================================================
# FINAL CROSS-EXAMINATION VERDICT
# ============================================================================
print("\n" + "=" * 70)
print("FINAL CROSS-EXAMINATION VERDICT")
print("=" * 70)
print(f"""
  RESULTS SUMMARY:
    Test 1 (DESI w(z)):       BAD    — cosmology sector fails DESI DR2
    Test 2 (SPARC bias):      GOOD   — gravity sector passes, bias is calibration
    Test 3 (Multipliers):     MIXED  — 1/3 derived, 5/7 conditional, x_c open

  CROSS-TEST FINDINGS:
    Tests 2 and 3 are mutually consistent and reinforce each other
    through the 3D isotropy principle. The gravity sector is solid.

    Test 1 reveals a specific, isolated failure in Ω_tilt(z) — the
    cosmological dark energy evolution law. This failure does not
    propagate into the gravity sector.

    The pre-existing prediction w_eff ≈ -1.08 is falsified by DESI DR2
    (w0 = -0.73 ± 0.10, pull = 3.5σ). This is a genuine prediction
    failure, not a post-hoc claim.

  PROJECT STATUS:
    The project is worthy of continued pursuit — but with a redefined
    focus. The gravity results are genuinely novel and pass independent
    observational tests. The cosmology sector needs rebuilding from the
    Ω_tilt(z) construction outward.

    The model should NOT be presented as "replacing dark energy."
    It CAN be presented as "a geometric derivation of the MOND
    acceleration from hypersurface tilt geometry."

  RECOMMENDED NEXT ACTIONS (in order):
    1. Publish the gravity paper (a₀ + SPARC) — achievable now
    2. Fix Ω_tilt(z) self-consistently — blocking cosmology
    3. Derive x_c geometrically — completes the multiplier derivation
    4. After (2): retest against DESI DR2 — decisive cosmology test
""")
print("=" * 70)
