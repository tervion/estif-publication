#!/usr/bin/env python3
"""
TASK ONE fixer -- v6.4.2 repair round, script-side corrections (D5 batch).
Run from the estif_publication repo root:  python3 task1_apply_corrections.py

Applies, in order:
  F1  C3  -- clean the duplicated draft "[3] MANDATORY FILTER" block out of
             audit/Jul 9 Session/estif_jwst_growth_spec.py and install the
             cleaned file at tests/estif_jwst_growth_spec.py (overwrites the
             pre-C3 repo copy; the audit original is left untouched).
  F2  C2  -- merge the two honest closing lines from the Jul 9 audit copy
             into tests/estif_pathone_cosmology.py (repo base kept).
  F3  C5  -- replace the E1 explanatory print block in
             tests/estif_task6_eddy_eos.py with the CORRECTIONS_v6.3.1 text.
  F4      -- relativize the four hardcoded /home/claude/ figure paths in
             tests/estif_flow_sim.py (Mac-mini portability, P2C-2).
  F5      -- supersession header above best_a0() in tests/btfr_lensing.py
             (Jul 11 finding A-1).
  F6      -- draft-duplicate docstrings prepended to tests/test_UKN.py and
             tests/test_UKN2.py (Jul 12 batch).
  F7      -- requirements.txt: add sympy + colossus, restamp header v6.4.2,
             scope-note the numpy pin. The APPROVED-FORK marker line is
             preserved byte-for-byte and verified after patching.

Every fix asserts its anchor exists before touching anything and reports
ALREADY APPLIED (skip) if the new text is found, so re-running is safe.
Nothing is deleted from audit/. Exit 0 iff every fix is applied or skipped.
"""
import sys
from pathlib import Path

if not Path("tests").is_dir() or not Path("audit").is_dir():
    raise SystemExit("Run from the estif_publication repo root (tests/ and audit/ not found).")

report = []

def apply(path, old, new, label, count=1):
    p = Path(path)
    t = p.read_text(encoding="utf-8")
    if new in t:
        report.append((label, "ALREADY APPLIED (skipped)"))
        return
    if old not in t:
        raise SystemExit(f"ANCHOR NOT FOUND for {label} in {path} -- file has drifted; aborting with no changes to it.")
    if count == 1 and t.count(old) != 1:
        raise SystemExit(f"ANCHOR NOT UNIQUE for {label} in {path} ({t.count(old)} hits); aborting.")
    p.write_text(t.replace(old, new, count if count else -1), encoding="utf-8")
    report.append((label, "applied"))

# ---------------------------------------------------------------- F1 : C3
src = Path("audit/Jul 9 Session/estif_jwst_growth_spec.py")
dst = Path("tests/estif_jwst_growth_spec.py")
t = src.read_text(encoding="utf-8")
m_start = 'print("[3] MANDATORY FILTER'
m_end = 'print("[3] HARD FILTER'
if m_start in t:
    i, j = t.index(m_start), t.index(m_end)
    if not (0 <= i < j):
        raise SystemExit("F1: unexpected block order in the audit jwst copy; aborting.")
    cleaned = t[:i] + t[j:]
else:
    cleaned = t
if m_start in cleaned:
    raise SystemExit("F1: draft block still present after cleaning; aborting.")
if m_end not in cleaned:
    raise SystemExit("F1: HARD FILTER block missing from cleaned file; aborting.")
existing = dst.read_text(encoding="utf-8") if dst.exists() else ""
if existing == cleaned:
    report.append(("F1 C3 jwst clean+install", "ALREADY APPLIED (skipped)"))
else:
    dst.write_text(cleaned, encoding="utf-8")
    report.append(("F1 C3 jwst clean+install", f"applied (installed {dst}, {len(cleaned)} chars)"))

# ---------------------------------------------------------------- F2 : C2
apply("tests/estif_pathone_cosmology.py",
'''    print("  geometric Omega_m + a physical dark-energy interpretation, zero fitted")
    print("  cosmological parameters.")''',
'''    print("  a geometric Omega_m CONSISTENCY relation + a physical dark-energy")
    print("  interpretation, with no fitted DARK-ENERGY parameters (w=-1 fixed).")''',
"F2 C2 closing-lines merge")

# ---------------------------------------------------------------- F3 : C5
apply("tests/estif_task6_eddy_eos.py",
'''print("""  Physical picture: the cosmic eddy is bulk rotation. For a rotating
  shell of comoving radius scaling as a, angular momentum L = I omega with
  moment of inertia I ~ M a^2. Conserving L per comoving patch as a grows
  gives omega ~ a^-2. Rotational energy density:
      rho_rot = (1/2) I omega^2 / Volume ~ (M a^2)(a^-2)^2 / a^3 = M a^-5.
  Wait -- per unit PROPER volume (~a^3) and with the a^2 in I absorbed into
  the comoving mass, the invariant scaling of ROTATIONAL energy density is
  rho_rot ~ a^-6  (the classic 'stiff' spin-energy scaling: L=const, E~L^2/I,
  I~a^2 per patch, energy density ~ (a^-2)^2 * a^? ...). The engine-agnostic
  robust statement: conserved-L rotation is a STIFF component.""")''',
'''print("""  Physical picture: the cosmic eddy is bulk rotation. A rotating patch of
  comoving size ~a has moment of inertia I ~ M a^2; conserving angular momentum
  L = I*omega per comoving patch gives omega ~ a^-2. The rotational energy density
  rho_rot ~ I omega^2 / a^3 then scales as a steep NEGATIVE power of a -- a
  "stiff"/blueshifting component (w > 0). The exact exponent depends on how the
  comoving mass and volume factors are booked, but for ANY such conserved-L
  reduction the exponent is >= 5, i.e. w >= 2/3: the component GROWS toward the
  past and is negligible today. It behaves as extra early matter/stiff fluid,
  NOT as dark energy. We take the representative stiff case w = +1 below.""")''',
"F3 C5 E1 print-block")

# ---------------------------------------------------------------- F4 : flow_sim paths
for n, name in enumerate(["fig1_single_drain", "fig2_mass_independence",
                          "fig3_H_tracking", "fig4_magnitude_gap"], 1):
    apply("tests/estif_flow_sim.py",
          f'fig.savefig("/home/claude/{name}.png", dpi=110)',
          f'fig.savefig("{name}.png", dpi=110)',
          f"F4.{n} flow_sim path {name}")

# ---------------------------------------------------------------- F5 : btfr header
apply("tests/btfr_lensing.py",
"def best_a0(pts):",
'''# --- SUPERSESSION NOTE (v6.4.2, 19 Jul 2026) -------------------------------
# The inverted-a0 block below (best_a0) propagates the 0.1 dex M*/L systematic
# PER BIN and is SUPERSEDED by a0_tension_corrected.py (fully correlated
# treatment: +1.17 sigma kinematic, +1.52 sigma lensing, letter section 6.3
# item vii). The forward chi2/N test (report) above remains the current
# data-level receipt.
# ---------------------------------------------------------------------------
def best_a0(pts):''',
"F5 btfr supersession header")

# ---------------------------------------------------------------- F6 : UKN headers
apply("tests/test_UKN.py",
"import sympy as sp\nx,y,z,G,A,B,C,H = sp.symbols",
'''"""Draft duplicate (AST-identical) of phase2_a1prime/estif_p2_door3_swirl_ledger.py.
Kept as the filename cited in RHAC-004; the canonical documented copy supersedes. (v6.4.2)"""
import sympy as sp
x,y,z,G,A,B,C,H = sp.symbols''',
"F6.1 test_UKN header")

apply("tests/test_UKN2.py",
"import sympy as sp\n\na,H0,Om,G = sp.symbols",
'''"""Draft duplicate (AST-identical) of phase2_a1prime/estif_a1prime_deepen_exact.py.
Kept as the filename cited in RHAC-005/006; the canonical documented copy supersedes. (v6.4.2)"""
import sympy as sp

a,H0,Om,G = sp.symbols''',
"F6.2 test_UKN2 header")

# ---------------------------------------------------------------- F7 : requirements.txt
apply("requirements.txt",
"# ESTIF v6.1 \u2014 Emergent Spacetime from Inward Flow",
"# ESTIF v6.4.2 \u2014 Emergent Spacetime from Inward Flow",
"F7.1 requirements header restamp")

apply("requirements.txt",
'''numpy>=1.24.0,<2.0.0
scipy>=1.10.0,<2.0.0
matplotlib>=3.7.0,<4.0.0
astropy>=5.2.0,<7.0.0''',
'''numpy>=1.24.0,<2.0.0
scipy>=1.10.0,<2.0.0
matplotlib>=3.7.0,<4.0.0
astropy>=5.2.0,<7.0.0
sympy>=1.12,<2.0.0
colossus>=1.3.0''',
"F7.2 requirements add sympy+colossus")

apply("requirements.txt",
'''#   incompatible. That warning is harmless \u2014 opencv is not used by ESTIF.''',
'''#   incompatible. That warning is harmless \u2014 opencv is not used by ESTIF.
#   Scope: the pin protects src/ model code. tests/ receipts have also run
#   green on numpy 2.x (audit machines, Jul 2026); the pin stays for
#   reproducibility of src/.
#
# sympy / colossus (added v6.4.2)
#   sympy: required by the phase2_a1prime suite, test_UKN*.py, and
#   estif_C15_gw_sector.py. colossus: required by estif_front1_*,
#   estif_front2_*, and estif_jwst_growth_spec.py (mass-function receipts).''',
"F7.3 requirements numpy scope note")

marker = "#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-16-10-25-V-2"
if marker not in Path("requirements.txt").read_text(encoding="utf-8"):
    raise SystemExit("F7 GUARD FAILED: APPROVED-FORK marker missing after patch; restore from git.")
report.append(("F7.4 APPROVED-FORK marker preserved", "verified"))

print("TASK ONE -- script-side corrections (v6.4.2)")
print("-" * 60)
for label, state in report:
    print(f"  {label:<38} {state}")
print("-" * 60)
print("Done. Nothing in audit/ was modified.")
