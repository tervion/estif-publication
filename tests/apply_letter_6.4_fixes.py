#!/usr/bin/env python3
"""
D11 — Apply the four letter fixes + chi-squared pick, write docs/Letter/6.4/.
Run from the repo root: python3 apply_letter_6.4_fixes.py
"""
import shutil
from pathlib import Path

SRC = Path("docs/Letter/6.3/ESTIF_letter_v6.3.md")
DST_DIR = Path("docs/Letter/6.4")
DST = DST_DIR / "ESTIF_letter_v6.4.md"

if not SRC.exists():
    raise SystemExit(f"Not found: {SRC} -- run this from the estif_publication repo root.")

text = SRC.read_text(encoding="utf-8")
original_len = len(text)

# --- Fix 0: header stamp ---
text = text.replace(
    "March 2026 (v6.2) · updated July 2026 (v6.3)",
    "March 2026 (v6.2) · updated July 2026 (v6.3) · updated July 2026 (v6.4)"
)

# --- Fix 1: C1 conditional on the LISA line ---
old_lisa = (
    "and the predicted LISA GW delay (491 \u00b5s, S/N = 49\u03c3) with no free parameters "
    "after calibration. The strong-field results are independent of the MOND construction; "
    "we mention them only to establish that the framework is internally consistent at "
    "multiple scales."
)
new_lisa = (
    "and the predicted LISA GW delay (491 \u00b5s, S/N = 49\u03c3) with no free parameters "
    "after calibration. **The LISA deviation-from-GR figure is conditional on the eddy-stress "
    "sector remaining as currently specified (C1); the ESTIF vacuum solution is exactly "
    "Schwarzschild, so any deviation from GR enters only through that sector.** "
    "The strong-field results are independent of the MOND construction; "
    "we mention them only to establish that the framework is internally consistent at "
    "multiple scales."
)
assert old_lisa in text, "Fix 1 anchor not found -- letter text may have changed."
text = text.replace(old_lisa, new_lisa)

# --- Fix 2: limitation (v) -- honorable-null rewrite, chi^2/N = 1.965 ---
old_v = (
    "(v) The ESTIF cosmological dark energy sector was found inconsistent with DESI DR2 BAO "
    "data and is under revision; no cosmological claims are made here. *(v6.3 note: Path One "
    "frozen-eddy limit now ties $\\Lambda$CDM on DESI DR2, $\\chi^2/N = 1.92$; Path Two "
    "thawing derivation is in progress.)*"
)
new_v = (
    "(v) The ESTIF cosmological dark energy sector has been fully re-examined since v6.3. "
    "Path Two (a dynamical, thawing dark-energy sector) was pre-registered and closed: the "
    "derivation returns $w = -1$ exactly, an honorable null -- no dynamical dark energy is "
    "found, and none is claimed. Path One (frozen-eddy, constant-$\\Lambda$ limit) ties "
    "$\\Lambda$CDM on DESI DR2, $\\chi^2/N = 1.965$ (13-bin pipeline). No cosmological claims "
    "beyond this are made here. *(v6.4 update, supersedes the v6.3 note.)*"
)
assert old_v in text, "Fix 2 anchor not found -- letter text may have changed."
text = text.replace(old_v, new_v)

# --- Fix 3: state which a0 the (vii) tensions are measured against ---
old_vii_anchor = "primary sample $R < 1000\\,\\mathrm{kpc}$), both fully accounted for"
new_vii_anchor = (
    "primary sample $R < 1000\\,\\mathrm{kpc}$), both measured against the bootstrap value "
    "$a_0 = 1.192\\times10^{-10}\\,\\mathrm{m/s^2}$ used internally by these two scripts -- "
    "not the paper's headline derived value $a_0 = 1.179\\times10^{-10}\\,\\mathrm{m/s^2}$ "
    "(Section 3); relative to the headline value the tensions are marginally larger "
    "(kinematic \u2248 +1.21\u03c3). Both figures are fully accounted for"
)
assert old_vii_anchor in text, "Fix 3 anchor not found -- letter text may have changed."
text = text.replace(old_vii_anchor, new_vii_anchor)

# --- Fix 4: Bullet Cluster -- HONEST placeholder, not a fabricated answer ---
old_vi_end = (
    "(vi) ESTIF reproduces the normalization of $a_0$ but inherits MOND's known cluster-scale "
    "missing-mass problem; this is not addressed here."
)
new_vi_plus_viii = (
    "(vi) ESTIF reproduces the normalization of $a_0$ but inherits MOND's known cluster-scale "
    "missing-mass problem; this is not addressed here.\n\n"
    "(viii) **Bullet Cluster / cluster-scale lensing-baryon offset.** *(v6.4 addition, flagged "
    "open -- not resolved)* The Bullet Cluster is the sharpest instance of the cluster-scale "
    "problem in (vi). No ESTIF-specific analysis of this system has been carried out; this is "
    "flagged here as an explicit open item rather than left unaddressed by omission."
)
assert old_vi_end in text, "Fix 4 anchor not found -- letter text may have changed."
text = text.replace(old_vi_end, new_vi_plus_viii)

# --- Changelog entry ---
old_changelog_head = "## Changelog\n\n**v6.3 (July 2026):**"
new_changelog_head = (
    "## Changelog\n\n"
    "**v6.4 (July 2026):**\n"
    "- \u00a72.1 box: LISA 491 \u00b5s line made explicitly conditional on C1 (open eddy-stress sector).\n"
    "- Limitation (v) rewritten: Path Two closed as an honorable null ($w=-1$ exactly); "
    "Path One figure corrected to $\\chi^2/N = 1.965$ (13-bin pipeline), replacing the "
    "task6 constant-$\\Lambda$-limit figure of 1.92.\n"
    "- Limitation (vii): added explicit statement that the +1.17\u03c3/+1.52\u03c3 tensions use the "
    "bootstrap $a_0$, not the paper's headline derived value.\n"
    "- New limitation (viii): Bullet Cluster flagged as an open, unresolved item.\n\n"
    "**v6.3 (July 2026):**"
)
assert old_changelog_head in text, "Changelog anchor not found -- letter text may have changed."
text = text.replace(old_changelog_head, new_changelog_head)

DST_DIR.mkdir(parents=True, exist_ok=True)
DST.write_text(text, encoding="utf-8")

print(f"Wrote {DST}  ({len(text)} chars, was {original_len} in v6.3)")
print("All 4 fixes + chi^2/N=1.965 applied. Bullet Cluster is a FLAG, not an answer -- see limitation (viii).")
print()
print("Optional PDF export (only if you want one -- mirrors the docs/Letter/6.2/ structure):")
print("  pandoc docs/Letter/6.4/ESTIF_letter_v6.4.md -o docs/Letter/6.4/PDF/ESTIF_letter_6.4.pdf")
