#!/usr/bin/env python3
"""Reorders limitations (vi)(viii)(vii) -> (vi)(vii)(viii). Run from repo root."""
from pathlib import Path

F = Path("docs/Letter/6.4/ESTIF_letter_v6.4.md")
text = F.read_text(encoding="utf-8")

viii_block = (
    "(viii) **Bullet Cluster / cluster-scale lensing-baryon offset.** "
    "*(v6.4 addition, flagged open -- not resolved)* The Bullet Cluster is the sharpest "
    "instance of the cluster-scale problem in (vi). No ESTIF-specific analysis of this "
    "system has been carried out; this is flagged here as an explicit open item rather "
    "than left unaddressed by omission."
)
vii_block_start = "(vii) **Weak-lensing BTFR consistency"

assert viii_block in text, "viii block not found -- check current file contents first."
assert vii_block_start in text, "vii block not found -- check current file contents first."

# Remove (viii) from its current spot (it's sandwiched between two blank lines after (vi))
text = text.replace("\n\n" + viii_block, "", 1)

# Re-insert it directly after the full (vii) paragraph, before the "---" section break
vii_para_end_marker = "and remains unresolved."
idx = text.index(vii_para_end_marker) + len(vii_para_end_marker)
text = text[:idx] + "\n\n" + viii_block + text[idx:]

F.write_text(text, encoding="utf-8")
print("Reordered: limitations now read (v)(vi)(vii)(viii).")
