"""
ESTIF fidelity audit — Step 1 evidence extractor
=================================================
The vacuum-flow derivation (estif_converse_flow_law.py) rests on three
axioms. The question of Step 1 is NOT whether those axioms give the right
physics -- they provably do -- but whether ESTIF v6.2, AS WRITTEN, actually
asserts them. That is a reading task. This script makes the reading
exhaustive and reproducible: it pulls every sentence in the corpus that
bears on each axiom, so the determination is made against the actual text
(including files a human reviewer may not have re-read) rather than memory.

The three axioms under audit
  A1  FLOW, NOT STRETCH: 3D space is a flat hypersurface carried through
      the bulk (Painleve-Gullstrand: flat 3-slices + shift). The competing
      picture -- a uniformly SHRINKING/SCALING space (FRW scale factor) --
      is a DIFFERENT mechanism and counts as a contradiction, not support.
  A2  UNIVERSAL SPEED c: every object moves through the 4D bulk at speed c,
      per-object, everywhere -- not merely "recession = c at the Hubble
      radius". A characteristic flow speed that is a FRACTION of c (e.g.
      cx0 ~ 0.31c) presented as THE flow speed counts against A2 unless
      explicitly framed as a 3-space projection of a total-c motion.
  A3  VACUUM CARRIES NO EFFECTIVE ENERGY: rho_eff = 0 where there is no
      matter (the condition that FORCES v^2 = 2A/r). Statements that T_mu_nu
      is unspecified / future work are the ABSENCE of A3, not its assertion.

What the script does
  - Walks the given paths (.md, .py, .tex, .txt), splits into sentences.
  - For each axiom, matches SUPPORT patterns and CONTRADICTION patterns.
  - Prints every hit with file:line and the sentence, grouped by axiom,
    support vs contradiction separated.
  - Prints a mechanical tally. The tally does NOT decide the audit -- it
    surfaces the evidence; the human makes the call on wording intent.

What the script explicitly does NOT do
  - It does not judge whether an ABSENT axiom is a fatal gap or a natural
    extension of the core principle. That is a physics/intent judgment,
    made in the companion written audit, not by keyword counting.

Run:  python3 estif_fidelity_audit.py [path ...]
      (defaults to the estif_publication tree if no path given)
Deps: standard library only.
"""

import os
import re
import sys

# ----------------------------------------------------------------------
# Axiom probes. Patterns are regexes matched case-insensitively against
# individual sentences. Keep them broad on SUPPORT (we would rather over-
# collect and let the reader discard) and specific on CONTRADICTION.
# ----------------------------------------------------------------------
AXIOMS = {
    "A1_FLOW_NOT_STRETCH": {
        "title": "A1  Flow, not stretch (flat 3-slices carried through bulk)",
        "support": [
            r"\bflows?\b(?!\s*chart)", r"\binward flow\b", r"\bmov(e|es|ing)\s+through\b",
            r"\bhypersurface\b", r"\bflat\s+(3-?slices?|space|hypersurface)\b",
            r"\bthrough (?:the )?4(?:th|-?d| ?dimension)\b",
            r"\borthogonal to (?:its|the) surface\b", r"\bpainlev", r"\bshift vector\b",
        ],
        "contradiction": [
            r"\bshrink(s|ing|age)?\b", r"\bscal(e|ing) factor\b",
            r"\beverything .* (?:shrink|scal)", r"\bruler shrinks\b",
            r"\bstretch(es|ing)?\b", r"\bexpand(s|ing) outward\b",
            r"\bS\(t\)\s*=\s*exp", r"\bshrinking together\b",
        ],
    },
    "A2_UNIVERSAL_C": {
        "title": "A2  Universal speed c (every object moves through bulk at c)",
        "support": [
            r"\beverything moves\b.*\bc\b", r"\bmov(e|es|ing) (?:through .*)?at (?:the )?speed of light\b",
            r"\bat speed c\b", r"\binfall velocity (?:is )?c\b",
            r"\bd\s*w\s*\^?2?\s*\+\s*d\s*sigma", r"\buniversal(?:ly)?\b.*\bc\b",
            r"\bmoves? through the (?:4d )?bulk at\b", r"\bspeed[- ]c\b",
        ],
        "contradiction": [
            r"\bv_?flow\s*=\s*c\s*[x*]\s*_?0", r"\bc\s*[x*]\s*_?0\b",
            r"\brecession velocity equals c\b", r"\b0?\.31\s*c\b",
            r"\bsubluminal\b", r"\bfraction .* of (?:the )?(?:hubble )?flow\b",
        ],
    },
    "A3_VACUUM_ZERO_ENERGY": {
        "title": "A3  Vacuum carries no effective energy (rho_eff = 0 in vacuum)",
        "support": [
            r"\bempty space\b", r"\bvacuum\b", r"\bwhere there is no (?:matter|mass)\b",
            r"\bno (?:matter|mass) .* no (?:gravity|energy)\b",
            r"\brho_?eff\s*=\s*0\b", r"\beffective energy .* (?:zero|vanish)",
            r"\bin the absence of (?:matter|mass)\b",
        ],
        "contradiction": [
            r"\bT_?mu_?nu\b.*\b(?:not|future|incomplete|unspecified|postulat)",
            r"\bstress[- ]energy .* (?:not|future|incomplete|projection has not)",
            r"\bnot yet\b.*\bT_?mu",
        ],
    },
}

SENT_SPLIT = re.compile(r'(?<=[.!?])\s+|\n{2,}')
EXT = (".md", ".py", ".tex", ".txt")


def gather_files(paths):
    out = []
    for p in paths:
        if os.path.isfile(p) and p.endswith(EXT):
            out.append(p)
        elif os.path.isdir(p):
            for root, _, files in os.walk(p):
                if any(s in root for s in ("/.git", "__pycache__", "/archive")):
                    continue
                for f in files:
                    if f.endswith(EXT):
                        out.append(os.path.join(root, f))
    return sorted(set(out))


def line_of(offset, text):
    return text.count("\n", 0, offset) + 1


def scan(files):
    results = {k: {"support": [], "contradiction": []} for k in AXIOMS}
    compiled = {
        k: {
            "support": [re.compile(p, re.I) for p in v["support"]],
            "contradiction": [re.compile(p, re.I) for p in v["contradiction"]],
        }
        for k, v in AXIOMS.items()
    }
    for path in files:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except Exception as e:
            print(f"  (skip {path}: {e})")
            continue
        pos = 0
        for sent in SENT_SPLIT.split(text):
            s = sent.strip()
            if not s:
                pos += len(sent) + 1
                continue
            off = text.find(s, pos)
            pos = off + len(s) if off >= 0 else pos + len(sent)
            ln = line_of(off if off >= 0 else pos, text)
            flat = re.sub(r"\s+", " ", s)[:220]
            for k in AXIOMS:
                for rx in compiled[k]["contradiction"]:
                    if rx.search(s):
                        results[k]["contradiction"].append((path, ln, flat, rx.pattern))
                        break
                for rx in compiled[k]["support"]:
                    if rx.search(s):
                        results[k]["support"].append((path, ln, flat, rx.pattern))
                        break
    return results


def rel(path):
    marker = "estif_publication/"
    i = path.find(marker)
    return path[i + len(marker):] if i >= 0 else path


def report(results):
    print("=" * 74)
    print("ESTIF FIDELITY AUDIT — textual evidence for the three axioms")
    print("=" * 74)
    for k, v in AXIOMS.items():
        sup = results[k]["support"]
        con = results[k]["contradiction"]
        print()
        print("-" * 74)
        print(v["title"])
        print("-" * 74)
        print(f"  SUPPORT sentences: {len(sup)}    CONTRADICTION sentences: {len(con)}")
        if con:
            print("\n  >>> CONTRADICTION evidence (competing mechanism / opposes axiom):")
            for path, ln, s, pat in con[:40]:
                print(f"    [{rel(path)}:{ln}]  {s}")
        if sup:
            print("\n  --- SUPPORT evidence (sentences touching the axiom):")
            for path, ln, s, pat in sup[:40]:
                print(f"    [{rel(path)}:{ln}]  {s}")
        if not sup and not con:
            print("\n  (no sentences in the corpus match this axiom either way)")
        # Mechanical flag only -- NOT the verdict.
        if not sup and not con:
            flag = "ABSENT (no textual basis found)"
        elif con and not sup:
            flag = "CONTRADICTED (only opposing evidence found)"
        elif con and sup:
            flag = "CONTESTED (both supporting and opposing evidence present)"
        else:
            flag = "TOUCHED (supporting sentences exist; judge if they ASSERT it)"
        print(f"\n  MECHANICAL FLAG: {flag}")
    print()
    print("=" * 74)
    print("READING THE OUTPUT")
    print("=" * 74)
    print("""  The flags are mechanical, not the audit's conclusion. In particular:
    - "TOUCHED" means sentences mention the concept; whether they ASSERT
      the axiom as a load-bearing premise is a judgment for the reader.
    - "CONTESTED" (esp. A1) means a competing mechanism appears in the
      same corpus -- the two must be reconciled before either can be a
      clean premise.
    - Absence of A3 support with T_mu_nu-is-future-work contradictions is
      the expected signature of a framework that POSTULATES its gravity
      rather than deriving it from a vacuum condition.
  The written audit interprets these; keyword counts do not close it.""")


if __name__ == "__main__":
    default = "/Users/peterangelov/estif_publication"
    paths = sys.argv[1:] if len(sys.argv) > 1 else [default]
    files = gather_files(paths)
    print(f"Scanning {len(files)} files under: {', '.join(paths)}\n")
    report(scan(files))
