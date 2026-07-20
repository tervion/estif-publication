#!/usr/bin/env python3
"""
ESTIF path fixer -- v6.4.2 -- session 20260721 v2 (verified against a clone of
tervion/estif-publication @ 870aaab). Run AFTER reorg.sh. Targets tests/scripts/**.py
in the FINAL layout. Idempotent: second run makes zero edits.

Rules (all depth-aware: ../docs and ../src from scripts/, ../../ from phase2_a1prime/):
  1. CACHE_DIR      -> every CACHE_DIR under scripts/ points at tests/docs
                       (uniform policy: docs holds ALL cached data, scripts stays
                       pure .py, no duplicate downloads beside scripts).
                       Stale "write beside this script" comments are corrected.
  2. sys.path hops  -> join forms (including calls WRAPPED across two lines and
                       single-arg '../src' / '../../estif_publication/src') are
                       collapsed to one canonical __file__-based line at the
                       correct depth. Plain-string forms are upgraded. A bare
                       variable arg (src_path) is left; its DEFINITION is fixed.
  3. IDENT = join(dirname(abspath(__file__)), ...) assignments:
                       CACHE_DIR -> docs; tail 'src' (src_path) -> src hop;
                       SPARC_LOCAL -> docs cache. Other idents untouched.
  4. join(HERE, "<moved cache>") -> docs splice (HERE definition verified first).
  5. '/home/claude/...' container leftovers -> repo docs paths. The four DESI
     entries in test_desi_wz_consistency map onto the caches that already live
     in docs (dr1_mean.txt family) -- shared with estif_task5_desi, no duplicates.
     Both SPARC scripts converge on docs/SPARC_Lelli2016c.xml (the VOTable that
     test_sparc_tully_fisher downloads); the dead sparc_vizier.tsv dependency is
     retired as a first choice but kept as a docs-relative fallback in tully.
  6. TEST_INDEX.md quoted path literals -> ../docs/TEST_INDEX.md (depth-correct).
     Prose/docstring mentions are never touched.

TWO-PHASE SAFETY: phase A scans everything and writes nothing. Any unrecognized
form aborts with a full drift report and exit 1 -- zero files touched. Phase B
applies edits bottom-up with anchor asserts on every replaced line.
  --dry-run : print the full plan, write nothing.
"""

import os
import re
import sys
from pathlib import Path

MOVED_PREFIXES = ("dr1_mean", "dr1_cov", "dr2_mean", "dr2_cov", "dr2b_",
                  "pathone_dr1", "pathone_dr2", "pathtwo_target_lock_output", "dr2t_")

HOME_CLAUDE_MAP = {
    "desi_dr1_mean.txt": "dr1_mean.txt",
    "desi_dr1_cov.txt": "dr1_cov.txt",
    "desi_dr2_mean.txt": "dr2_mean.txt",
    "desi_dr2_cov.txt": "dr2_cov.txt",
    "sparc_vizier.tsv": "SPARC_Lelli2016c.xml",
    "SPARC_Lelli2016c.mrt": "SPARC_Lelli2016c.mrt",
}
TULLY_FALLBACK_KEEP = {"sparc_vizier.tsv": "sparc_vizier.tsv",
                       "SPARC_Lelli2016c.mrt": "SPARC_Lelli2016c.mrt"}

DIRNAME = r"os\.path\.dirname\(\s*(?:os\.path\.abspath\(\s*__file__\s*\)|__file__)\s*\)"

RE_ASSIGN_DIRNAME = re.compile(
    r"^(?P<ind>\s*)(?P<ident>[A-Za-z_]\w*)\s*=\s*" + DIRNAME + r"\s*(?P<cmt>#.*)?$")
RE_ASSIGN_JOIN = re.compile(
    r"^(?P<ind>\s*)(?P<ident>[A-Za-z_]\w*)\s*=\s*os\.path\.join\(\s*" + DIRNAME +
    r"\s*,(?P<args>.+)\)\s*(?P<cmt>#.*)?$")
RE_SP_ANY = re.compile(r"sys\.path\.(append|insert)")
RE_SP_JOIN = re.compile(
    r"^(?P<ind>\s*)sys\.path\.(?P<meth>append|insert)\(\s*(?P<zero>0\s*,\s*)?"
    r"os\.path\.join\(\s*" + DIRNAME + r"\s*,(?P<args>.+)\)\s*\)\s*(?P<cmt>#.*)?$")
RE_SP_STR = re.compile(
    r"^(?P<ind>\s*)sys\.path\.(?P<meth>append|insert)\(\s*(?P<zero>0\s*,\s*)?"
    r"(?P<q>['\"])(?P<path>[^'\"]+)(?P=q)\s*\)\s*(?P<cmt>#.*)?$")
RE_SP_IDENT = re.compile(
    r"^\s*sys\.path\.(append|insert)\(\s*(?:0\s*,\s*)?[A-Za-z_]\w*\s*\)\s*(#.*)?$")
RE_QUOTED = re.compile(r"['\"]([^'\"]+)['\"]")
RE_IDENT_JOIN = re.compile(
    r"os\.path\.join\(\s*(?P<ident>[A-Za-z_]\w*)\s*,\s*(?P<inner>(?:['\"][^'\"]+['\"]\s*,\s*)*)"
    r"(?P<q>['\"])(?P<name>[^'\"]+)(?P=q)\s*\)")
RE_HOME = re.compile(r"(?P<q>['\"])/home/claude/(?P<name>[^'\"]+)(?P=q)")
RE_TI_ANY = re.compile(r"(['\"])((?:\.\./)*(?:docs_tests/|docs/)?TEST_INDEX\.md)\1")
RE_SRC_IMPORT = re.compile(r"^\s*(from\s+src[\s.]|import\s+estif_ec|from\s+estif_ec)")
STALE_CMT = "beside this script"


def norm_args(argtext):
    lits = RE_QUOTED.findall(argtext)
    only_lits = not re.sub(r"['\"][^'\"]*['\"]|[\s,]", "", argtext)
    if not only_lits:
        return None
    out = []
    for a in lits:
        out.extend(seg for seg in a.split("/") if seg not in ("", "."))
    return out


def qargs(parts):
    return ", ".join(f'"{p}"' for p in parts)


def joinexpr(parts, tail=None):
    ps = list(parts) + ([tail] if tail else [])
    return f'os.path.join(os.path.dirname(os.path.abspath(__file__)), {qargs(ps)})'


def cmt(c):
    return f"  {c}" if c else ""


def sp_line(ind, meth, zero, parts, c):
    z = "0, " if zero else ""
    return f"{ind}sys.path.{meth}({z}{joinexpr(parts)}){cmt(c)}"


def logical_units(lines):
    i = 0
    while i < len(lines):
        j, buf = i, lines[i]
        while buf.count("(") > buf.count(")") and j + 1 < len(lines) and j - i < 4:
            j += 1
            buf = buf.rstrip() + " " + lines[j].strip()
        yield i, j, buf
        i = j + 1


def plan_file(py: Path, root: Path):
    rel = py.relative_to(root)
    d = py.parent
    parts_docs = os.path.relpath(str(root / "tests" / "docs"), str(d)).split(os.sep)
    parts_src = os.path.relpath(str(root / "src"), str(d)).split(os.sep)
    parts_root = os.path.relpath(str(root), str(d)).split(os.sep)
    ti_final = "/".join(parts_docs) + "/TEST_INDEX.md"

    raw = py.read_text(encoding="utf-8")
    lines = raw.splitlines()
    text = raw
    filedir_idents = set(m.group("ident") for m in
                         (RE_ASSIGN_DIRNAME.match(l) for l in lines) if m)
    edits, notes, drifts = [], [], []   # edits: (start, end, [old lines], new line)

    for i, j, unit in logical_units(lines):
        n = i + 1
        span = lines[i:j + 1]

        m = RE_ASSIGN_DIRNAME.match(unit)
        if m:
            if m.group("ident") == "CACHE_DIR":
                c = m.group("cmt")
                if c and STALE_CMT in c:
                    c = "# caches live in tests/docs"
                new = f'{m.group("ind")}CACHE_DIR = {joinexpr(parts_docs)}{cmt(c)}'
                edits.append((i, j, span, new, f"CACHE_DIR -> {parts_docs} (uniform docs policy)"))
            continue

        m = RE_ASSIGN_JOIN.match(unit)
        if m:
            ident, na = m.group("ident"), norm_args(m.group("args"))
            if na is None:
                if ident in ("CACHE_DIR", "SPARC_LOCAL", "src_path"):
                    drifts.append((n, unit, f"{ident} join uses non-literal args"))
                continue
            if ident == "CACHE_DIR":
                if na == parts_docs:
                    notes.append((n, "ALREADY-FINAL: CACHE_DIR"))
                else:
                    c = m.group("cmt")
                    if c and STALE_CMT in c:
                        c = "# caches live in tests/docs"
                    edits.append((i, j, span,
                                  f'{m.group("ind")}CACHE_DIR = {joinexpr(parts_docs)}{cmt(c)}',
                                  f"CACHE_DIR {na} -> {parts_docs}"))
            elif ident == "SPARC_LOCAL":
                want = parts_docs + ["SPARC_Lelli2016c.xml"]
                if na == want:
                    notes.append((n, "ALREADY-FINAL: SPARC_LOCAL"))
                else:
                    edits.append((i, j, span,
                                  f'{m.group("ind")}SPARC_LOCAL = {joinexpr(parts_docs, "SPARC_Lelli2016c.xml")}{cmt(m.group("cmt"))}',
                                  f"SPARC_LOCAL -> docs cache"))
            elif na and na[-1] == "src":
                if na == parts_src:
                    notes.append((n, f"ALREADY-FINAL: {ident} -> src"))
                else:
                    edits.append((i, j, span,
                                  f'{m.group("ind")}{ident} = {joinexpr(parts_src)}{cmt(m.group("cmt"))}',
                                  f"{ident} src hop {na} -> {parts_src}"))
            continue

        if RE_SP_ANY.search(unit):
            m = RE_SP_JOIN.match(unit)
            s = RE_SP_STR.match(unit)
            if m:
                na = norm_args(m.group("args"))
                if na is None:
                    drifts.append((n, unit, "sys.path join uses non-literal args"))
                elif na and na[-1] == "src":
                    if na == parts_src:
                        notes.append((n, "ALREADY-FINAL: sys.path -> src"))
                    else:
                        edits.append((i, j, span,
                                      sp_line(m.group("ind"), m.group("meth"), m.group("zero"),
                                              parts_src, m.group("cmt")),
                                      f"sys.path src hop {na} -> {parts_src}"))
                elif na and all(a == ".." for a in na):
                    if na == parts_root:
                        notes.append((n, "ALREADY-FINAL: sys.path -> repo root"))
                    else:
                        edits.append((i, j, span,
                                      sp_line(m.group("ind"), m.group("meth"), m.group("zero"),
                                              parts_root, m.group("cmt")),
                                      f"sys.path root hop {na} -> {parts_root}"))
                else:
                    drifts.append((n, unit, f"sys.path join args {na} unrecognized"))
            elif s:
                p = s.group("path").rstrip("/")
                segs = [x for x in p.split("/") if x not in ("", ".")]
                if segs and segs[-1] == "src":
                    edits.append((i, j, span,
                                  sp_line(s.group("ind"), s.group("meth"), s.group("zero"),
                                          parts_src, s.group("cmt")),
                                  f"sys.path string '{p}' UPGRADED to __file__-based src hop"))
                elif segs and all(x == ".." for x in segs):
                    edits.append((i, j, span,
                                  sp_line(s.group("ind"), s.group("meth"), s.group("zero"),
                                          parts_root, s.group("cmt")),
                                  f"sys.path string '{p}' UPGRADED to __file__-based root hop"))
                else:
                    drifts.append((n, unit, f"sys.path string path '{p}' unrecognized"))
            elif RE_SP_IDENT.match(unit):
                notes.append((n, "sys.path uses a variable -- its definition is handled above"))
            else:
                if "sys.path" in unit.split("#")[0]:
                    drifts.append((n, unit, "sys.path manipulation in unrecognized form"))
            continue

        if RE_HOME.search(unit):
            if i != j:
                drifts.append((n, unit, "/home/claude path on a wrapped line -- extend fixer"))
                continue
            line = lines[i]
            ok = True

            def home_sub(mm):
                nonlocal ok
                name = mm.group("name")
                table = TULLY_FALLBACK_KEEP if py.name == "test_sparc_tully_fisher.py" else HOME_CLAUDE_MAP
                if name not in table:
                    ok = False
                    return mm.group(0)
                return joinexpr(parts_docs, table[name])

            new = RE_HOME.sub(home_sub, line)
            if not ok:
                drifts.append((n, line, "unmapped /home/claude path -- extend HOME_CLAUDE_MAP"))
            elif new != line:
                edits.append((i, i, [line], new, "/home/claude container path -> tests/docs"))
            continue

        hit = False
        for mm in RE_IDENT_JOIN.finditer(unit):
            ident, name, inner = mm.group("ident"), mm.group("name"), mm.group("inner")
            if ident == "CACHE_DIR" or not name.startswith(MOVED_PREFIXES):
                continue
            inner_parts = RE_QUOTED.findall(inner)
            if inner_parts == parts_docs:
                notes.append((n, f"ALREADY-FINAL: join({ident}, docs, {name})"))
                continue
            if inner_parts:
                drifts.append((n, unit, f"join({ident}, ...) has unexpected parts {inner_parts}"))
                hit = True
                continue
            if ident not in filedir_idents:
                drifts.append((n, unit,
                               f"join({ident}, '{name}') -- {ident} is not a verified "
                               f"__file__-dirname variable in this file"))
                hit = True
                continue
            if i != j:
                drifts.append((n, unit, "moved-cache join on a wrapped line -- extend fixer"))
                hit = True
                continue
            line = lines[i]
            new = line.replace(f"{ident}, {mm.group('q')}{name}{mm.group('q')}",
                               f"{ident}, {qargs(parts_docs)}, {mm.group('q')}{name}{mm.group('q')}")
            if new != line:
                edits.append((i, i, [line], new, f"join({ident}, '{name}') spliced -> docs"))
                hit = True
        if hit:
            continue

        if "TEST_INDEX" in unit and i == j:
            line = lines[i]
            def ti_sub(mm):
                return mm.group(1) + ti_final + mm.group(1)
            new = RE_TI_ANY.sub(ti_sub, line)
            if new != line and RE_TI_ANY.search(line).group(2) != ti_final:
                edits.append((i, i, [line], new, f"TEST_INDEX literal -> '{ti_final}'"))
            elif RE_TI_ANY.search(line) and RE_TI_ANY.search(line).group(2) == ti_final:
                notes.append((n, "ALREADY-FINAL: TEST_INDEX literal"))

    if RE_SRC_IMPORT.search(text) and "sys.path" not in text:
        notes.append((0, "WARNING: imports src/estif_ec but has no sys.path hop -- relies on cwd"))

    return edits, notes, drifts


def main():
    argv = [a for a in sys.argv[1:]]
    dry = "--dry-run" in argv
    argv = [a for a in argv if a != "--dry-run"]
    here = Path(__file__).resolve()
    root = Path(argv[0]).resolve() if argv else here.parent.parent
    scripts = root / "tests" / "scripts"
    if not scripts.is_dir():
        print(f"ABORT: {scripts} not found -- run reorg.sh first (move, THEN fix).")
        sys.exit(2)

    targets = sorted(p for p in scripts.rglob("*.py") if "__pycache__" not in p.parts)
    print("=" * 78)
    print(f"ESTIF PATH FIXER -- v6.4.2 -- session 20260721 v2 -- {len(targets)} script(s)")
    print("=" * 78)

    plans, all_drifts = {}, []
    for py in targets:
        edits, notes, drifts = plan_file(py, root)
        plans[py] = (edits, notes)
        for n, line, why in drifts:
            all_drifts.append((py, n, line, why))

    if all_drifts:
        print(f"\nDRIFT DETECTED in {len(set(p for p, *_ in all_drifts))} file(s) -- "
              "NOTHING WAS MODIFIED. Full report:\n")
        for py, n, line, why in all_drifts:
            print(f"  {py.relative_to(root)}:{n}")
            print(f"    line : {line.strip()[:110]}")
            print(f"    why  : {why}\n")
        print("Paste this report back to Claude; the fixer will be extended for these "
              "exact forms, then re-run. Zero files were touched.")
        sys.exit(1)

    changed_files, total_edits = 0, 0
    for py in targets:
        edits, notes = plans[py]
        if not edits and not notes:
            continue
        rel = py.relative_to(root)
        print(f"\n  {rel}")
        for n, msg in notes:
            print(f"    note  line {n}: {msg}" if n else f"    {msg}")
        if not edits:
            continue
        lines = py.read_text(encoding="utf-8").splitlines()
        for start, end, old, new, why in sorted(edits, key=lambda e: -e[0]):
            assert lines[start:end + 1] == old, f"anchor moved in {rel}:{start + 1} -- aborting"
            lines[start:end + 1] = [new]
            total_edits += 1
            tag = "plan " if dry else "FIXED"
            print(f"    {tag} line {start + 1}: {why}")
        if not dry:
            py.write_text("\n".join(lines) + "\n", encoding="utf-8")
        changed_files += 1

    print("\n" + "=" * 78)
    print(f"{'DRY RUN -- would touch' if dry else 'DONE --'} {changed_files} file(s), "
          f"{total_edits} edit(s). Second run must report 0 edits (idempotency check).")
    leftovers = []
    for py in targets:
        t = py.read_text(encoding="utf-8")
        if "/home/claude" in t:
            leftovers.append(f"{py.relative_to(root)} still contains /home/claude")
    for x in leftovers:
        print(f"  RESIDUAL: {x}")
    if not leftovers:
        print("  residual check: no /home/claude paths remain under scripts/.")


if __name__ == "__main__":
    main()
