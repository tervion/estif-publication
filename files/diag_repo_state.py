#!/usr/bin/env python3
"""
ESTIF repo diagnostic -- v6.4.2 -- session 20260721 v2 (read-only).
  BEFORE reorg -> before-picture of the layout and every path construct.
  AFTER reorg + fix_paths.py -> end-state receipt (all six checks PASS = done).

v2 changes vs v1: .DS_Store is macOS junk, counted separately and never fails a
check (reorg.sh deletes it); src-import detection ignores comment lines; CACHE
policy is uniform (every CACHE_DIR must point at tests/docs -- scripts/ stays
pure .py); dynamic f-string cache families (pathone_/dr2b_/dr2t_/aicbic_/
closure_) are detected; /home/claude container leftovers and wrapped two-line
sys.path calls are recognized.

Read-only. Nothing is modified. Optional argv: repo root (default: parent of
this script's folder).
"""

import os
import re
import subprocess
import sys
from pathlib import Path

MOVED_PREFIXES = ("dr1_mean", "dr1_cov", "dr2_mean", "dr2_cov", "dr2b_",
                  "pathone_dr1", "pathone_dr2", "pathtwo_target_lock_output", "dr2t_")
DYNAMIC_FAMILIES = ('f"pathone_{', 'f"dr2b_{', 'f"dr2t_{', 'f"aicbic_{',
                    'f"closure_{', "f'pathone_{", "f'dr2b_{", "f'dr2t_{",
                    "f'aicbic_{", "f'closure_{", 'f"{key}.txt"', "f'{key}.txt'")

DIRNAME = r"os\.path\.dirname\(\s*(?:os\.path\.abspath\(\s*__file__\s*\)|__file__)\s*\)"
RE_CACHE_BARE = re.compile(r"^\s*CACHE_DIR\s*=\s*" + DIRNAME + r"\s*(#.*)?$")
RE_CACHE_JOIN = re.compile(r"^\s*CACHE_DIR\s*=\s*os\.path\.join\(\s*" + DIRNAME +
                           r"\s*,(?P<args>.+)\)\s*(#.*)?$")
RE_CACHE_ANY = re.compile(r"^\s*CACHE_DIR\s*=")
RE_SP_ANY = re.compile(r"sys\.path\.(append|insert)")
RE_SP_JOIN = re.compile(r"^\s*sys\.path\.(append|insert)\(\s*(?:0\s*,\s*)?"
                        r"os\.path\.join\(\s*" + DIRNAME + r"\s*,(?P<args>.+)\)\s*\)\s*(#.*)?$")
RE_SP_IDENT = re.compile(r"^\s*sys\.path\.(append|insert)\(\s*(?:0\s*,\s*)?[A-Za-z_]\w*\s*\)\s*(#.*)?$")
RE_QUOTED = re.compile(r"['\"]([^'\"]+)['\"]")
RE_SRC_IMPORT = re.compile(r"^\s*(from\s+src[\s.]|import\s+estif_ec|from\s+estif_ec)")
RE_TI = re.compile(r"TEST_INDEX")
RE_OUT = re.compile(r"savefig|\.png")
RE_IDENT_JOIN = re.compile(r"os\.path\.join\(\s*([A-Za-z_]\w*)\s*,\s*"
                           r"((?:['\"][^'\"]+['\"]\s*,\s*)*)(['\"])([^'\"]+)\3\s*\)")


def norm_args(argtext):
    lits = RE_QUOTED.findall(argtext)
    only_lits = not re.sub(r"['\"][^'\"]*['\"]|[\s,]", "", argtext)
    if not only_lits:
        return None
    out = []
    for a in lits:
        out.extend(seg for seg in a.split("/") if seg not in ("", "."))
    return out


def logical_units(lines):
    i = 0
    while i < len(lines):
        j, buf = i, lines[i]
        while buf.count("(") > buf.count(")") and j + 1 < len(lines) and j - i < 4:
            j += 1
            buf = buf.rstrip() + " " + lines[j].strip()
        yield i, j, buf
        i = j + 1


def git(root, *args):
    try:
        r = subprocess.run(["git", "-C", str(root)] + list(args),
                           capture_output=True, text=True, timeout=30)
        return r.stdout.strip()
    except Exception:
        return ""


def main():
    here = Path(__file__).resolve()
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else here.parent.parent
    tests = root / "tests"
    print("=" * 78)
    print(f"ESTIF DIAGNOSTIC (read-only) -- v6.4.2 -- session 20260721 v2 -- repo root: {root}")
    print("=" * 78)

    head = git(root, "rev-parse", "--short", "HEAD")
    dirty = [l for l in git(root, "status", "--porcelain").splitlines() if l.strip()]
    print(f"\n[GIT] HEAD {head or 'n/a'} | uncommitted changes: {len(dirty)}")
    if dirty:
        print("      WARNING: tree not clean -- commit before/after each reorg step "
              "so every move is its own recoverable commit.")

    if not tests.is_dir():
        print(f"\nABORT: {tests} not found.")
        sys.exit(2)

    junk = [p for p in tests.iterdir() if p.is_file() and p.name == ".DS_Store"]
    loose = sorted(p for p in tests.iterdir() if p.is_file() and p.name != ".DS_Store")
    folders = sorted(p.name for p in tests.iterdir() if p.is_dir())
    by_ext = {}
    for p in loose:
        by_ext.setdefault(p.suffix or "(none)", []).append(p.name)
    print(f"\n[TESTS/ TOP LEVEL] {len(loose)} loose file(s), folders: {folders}")
    if junk:
        print(f"  plus {len(junk)} .DS_Store (macOS junk -- ignored here; reorg.sh deletes it)")
    for ext, names in sorted(by_ext.items()):
        head6 = ", ".join(names[:6]) + (" ..." if len(names) > 6 else "")
        print(f"  {ext:8s} x {len(names):2d}  {head6}")

    flags = [
        ("tests/docs_tests exists", (tests / "docs_tests").is_dir()),
        ("tests/docs exists", (tests / "docs").is_dir()),
        ("tests/scripts exists", (tests / "scripts").is_dir()),
        ("tests/plots exists", (tests / "plots").is_dir()),
        ("tests/_index_view exists", (tests / "_index_view").is_dir()),
        ("tests/phase2_a1prime (top level)", (tests / "phase2_a1prime").is_dir()),
        ("src/files.zip exists", (root / "src" / "files.zip").is_file()),
    ]
    print("\n[FLAGS]")
    for label, val in flags:
        print(f"  {'YES' if val else 'no '}  {label}")
    iv = tests / "_index_view"
    if iv.is_dir():
        n = sum(1 for _ in iv.rglob("*") if _.is_file())
        print(f"       _index_view holds {n} file(s) -- redundancy ruling: DELETE (reorg.sh does it)")

    scan_root_dirs = [tests]
    skip = {"_index_view", "audit", "__pycache__"}
    pys = sorted(p for p in tests.rglob("*.py")
                 if not (set(p.parts) & skip) and not p.name.startswith("._"))

    print("\n[PER-SCRIPT SCAN] (tests/ recursive; _index_view/ and audit/ skipped)")
    n_needfix, n_finalok, n_home, cache_ref_files = 0, 0, 0, 0
    for py in pys:
        d = py.parent
        parts_docs = os.path.relpath(str(tests / "docs"), str(d)).split(os.sep)
        parts_src = os.path.relpath(str(root / "src"), str(d)).split(os.sep)
        parts_root = os.path.relpath(str(root), str(d)).split(os.sep)
        try:
            raw = py.read_text(encoding="utf-8")
        except Exception as e:
            print(f"\n  {py.relative_to(tests.parent)}  (unreadable: {e})")
            continue
        lines = raw.splitlines()
        findings = []
        needs = []

        refs = sorted({b for b in MOVED_PREFIXES if b in raw} |
                      {fam.split('"')[1].split("'")[0].rstrip("{") for fam in DYNAMIC_FAMILIES
                       if fam in raw and "{key}" not in fam})
        dyn = any(fam in raw for fam in DYNAMIC_FAMILIES)
        if refs or dyn:
            cache_ref_files += 1
            tag = " + dynamic f-string cache names" if dyn else ""
            findings.append(("cache use", 0, f"{refs if refs else '(dynamic only)'}{tag}"))

        if "/home/claude" in raw:
            n_home += 1
            for k, l in enumerate(raw.splitlines(), 1):
                if "/home/claude" in l:
                    findings.append(("HOME-PATH", k, l.strip()[:100]))
                    needs.append(f"line {k}: /home/claude container leftover")

        for i, j, unit in logical_units(lines):
            n = i + 1
            if RE_CACHE_ANY.match(unit):
                if RE_CACHE_BARE.match(unit):
                    findings.append(("CACHE_DIR", n, unit.strip()[:100]))
                    findings.append(("", 0, f"   -> points at own folder -- NEEDS FIX "
                                            f"(uniform policy: all caches in tests/docs)"))
                    needs.append(f"line {n}: CACHE_DIR not final")
                else:
                    m = RE_CACHE_JOIN.match(unit)
                    na = norm_args(m.group("args")) if m else None
                    if na == parts_docs:
                        findings.append(("CACHE_DIR", n, "FINAL-OK (join -> tests/docs)"))
                        n_finalok += 1
                    else:
                        findings.append(("CACHE_DIR", n, unit.strip()[:100]))
                        findings.append(("", 0, f"   -> form/args {na} -- NEEDS FIX"))
                        needs.append(f"line {n}: CACHE_DIR not final")
                continue
            if RE_SP_ANY.search(unit):
                m = RE_SP_JOIN.match(unit)
                if m:
                    na = norm_args(m.group("args"))
                    ok = (na == parts_src) or (na == parts_root)
                    hint = "FINAL-OK" if ok else f"args {na} -- NEEDS FIX after move"
                    findings.append(("sys.path", n, f"{unit.strip()[:88]}  [{hint}]"))
                    if not ok:
                        needs.append(f"line {n}: sys.path hop not final")
                elif RE_SP_IDENT.match(unit):
                    findings.append(("sys.path", n, unit.strip()[:100] + "  [via variable]"))
                else:
                    findings.append(("sys.path", n, unit.strip()[:100] + "  [unrecognized form]"))
                    needs.append(f"line {n}: sys.path unrecognized")
                continue
            for mm in RE_IDENT_JOIN.finditer(unit):
                ident, inner, name = mm.group(1), mm.group(2), mm.group(4)
                if ident != "CACHE_DIR" and name.startswith(MOVED_PREFIXES):
                    inner_parts = RE_QUOTED.findall(inner)
                    ok = inner_parts == parts_docs
                    findings.append(("cache join", n,
                                     f"join({ident}, ..., '{name}')  "
                                     f"[{'FINAL-OK' if ok else 'NEEDS docs splice'}]"))
                    if not ok:
                        needs.append(f"line {n}: join({ident}, '{name}') not spliced to docs")

        for k, l in enumerate(lines, 1):
            if RE_SRC_IMPORT.match(l):
                findings.append(("src import", k, l.strip()[:100]))
            elif RE_TI.search(l) and RE_QUOTED.search(l) and "TEST_INDEX.md" in "".join(RE_QUOTED.findall(l)):
                findings.append(("TEST_INDEX", k, l.strip()[:100]))
            elif RE_OUT.search(l) and "savefig" in l:
                findings.append(("output", k, l.strip()[:100]))

        has_import = any(RE_SRC_IMPORT.match(l) for l in lines)
        if has_import and "sys.path" not in raw:
            findings.append(("WARNING", 0, "imports src/estif_ec but has NO sys.path hop -- relies on cwd"))

        if needs:
            n_needfix += 1
        if findings:
            print(f"\n  {py.relative_to(tests.parent)}")
            for tag, ln, msg in findings:
                if tag == "":
                    print(f"              {msg}")
                elif ln:
                    print(f"    {tag:11s}: line {ln}: {msg}")
                else:
                    print(f"    {tag:11s}: {msg}")

    print(f"\n[SCAN TOTALS] {len(pys)} .py scanned | {cache_ref_files} touch cached data | "
          f"{n_finalok} CACHE_DIR final | {n_home} with /home/claude leftovers | "
          f"NEED FIXING: {n_needfix}")

    print("\n[ARCHIVE DIAGNOSTICS] folders under archive/ with 'diagnostic' in the name")
    arch = root / "archive"
    live_dirs = [root / "src", tests, root / "docs", root / "results", root / "data"]
    if arch.is_dir():
        for d in sorted(arch.iterdir()):
            if d.is_dir() and "diagnostic" in d.name.lower():
                files = sorted(p.name for p in d.rglob("*") if p.is_file())
                print(f"  {d.relative_to(root)}/ -- {len(files)} file(s) (snapshot; audit/ lives inside -- untouched)")
    else:
        print("  (no archive/ folder)")

    print("\n[NAMED FILES] UKN / joint_calibration presence + last commit date")
    for pat in ("test_UKN.py", "test_UKN2.py", "test_joint_calibration.py",
                "test_joint_calibration_derived.py", "joint_calibration_results.png"):
        hits = [p for p in root.rglob(pat) if ".git" not in p.parts]
        if not hits:
            print(f"  {pat}: not present outside .git")
        for p in hits:
            when = git(root, "log", "-1", "--format=%as", "--", str(p.relative_to(root)))
            print(f"  {p.relative_to(root)}  last commit: {when or 'NO GIT HISTORY'}")

    docs = tests / "docs"
    checks = [
        ("zero loose files in tests/ (.DS_Store excluded)", len(loose) == 0),
        ("tests/ has scripts/ docs/ plots/ and nothing else",
         (tests / "scripts").is_dir() and docs.is_dir() and (tests / "plots").is_dir()
         and set(folders) <= {"scripts", "docs", "plots"} or
         {"scripts", "docs", "plots"} == set(p.name for p in tests.iterdir() if p.is_dir())),
        ("docs_tests/ gone", not (tests / "docs_tests").exists()),
        ("_index_view/ gone", not (tests / "_index_view").exists()),
        ("src/files.zip gone", not (root / "src" / "files.zip").exists()),
        ("no path construct left to fix (CACHE_DIR final, hops final, "
         "no /home/claude, no unspliced cache joins)", n_needfix == 0),
    ]
    print("\n[END-STATE CHECK] (meaningful after reorg + fix; PASS on all = done)")
    for label, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {label}")

    print("\nDiagnostic complete. Read-only -- nothing was modified.")


if __name__ == "__main__":
    main()
