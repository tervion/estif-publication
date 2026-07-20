echo "ESTIF reorg - v6.4.2 - session 20260721 v2 - move ONCE, then run fix_paths.py"
cd "${1:-/Users/peterangelov/estif_publication}"
echo "working in: $(pwd)"
test -d tests || echo "ABORT: no tests/ folder here - wrong directory, stopping now"
test -d tests || exit 1
test -d .git || echo "note: no .git found - will use plain mv instead of git mv"
echo "step 1 - docs_tests becomes docs (single data+index home)"
if test -d tests/docs_tests -a -d tests/docs; then echo "both exist - merging docs_tests contents into docs/"; for f in tests/docs_tests/*(N); do git mv "$f" tests/docs/ 2>/dev/null || mv "$f" tests/docs/; done; rmdir tests/docs_tests 2>/dev/null; fi
if test -d tests/docs_tests; then git mv tests/docs_tests tests/docs 2>/dev/null || mv tests/docs_tests tests/docs; fi
mkdir -p tests/docs tests/scripts tests/plots
echo "step 2 - every loose .py to scripts/"
for f in tests/*.py(N); do git mv "$f" tests/scripts/ 2>/dev/null || mv "$f" tests/scripts/; done
echo "step 3 - phase2_a1prime to scripts/phase2_a1prime (kept whole as a coherent receipt suite - confirmed 21 Jul)"
if test -d tests/phase2_a1prime; then git mv tests/phase2_a1prime tests/scripts/phase2_a1prime 2>/dev/null || mv tests/phase2_a1prime tests/scripts/phase2_a1prime; fi
echo "step 4 - loose images to plots/"
for f in tests/*.png(N) tests/*.jpg(N) tests/*.jpeg(N); do git mv "$f" tests/plots/ 2>/dev/null || mv "$f" tests/plots/; done
echo "step 5 - loose .md and .txt to docs/"
for f in tests/*.md(N) tests/*.txt(N); do git mv "$f" tests/docs/ 2>/dev/null || mv "$f" tests/docs/; done
echo "step 6 - delete _index_view (redundancy ruling - TEST_INDEX.md in docs/ is the single source of truth)"
if test -d tests/_index_view; then git rm -r -q --ignore-unmatch tests/_index_view 2>/dev/null; rm -rf tests/_index_view; fi
echo "step 7 - delete src/files.zip (decision D6)"
if test -f src/files.zip; then git rm -q --ignore-unmatch src/files.zip 2>/dev/null; rm -f src/files.zip; fi
echo "step 8 - remove macOS Finder junk under tests (.DS_Store)"
find tests -name ".DS_Store" -type f -delete
echo "final check - loose files left in tests/ (this line must be followed by NOTHING):"
find tests -maxdepth 1 -type f
echo "tests/ now contains:"
ls tests
echo "reorg done - NEXT: python3 _session_20260721/fix_paths.py - then commit this as its own sequential commit"
