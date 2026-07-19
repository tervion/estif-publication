# AUDIT NOTES — "The hidden catch in P2" session vs repo state
**Audit date:** 19 July 2026 · **Scope:** `audit/The hidden catch in P2/` (18 files + .DS_Store) vs `tests/`, `docs/plan/`, `docs/latex/`, the four session-folder passes, and the seven-document pass · **Method:** full read of every text artifact via Filesystem MCP; size-level diff against repo copies (byte-diff pending only if desired — sizes match exactly, the Jul 11/12 pattern); past-session retrieval for context. Written because this folder sat outside the "Jul N Session" naming convention and was missed by the entire audit corpus — the author flagged the omission on 19 Jul.

**Overall verdict: THE MISSING 13 JULY SESSION FOLDER, FOUND UNDER ITS SESSION TITLE — AND EVERYTHING IN IT IS ALREADY COMMITTED.** This folder is the origin record of RHAC-010 (the a₀-horizon doctrine). All nine text artifacts match their repo copies size-for-size; nothing is lost and no commit is missing. What the folder adds: session provenance for the five scripts previously flagged "no session copy exists anywhere" (a₀-horizon four + `ripple_speed.py`), the full statement of the P-2 catch, one **letter blocker absent from the entire audit corpus** (Bullet Cluster answer), one **Mac-mini portability bug** in `estif_flow_sim.py`, one open receipt slot (`ripple_speed.py` EoM placeholder), and two minor wording/bookkeeping items.

---

## A. What "P-2" is, and the catch itself

**P-2 = the second run at principle P**: the attempt to rebuild the a₀ story on √Λ alone — "matter genuinely gone." The receipts self-label it: `a0_horizon_test.py` — "reproducible receipt for the P-2 sqrt(Lambda) test"; `a0_prefactor_derivation.py` — "Step 1 of the P-2 test"; guide §2 — "(This is the P-2 rewrite, now confirmed.)"

**The catch, three layers (each on a runnable receipt):**

1. **The matter cannot be fired** (`a0_horizon_test.py`, guide §2). The identity c²√Λ = √3·√Ω_Λ·cH₀ proves any Λ-native form of a₀ is algebraically the Hubble form up to O(1); the matter dependence hides inside H₀ via Friedmann and is not removed by switching horizons. Consequences: a₀ anchors to c·H (de Sitter surface gravity); **P is demoted from foundation to "downstream curiosity"** — a stronger demotion than RHAC-009's closure (009 closed P-as-derivable; the guide removes P even as a premise for a₀); and the √3-from-isotropy story collapses — the √3 in the identity is the Friedmann conversion factor, the same number wearing two stories. This is the birth of RHAC-010's "working form, not derived."
2. **The number resists** (`a0_prefactor_derivation.py`, guide §3). Target is a **band** [0.118, 0.133], not the point 0.128 (that was one corner of the parameter grid); 1/(2π) = 0.159 and 1/6 both overshoot; the numerology floor is 1.3% (≈1 in 80 random O(1) expressions hits the band by luck — an in-band match alone is not evidence); the river-model edge gives 0.26–0.28·cH_Λ, ~20% high, and that 20% is exactly the ambiguity between acceleration definitions (proper vs surface-gravity vs flow-gradient). Inverse form: the edge law a₀ = cH_Λ·x_c needs x_c ≈ 0.221 vs the 0.272 in use → **the a₀-edge and the mode-crossover may be different surfaces** (guide §6 double-duty flag — the origin of the item CHANGELOG's [6.4.1] cross-references).
3. **The river cannot do it locally** (`estif_flow_sim.py` → `estif_horizon.py`, guide §4–§5). Direct simulation, both combination rules: the local flow yields (H/2)·v_gal — correct MOND-like direction, but ∝ M^(1/3) (fails mass-independence) and ~1000× too small; c enters the flow only at the horizon (r = c/H). **a₀ is non-local.** Injecting only the horizon acceleration a_H = k·cH turns on the full phenomenology at once: flat rotation curves, BTFR v⁴ = G·M·a_H with correct normalization, a₀ ∝ M⁰, magnitude cH/2π = 1.0×10⁻¹⁰ vs observed 1.2×10⁻¹⁰. Open problem 4 (guide §8) is the surviving prize: the **local mechanism** — what makes cH act at galaxy radii. Verlinde and Padmanabhan assume it; deriving it from the 4D inflow "would be ESTIF's genuine contribution." `ripple_speed.py` is the prepared instrument for that question.

**Cross-sleeve link recorded in the session (guide §7):** the Phase-2 dead constant Λ is exactly what a₀ = cH needs — constant Λ ⇒ fixed de Sitter horizon ⇒ constant a₀. The dark-energy null is load-bearing for the a₀ sleeve.

---

## B. Inventory and repo integrity (size-level)

| Folder artifact | Repo location | Size match |
|---|---|---|
| `a0_horizon_test.py` (2.75 KB) | `tests/` | ✅ identical |
| `a0_prefactor_derivation.py` (3.81 KB) | `tests/` | ✅ identical |
| `estif_flow_sim.py` (10.62 KB) | `tests/` | ✅ identical (carries the P2C-2 bug in both copies) |
| `estif_horizon.py` (7.91 KB) | `tests/` | ✅ identical |
| `ripple_speed.py` (4.44 KB) | `tests/` | ✅ identical |
| `ESTIF_document_update_guide.md` (10.87 KB) | `docs/plan/` | ✅ identical |
| `GAZTANAGA_COMPARISON.md` (10.56 KB) | `docs/plan/` | ✅ identical |
| `PHASE2_DECLARATION.md` (7.46 KB) | `docs/plan/` | ✅ identical |
| `estif_field_dynamics.tex` (8.09 KB) + preview PDF | `docs/latex/` | ✅ identical |
| `fig1–fig4` (flow_sim run) + `h1–h4` (horizon run) PNGs | — (session outputs) | preserved here only |

No missing commits. The figure set preserves both simulation runs' outputs — the only place they exist.

---

## C. Findings

### P2C-1 — HIGH · A letter blocker absent from the entire audit corpus: the Bullet Cluster answer
Guide §9 gates the McGaugh letter on three items: "(a) relativistic-consistency section written, (b) **Bullet Cluster answer prepared**, (c) Gaztañaga overlap checked." (a) is tracked corpus-wide; (c) is delivered (the memo, with its no-novelty verdict); **(b) appears nowhere else** — not in the five-document letter formulation, not in README's open problems, not in ROADMAP. The agreed letter-gating sentence must grow one clause, and a home for the answer must be chosen (letter appendix vs README FAQ vs a VALIDATION note). Rides D11.

### P2C-2 — MEDIUM-HIGH · Mac-mini portability bug in `estif_flow_sim.py` (both copies)
The script hardcodes its four figure outputs to `/home/claude/fig*.png` — a container path that does not exist on macOS, so the script **crashes on the Mac mini** as written (`estif_horizon.py`, by contrast, writes to the working directory and is fine). Since the repo copy is size-identical to the session copy, the bug is committed. Consequence: this RHAC-010 receipt cannot yet have passed a genuine mini run. Fix: four one-line path changes (relative paths, `estif_horizon.py` style). Should land before or with the D10 mini batch so the a₀-horizon four can be mini-verified as a set.

### P2C-3 — MEDIUM · `ripple_speed.py` is an instrument, not yet a receipt
Its machinery is validated (wave equation, Klein–Gordon, the acoustic/Unruh flow template all print correctly), but `estif_eom` is still the wave-equation **placeholder** — the ESTIF flow law's equation of motion has never been inserted. This is the concrete open task feeding open problem 4 (local mechanism) and the "stiffness = c²" grounding, and it connects to C-15 (the script's own comment: "these ripples ARE gravitational waves — the same object as your LISA prediction"). Even a clean speed-of-c result would not fix the prefactor (the script says so itself). TEST_INDEX role wording should read "instrument/template; EoM slot open," not imply a completed test. This is a research item for the v6.4 roadmap section (D1/D2 content branch), not doc repair.

### P2C-4 — MEDIUM · Session-provenance corrections to prior passes (good direction)
AUDIT_NOTES_JUL11 §C-3 and AUDIT_NOTES_JUL13 D-4 stated "no session-folder copy exists anywhere" for `ripple_speed.py` and the a₀-horizon four, and that `audit/Jul 13 Session/files/` being empty closed the session-folder route. **Falsified:** the copies exist here; the Jul 13 material lived under the session title. The git-only-dating list shrinks from eight items to three: `test_UKN.py`, `test_UKN2.py`, `test_joint_calibration_derived.py`. TEST_INDEX rows for the five scripts gain: "session copy: `audit/The hidden catch in P2/` (13 Jul), size-identical."

### P2C-5 — LOW · "Five-sector" vs four doors, in both copies of the declaration
PHASE2_DECLARATION §6 says "the five-sector census was re-run under A1′" while its own §3 defines four doors (shape, rate, slosh, swirl); RHAC-006 carries the same "five sectors" phrase. Decide the canonical count wording (four doors, with the flow split making slosh+swirl two of them — or name the fifth explicitly) and fix it in the same touch as the D4/UKN batch, since those exact lines are already being edited.

### P2C-6 — LOW · Gaztañaga memo item 4 remains open (research, sleeve 5)
The memo's required-action 4 — check numerically whether Λ = 4πG⟨ρ+3p⟩ reduces to the cH₀/2 pull condition — was explicitly left open ("verdict delivered; numerical cross-check item 4 remains"). If it reduces, even the "different mechanism" claim weakens; if it diverges, the divergence is the only place original Ωm-side content could live. Belongs on the roadmap, not the repair batch.

### P2C-7 — LOW · Folder naming caused the audit miss
`The hidden catch in P2` sits outside the `Jul N Session` convention, which is why four session passes and seven document passes never opened it. Options (author's call, per the intentional-naming rule): rename to `Jul 13 Session — The hidden catch in P2`, or drop a one-line pointer README inside `audit/Jul 13 Session/`. Either prevents a repeat.

---

## D. What this folder settles (for the record)

The RHAC-010 doctrine is fully receipted and the receipts are committed; the doctrine's three-part structure (target band / numerology floor / river edge) traces line-by-line to `a0_prefactor_derivation.py`; the §6 x_c double-duty flag and the §8 open-problems register are the current canonical statement of the gravity sleeve's challenges; the guide's §7 "do not touch" list confirms the empirical fits and the black-hole exterior were never in question. The catch is a discovery, not a failure: P-2 was run as an honest test, returned "scale yes, number no, locality no," and RHAC-010 recorded exactly that.

---

## E. Proposed actions — for decision, nothing executed

1. **P2C-1:** add the Bullet Cluster clause to the single letter formulation (now six documents + the letter itself) and choose the answer's home. Rides D11.
2. **P2C-2:** relativize the four `savefig` paths in `tests/estif_flow_sim.py`; then the a₀-horizon four join the D10 mini batch as a set.
3. **P2C-3:** schedule the `ripple_speed.py` EoM insertion as a sleeve-1 research item (v6.4 roadmap section, if the D1/D2 content branch is taken); reword its TEST_INDEX role to "instrument; EoM slot open."
4. **P2C-4:** fold the provenance corrections into the standing TEST_INDEX batch; strike the two "no session copy" claims in the Jul 11/13 notes with a dated annotation.
5. **P2C-5:** four-vs-five wording fix rides the D4/UKN batch.
6. **P2C-6:** Gaztañaga item 4 onto the roadmap (sleeve 5, low priority).
7. **P2C-7:** folder rename or pointer README — author's call.
8. **File this document** into `audit/` beside the other AUDIT_NOTES so the corpus is complete.

---

**Notes version:** 1.0 · 19 July 2026
