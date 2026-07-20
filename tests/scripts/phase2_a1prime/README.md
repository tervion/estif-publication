# Phase-2 & A1′ receipt suite

Constraint-algebra receipts for the residual-sector census (Phase 2)
and the strict-A1 growth no-go theorem plus its A1′ resolution.

Run with: `python3 <filename>.py` (requires `sympy`; script 6 also
uses `math`, standard library).

| # | File | Verdict | Weight |
|---|---|---|---|
| 1 | `estif_p2_door1_rate_dial.py` | Door 1 (rate) forbidden by A2; costs nothing | consistency check |
| 2 | `estif_p2_door2_slosh_divergence.py` | Door 2 (slosh) ledger = 0, cosmic average | theorem (strongest receipt) |
| 3 | `estif_p2_door3_swirl_ledger.py` | Door 3 (swirl) forbidden free / negative sourced | theorem for solenoidal class, one instance |
| 4 | `estif_growth_nogo_law.py` | strict-A1: δ=H/H₀ solves growth law; f=−0.911 vs DESI +0.76 | verifies solution + sign flip |
| 5 | `estif_growth_nogo_audit.py` | momentum gate + Hamiltonian slaving + closure | 4 of 5 audit attacks receipt-backed |
| 6 | `estif_a1prime_growth_restored.py` | fade mode confirmed; deepen mode f=0.76 (numeric) | numeric corroboration |
| 7 | `estif_a1prime_deepen_exact.py` | deepen mode exact symbolic proof | authoritative receipt |

## Known flags (do not overstate these results)

- **Script 1**: the fact that the ADM lapse carries no stress is a
  standard GR fact, *asserted* here, not derived by the script.
- **Script 3**: the ~10¹¹ suppression and dilution claims discussed
  around this result are order-of-magnitude estimates, not receipts.
- **Script 5**: the LTB (Lemaître–Tolman–Bondi) correspondence
  mentioned in discussion is literature-backed, not machine-verified.
- **Script 6**: does NOT check gravitational-wave propagation speed.
  The GW-speed-at-c claim is a separate, open derivation (item C-15),
  only *expected* from axiom A2, not derived in this suite.

## Expected outputs (for regression-checking future re-runs)

```
1: sqrt(-2*G*M/r + 1)  /  3*H**2/(8*pi*G)  /  3*H
2: 0  /  [0, 0, 0]
3: 0  /  Matrix([[0],[0],[0]])  /  Matrix([[A**2/2+C**2/2,0,0],[0,A**2/2+B**2/2,0],[0,0,B**2/2+C**2/2]])  /  -(A**2+B**2+C**2)/(32*pi*G)
4: 0  /  -0.911
5: [0, 0, 0]  /  0  /  0  /  0
6: 0  /  0.76
7: 0
```
