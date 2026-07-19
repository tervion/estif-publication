"""
Phase 2 -- Door 2 (curl-free slosh perturbation) receipt.
Verdict: ZERO. For ANY symbolic scalar field phi (arbitrary curl-free
perturbation), the full Hamiltonian-constraint ledger entry (linear +
quadratic terms together) is identically a pure divergence -- so its
cosmic average is exactly zero. Genuine theorem for the whole scalar
class, not an example. Caveat: surface terms vanish for periodic or
decaying fields, which is the correct setting for a cosmic average.

Also verified: the momentum constraint is satisfied identically for
ANY curl-free perturbation (slosh always passes the entry gate).

Status: theorem, strongest receipt in the suite.
Expected output:
    0
    [0, 0, 0]
"""
import sympy as sp

x1, x2, x3, H = sp.symbols('x1 x2 x3 H', real=True)
X = (x1, x2, x3)
phi = sp.Function('phi')(*X)                              # ARBITRARY slosh field
S = sp.Matrix(3, 3, lambda i, j: sp.diff(phi, X[i], X[j]))  # its K-perturbation
lap = S.trace()

# Door-2 ledger entry: total (background+slosh) minus background
resid = (3*H + lap)**2 - (3*H**2 + 2*H*lap +
                           sum(S[i, j]**2 for i in range(3) for j in range(3))) - 6*H**2

# claim: the WHOLE entry is a divergence  ->  cosmic average = 0
F = [4*H*sp.diff(phi, X[i]) + sp.diff(phi, X[i])*lap
     - sum(sp.diff(phi, X[j])*S[i, j] for j in range(3)) for i in range(3)]
print(sp.simplify(resid - sum(sp.diff(F[i], X[i]) for i in range(3))))    # -> 0

# momentum constraint: ANY slosh passes, identically
print([sp.simplify(sum(sp.diff(S[i, j] - sp.eye(3)[i, j]*lap, X[j]) for j in range(3)))
       for i in range(3)])  # -> [0, 0, 0]
