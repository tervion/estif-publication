# early_term_diagnostic.py

import estif_ec_fd_model as estif
import numpy as np

# Test S(t) at early times
times = np.logspace(10, 17, 20)  # From 10^10 to 10^17 seconds
for t in times:
    S = estif.global_S(t)
    H = estif.H_variable(t)
    print(f"t={t:.2e} s: S(t)={S:.6f}, H(t)={H:.3e}")