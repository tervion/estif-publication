# debug_high_z.py

import estif_ec_fd_model as estif
import numpy as np

z_problem = 1100
t_obs = 4.35e17

print("Debugging z=2.0 failure:")
t_emit = estif.t_from_z(z_problem, t_obs)
print(f"t_emit = {t_emit:.3e} s")

S_emit = estif.global_S(t_emit)
S_obs = estif.global_S(t_obs)
print(f"S_emit = {S_emit:.6f}, S_obs = {S_obs:.6f}")

z_check = S_emit/S_obs - 1
print(f"Recovered z = {z_check:.4f} (expected {z_problem})")

# Check if H(t) is reasonable
H_emit = estif.H_variable(t_emit)
print(f"H(t_emit) = {H_emit:.3e} s⁻¹")
print(f"H₀ = {estif.const.H_0:.3e} s⁻¹")
print(f"Ratio H(t_emit)/H₀ = {H_emit/estif.const.H_0:.2e}")