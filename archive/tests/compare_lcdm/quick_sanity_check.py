# quick_sanity_check.py

import estif_ec_fd_model as estif

z_test = 0.5
mu = estif.distance_modulus_estif_numerical(z_test)
print(f"At z={z_test}: μ = {mu:.2f} mag")

# Expected: ~42 mag (reasonable for z=0.5)
if 35 < mu < 50:
    print("✅ Sanity check PASSED - value in reasonable range")
else:
    print("❌ RED FLAG - value outside expected range!")