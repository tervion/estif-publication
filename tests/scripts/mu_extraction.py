import numpy as np
NMAX, B = 33.265, 15.429                       # joint EHT+Lambda calibration
n  = lambda x: NMAX*np.exp(-B*x)               # dynamic n(x)
F  = lambda x: x**(2*n(x))                     # (omega/H0)^2 = x^{2n(x)}

# Force law a = -(c^2/2) dF/dr, x = Rs/r  =>  g_obs/g_N = |dF/dx|
# READING A (n local, differentiate through):
#   dF/dx = x^{2n-1} * 2n * (1 - B*x*ln x)
nuA = lambda x: 2*n(x)*x**(2*n(x)-1)*(1 - B*x*np.log(x))
# READING B (n frozen at 1/2): dF/dx = 1  ->  pure Newton
# numeric derivative check of Reading A:
num = lambda x,h=1e-9: (F(x*(1+h))-F(x*(1-h)))/(2*x*h)

print("x = Rs/r      n(x)      2n     g_obs/g_N (A)   numeric check")
for x in [0.35, 0.30, 0.272, 0.25, 0.20, 0.15, 0.10, 0.05, 1e-2, 1e-3, 1e-4, 1e-6]:
    print(f"{x:9.3g} {n(x):9.4f} {2*n(x):7.3f} {nuA(x):15.6g} {num(x):15.6g}")

print()
print("MOND requires g_obs/g_N >= 1, RISING as field weakens (x -> 0).")
print("Galactic disks live at x ~ 1e-6 to 1e-8.")
