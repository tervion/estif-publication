# estif_ec_gr_constants.py

"""
ESTIF Physical Constants (v6.4.1 — July 2026)

Contains all fundamental constants used by the ESTIF model.
Grouped by type: physical, cosmological, ESTIF-derived.
All numerical values are unchanged since v6.2 — this file is part of the
frozen v6.2 analytical receipt suite. Claim-status labels updated to v6.4.1.

Key derived constants:
    x₀ = R_H / R_UNIVERSE_0 = 0.3107 ≈ Ωm  (consistency relation, C2 — not an
        independent prediction: R_UNIVERSE_0 is the ΛCDM particle horizon)
    a₀ = H_0 × c × x₀ / √3 = 1.179×10⁻¹⁰ m/s²  (horizon-scale acceleration,
        a₀ ≈ c·H; the √3 prefactor is a working form, NOT derived — 13 Jul 2026)
    N_MAX ≈ 5/7 × ln(r_e / l_P)  (electron radius connection; 5/7 conditional
        on the open x_c = 0.272 derivation)
"""

import numpy as np

# ============================================================================
# FUNDAMENTAL PHYSICAL CONSTANTS
# ============================================================================

G: float = 6.67430e-11         # Gravitational constant [m³ kg⁻¹ s⁻²]
c:  float = 2.99792458e8       # Speed of light [m/s]

# ============================================================================
# MASS CONSTANTS
# ============================================================================

M_sun   = 1.989e30             # Solar mass [kg]
M_earth = 5.972e24             # Earth mass [kg]
M_proton= 1.67262192369e-27    # Proton mass [kg]
M_electron = 9.10938e-31       # Electron mass [kg]

# ============================================================================
# DISTANCE CONSTANTS
# ============================================================================

R_sun   = 6.96e8               # Solar radius [m]
R_earth = 6.371e6              # Earth radius [m]
AU      = 1.496e11             # Astronomical unit [m]
KPC     = 3.085677581e19       # 1 kiloparsec [m]
MPC     = 3.085677581e22       # 1 megaparsec [m]

# Classical electron radius — the natural scale for N_MAX and B
# r_e = e²/(4πε₀ m_e c²) — where EM self-energy equals rest mass energy
# This is the boundary between electromagnetism and gravity
R_ELECTRON = 2.8179403227e-15  # Classical electron radius [m]

# Planck length — the quantum gravity scale
L_PLANCK = 1.616255e-35        # Planck length [m]

# ============================================================================
# COSMOLOGICAL CONSTANTS
# ============================================================================

# Hubble constant — Planck 2018 (67.66 km/s/Mpc)
H_0: float = 2.1927e-18       # Hubble constant [s⁻¹]

# Observable universe radius today — the ΛCDM particle horizon (an import
# that itself contains Ωm; hence Ωm = x₀ is a consistency relation, C2)
R_UNIVERSE_0 = 4.4e26          # Observable universe radius [m]

# ============================================================================
# ESTIF FORMULA PARAMETERS
# ============================================================================

# Calibrated parameters — do NOT modify without re-running test_joint_calibration.py
# All three (EHT, Planck Λ, LISA) must pass simultaneously after any change.
N_MAX_COMBINED = 33.265        # Tilt exponent in flat space
B_COMBINED     = 15.429        # Exponential decay rate

# Connection to electron radius (confirmed to 0.08% and 0.69%):
# N_MAX ≈ 5/7 × ln(R_ELECTRON / L_PLANCK)
# B     ≈ 1/3 × ln(R_ELECTRON / L_PLANCK)
# ln(r_e / l_P) ≈ 46.608

# ============================================================================
# ESTIF COSMOLOGICAL PARAMETERS
# ============================================================================

# Planck 2018 matter density — equals x₀ to 0.12%
OMEGA_M      = 0.3111          # Total matter density (Ωb + Ωdm)
OMEGA_LAMBDA = 0.6889          # Dark energy density (imported Λ; the Ω_tilt replacement was retired in v6.3)
OMEGA_B      = 0.049           # Baryonic matter density (BBN measured)
OMEGA_DM     = 0.262           # Dark matter density (= x₀ − Ωb to 0.10%)

# Tilt exponent for dark energy evolution — RETIRED sector (v6.3): the
# Ω_tilt(z) cosmology claim is withdrawn; the constant is kept for the
# historical receipts.
# Exact formula: x(z) = x₀ × (1+z) × H₀/H_ΛCDM(z) gives α ≈ 0.077–0.089
ALPHA_COSMO  = 0.1036          # Best-fit tilt exponent (within 2σ of geometric range)

# ============================================================================
# ESTIF DERIVED CONSTANTS
# ============================================================================

# Cosmological curvature ratio x₀ = R_H / R_universe
# Equals Ωm to 0.12% — a consistency relation (C2), not an independent prediction
# x₀ = (c/H₀) / R_UNIVERSE_0
X_0 = (c / H_0) / R_UNIVERSE_0   # ≈ 0.3107

# Critical density today
RHO_CRIT_0 = 3 * H_0**2 / (8 * np.pi * G)   # ≈ 8.60e-27 kg/m³

# Eddy background density = x₀ × ρ_crit (equals Ωm × ρ_crit)
RHO_EDDY_0 = X_0 * RHO_CRIT_0               # ≈ 2.67e-27 kg/m³

# MOND acceleration constant — horizon-scale form (doctrine of 13 Jul 2026):
# a₀ is a horizon quantity, a₀ ≈ c·H (de Sitter surface-gravity scale). Flat
# rotation curves, mass-independence, the BTFR and the magnitude follow from
# that alone. The exact O(1) prefactor (written here as x₀/√3) is a working
# form, NOT derived. Numerical value unchanged: 1.72% from MOND's 1.2×10⁻¹⁰.
A0_MOND_ESTIF = H_0 * c * X_0 / np.sqrt(3)  # ≈ 1.179×10⁻¹⁰ m/s²

# ============================================================================
# BETA_DRAG (Legacy — For Compatibility Only)
# ============================================================================
# In v4.0+, tilt suppression is purely geometric via beta_combined(x).
# BETA_DRAG = 1.0 represents the unsuppressed limit.
# The actual observable correction is position-dependent: √β(x).
BETA_DRAG = 1.0

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def schwarzschild_radius(mass):
    """
    Schwarzschild radius: Rs = 2GM/c²

    Args:
        mass: Mass [kg]
    Returns:
        Schwarzschild radius [m]
    """
    return 2 * G * mass / c**2


def escape_velocity(mass, radius):
    """
    Escape velocity: v_esc = √(2GM/r)

    Args:
        mass: Mass [kg]
        radius: Radius [m]
    Returns:
        Escape velocity [m/s]
    """
    return np.sqrt(2 * G * mass / radius)


def hubble_radius():
    """
    Current Hubble radius: R_H = c/H₀

    Returns:
        Hubble radius [m]
    """
    return c / H_0


def x0_value():
    """
    Cosmological curvature ratio x₀ = R_H / R_universe.
    Equals Ωm to 0.12% — a consistency relation (C2), not an
    independent prediction.

    Returns:
        x₀ (dimensionless)
    """
    return X_0


# ============================================================================
# VALIDATION — print key derived values on import (if verbose)
# ============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("ESTIF Constants — Derived Values")
    print("=" * 60)
    print(f"  x₀ = R_H/R_universe  = {X_0:.6f}")
    print(f"  Ωm (Planck 2018)     = {OMEGA_M:.6f}")
    print(f"  Agreement:             {abs(X_0-OMEGA_M)/OMEGA_M*100:.4f}%")
    print(f"  x₀ − Ωb              = {X_0 - OMEGA_B:.6f}")
    print(f"  Ωdm (Planck 2018)    = {OMEGA_DM:.6f}")
    print(f"  Agreement:             {abs(X_0-OMEGA_B-OMEGA_DM)/OMEGA_DM*100:.4f}%")
    print()
    print(f"  a₀ (ESTIF/√3)        = {A0_MOND_ESTIF:.4e} m/s²")
    print(f"  a₀ (MOND empirical)  = 1.2000e-10 m/s²")
    print(f"  Agreement:             {abs(A0_MOND_ESTIF/1.2e-10-1)*100:.2f}%")
    print()
    ln_re_lP = np.log(R_ELECTRON / L_PLANCK)
    print(f"  ln(r_e/l_P)          = {ln_re_lP:.3f}")
    print(f"  5/7 × ln(r_e/l_P)   = {5/7*ln_re_lP:.3f}  (N_MAX = {N_MAX_COMBINED})")
    print(f"  1/3 × ln(r_e/l_P)   = {1/3*ln_re_lP:.3f}  (B = {B_COMBINED})")
