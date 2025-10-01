# estif_ec_fd_constants.py

import numpy as np

# Physical constants
G: float = 6.67430e-11     # Gravitational constant, m³ kg⁻¹ s⁻²
c = 2.99792458e8    # Speed of light, m/s

# Mass constants
M_sun = 1.989e30    # Solar mass, kg
M_earth = 5.972e24  # Earth mass, kg
M_proton = 1.67262192369e-27  # Proton mass, kg

# Distance constants
R_sun = 6.96e8      # Solar radius, m
R_earth = 6.371e6   # Earth radius, m
AU = 1.496e11       # Astronomical unit, m

# Cosmological constants
H_0: float = 2.1927e-18    # Hubble constant in s⁻¹ (67.66 km/s/Mpc from Planck 2018)

# --- ESTIF Specific Constants ---

# Friction Drag Constant (mass-based resistance in inward flow)
# Fitted to: BBN helium Y_p ~ 0.245, supernova χ² ~ 1.1, weak-field GR limits <1%
# Physical interpretation: Drag coefficient for mass resistance in 4D flow
BETA_DRAG = 0.05  # Dimensionless

# Early Surge Constant (controls rapid early-universe flow)
# Fitted to: BBN helium abundance and cosmological data
# Naturally produces CMB age ~ 377,000 years (matching ΛCDM ~380,000 years)
A_DEFAULT = 0.0005  # Units: dimensionless (multiplies 1/t^0.75 term)

# Physical interpretation of parameters:
# - H_0: Baseline inward flow rate (late universe, steady state)
# - A_DEFAULT: Early surge strength (rapid inward flow after Big Bang)
# - BETA_DRAG: Mass resistance coefficient (larger mass → larger eddy)

# Note: No Planck constant (ℏ) needed—ESTIF is purely classical

def schwarzschild_radius(mass):
    """Calculate Schwarzschild radius for given mass."""
    return 2 * G * mass / c**2

def escape_velocity(mass, radius):
    """Calculate escape velocity from surface."""
    return np.sqrt(2 * G * mass / radius)

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-01-10-25-V-2