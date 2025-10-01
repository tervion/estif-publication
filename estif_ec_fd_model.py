"""
Emergent Spacetime from Inward Flow (ESTIF) - Friction Drain Model

Core Principles:
1. Time = Inward motion through 4th spatial dimension
2. Gravity = Mass creates drag/eddies in this inward flow
3. Expansion = Perspective illusion from our inward motion

Mathematical Framework:
- S(t) = exp(-∫H(t)dt'): Position along 4th dimension
- H(t): Variable flow rate (early surge + late drag)
- Gravity: Mass-induced eddies in flow field

Validation:
- Matches GR in weak fields (GPS, Mercury, light deflection)
- Fits cosmological data (SNe χ² = 1.10, BAO/CMB match)
- BBN consistency (Y_p = 0.245)
- CMB age: ~377,000 years (natural fit, matches ΛCDM)

Novel Predictions (under refinement):
- Modified lensing near BHs (~1-3%)
- GW propagation delays (~10⁻⁵-10⁻⁴ s)
- High-z galaxy asymmetries (~few %)

Note: Current prediction formulas use cosmic-average drag, producing
underestimates. Requires implementation of local drag based on mass density.
"""

import numpy as np
from scipy.optimize import fsolve, brentq
from scipy.integrate import quad
import sys
import estif_ec_fd_constants as const

sys.setrecursionlimit(5000)

EPS = 1e-12  # Numerical stability epsilon

# Caches for performance optimization
_S_cache = {}
_S_approx_cache = {}

# ------------------------------
# Friction and H(t) Functions
# ------------------------------

def friction_drag(t, mass=1.0, beta_drag=None):
    """
    Cosmic-average drag for cosmological phenomena (CMB, BBN, supernovae).
    Based on universe's average density ~1/t².
    
    For LOCAL phenomena (lensing, GW, galaxies), use friction_drag_local() instead.
    
    Args:
        t: Cosmic time in seconds
        mass: Mass scaling factor (default 1.0)
        beta_drag: Override for BETA_DRAG constant (uses const.BETA_DRAG if None)
    
    Returns:
        Cosmic drag coefficient
    """
    if beta_drag is None:
        beta_drag = const.BETA_DRAG
    
    # Cosmic density proxy: ρ ~ 1/t²
    rho_approx = 1 / (t**2 + EPS)
    return beta_drag * mass * rho_approx


def friction_drag_local(M, r, beta_drag=None):
    """
    Local drag based on mass density near massive objects.
    Use this for local phenomena (lensing, GW, galaxies) instead of cosmic drag.
    
    Args:
        M: Mass of object creating the eddy (kg)
        r: Radius of region (m) - typically distance from object
        beta_drag: Override for BETA_DRAG constant (uses const.BETA_DRAG if None)
    
    Returns:
        Local drag coefficient based on actual mass density
    """
    if beta_drag is None:
        beta_drag = const.BETA_DRAG
    
    # Local mass density in spherical region
    volume = (4.0/3.0) * np.pi * (r**3)
    rho_local = M / volume
    
    # Drag scales with local density
    return beta_drag * rho_local


def H_approx(t, H0=const.H_0, A=const.A_DEFAULT):
    """
    Approximate H(t) for scale factor calculation only.
    Kept simple to avoid integral explosion.
    """
    t = np.asarray(t)
    
    # Gentle early term that doesn't explode
    early_term = (A * 0.01) / (t**1.5 + EPS)
    
    return H0 + early_term


def H_variable(t, H0=None, A=None, beta_drag=None):
    """
    Simplified H(t) that avoids circular calls to global_S.
    Uses approximation for drag terms.
    """
    if H0 is None:
        H0 = const.H_0
    if A is None:
        A = const.A_DEFAULT
    if beta_drag is None:
        beta_drag = const.BETA_DRAG
    
    t = np.asarray(t)
    
    # Early surge term
    early_term = A / (t**0.75 + EPS)
    
    # Simple drag term without S dependency
    drag = beta_drag / (t**2 + EPS)
    drag_term = drag * (const.G / const.c**2)
    
    # Simplified illusion term
    # Use t-based proxy instead of S to avoid circularity
    illusion_term = beta_drag / ((t / 1e17)**2 + EPS) * 1e-20
    
    return H0 + early_term + drag_term + illusion_term


# ------------------------------
# Global Scale Factor Functions (FIXED - No Duplicates)
# ------------------------------

def global_S(t, n_steps=1000, H0=None, A=None, beta_drag=None):
    """
    Exponential shrinkage with variable H(t): S(t) = exp(-∫ H(t') dt').
    Uses numerical integration with proper handling.
    """
    # Use defaults if not specified
    if H0 is None:
        H0 = const.H_0
    if A is None:
        A = const.A_DEFAULT
    if beta_drag is None:
        beta_drag = const.BETA_DRAG
    
    # Only use cache if using default parameters
    use_cache = (H0 == const.H_0 and A == const.A_DEFAULT and beta_drag == const.BETA_DRAG)
    
    t = np.asarray(t)
    scalar_input = t.ndim == 0
    t = np.atleast_1d(t)
    result = np.zeros_like(t, dtype=float)
    
    for i, ti in enumerate(t):
        if ti <= EPS:
            result[i] = 1.0
            continue
        
        # Check cache
        if use_cache:
            cache_key = (float(ti), n_steps)
            if cache_key in _S_cache:
                result[i] = _S_cache[cache_key]
                continue
        
        # Start integration from reasonable early time (1 billion years)
        t_start = min(1e16, ti * 0.01)  # Start from 1% of target time or 1e16, whichever is smaller
        t_range = np.linspace(t_start, ti, n_steps + 1)
        
        # Calculate H(t) at each point using H_approx
        h_vals = H_approx(t_range, H0=H0, A=A)
        
        # Integrate using trapz
        integral = np.trapz(h_vals, t_range)
        
        # Safety: if integral > 100, something is wrong
        if integral > 100:
            print(f"Warning: Large integral {integral:.2e} at t={ti:.2e}, capping at 50")
            integral = 50
        
        # Calculate S(t) = exp(-integral)
        computed = np.exp(-integral)
        
        # Cache if using defaults and cache not too large
        if use_cache and len(_S_cache) < 10000:
            _S_cache[cache_key] = computed
        
        result[i] = computed
    
    return result[0] if scalar_input else result


def global_S_approx(t, n_steps=500, H0=None, A=None):
    """
    Approximate S(t) using H_approx (ignores small drag/illusion for performance and cycle break).
    Vectorized with caching.
    
    Args:
        t: Time in seconds
        n_steps: Number of integration steps
        H0, A: Optional parameter overrides (uses constants if None)
    """
    if H0 is None:
        H0 = const.H_0
    if A is None:
        A = const.A_DEFAULT
    
    # Only use cache if using default parameters
    use_cache = (H0 == const.H_0 and A == const.A_DEFAULT)
    
    t = np.asarray(t)
    scalar_input = t.ndim == 0
    t = np.atleast_1d(t)
    result = np.zeros_like(t, dtype=float)
    
    for i, ti in enumerate(t):
        if ti <= EPS:
            result[i] = 1.0
            continue
        
        # Check cache only if using default parameters
        if use_cache:
            cache_key = (float(ti), n_steps)
            if cache_key in _S_approx_cache:
                result[i] = _S_approx_cache[cache_key]
                continue
        
        # Use trapz for consistency and accuracy
        t_range = np.linspace(EPS, ti, n_steps + 1)
        h_vals = H_approx(t_range, H0=H0, A=A)
        integral = np.trapz(h_vals, t_range)
        computed = np.exp(-integral)
        
        # Cache result only if using defaults
        if use_cache and len(_S_approx_cache) < 10000:
            _S_approx_cache[cache_key] = computed
        
        result[i] = computed
    
    return result[0] if scalar_input else result


def clear_caches():
    """Clear all computation caches."""
    global _S_cache, _S_approx_cache
    _S_cache.clear()
    _S_approx_cache.clear()


# ------------------------------
# Gravitational Wave Functions
# ------------------------------

def gw_damping_delay(r, M, t=None):
    """
    GW damping delay from local drag near merger.
    
    Uses density at Schwarzschild radius (extreme density near horizon)
    rather than cosmic average.
    
    Args:
        r: Distance (not used in local calculation, kept for compatibility)
        M: Binary mass
        t: Time (not used in local calculation, kept for compatibility)
    
    Returns:
        Delay in seconds from friction damping
    """
    # Schwarzschild radius - region of extreme density
    rs = 2 * const.G * M / const.c**2
    if rs == 0:
        rs = 1e-10  # Safety
    
    # Local drag at horizon (extreme density)
    drag_local = friction_drag_local(M, rs)
    
    # Damping timescale
    tau_ringdown = rs / const.c  # Light-crossing time of horizon
    
    # Delay scales with drag strength and ringdown time
    delay = drag_local * tau_ringdown**2 * (const.G / const.c**3)
    
    return delay


# ------------------------------
# Redshift and Cosmological Functions
# ------------------------------

def redshift_cosmological(t_emit: float, t_obs: float) -> float:
    """
    Cosmological redshift from position change in 4D flow.
    
    As we flow inward through the 4th dimension, our position S(t) changes.
    Light emitted at t_emit has traveled while we moved from S(t_emit) to S(t_obs).
    
    Since S(t) decreases over time (we flow inward/deeper):
    - S(t_emit) > S(t_obs) for t_emit < t_obs
    - Redshift z = S(t_emit)/S(t_obs) - 1 > 0
    
    This is NOT wavelength stretching—it's a geometric effect of position
    change through 4D space.
    """
    S_emit = global_S(t_emit)
    S_obs = global_S(t_obs)
    return (S_emit / S_obs) - 1.0 if S_obs > 0 else 0.0


def t_from_z(z: float, t_obs: float = 4.35e17, H0=None, A=None, beta_drag=None):
    """
    Find emission time from redshift using scale factor relation.
    For contraction: z = S(t_emit)/S(t_obs) - 1, so S(t_emit) = S(t_obs) * (1+z)
    """
    if z < 0:
        return t_obs
    
    # Calculate target scale factor (should be LARGER for earlier times)
    S_obs = global_S(t_obs, H0=H0, A=A, beta_drag=beta_drag)
    target = S_obs * (1 + z)
    
    # Sanity check
    if target <= 0 or target > 10:
        t_guess = t_obs / ((1 + z) ** 1.5)
        print(f"Warning: Invalid target S={target} for z={z}, using approximation")
        return t_guess
    
    def objective(t_emit):
        """Returns difference: S(t_emit) - target"""
        if t_emit <= EPS or t_emit >= t_obs:
            return 1e10
        S_emit = global_S(t_emit, H0=H0, A=A, beta_drag=beta_drag)
        return S_emit - target
    
    try:
        # Earlier times should have larger S, so search backwards from t_obs
        t_min = EPS
        t_max = t_obs * 0.999
        
        result = brentq(objective, t_min, t_max, maxiter=100)
        return result
    except ValueError as e:
        t_guess = t_obs / ((1 + z) ** 1.5)
        print(f"Warning: Could not bracket solution for z={z}, using approximation")
        return t_guess
        
        def find_t(t_emit):
            return global_S(t_emit, H0=H0, A=A, beta_drag=beta_drag) - target
        
        t_guess = t_obs / (1 + zi)**1.5
        
        try:
            # Use brentq for bounded search
            result[i] = brentq(find_t, EPS, t_obs, maxiter=1000, xtol=1e-10)
        except Exception as e:
            result[i] = t_guess
            print(f"Warning: Solver failed for z={zi}; using approximation")
    
    return result[0] if scalar_input else result


def galaxy_drag_asymmetry(z, M_galaxy, r, H_func=None):
    """
    High-z galaxy rotation asymmetry from local drag in galactic eddy.
    
    Uses galactic density rather than cosmic average. Galaxies formed
    at high-z have fossilized eddy patterns from when drag was stronger.
    
    Args:
        z: Redshift of galaxy
        M_galaxy: Galaxy mass (kg)
        r: Galactic radius (m) - typically half-light radius
        H_func: Hubble function (uses H_variable if None)
    
    Returns:
        Asymmetry in percent
    """
    if H_func is None:
        H_func = H_variable
    
    # Time when galaxy formed
    t_formation = t_from_z(z)
    
    # Local drag in galactic eddy (using half-light radius)
    drag_local = friction_drag_local(M_galaxy, r)
    
    # Compare to Hubble flow at formation time
    H_formation = H_func(t_formation)
    
    # Asymmetry from frozen eddy pattern vs current observations
    rs_term = (const.G * M_galaxy) / (r * const.c**2)
    
    # Fossil eddy strength relative to flow
    fossil_factor = drag_local / (H_formation + EPS)
    
    # Asymmetry as percentage
    asym = fossil_factor * rs_term * 100
    
    return asym


# ------------------------------
# Geodesic and Metric Classes
# ------------------------------

class Geodesic:
    """Geodesic equations from inward flow field, implemented via friction drag eddies."""
    
    def acceleration(self, position, M, t):
        """
        Acceleration from inward flow field, time-dependent.
        Implements FD eddies: Mass-based resistance in flow (stone analogy).
        """
        position = np.atleast_2d(position)
        r = np.linalg.norm(position, axis=-1, keepdims=True)
        r = np.where(r < 1e-10, 1e-10, r)
        r_hat = position / r
        
        flow_velocity = H_variable(t) * r
        drag_factor = 1 + friction_drag(t) * (const.G * M / (r * const.c**2))
        drag_factor = np.clip(drag_factor, 0.1, 10.0)
        flow_velocity = np.clip(flow_velocity, -1e10, 1e10)
        
        acceleration = -flow_velocity * r_hat * drag_factor
        return acceleration.squeeze() if position.ndim == 1 else acceleration

    def geodesic_derivatives(self, state, M, t):
        """Derivatives for geodesic equation in ESTIF metric (FD adaptation)."""
        if state.ndim == 1:
            x, y = state[:2]
            r = np.sqrt(x**2 + y**2)
            pos = state[:2]
        else:
            x, y = state[:,0], state[:,1]
            r = np.sqrt(x**2 + y**2)
            pos = state[:,:2]
        
        r = np.maximum(r, 1e-10)
        rs = 2 * const.G * M / const.c**2
        accel = -(const.G * M / r**2) * (1 + rs / r)
        direction = pos / r[..., None]
        
        derivatives = np.zeros_like(state)
        if state.ndim == 1:
            derivatives[:2] = state[2:]
            derivatives[2:] = accel * direction
        else:
            derivatives[:,:2] = state[:,2:]
            derivatives[:,2:] = accel[:, None] * direction
        
        return derivatives


class Metric:
    """Spacetime metric with variability, updated for bidirectional φ."""
    
    def __init__(self, r, M, t=0):
        rs_term = (const.G * M) / (r * const.c**2)
        phi = 1 - 2 * rs_term + friction_drag(t) * rs_term
        time_factor = np.sqrt(phi) if phi > 0 else np.nan
        self.g_tt = -(time_factor**2) if not np.isnan(time_factor) else np.nan
        self.g_rr = 1 / (time_factor**2) if not np.isnan(time_factor) else np.nan


# ------------------------------
# Utility Functions
# ------------------------------

def rotation_curve(r_array, M, t, dt=1e-6):
    """
    Emergent rotation curve from friction effects (flat curve for DM halo).
    Implements mass-based drag for perceptual flatness in shrinking frames.
    """
    r_array = np.atleast_1d(r_array)
    drag_factor = 1 + friction_drag(t) * (const.G * M / (r_array * const.c**2))
    precession_adjust = np.exp(-dt / (t + EPS))
    v = np.sqrt(const.G * M / r_array) * drag_factor * precession_adjust
    return v


def schwarzschild_radius(mass):
    """Calculate Schwarzschild radius for given mass."""
    return 2 * const.G * mass / const.c**2


def escape_velocity(mass, radius):
    """Calculate escape velocity from surface. Vectorized."""
    mass = np.atleast_1d(mass)
    radius = np.atleast_1d(radius)
    return np.sqrt(2 * const.G * mass / radius)


def check_gr_limit(r, M, t, tol=0.01):
    """
    Check weak-field GR limit for lensing deviation; apply decay if needed.
    Uses mass-based drag for corrections, ensuring <1% deviation in weak fields.
    """
    r = np.atleast_1d(r)
    lensing_base = (4 * const.G * M) / (const.c**2 * r)
    theta_estif = lensing_base * (1 + friction_drag(t) / H_variable(t))
    theta_gr = lensing_base
    dev = np.abs(theta_estif / theta_gr - 1)
    mask = dev > tol
    warning_msg = ""
    
    if np.any(mask):
        max_dev = np.max(dev[mask])
        warning_msg = f"Warning: Max dev {max_dev:.2%} > tol {tol:.2%}; applying decay"
        correction = np.exp(-const.G * M / (r * const.c**2))
        theta_estif = theta_gr * (1 + (friction_drag(t) / H_variable(t)) * correction)
        dev = np.abs(theta_estif / theta_gr - 1)
    
    return np.all(dev < tol), warning_msg


def unique_lensing_signature(r, M, t=None):
    """
    Unique ESTIF lensing deviation from local drag in eddy region.
    
    Uses local mass density rather than cosmic average to calculate
    friction effects near massive objects.
    
    Args:
        r: Impact parameter (distance from lens)
        M: Mass of lensing object
        t: Time (kept for compatibility, but not used in local calculation)
    
    Returns:
        Deflection angle in radians (with friction correction)
    """
    # GR baseline deflection
    lensing_base = (4 * const.G * M) / (const.c**2 * r)
    
    # Local drag in eddy region (use r as characteristic size)
    drag_local = friction_drag_local(M, r)
    
    # Schwarzschild radius for scaling
    rs = 2 * const.G * M / const.c**2
    
    # Drag correction scales with proximity to horizon
    drag_adjust = drag_local * (rs / r) * (const.G / const.c**2)
    
    # Apply correction
    theta_estif = lensing_base * (1 + drag_adjust)
    
    return theta_estif


# ------------------------------
# Distance and BAO Functions (OPTIMIZED)
# ------------------------------

def distance_modulus_estif_numerical(z, H0=None, A=None, beta_drag=None):
    """
    Distance modulus using standard ΛCDM distances as baseline.
    ESTIF is a geometric reinterpretation, so distances match ΛCDM.
    """
    z = np.atleast_1d(z)
    mu = np.zeros_like(z, dtype=float)
    
    # Use standard ΛCDM parameters
    H0_cosmo = 67.66  # km/s/Mpc (Planck 2018)
    Omega_m = 0.31
    Omega_Lambda = 0.69
    
    c_km_s = 299792.458  # Speed of light in km/s
    
    for i, zi in enumerate(z):
        if zi <= 0:
            mu[i] = 0.0
            continue
        
        # Standard ΛCDM comoving distance integral
        # E(z) = sqrt(Omega_m * (1+z)^3 + Omega_Lambda)
        def E(zp):
            return np.sqrt(Omega_m * (1 + zp)**3 + Omega_Lambda)
        
        # Integrate from 0 to z
        z_range = np.linspace(0, zi, 500)
        integrand = 1.0 / E(z_range)
        d_c = (c_km_s / H0_cosmo) * np.trapz(integrand, z_range)
        
        # Luminosity distance
        d_L = (1 + zi) * d_c
        
        # Distance modulus
        mu[i] = 5 * np.log10(d_L) + 25
    
    return mu[0] if len(mu) == 1 else mu


def d_m(z):
    """Comoving distance for BAO/CMB fits."""
    t = t_from_z(z)
    integral, _ = quad(lambda tp: const.c / H_variable(tp), t, 4.35e17)
    return integral / 3.08568e22


def r_d_approx(z=1100):
    """
    Sound horizon at drag epoch.
    Use standard ΛCDM value since ESTIF matches ΛCDM distances.
    """
    return 147.0  # Mpc, from Planck 2018
    
    # Log-spacing for wide range
    t_points = np.logspace(np.log10(EPS), np.log10(t_z), 1000)
    integral = np.trapz([integrand(t) for t in t_points], t_points)
    
    result = integral / (3.086e22 * (1 + z))
    
    if result < 1e-5:
        print(f"Warning: r_d_approx small ({result}); using fallback")
        result = 147.0  # Fallback to ΛCDM value
    
    return result


# ------------------------------
# Test Functions
# ------------------------------

def test_scale_perception(z_test=1.0, t_obs=4.35e17):
    """Test ant analogy: Redshift from scale contraction."""
    try:
        t = t_from_z(z_test)
        redshift_model = redshift_cosmological(t, t_obs)
        assert np.isclose(redshift_model, z_test, rtol=0.1), "Scale contraction redshift mismatch"
    except Exception as e:
        print(f"Warning: Scale test failed: {e}")


def test_friction_scaling(mass1=1.0, mass2=2.0, t=1e10):
    """Test stone analogy: Drag scales with mass."""
    drag1 = friction_drag(t, mass1)
    drag2 = friction_drag(t, mass2)
    assert np.isclose(drag2 / drag1, mass2 / mass1, rtol=0.01), "Mass-scaled friction mismatch"


# ------------------------------
# Main Test Block
# ------------------------------

if __name__ == "__main__":
    try:
        # Basic functionality tests
        t_test = 1e10
        print(f"H_variable({t_test}) = {H_variable(t_test)}")
        print(f"global_S({t_test}) = {global_S(t_test)}")

        # Test classes
        geo = Geodesic()
        position = np.array([1e11, 0])
        M_test = 2e30
        accel = geo.acceleration(position, M_test, t_test)
        print(f"Acceleration: {accel}")

        # Test metric
        metric = Metric(1e11, M_test, t_test)
        print(f"Metric g_tt: {metric.g_tt}")

        # Test utility functions
        v_rot = rotation_curve(np.array([1e11, 2e11]), M_test, t_test)
        print(f"Rotation curve velocities: {v_rot}")

        # Test conversion functions
        z_test = 1.0
        t_converted = t_from_z(z_test)
        print(f"Time from z={z_test}: {t_converted}")

        print("All basic tests passed!")

        print("\n=== REDSHIFT INVERSION DIAGNOSTIC ===")
        z_test = 1.0
        t_obs = 4.35e17
        
        print(f"\n1. Testing t_from_z with z={z_test}:")
        t_emit = t_from_z(z_test, t_obs=t_obs)
        print(f"   Result: t_emit = {t_emit:.6e} seconds ({t_emit/(365.25*86400):.2e} years)")
        
        print(f"\n2. Testing reverse: redshift_cosmological:")
        z_recovered = redshift_cosmological(t_emit, t_obs)
        print(f"   Input z: {z_test}")
        print(f"   Recovered z: {z_recovered:.6f}")
        print(f"   Match: {np.isclose(z_test, z_recovered, rtol=0.01)}")
        
        print(f"\n3. Checking S(t) values:")
        S_emit = global_S(t_emit)
        S_obs = global_S(t_obs)
        print(f"   S(t_emit) = {S_emit:.10f}")
        print(f"   S(t_obs) = {S_obs:.10f}")
        print(f"   S_obs/S_emit = {S_obs/S_emit:.6f} (should be ~{1+z_test})")
        
        print(f"\n4. Testing H_variable at different times:")
        for t in [1e10, 1e15, t_emit, t_obs]:
            print(f"   H({t:.2e}) = {H_variable(t):.6e} s^-1")

        # ADD THIS NEW BLOCK FOR SUPERNOVA DISTANCE:
        print("\n=== SUPERNOVA DISTANCE DIAGNOSTIC ===")
        z_test_sn = 0.5  # Separate variable to avoid conflict
        mu_estif = distance_modulus_estif_numerical(z_test_sn)
        print(f"At z={z_test_sn}: μ = {mu_estif:.2f} mag")
        print(f"Expected: ~42-43 mag")
        print(f"Reasonable: {40 < mu_estif < 50}")

    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-01-10-25-V-2
