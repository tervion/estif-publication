# compare_eht_m87.py

"""
Compare ESTIF-Gravity lensing prediction to EHT M87* observations.

Data Source: Event Horizon Telescope Collaboration (2019)
Paper: "First M87 Event Horizon Telescope Results. I. The Shadow of the 
        Supermassive Black Hole"
ApJL 875, L1 (2019)

Measurements:
- M87* black hole mass: 6.5 ± 0.7 × 10⁹ M_sun
- Shadow diameter: 42 ± 3 μas (microarcseconds)
- Distance: 16.8 ± 0.8 Mpc

Physics:
- GR predicts shadow at r_shadow ≈ √27/2 × R_s (photon sphere)
- ESTIF modifies light bending via friction_drag_local(M, r)
- Test: Does friction produce observable shadow size change?

Status: EHT precision ~7% (2019), improving to ~1% by 2030
"""

import sys
import os
# Add src directory to Python path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import numpy as np
import estif_ec_gr_model as estif
import estif_ec_gr_constants as const
import matplotlib.pyplot as plt

def main():
    print("="*80)
    print("EHT M87* BLACK HOLE SHADOW: ESTIF-GRAVITY vs OBSERVATION")
    print("="*80)
    print("\nData: Event Horizon Telescope Collaboration (2019)")
    print("ApJL 875, L1 - First M87 Event Horizon Telescope Results")
    print("="*80)
    
    # =========================================================================
    # Observational Data (EHT 2019)
    # =========================================================================
    
    # M87* parameters with uncertainties
    M_m87_central = 6.5e9 * const.M_sun
    M_m87_error = 0.7e9 * const.M_sun
    
    # Shadow angular diameter (microarcseconds)
    shadow_diameter_uas_observed = 42.0
    shadow_diameter_error_uas = 3.0
    
    # Distance to M87
    distance_mpc = 16.8
    distance_error_mpc = 0.8
    
    print(f"\n📊 OBSERVATIONAL DATA:")
    print(f"   Black hole mass: {M_m87_central/const.M_sun:.2e} ± {M_m87_error/const.M_sun:.2e} M_sun")
    print(f"   Shadow diameter: {shadow_diameter_uas_observed:.1f} ± {shadow_diameter_error_uas:.1f} μas")
    print(f"   Distance: {distance_mpc:.1f} ± {distance_error_mpc:.1f} Mpc")
    
    # Convert to physical units
    distance_m = distance_mpc * 3.086e22  # Mpc to meters
    R_s = estif.schwarzschild_radius(M_m87_central)
    
    # Shadow angular size to physical radius
    shadow_angle_rad = shadow_diameter_uas_observed * 1e-6 / 206265
    R_shadow_observed = shadow_angle_rad * distance_m
    
    print(f"   Schwarzschild radius: {R_s:.3e} m")
    print(f"   Observed shadow: {R_shadow_observed/R_s:.2f} R_s")
    
    # =========================================================================
    # General Relativity Prediction
    # =========================================================================
    
    print(f"\n🔵 GENERAL RELATIVITY PREDICTION:")
    
    # For Schwarzschild black hole, photon sphere at r = 3R_s/2
    # Shadow radius = sqrt(27) R_s / 2 ≈ 2.598 R_s
    r_photon_sphere = 1.5 * R_s
    R_shadow_gr = np.sqrt(27) * Rs
    
    # Angular size from Earth
    theta_shadow_gr_rad = R_shadow_gr / distance_m
    theta_shadow_gr_uas = theta_shadow_gr_rad * 206265 * 1e6  # Convert to μas
    
    print(f"   Photon sphere: {r_photon_sphere/R_s:.1f} R_s")
    print(f"   Shadow radius: {R_shadow_gr/R_s:.3f} R_s")
    print(f"   Shadow diameter: {theta_shadow_gr_uas:.2f} μas")
    
    # Compare to observation
    gr_deviation_sigma = (theta_shadow_gr_uas - shadow_diameter_uas_observed) / shadow_diameter_error_uas
    print(f"   Deviation from EHT: {gr_deviation_sigma:.2f}σ")
    
    if abs(gr_deviation_sigma) < 1:
        print(f"   ✅ GR is consistent with observation (within 1σ)")
    else:
        print(f"   ⚠️  GR shows {abs(gr_deviation_sigma):.1f}σ tension")
    
    # =========================================================================
    # ESTIF-Gravity Prediction
    # =========================================================================
    
    print(f"\n🔴 ESTIF-GRAVITY PREDICTION:")
    
    # Calculate lensing at photon sphere
    theta_deflection_gr = (4 * const.G * M_m87_central) / (const.c**2 * r_photon_sphere)
    theta_deflection_estif = estif.unique_lensing_signature(r_photon_sphere, M_m87_central)
    
    # Fractional change in deflection
    deflection_change = (theta_deflection_estif / theta_deflection_gr - 1)
    
    print(f"   Using friction_drag_local(M={M_m87_central/const.M_sun:.2e} M_sun, r={r_photon_sphere/R_s:.1f} R_s)")
    
    # Local drag calculation
    drag_local = estif.friction_drag_local(M_m87_central, r_photon_sphere)
    print(f"   Local drag coefficient: {drag_local:.3e}")
    print(f"   GR deflection: {theta_deflection_gr:.3e} rad")
    print(f"   ESTIF deflection: {theta_deflection_estif:.3e} rad")
    print(f"   Fractional change: {deflection_change*100:.6f}%")
    
    # Shadow size scales approximately with deflection angle
    # (First-order approximation: more accurate would require full ray tracing)
    theta_shadow_estif_uas = theta_shadow_gr_uas * (1 + deflection_change)
    
    print(f"   Predicted shadow diameter: {theta_shadow_estif_uas:.2f} μas")
    print(f"   Change from GR: {(theta_shadow_estif_uas - theta_shadow_gr_uas):.4f} μas")
    
    # Compare to observation
    estif_deviation_sigma = (theta_shadow_estif_uas - shadow_diameter_uas_observed) / shadow_diameter_error_uas
    print(f"   Deviation from EHT: {estif_deviation_sigma:.2f}σ")
    
    # =========================================================================
    # Statistical Assessment
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 STATISTICAL ASSESSMENT")
    print(f"{'='*80}")
    
    print(f"\nModel Comparison:")
    print(f"   GR:    {theta_shadow_gr_uas:.2f} μas ({gr_deviation_sigma:+.2f}σ from observation)")
    print(f"   ESTIF: {theta_shadow_estif_uas:.2f} μas ({estif_deviation_sigma:+.2f}σ from observation)")
    print(f"   Obs:   {shadow_diameter_uas_observed:.1f} ± {shadow_diameter_error_uas:.1f} μas")
    
    # Determine verdict
    print(f"\n{'='*80}")
    print("🎯 VERDICT")
    print(f"{'='*80}")
    
    if abs(estif_deviation_sigma) < 1:
        print("\n✅ ESTIF-GRAVITY is CONSISTENT with EHT observation (within 1σ)")
        status = "consistent"
    elif abs(estif_deviation_sigma) < 2:
        print("\n⚠️  ESTIF-GRAVITY shows MILD TENSION with EHT (1-2σ)")
        status = "marginal"
    elif abs(estif_deviation_sigma) < 3:
        print("\n⚠️  ESTIF-GRAVITY shows MODERATE TENSION with EHT (2-3σ)")
        status = "tension"
    else:
        print("\n❌ ESTIF-GRAVITY is RULED OUT by EHT observation (>3σ discrepancy)")
        status = "ruled_out"
    
    print(f"\nQuantitative deviation: {abs(estif_deviation_sigma):.3f}σ")
    print(f"Fractional difference from GR: {abs(deflection_change)*100:.6f}%")
    
    # Future prospects
    print(f"\n🔭 OBSERVATIONAL PROSPECTS:")
    print(f"   Current EHT precision: ~7% ({shadow_diameter_error_uas:.1f} μas)")
    print(f"   Required for 3σ detection: {3*shadow_diameter_error_uas:.1f} μas difference")
    print(f"   ESTIF prediction: {abs(theta_shadow_estif_uas - theta_shadow_gr_uas):.4f} μas difference")
    
    required_precision_uas = abs(theta_shadow_estif_uas - theta_shadow_gr_uas) / 3
    required_precision_percent = required_precision_uas / shadow_diameter_uas_observed * 100
    
    print(f"   Precision needed to test ESTIF: {required_precision_percent:.3f}% ({required_precision_uas:.4f} μas)")
    print(f"   EHT 2025-2030 goal: ~1% precision")
    
    if required_precision_percent < 1:
        print(f"   ✅ TESTABLE with next-generation EHT")
    elif required_precision_percent < 5:
        print(f"   ⚠️  MARGINALLY TESTABLE with improved EHT")
    else:
        print(f"   ❌ BELOW foreseeable EHT precision")
    
    # =========================================================================
    # Visualization
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 GENERATING COMPARISON PLOT")
    print(f"{'='*80}")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Shadow size comparison
    models = ['EHT 2019\n(Observed)', 'General\nRelativity', 'ESTIF-Gravity']
    diameters = [shadow_diameter_uas_observed, theta_shadow_gr_uas, theta_shadow_estif_uas]
    errors = [shadow_diameter_error_uas, 0, 0]
    colors = ['blue', 'green', 'red']
    
    bars = ax1.bar(models, diameters, yerr=[[errors[0], 0, 0], [errors[0], 0, 0]], 
                   capsize=10, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
    
    # Add observation band
    ax1.axhline(shadow_diameter_uas_observed, color='blue', linestyle='--', 
                alpha=0.3, linewidth=2)
    ax1.fill_between([-0.5, 2.5], 
                     shadow_diameter_uas_observed - shadow_diameter_error_uas,
                     shadow_diameter_uas_observed + shadow_diameter_error_uas,
                     alpha=0.2, color='blue', label='1σ EHT uncertainty')
    
    ax1.set_ylabel('Shadow Diameter (μas)', fontsize=13, fontweight='bold')
    ax1.set_title('M87* Black Hole Shadow Size', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper left')
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_ylim(38, 46)
    
    # Add value labels on bars
    for i, (bar, diam) in enumerate(zip(bars, diameters)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{diam:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Plot 2: Deviation in sigma
    sigma_values = [0, gr_deviation_sigma, estif_deviation_sigma]
    colors_sigma = ['blue', 'green', 'red']
    
    bars2 = ax2.bar(models, sigma_values, color=colors_sigma, alpha=0.7, 
                    edgecolor='black', linewidth=2)
    
    # Add significance bands
    ax2.axhspan(-1, 1, alpha=0.15, color='green', label='Consistent (<1σ)')
    ax2.axhspan(1, 2, alpha=0.1, color='yellow', label='Mild tension (1-2σ)')
    ax2.axhspan(2, 3, alpha=0.1, color='orange', label='Moderate tension (2-3σ)')
    ax2.axhspan(-2, -1, alpha=0.1, color='yellow')
    ax2.axhspan(-3, -2, alpha=0.1, color='orange')
    
    ax2.axhline(0, color='black', linestyle='-', linewidth=1)
    ax2.set_ylabel('Deviation (σ)', fontsize=13, fontweight='bold')
    ax2.set_title('Statistical Significance', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='upper right')
    ax2.grid(axis='y', alpha=0.3)
    
    # Add value labels
    for bar, sigma in zip(bars2, sigma_values):
        height = bar.get_height()
        if height != 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1*np.sign(height),
                    f'{sigma:.2f}σ', ha='center', va='bottom' if height > 0 else 'top',
                    fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('eht_m87_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("   ✓ Plot saved: eht_m87_comparison.png")
    
    # =========================================================================
    # Summary and Recommendations
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📋 SUMMARY AND RECOMMENDATIONS")
    print(f"{'='*80}")
    
    print(f"\n1. Current Status:")
    print(f"   • ESTIF predicts {abs(deflection_change)*100:.6f}% deviation from GR")
    print(f"   • Corresponds to {abs(theta_shadow_estif_uas - shadow_diameter_uas_observed):.4f} μas shift")
    print(f"   • Statistical significance: {abs(estif_deviation_sigma):.3f}σ")
    
    print(f"\n2. Observational Context:")
    print(f"   • EHT 2019 precision: ~7%")
    print(f"   • ESTIF signal: ~{abs(deflection_change)*100:.6f}%")
    print(f"   • Signal-to-noise: {abs(deflection_change)*100 / 7:.4f}")
    
    print(f"\n3. Future Prospects:")
    if required_precision_percent < 1:
        print(f"   ✅ Definitive test possible with EHT improvements by 2030")
        print(f"   • Need: {required_precision_percent:.3f}% precision")
        print(f"   • Goal: ~1% precision by 2030")
    else:
        print(f"   ⚠️  Signal below foreseeable EHT precision")
        print(f"   • Need: {required_precision_percent:.3f}% precision")
        print(f"   • EHT limit: ~1% precision")
    
    print(f"\n4. Recommendations:")
    if status == "consistent":
        print(f"   • ESTIF remains viable - continue to Phase 2 (LIGO)")
        print(f"   • Monitor EHT precision improvements")
        print(f"   • Consider parameter space exploration (β_drag)")
    elif status == "marginal":
        print(f"   • ESTIF marginally consistent - proceed with caution")
        print(f"   • Priority: Examine other predictions (GW, galaxies)")
        print(f"   • Consider parameter adjustments if needed")
    elif status == "tension":
        print(f"   • ESTIF shows tension - examine assumptions")
        print(f"   • Check: Is photon sphere approximation valid?")
        print(f"   • Consider: Are there systematic errors in EHT data?")
    else:
        print(f"   • ESTIF ruled out at >3σ level")
        print(f"   • Recommend: Re-examine fundamental assumptions")
        print(f"   • Consider: Model may need substantial revision")
    
    print(f"\n{'='*80}")
    print("END OF EHT M87* COMPARISON")
    print(f"{'='*80}\n")
    
    return {
        'status': status,
        'sigma_deviation': estif_deviation_sigma,
        'fractional_change': deflection_change,
        'testable': required_precision_percent < 1
    }


if __name__ == "__main__":
    results = main()


#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2


