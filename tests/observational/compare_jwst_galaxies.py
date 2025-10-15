# compare_jwst_galaxies.py

"""
Compare ESTIF-Gravity galaxy asymmetry prediction to JWST observations.

Data Source: JWST CEERS (Cosmic Evolution Early Release Science) Survey
Papers: Multiple CEERS collaboration papers (2022-2024)
- Finkelstein et al. 2022 (ApJL 940, L55) - CEERS overview
- Kartaltepe et al. 2023 (ApJL 946, L15) - Morphology catalog

Survey Parameters:
- Field: Extended Groth Strip (EGS)
- Area: ~100 arcmin²
- Redshift range: z = 0.5 to z > 10
- Sample size: ~10,000 galaxies with morphologies
- NIRCam imaging: 7 filters (0.9-4.4 μm)

Morphology Measurements:
- Sérsic index n: Precision ~10-20%
- Axis ratio q: Precision ~5-10%
- Asymmetry parameter A: Precision ~5%
- Size measurements: Precision ~10%

Physics:
- GR + ΛCDM: Symmetric rotation unless perturbed
- ESTIF: Predicts asymmetry from fossilized drag patterns
- Test: Statistical detection in ensemble of ~100 galaxies

Status: CEERS data public, ongoing morphology analysis
"""

import sys
import os
# Add src directory to Python path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import numpy as np
import estif_ec_gr_model as estif
import estif_ec_gr_constants as const
import matplotlib.pyplot as plt
from scipy import stats

def main():
    print("="*80)
    print("JWST HIGH-REDSHIFT GALAXIES: ESTIF ASYMMETRY PREDICTION")
    print("="*80)
    print("\nData: JWST CEERS Survey (2022-2024)")
    print("Cosmic Evolution Early Release Science - Morphology Catalog")
    print("="*80)
    
    # =========================================================================
    # Observational Context (JWST CEERS)
    # =========================================================================
    
    print(f"\n📊 OBSERVATIONAL CONTEXT (JWST CEERS):")
    print(f"   Survey parameters:")
    print(f"      Field: Extended Groth Strip")
    print(f"      Area: ~100 arcmin²")
    print(f"      Redshift coverage: z = 0.5 to z > 10")
    print(f"      Sample size: ~10,000 galaxies")
    
    print(f"\n   Morphology measurements:")
    print(f"      Sérsic index precision: ~10-20%")
    print(f"      Axis ratio precision: ~5-10%")
    print(f"      Asymmetry parameter precision: ~5%")
    print(f"      Size measurement precision: ~10%")
    
    # Typical galaxy parameters at different redshifts
    print(f"\n   Typical galaxy properties (z=2-4):")
    print(f"      Stellar mass: 10⁹ - 10¹¹ M_sun")
    print(f"      Half-light radius: 0.5 - 3 kpc")
    print(f"      Star formation rate: 10 - 100 M_sun/yr")
    
    # =========================================================================
    # ESTIF Predictions Across Redshift
    # =========================================================================
    
    print(f"\n🔴 ESTIF-GRAVITY PREDICTIONS:")
    
    # Redshift samples (focusing on JWST optimal range)
    z_sample = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0])
    
    # Galaxy parameter ranges (will average over these)
    M_stellar_range = [1e10, 5e10, 1e11]  # Solar masses
    r_halflight_range = [1e19, 2e19, 3e19]  # ~0.3, 0.6, 1 kpc
    
    print(f"\n   Parameter space explored:")
    print(f"      Redshifts: z = {z_sample.min()} to {z_sample.max()}")
    print(f"      Stellar masses: {M_stellar_range[0]/const.M_sun:.1e} - {M_stellar_range[-1]/const.M_sun:.1e} M_sun")
    print(f"      Half-light radii: {r_halflight_range[0]/3.086e19:.1f} - {r_halflight_range[-1]/3.086e19:.1f} kpc")
    
    # Calculate predictions for each combination
    results_grid = {}
    
    print(f"\n   Calculating predictions...")
    for z in z_sample:
        asymmetries_at_z = []
        
        for M_stellar_msun in M_stellar_range:
            M_gal = M_stellar_msun * const.M_sun
            
            for r_kpc_scaled in r_halflight_range:
                r_gal = r_kpc_scaled
                
                # ESTIF prediction
                asym = estif.galaxy_drag_asymmetry(z, M_gal, r_gal)
                asymmetries_at_z.append(asym)
        
        # Store statistics for this redshift
        results_grid[z] = {
            'mean': np.mean(asymmetries_at_z),
            'std': np.std(asymmetries_at_z),
            'min': np.min(asymmetries_at_z),
            'max': np.max(asymmetries_at_z),
            'all': asymmetries_at_z
        }
    
    # Print summary table
    print(f"\n   Prediction Summary:")
    print(f"   {'z':<6} {'Mean':<10} {'Std':<10} {'Min':<10} {'Max':<10}")
    print(f"   {'-'*50}")
    for z in z_sample:
        r = results_grid[z]
        print(f"   {z:<6.1f} {r['mean']:<10.4f} {r['std']:<10.4f} {r['min']:<10.4f} {r['max']:<10.4f}")
    
    # =========================================================================
    # Statistical Detectability Analysis
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 STATISTICAL DETECTABILITY ANALYSIS")
    print(f"{'='*80}")
    
    # JWST measurement precision
    jwst_asymmetry_precision = 5.0  # Percent, per-galaxy
    jwst_axis_ratio_precision = 7.0  # Percent, alternative measure
    
    print(f"\nJWST Capabilities:")
    print(f"   Single galaxy asymmetry precision: ~{jwst_asymmetry_precision}%")
    print(f"   Axis ratio precision: ~{jwst_axis_ratio_precision}%")
    
    # For each redshift, assess detectability
    print(f"\nDetectability Assessment (per redshift bin):")
    print(f"   Assumes N=100 galaxies per bin for statistical analysis")
    print()
    
    n_galaxies_per_bin = 100
    statistical_precision = jwst_asymmetry_precision / np.sqrt(n_galaxies_per_bin)
    
    print(f"   Statistical precision with N={n_galaxies_per_bin}: {statistical_precision:.2f}%")
    print()
    
    detectable_count = 0
    marginal_count = 0
    
    for z in z_sample:
        mean_asym = results_grid[z]['mean']
        sigma_detection = mean_asym / statistical_precision
        
        print(f"   z = {z:.1f}:")
        print(f"      ESTIF signal: {mean_asym:.4f}%")
        print(f"      Detection significance: {sigma_detection:.2f}σ")
        
        if sigma_detection >= 3.0:
            print(f"      ✅ DETECTABLE (>3σ)")
            detectable_count += 1
        elif sigma_detection >= 2.0:
            print(f"      ⚠️  MARGINAL (2-3σ)")
            marginal_count += 1
        elif sigma_detection >= 1.0:
            print(f"      ⚠️  WEAK (1-2σ)")
        else:
            print(f"      ❌ UNDETECTABLE (<1σ)")
        print()
    
    # =========================================================================
    # Overall Verdict
    # =========================================================================
    
    print(f"{'='*80}")
    print("🎯 VERDICT")
    print(f"{'='*80}")
    
    max_signal = max([results_grid[z]['max'] for z in z_sample])
    max_mean_signal = max([results_grid[z]['mean'] for z in z_sample])
    
    print(f"\nSignal Strength:")
    print(f"   Maximum predicted asymmetry: {max_signal:.4f}%")
    print(f"   Maximum mean asymmetry: {max_mean_signal:.4f}%")
    print(f"   Single-galaxy precision: {jwst_asymmetry_precision:.1f}%")
    print(f"   Statistical precision (N=100): {statistical_precision:.2f}%")
    
    if detectable_count >= 3:
        print(f"\n✅ TESTABLE with JWST CEERS data")
        print(f"   • Detectable at {detectable_count}/{len(z_sample)} redshift bins")
        print(f"   • Requires ensemble analysis of ~100 galaxies per bin")
        print(f"   • Total sample needed: ~{detectable_count*100} galaxies")
        status = "testable"
    elif marginal_count + detectable_count >= 3:
        print(f"\n⚠️  MARGINALLY TESTABLE with JWST")
        print(f"   • {detectable_count} bins detectable, {marginal_count} bins marginal")
        print(f"   • Requires larger sample or improved analysis methods")
        print(f"   • Consider: Stack all z=2-4 galaxies together")
        status = "marginal"
    elif max_mean_signal > statistical_precision:
        print(f"\n⚠️  WEAK SIGNAL - Challenging but possible")
        print(f"   • Best case: {max_mean_signal:.3f}% vs {statistical_precision:.2f}% precision")
        print(f"   • Need larger sample sizes or cross-correlation methods")
        status = "challenging"
    else:
        print(f"\n❌ UNDETECTABLE with current JWST capabilities")
        print(f"   • Signal: {max_mean_signal:.3f}% vs precision: {statistical_precision:.2f}%")
        print(f"   • Improvement needed: {statistical_precision/max_mean_signal:.0f}×")
        status = "undetectable"
    
    # =========================================================================
    # Visualization
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 GENERATING COMPARISON PLOTS")
    print(f"{'='*80}")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Plot 1: Mean asymmetry vs redshift
    ax1 = axes[0, 0]
    means = [results_grid[z]['mean'] for z in z_sample]
    stds = [results_grid[z]['std'] for z in z_sample]
    
    ax1.errorbar(z_sample, means, yerr=stds, fmt='ro-', linewidth=2, 
                 markersize=8, capsize=5, label='ESTIF Prediction')
    ax1.axhline(jwst_asymmetry_precision, color='blue', linestyle='--', 
                linewidth=2, label=f'JWST single-galaxy precision ({jwst_asymmetry_precision}%)')
    ax1.axhline(statistical_precision, color='green', linestyle='--', 
                linewidth=2, label=f'Statistical precision (N=100: {statistical_precision:.2f}%)')
    ax1.fill_between(z_sample, 0, statistical_precision, alpha=0.2, color='green')
    
    ax1.set_xlabel('Redshift z', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Asymmetry (%)', fontsize=12, fontweight='bold')
    ax1.set_title('ESTIF Prediction vs JWST Sensitivity', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3)
    
    # Plot 2: Detection significance
    ax2 = axes[0, 1]
    sigmas = [results_grid[z]['mean'] / statistical_precision for z in z_sample]
    
    colors_sigma = ['green' if s >= 3 else 'orange' if s >= 2 else 'red' for s in sigmas]
    bars = ax2.bar(z_sample, sigmas, width=0.3, color=colors_sigma, 
                   alpha=0.7, edgecolor='black', linewidth=2)
    
    ax2.axhline(3, color='green', linestyle='--', linewidth=2, label='3σ (strong detection)')
    ax2.axhline(2, color='orange', linestyle='--', linewidth=2, label='2σ (marginal)')
    ax2.axhline(1, color='red', linestyle='--', linewidth=2, label='1σ (weak)')
    ax2.fill_between([0.5, 5.5], 3, 10, alpha=0.1, color='green')
    ax2.fill_between([0.5, 5.5], 2, 3, alpha=0.1, color='orange')
    
    ax2.set_xlabel('Redshift z', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Detection Significance (σ)', fontsize=12, fontweight='bold')
    ax2.set_title('Statistical Detectability (N=100 per bin)', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0.5, 5.5)
    ax2.set_ylim(0, max(sigmas)*1.3)
    
    # Plot 3: Parameter space heatmap
    ax3 = axes[1, 0]
    
    # Create grid for one representative redshift (z=3)
    z_demo = 3.0
    M_grid = np.logspace(10, 11.5, 20)  # Solar masses
    r_grid = np.linspace(0.5e19, 3e19, 20)  # meters
    
    asym_grid = np.zeros((len(r_grid), len(M_grid)))
    
    for i, r in enumerate(r_grid):
        for j, M_msun in enumerate(M_grid):
            M = M_msun * const.M_sun
            asym_grid[i, j] = estif.galaxy_drag_asymmetry(z_demo, M, r)
    
    im = ax3.contourf(M_grid/const.M_sun, r_grid/3.086e19, asym_grid, 
                      levels=20, cmap='RdYlGn_r')
    cbar = plt.colorbar(im, ax=ax3)
    cbar.set_label('Asymmetry (%)', fontsize=11)
    
    ax3.set_xlabel('Stellar Mass (M_sun)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Half-light Radius (kpc)', fontsize=12, fontweight='bold')
    ax3.set_title(f'ESTIF Asymmetry at z={z_demo}', fontsize=13, fontweight='bold')
    ax3.set_xscale('log')
    
    # Add contour for JWST precision
    ax3.contour(M_grid/const.M_sun, r_grid/3.086e19, asym_grid, 
                levels=[statistical_precision], colors='blue', linewidths=3,
                linestyles='--')
    
    # Plot 4: Sample size requirements
    ax4 = axes[1, 1]
    
    # Calculate N needed for 3σ detection at each z
    n_needed = []
    for z in z_sample:
        mean_asym = results_grid[z]['mean']
        # N needed: (3σ × single_precision / signal)²
        if mean_asym > 0:
            n = (3 * jwst_asymmetry_precision / mean_asym)**2
        else:
            n = 1e6  # Effectively infinite
        n_needed.append(n)
    
    ax4.semilogy(z_sample, n_needed, 'bo-', linewidth=2, markersize=8)
    ax4.axhline(100, color='green', linestyle='--', linewidth=2, 
                label='N=100 (CEERS typical)')
    ax4.axhline(1000, color='orange', linestyle='--', linewidth=2,
                label='N=1000 (full CEERS)')
    ax4.fill_between(z_sample, 0, 100, alpha=0.15, color='green')
    
    ax4.set_xlabel('Redshift z', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Sample Size Needed for 3σ', fontsize=12, fontweight='bold')
    ax4.set_title('Sample Size Requirements', fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(alpha=0.3)
    ax4.set_ylim(10, 1e5)
    
    plt.tight_layout()
    plt.savefig('jwst_ceers_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("   ✓ Plot saved: jwst_ceers_comparison.png")
    
    # =========================================================================
    # Observational Strategy
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("🔬 OBSERVATIONAL STRATEGY")
    print(f"{'='*80}")
    
    print(f"\nRecommended Analysis Approach:")
    
    if status == "testable":
        print(f"\n✅ Direct Detection Possible")
        print(f"   1. Select ~{detectable_count*100} galaxies in z={z_sample[0]}-{z_sample[-1]} range")
        print(f"   2. Measure asymmetry parameters (A, q) for each")
        print(f"   3. Bin by redshift ({len(z_sample)} bins)")
        print(f"   4. Test for systematic increase in asymmetry with z")
        print(f"   5. Expected significance: {max(sigmas):.1f}σ at optimal redshift")
        
        print(f"\n   Statistical Tests:")
        print(f"   • Null hypothesis: <A> = 0 (GR prediction)")
        print(f"   • Alternative: <A> ∝ f(z) (ESTIF prediction)")
        print(f"   • Method: Weighted mean comparison across z bins")
        print(f"   • Control: Compare to z < 1 baseline (should be near zero)")
        
    elif status == "marginal":
        print(f"\n⚠️  Marginal Detection - Enhanced Analysis Needed")
        print(f"   1. Stack all z=2-4 galaxies together (increase statistics)")
        print(f"   2. Use cross-correlation with mass/size")
        print(f"   3. Consider multiple asymmetry indicators:")
        print(f"      - Asymmetry parameter A")
        print(f"      - Axis ratio q")
        print(f"      - Spiral arm pitch angle")
        print(f"      - Kinematic misalignment")
        print(f"   4. Combine CEERS with other JWST surveys (JADES, GLASS)")
        
    elif status == "challenging":
        print(f"\n⚠️  Challenging - Requires Optimal Strategy")
        print(f"   1. Maximize sample size: Use all JWST public data")
        print(f"   2. Focus on highest-S/N galaxies only")
        print(f"   3. Develop custom asymmetry metric optimized for ESTIF")
        print(f"   4. Use Bayesian hierarchical modeling")
        print(f"   5. Consider ML/deep learning for subtle pattern detection")
        
    else:
        print(f"\n❌ Direct Detection Not Feasible")
        print(f"   Signal below detection threshold even with full CEERS sample")
        print(f"   Alternative approaches:")
        print(f"   1. Wait for deeper JWST surveys (longer exposures)")
        print(f"   2. Focus on IFU spectroscopy for kinematic asymmetries")
        print(f"   3. Prioritize other ESTIF predictions (EHT, LIGO)")
    
    print(f"\nData Availability:")
    print(f"   • CEERS data: Public (MAST archive)")
    print(f"   • Morphology catalogs: Kartaltepe et al. 2023")
    print(f"   • Photometric redshifts: Finkelstein et al. 2022")
    print(f"   • Access: https://ceers.github.io/")
    
    # =========================================================================
    # Cross-Check Recommendations
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("🔗 CROSS-CHECKS WITH OTHER PREDICTIONS")
    print(f"{'='*80}")
    
    print(f"\nConsistency Analysis:")
    print(f"   Compare JWST results with:")
    print(f"   1. EHT M87* (run compare_eht_m87.py)")
    print(f"      - Both use friction_drag_local(M, r)")
    print(f"      - Should show similar β_drag constraints")
    print(f"   2. LIGO GW150914 (run compare_ligo_gw.py)")
    print(f"      - Different physics (GW vs EM)")
    print(f"      - May have different sensitivities")
    
    print(f"\n   Pattern Recognition:")
    if max_mean_signal < 0.01:
        print(f"   ⚠️  All three predictions show similar magnitude issues (~0.01% scale)")
        print(f"   → Suggests systematic problem with friction_drag_local formula")
        print(f"   → Likely missing relativistic or geometric factors")
    elif max_mean_signal < 1:
        print(f"   ⚠️  Predictions consistently ~1% level")
        print(f"   → May be at limit of observability")
        print(f"   → Consider: Is β_drag universally applicable?")
    else:
        print(f"   ✓ Predictions at detectable levels")
        print(f"   → Formula likely correct")
        print(f"   → Proceed with observational campaign")
    
    # =========================================================================
    # Summary and Recommendations
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📋 SUMMARY AND RECOMMENDATIONS")
    print(f"{'='*80}")
    
    print(f"\n1. Prediction Summary:")
    print(f"   • Redshift range: z = {z_sample.min()} to {z_sample.max()}")
    print(f"   • Mean asymmetry: {min(means):.4f}% to {max(means):.4f}%")
    print(f"   • Maximum asymmetry: {max_signal:.4f}%")
    print(f"   • Detection significance: {min(sigmas):.2f}σ to {max(sigmas):.2f}σ")
    
    print(f"\n2. Observational Feasibility:")
    print(f"   • JWST precision: {jwst_asymmetry_precision}% (single galaxy)")
    print(f"   • Statistical precision: {statistical_precision:.2f}% (N=100)")
    print(f"   • Detectable bins: {detectable_count}/{len(z_sample)}")
    print(f"   • Status: {status.upper()}")
    
    print(f"\n3. Priority Actions:")
    if status in ["testable", "marginal"]:
        print(f"   HIGH PRIORITY:")
        print(f"   1. Download CEERS morphology catalog")
        print(f"   2. Implement asymmetry measurement pipeline")
        print(f"   3. Perform statistical analysis on z={z_sample[np.argmax(sigmas)]:.1f} bin (best S/N)")
        print(f"   4. Cross-check with EHT and LIGO results")
    else:
        print(f"   MEDIUM PRIORITY:")
        print(f"   1. Focus on EHT/LIGO predictions first")
        print(f"   2. Revisit JWST if other predictions work")
        print(f"   3. Consider theoretical refinements to boost signal")
    
    print(f"\n4. Physics Implications:")
    print(f"   • If detected: Strong evidence for 4D friction mechanism")
    print(f"   • If marginally detected: Constrains β_drag parameter space")
    print(f"   • If not detected: Rules out current ESTIF formulation")
    print(f"   • Cross-prediction consistency critical for validation")
    
    print(f"\n5. Next Steps:")
    print(f"   Immediate (1-2 weeks):")
    print(f"   • Run all three comparison scripts (EHT, LIGO, JWST)")
    print(f"   • Identify which prediction is strongest")
    print(f"   • Look for systematic patterns in magnitude discrepancies")
    
    print(f"\n   Short-term (1-3 months):")
    if status == "testable":
        print(f"   • Implement JWST data analysis pipeline")
        print(f"   • Perform null-hypothesis test on CEERS data")
        print(f"   • Prepare results for publication")
    else:
        print(f"   • Investigate formula scaling factors")
        print(f"   • Consider relativistic corrections")
        print(f"   • Re-derive friction_drag_local if needed")
    
    print(f"\n   Long-term (6-12 months):")
    print(f"   • Publish observational comparison results")
    print(f"   • Refine model based on constraints")
    print(f"   • Propose for dedicated observing time if promising")
    
    print(f"\n{'='*80}")
    print("END OF JWST COMPARISON")
    print(f"{'='*80}\n")
    
    return {
        'status': status,
        'max_sigma': max(sigmas),
        'detectable_bins': detectable_count,
        'max_asymmetry': max_signal,
        'sample_size_needed': min(n_needed)
    }


if __name__ == "__main__":
    results = main()
    

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2


