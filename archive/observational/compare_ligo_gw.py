# compare_ligo_gw.py

"""
Compare ESTIF-Gravity GW damping prediction to LIGO/Virgo observations.

Data Source: LIGO Scientific Collaboration (2016)
Paper: "Observation of Gravitational Waves from a Binary Black Hole Merger"
PRL 116, 061102 (2016) - GW150914

Event Parameters:
- Detection date: September 14, 2015
- Component masses: 36 (+5/-4) M_sun and 29 (+4/-4) M_sun
- Final mass: 62 (+4/-4) M_sun
- Radiated energy: 3.0 (+0.5/-0.5) M_sun c²
- Peak luminosity: 3.6×10⁵⁶ erg/s
- Distance: 410 (+160/-180) Mpc

Timing Precision:
- LIGO H1/L1 arrival time difference: 6.9 (+0.5/-0.4) ms
- Waveform match precision: ~1 ms
- Phase coherence: ~0.1 ms

Physics:
- GR predicts specific chirp waveform from inspiral-merger-ringdown
- ESTIF predicts damping from friction_drag_local at merger
- Test: Does friction produce measurable timing delay?

Status: LIGO/Virgo ongoing, ~90 events detected (O1-O3)
        LISA launch ~2037, precision improvement ~100×
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
    print("LIGO GW150914: ESTIF-GRAVITY DAMPING PREDICTION")
    print("="*80)
    print("\nData: LIGO Scientific Collaboration (2016)")
    print("PRL 116, 061102 - Observation of Gravitational Waves")
    print("="*80)
    
    # =========================================================================
    # Observational Data (LIGO GW150914)
    # =========================================================================
    
    # Binary component masses (solar masses)
    m1_msun = 36.0
    m1_error_plus = 5.0
    m1_error_minus = 4.0
    
    m2_msun = 29.0
    m2_error_plus = 4.0
    m2_error_minus = 4.0
    
    # Convert to SI
    m1 = m1_msun * const.M_sun
    m2 = m2_msun * const.M_sun
    
    # Total and chirp mass
    M_total = m1 + m2
    M_chirp = ((m1 * m2)**(3/5)) / ((m1 + m2)**(1/5))
    
    # Final mass and radiated energy
    M_final_msun = 62.0
    M_final_error = 4.0
    M_final = M_final_msun * const.M_sun
    
    M_radiated_msun = (m1_msun + m2_msun) - M_final_msun
    M_radiated = M_radiated_msun * const.M_sun
    
    # Distance
    distance_mpc = 410.0
    distance_error_plus = 160.0
    distance_error_minus = 180.0
    
    # Timing precision
    ligo_timing_precision_ms = 1.0  # Waveform match precision
    ligo_phase_precision_ms = 0.1   # Phase coherence precision
    
    print(f"\n📊 OBSERVATIONAL DATA (GW150914):")
    print(f"   Component masses:")
    print(f"      m1 = {m1_msun:.0f} (+{m1_error_plus:.0f}/-{m1_error_minus:.0f}) M_sun")
    print(f"      m2 = {m2_msun:.0f} (+{m2_error_plus:.0f}/-{m2_error_minus:.0f}) M_sun")
    print(f"   System properties:")
    print(f"      Total mass: {M_total/const.M_sun:.0f} M_sun")
    print(f"      Chirp mass: {M_chirp/const.M_sun:.1f} M_sun")
    print(f"      Final mass: {M_final_msun:.0f} ± {M_final_error:.0f} M_sun")
    print(f"      Radiated energy: {M_radiated_msun:.1f} M_sun c²")
    print(f"   Distance: {distance_mpc:.0f} (+{distance_error_plus:.0f}/-{distance_error_minus:.0f}) Mpc")
    print(f"   LIGO timing precision:")
    print(f"      Waveform match: ~{ligo_timing_precision_ms:.1f} ms")
    print(f"      Phase coherence: ~{ligo_phase_precision_ms:.1f} ms")
    
    # =========================================================================
    # General Relativity Expectation
    # =========================================================================
    
    print(f"\n🔵 GENERAL RELATIVITY PREDICTION:")
    
    # GR predicts specific waveform with no additional delays
    # Ringdown timescale (damping time of final black hole)
    R_s_final = estif.schwarzschild_radius(M_final)
    tau_ringdown_gr = R_s_final / const.c  # Light-crossing time
    
    # Frequency of fundamental quasi-normal mode
    f_qnm = const.c**3 / (2 * np.pi * const.G * M_final) / 10  # Approximate for a=0
    T_qnm = 1 / f_qnm
    
    print(f"   Final BH Schwarzschild radius: {R_s_final:.3e} m ({R_s_final/1000:.1f} km)")
    print(f"   Ringdown timescale: {tau_ringdown_gr*1000:.3f} ms")
    print(f"   QNM frequency: {f_qnm:.1f} Hz")
    print(f"   QNM period: {T_qnm*1000:.3f} ms")
    print(f"   Expected GW signature: Inspiral-merger-ringdown with no delays")
    
    # =========================================================================
    # ESTIF-Gravity Prediction
    # =========================================================================
    
    print(f"\n🔴 ESTIF-GRAVITY PREDICTION:")
    
    # Calculate damping delay using ESTIF formula
    delay_estif = estif.gw_damping_delay(0, M_total)
    
    print(f"   Using friction_drag_local at merger:")
    print(f"      Total mass: {M_total/const.M_sun:.0f} M_sun")
    print(f"      Schwarzschild radius: {estif.schwarzschild_radius(M_total):.3e} m")
    
    # Show friction calculation details
    R_s_total = estif.schwarzschild_radius(M_total)
    drag_local = estif.friction_drag_local(M_total, R_s_total)
    
    print(f"   Friction drag coefficient: {drag_local:.3e}")
    print(f"   Predicted damping delay: {delay_estif:.3e} s")
    print(f"   Delay in milliseconds: {delay_estif*1000:.3e} ms")
    
    # Compare to observational precision
    print(f"\n   Comparison to detector sensitivity:")
    
    # Signal-to-noise ratios
    snr_waveform = delay_estif / (ligo_timing_precision_ms * 1e-3)
    snr_phase = delay_estif / (ligo_phase_precision_ms * 1e-3)
    
    print(f"      Delay / Waveform precision = {snr_waveform:.3e}")
    print(f"      Delay / Phase precision = {snr_phase:.3e}")
    
    # =========================================================================
    # LISA Comparison (Future)
    # =========================================================================
    
    print(f"\n🛰️  LISA PROSPECTS (Launch ~2037):")
    
    lisa_timing_precision_ms = 0.01  # ~10 μs precision expected
    snr_lisa = delay_estif / (lisa_timing_precision_ms * 1e-3)
    
    print(f"   Expected timing precision: ~{lisa_timing_precision_ms:.3f} ms")
    print(f"   ESTIF delay / LISA precision = {snr_lisa:.3e}")
    
    # =========================================================================
    # Statistical Assessment
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 STATISTICAL ASSESSMENT")
    print(f"{'='*80}")
    
    print(f"\nSignal Strength:")
    print(f"   ESTIF prediction: {delay_estif:.3e} s ({delay_estif*1000:.3e} ms)")
    print(f"   LIGO precision:   {ligo_timing_precision_ms:.3e} ms")
    print(f"   LISA precision:   {lisa_timing_precision_ms:.3e} ms (future)")
    
    print(f"\nSignal-to-Noise Ratios:")
    print(f"   LIGO waveform: {snr_waveform:.3e} (need >1 for detection)")
    print(f"   LIGO phase:    {snr_phase:.3e}")
    print(f"   LISA:          {snr_lisa:.3e} (future)")
    
    # Determine verdict
    print(f"\n{'='*80}")
    print("🎯 VERDICT")
    print(f"{'='*80}")
    
    if snr_waveform >= 1.0:
        print("\n✅ POTENTIALLY DETECTABLE by LIGO/Virgo")
        print(f"   Signal exceeds waveform precision by {snr_waveform:.1f}×")
        status = "detectable_ligo"
    elif snr_lisa >= 1.0:
        print("\n⚠️  NOT DETECTABLE by LIGO, but LISA-ACCESSIBLE")
        print(f"   LIGO S/N: {snr_waveform:.3e} (below threshold)")
        print(f"   LISA S/N: {snr_lisa:.3e} (above threshold)")
        status = "detectable_lisa"
    elif snr_lisa >= 0.1:
        print("\n⚠️  MARGINALLY DETECTABLE with future LISA")
        print(f"   LISA S/N: {snr_lisa:.3e}")
        print(f"   May require statistical analysis of multiple events")
        status = "marginal_lisa"
    else:
        print("\n❌ BELOW DETECTION THRESHOLD for foreseeable detectors")
        print(f"   LIGO S/N: {snr_waveform:.3e}")
        print(f"   LISA S/N: {snr_lisa:.3e}")
        print(f"   Improvement needed: {1/snr_lisa:.0e}×")
        status = "undetectable"
    
    # =========================================================================
    # Parameter Space Exploration
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("⚙️  PARAMETER SPACE ANALYSIS")
    print(f"{'='*80}")
    
    # What β_drag would be needed for LIGO detection?
    required_delay_ligo = ligo_timing_precision_ms * 1e-3  # 1 ms in seconds
    current_delay = delay_estif
    beta_needed_ligo = const.BETA_DRAG * (required_delay_ligo / current_delay)
    
    # What about LISA?
    required_delay_lisa = lisa_timing_precision_ms * 1e-3
    beta_needed_lisa = const.BETA_DRAG * (required_delay_lisa / current_delay)
    
    print(f"\nCurrent model:")
    print(f"   β_drag = {const.BETA_DRAG}")
    print(f"   Predicted delay = {delay_estif:.3e} s")
    
    print(f"\nRequired for LIGO detection (1 ms threshold):")
    print(f"   β_drag needed: {beta_needed_ligo:.3e}")
    print(f"   Adjustment factor: {beta_needed_ligo/const.BETA_DRAG:.0e}×")
    
    if beta_needed_ligo / const.BETA_DRAG < 10:
        print(f"   ✓ Within reasonable parameter space")
    elif beta_needed_ligo / const.BETA_DRAG < 100:
        print(f"   ⚠️  Requires moderate adjustment")
    else:
        print(f"   ⚠️  Requires large adjustment - may indicate model issue")
    
    print(f"\nRequired for LISA detection (10 μs threshold):")
    print(f"   β_drag needed: {beta_needed_lisa:.3e}")
    print(f"   Adjustment factor: {beta_needed_lisa/const.BETA_DRAG:.0e}×")
    
    if beta_needed_lisa / const.BETA_DRAG < 10:
        print(f"   ✓ Within reasonable parameter space")
    elif beta_needed_lisa / const.BETA_DRAG < 100:
        print(f"   ⚠️  Requires moderate adjustment")
    else:
        print(f"   ⚠️  Requires large adjustment - may indicate model issue")
    
    # =========================================================================
    # Visualization
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📊 GENERATING COMPARISON PLOTS")
    print(f"{'='*80}")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Delay comparison on log scale
    detectors = ['ESTIF\nPrediction', 'LIGO\nPrecision', 'LISA\nPrecision\n(future)']
    delays_ms = [delay_estif*1000, ligo_timing_precision_ms, lisa_timing_precision_ms]
    colors = ['red', 'blue', 'green']
    
    bars1 = ax1.bar(detectors, delays_ms, color=colors, alpha=0.7, 
                    edgecolor='black', linewidth=2)
    ax1.set_yscale('log')
    ax1.set_ylabel('Time Delay / Precision (ms)', fontsize=13, fontweight='bold')
    ax1.set_title('GW Timing: ESTIF Prediction vs Detector Sensitivity', 
                  fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add detection threshold line
    ax1.axhline(ligo_timing_precision_ms, color='blue', linestyle='--', 
                alpha=0.5, linewidth=2, label='LIGO threshold')
    ax1.axhline(lisa_timing_precision_ms, color='green', linestyle='--', 
                alpha=0.5, linewidth=2, label='LISA threshold')
    
    # Add value labels
    for bar, delay in zip(bars1, delays_ms):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height*1.5,
                f'{delay:.2e}', ha='center', va='bottom', 
                fontsize=10, fontweight='bold')
    
    ax1.legend(fontsize=11)
    
    # Plot 2: Signal-to-noise ratios
    snr_labels = ['LIGO\nWaveform', 'LIGO\nPhase', 'LISA\n(future)']
    snr_values = [snr_waveform, snr_phase, snr_lisa]
    colors_snr = ['blue', 'cyan', 'green']
    
    bars2 = ax2.bar(snr_labels, snr_values, color=colors_snr, alpha=0.7,
                    edgecolor='black', linewidth=2)
    ax2.set_yscale('log')
    ax2.axhline(1.0, color='red', linestyle='--', linewidth=2, 
                label='Detection threshold (S/N = 1)')
    ax2.fill_between([-0.5, 2.5], 1, 1000, alpha=0.15, color='green', 
                     label='Detectable region')
    ax2.set_ylabel('Signal-to-Noise Ratio', fontsize=13, fontweight='bold')
    ax2.set_title('ESTIF Detectability', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(axis='y', alpha=0.3)
    
    # Add value labels
    for bar, snr in zip(bars2, snr_values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height*1.5,
                f'{snr:.2e}', ha='center', va='bottom',
                fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('ligo_gw150914_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("   ✓ Plot saved: ligo_gw150914_comparison.png")
    
    # Additional plot: Mass dependence
    print("\n   Generating mass-dependent prediction plot...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    masses_msun = np.logspace(0, 2.5, 50)  # 1 to ~300 M_sun
    delays = [estif.gw_damping_delay(0, m * const.M_sun) for m in masses_msun]
    
    ax.loglog(masses_msun, np.array(delays)*1000, 'r-', linewidth=2, 
              label='ESTIF Prediction')
    ax.axhline(ligo_timing_precision_ms, color='blue', linestyle='--', 
               linewidth=2, label='LIGO Precision (~1 ms)')
    ax.axhline(lisa_timing_precision_ms, color='green', linestyle='--', 
               linewidth=2, label='LISA Precision (~10 μs)')
    
    # Mark GW150914
    ax.plot(M_total/const.M_sun, delay_estif*1000, 'ro', markersize=12, 
            label='GW150914', zorder=5)
    
    ax.fill_between(masses_msun, ligo_timing_precision_ms, 1000, 
                    alpha=0.15, color='blue', label='LIGO detectable')
    ax.fill_between(masses_msun, lisa_timing_precision_ms, ligo_timing_precision_ms, 
                    alpha=0.15, color='green', label='LISA-only detectable')
    
    ax.set_xlabel('Binary Total Mass (M_sun)', fontsize=13, fontweight='bold')
    ax.set_ylabel('ESTIF Damping Delay (ms)', fontsize=13, fontweight='bold')
    ax.set_title('GW Damping vs Binary Mass', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1e-50, 1e3)
    
    plt.tight_layout()
    plt.savefig('gw_mass_dependence.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print("   ✓ Plot saved: gw_mass_dependence.png")
    
    # =========================================================================
    # Summary and Recommendations
    # =========================================================================
    
    print(f"\n{'='*80}")
    print("📋 SUMMARY AND RECOMMENDATIONS")
    print(f"{'='*80}")
    
    print(f"\n1. Current Status:")
    print(f"   • ESTIF predicts {delay_estif:.3e} s damping delay")
    print(f"   • LIGO S/N ratio: {snr_waveform:.3e}")
    print(f"   • LISA S/N ratio: {snr_lisa:.3e} (future)")
    
    print(f"\n2. Physical Interpretation:")
    print(f"   • Friction drag at R_s: {drag_local:.3e}")
    print(f"   • Delay scales with (R_s/c)² × drag")
    print(f"   • Comparable to ringdown time: {(delay_estif/tau_ringdown_gr):.3e}")
    
    print(f"\n3. Observational Context:")
    if status == "detectable_ligo":
        print(f"   ✅ Signal detectable with current LIGO/Virgo")
        print(f"   • Recommend: Search LIGO O1-O3 data for systematic delays")
        print(f"   • Method: Stack multiple events for statistical significance")
    elif status == "detectable_lisa":
        print(f"   ⚠️  Signal below LIGO threshold, wait for LISA")
        print(f"   • Current facilities: Cannot test")
        print(f"   • LISA (2037+): Should detect {snr_lisa:.1f}σ signal")
    elif status == "marginal_lisa":
        print(f"   ⚠️  Marginally accessible with LISA")
        print(f"   • Single event: Difficult (~{snr_lisa:.1f}σ)")
        print(f"   • Statistical sample: May work with O(100) events")
    else:
        print(f"   ❌ Signal below foreseeable detector sensitivity")
        print(f"   • LIGO: {1/snr_waveform:.0e}× too weak")
        print(f"   • LISA: {1/snr_lisa:.0e}× too weak")
    
    print(f"\n4. Model Implications:")
    if beta_needed_ligo / const.BETA_DRAG > 1e3:
        print(f"   ⚠️  CRITICAL: β_drag adjustment factor {beta_needed_ligo/const.BETA_DRAG:.0e}×")
        print(f"   • This large factor suggests:")
        print(f"     - Formula may be missing scaling factors")
        print(f"     - Physical mechanism may be different")
        print(f"     - Need to revisit friction_drag_local derivation")
    elif beta_needed_ligo / const.BETA_DRAG > 100:
        print(f"   ⚠️  Moderate adjustment needed: {beta_needed_ligo/const.BETA_DRAG:.0e}×")
        print(f"   • Consider:")
        print(f"     - Is β_drag truly constant across scales?")
        print(f"     - Are there suppression factors from GW propagation?")
    else:
        print(f"   ✓ Current β_drag in reasonable range")
        print(f"   • Adjustment factor: {beta_needed_ligo/const.BETA_DRAG:.1f}×")
    
    print(f"\n5. Recommendations:")
    print(f"\n   Immediate Actions:")
    print(f"   • Cross-check with EHT M87* results (run compare_eht_m87.py)")
    print(f"   • If EHT shows similar magnitude issues → systematic problem")
    print(f"   • If EHT works but GW doesn't → GW-specific physics")
    
    print(f"\n   Research Priorities:")
    if status == "detectable_ligo":
        print(f"   1. Analyze LIGO public data for timing residuals")
        print(f"   2. Develop waveform template with friction corrections")
        print(f"   3. Propose targeted analysis to LIGO/Virgo collaboration")
    else:
        print(f"   1. Investigate formula scaling (may need relativistic corrections)")
        print(f"   2. Check if GW propagation through 4D adds suppression factors")
        print(f"   3. Focus on other predictions (EHT, JWST) for near-term tests")
    
    print(f"\n   Long-term Strategy:")
    print(f"   • If other predictions fail similarly → model needs major revision")
    print(f"   • If other predictions work → GW physics needs deeper examination")
    print(f"   • Consider: Is friction effect stronger in EM (lensing) vs GW?")
    
    print(f"\n{'='*80}")
    print("END OF LIGO GW150914 COMPARISON")
    print(f"{'='*80}\n")
    
    return {
        'status': status,
        'snr_ligo': snr_waveform,
        'snr_lisa': snr_lisa,
        'delay': delay_estif,
        'beta_adjustment_needed': beta_needed_ligo / const.BETA_DRAG
    }


if __name__ == "__main__":
    results = main()
    

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2


