# run_all_comparisons.py

"""
Master script to run all ESTIF-Gravity observational comparisons.

Executes:
1. compare_eht_m87.py - Black hole shadow lensing
2. compare_ligo_gw.py - Gravitational wave damping
3. compare_jwst_galaxies.py - High-z galaxy asymmetries

Generates comprehensive assessment of ESTIF testability.
"""

import sys
import os
import subprocess
import json
from datetime import datetime

# Add src directory to Python path for subprocess calls
src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src')
sys.path.insert(0, src_path)

def run_comparison(script_name, description):
    """Run a comparison script and capture results."""
    print(f"\n{'='*80}")
    print(f"RUNNING: {description}")
    print(f"Script: {script_name}")
    print(f"{'='*80}\n")
    
    # Determine the full path to the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, script_name)
    
    try:
        result = subprocess.run(
            ['python3', script_path],
            capture_output=False,
            text=True
        )
        
        if result.returncode == 0:
            print(f"\n✓ {script_name} completed successfully")
            return True
        else:
            print(f"\n❌ ERROR running {script_name}")
            print(f"Return code: {result.returncode}")
            return False
            
    except Exception as e:
        print(f"\n❌ ERROR running {script_name}")
        print(f"Exception: {e}")
        return False

def main():
    """Run all three comparison scripts."""
    print("="*80)
    print("ESTIF-GRAVITY: COMPREHENSIVE OBSERVATIONAL VALIDATION")
    print("="*80)
    print(f"\nDate: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Running all three observational comparison scripts...")
    print("="*80)
    
    results = {}
    
    # Run EHT M87* comparison
    results['eht'] = run_comparison('observational/compare_eht_m87.py', 'EHT M87* Black Hole Shadow')
    
    # Run LIGO GW comparison
    results['ligo'] = run_comparison('observational/compare_ligo_gw.py', 'LIGO GW150914 Gravitational Waves')
    
    # Run JWST galaxy comparison
    results['jwst'] = run_comparison('observational/compare_jwst_galaxies.py', 'JWST High-Redshift Galaxies')
    
    # Summary
    print("\n" + "="*80)
    print("COMPREHENSIVE SUMMARY")
    print("="*80)
    
    successful = sum([1 for v in results.values() if v is True])
    total = 3
    
    print(f"\nExecution Summary:")
    print(f"   Scripts completed: {successful}/{total}")
    print(f"   EHT M87*: {'✓' if results['eht'] else '✗'}")
    print(f"   LIGO GW: {'✓' if results['ligo'] else '✗'}")
    print(f"   JWST:    {'✓' if results['jwst'] else '✗'}")
    
    print(f"\n📊 Generated Files:")
    print(f"   • eht_m87_comparison.png")
    print(f"   • ligo_gw150914_comparison.png")
    print(f"   • gw_mass_dependence.png")
    print(f"   • jwst_ceers_comparison.png")
    
    print(f"\n📋 Next Steps:")
    print(f"   1. Review all three comparison reports above")
    print(f"   2. Identify which prediction has strongest signal")
    print(f"   3. Look for systematic patterns:")
    print(f"      - If all ~0.01%: Formula needs major revision")
    print(f"      - If all ~1%: At edge of detectability")
    print(f"      - If mixed: Physics-dependent effects")
    print(f"   4. Prioritize observational follow-up on strongest prediction")
    
    print(f"\n🔬 Research Strategy:")
    print(f"   • If 0/3 detectable → Model needs fundamental revision")
    print(f"   • If 1/3 detectable → Focus on that channel, investigate others")
    print(f"   • If 2/3 detectable → Strong case for ESTIF, investigate failure")
    print(f"   • If 3/3 detectable → Proceed with observational campaign")
    
    print(f"\n{'='*80}")
    print("END OF COMPREHENSIVE VALIDATION")
    print(f"{'='*80}\n")
    
    # Save summary
    with open('validation_summary.txt', 'w') as f:
        f.write(f"ESTIF-Gravity Validation Summary\n")
        f.write(f"{'='*80}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Scripts executed:\n")
        f.write(f"  EHT M87*: {'SUCCESS' if results['eht'] else 'FAILED'}\n")
        f.write(f"  LIGO GW:  {'SUCCESS' if results['ligo'] else 'FAILED'}\n")
        f.write(f"  JWST:     {'SUCCESS' if results['jwst'] else 'FAILED'}\n")
        f.write(f"\nReview individual script outputs for detailed results.\n")
    
    print("Summary saved to: validation_summary.txt\n")

if __name__ == "__main__":
    main()
    

#APPROVED-FORK-CONVERSION-SYNTAX-PROVEN-15-10-25-V-2


