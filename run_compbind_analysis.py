#!/usr/bin/env python
"""
Run CompBind analysis for DIPC/DOPS EGFR system
"""

import sys
import os

# Add compbind to path
compbind_path = '/Users/takeshi/Desktop/EGFR_MD_Exp/CG/250423EGFRdimer_DIPCDIPSCHOLDPSMGM3/compbind-main'
sys.path.insert(0, compbind_path)

from config import AnalysisConfig
from io_utils import load_universe, save_results
from core import select_proteins_and_regions
from trajectory_analysis import analyze_trajectory_detailed, collect_time_series_data
from stats_analysis import process_results_with_statistics
from visualization import (
    plot_time_series,
    plot_free_energy_with_errors,
    plot_dimerization_free_energy
)
from bayesian import bayesian_analysis_multi_regions, plot_bayesian_results
from residue_analysis import (
    analyze_residue_contacts_detailed,
    plot_residue_contact_heatmap,
    calculate_residue_importance_scores,
    plot_residue_importance
)
import pickle
import traceback


def main():
    """Run full analysis for DIPC/DOPS system"""

    print("\n" + "="*60)
    print("CompBind Analysis for DIPC/DOPS EGFR System")
    print("="*60)

    # Create configuration
    config = AnalysisConfig()

    # File paths
    config.psf_file = '../step5_assembly.psf'
    config.xtc_file = '../step7_production.xtc'
    config.output_dir = 'dimerization_mechanism'

    # Frame range
    config.start_frame = 20000
    config.stop_frame = 50000
    config.step_frame = 20

    # Analysis parameters (6Å cutoff - LIPAC compatible)
    config.protein_contact_cutoff = 6.0
    config.dimer_cutoff = 20.0

    # Temperature and constants
    config.temperature = 298.0
    config.boltzmann_const = 0.001987 * 4.184  # kJ/(mol⋅K)

    # Define binding regions
    config.add_region(
        name="TM_N_term",
        resid_start=65,
        resid_end=75,
        description="TM N-terminal interface",
        color='#2271B5'
    )

    config.add_region(
        name="JM",
        resid_start=95,
        resid_end=104,
        description="Juxtamembrane region",
        color='#F58700'
    )

    # Bayesian parameters
    config.bayesian_samples = 2000
    config.bayesian_tune = 1000

    # Output formats
    config.figure_format = ['png', 'pdf', 'ps', 'svg']
    config.dpi = 300

    os.makedirs(config.output_dir, exist_ok=True)

    print("\nConfiguration:")
    print(f"  PSF: {config.psf_file}")
    print(f"  XTC: {config.xtc_file}")
    print(f"  Output: {config.output_dir}")
    print(f"  Frames: {config.start_frame} to {config.stop_frame} (step {config.step_frame})")
    print(f"  Contact cutoff: {config.protein_contact_cutoff} Å (LIPAC-compatible)")
    for region in config.binding_regions:
        print(f"  {region.name}: residues {region.resid_start}-{region.resid_end}")
    print("="*60)

    # Validate configuration
    if not config.validate():
        print("Configuration validation failed!")
        return 1

    # Step 1: Load trajectory
    print("\n[1/7] Loading trajectory...")
    universe = load_universe(config.psf_file, config.xtc_file)
    if not universe:
        print("  Failed to load trajectory")
        return 1
    print(f"  Loaded {len(universe.trajectory)} frames")

    # Step 2: Select proteins and regions
    print("\n[2/7] Selecting proteins and regions...")

    # Create compatibility wrapper
    class CompatConfig:
        def __init__(self, config):
            self.REGION_1_NAME = config.binding_regions[0].name
            self.REGION_2_NAME = config.binding_regions[1].name
            self.REGION_1_RESID_START = config.binding_regions[0].resid_start
            self.REGION_1_RESID_END = config.binding_regions[0].resid_end
            self.REGION_2_RESID_START = config.binding_regions[1].resid_start
            self.REGION_2_RESID_END = config.binding_regions[1].resid_end

    compat_config = CompatConfig(config)
    proteins, region_1_selections, region_2_selections = select_proteins_and_regions(universe, compat_config)

    if not proteins:
        print("  No proteins found")
        return 1
    print(f"  Selected {len(proteins)} proteins")

    # Convert to new format
    region_selections = {
        config.binding_regions[0].name: region_1_selections,
        config.binding_regions[1].name: region_2_selections
    }

    # Step 3: Analyze trajectory
    print("\n[3/7] Analyzing trajectory...")
    all_results, detailed_results = analyze_trajectory_detailed(
        universe, proteins, region_selections, config, sample_frames=100
    )
    print(f"  Analyzed {len(all_results)} frames")

    # Step 4: Process results with statistics
    print("\n[4/7] Processing results with statistics...")
    stats = process_results_with_statistics(all_results, config)
    time_series = collect_time_series_data(all_results)
    print(f"  Processed {len(stats['dimer_probability'])} protein pairs")

    # Show dimers
    dimers = [(pair, prob) for pair, prob in stats['dimer_probability'].items() if prob > 0]
    if dimers:
        print("\n  Dimers found:")
        for pair, prob in dimers:
            print(f"    {pair}: {prob*100:.1f}%")
            print(f"      {config.binding_regions[0].name}: {stats['region_1_contact_prob'][pair]*100:.1f}%")
            print(f"      {config.binding_regions[1].name}: {stats['region_2_contact_prob'][pair]*100:.1f}%")

    # Step 5: Generate visualization plots
    print("\n[5/7] Generating plots...")

    print("  - Time series plots...")
    plot_time_series(time_series, config)

    print("  - Free energy plots...")
    plot_free_energy_with_errors(stats, config)

    print("  - Dimerization energy plot...")
    plot_dimerization_free_energy(stats, config)

    # Step 6: Residue analysis
    print("\n[6/7] Residue-level analysis...")
    residue_contact_maps = analyze_residue_contacts_detailed(
        universe, proteins, region_selections, config, sample_frames=50
    )

    if residue_contact_maps:
        plot_residue_contact_heatmap(residue_contact_maps, config)
        importance_scores = calculate_residue_importance_scores(residue_contact_maps, config)
        plot_residue_importance(importance_scores, config)

    # Step 7: Bayesian analysis
    print("\n[7/7] Bayesian analysis...")

    # Convert stats format for bayesian module
    stats['region_contact_frames'] = {}
    stats['region_free_energy'] = {}

    for pair in stats['dimer_probability']:
        stats['region_contact_frames'][pair] = {}
        stats['region_free_energy'][pair] = {}

        region1_key = f"{config.binding_regions[0].name}-{config.binding_regions[0].name}"
        stats['region_contact_frames'][pair][region1_key] = stats['region_1_contact_frames'][pair]
        stats['region_free_energy'][pair][region1_key] = stats['region_1_free_energy'][pair]

        region2_key = f"{config.binding_regions[1].name}-{config.binding_regions[1].name}"
        stats['region_contact_frames'][pair][region2_key] = stats['region_2_contact_frames'][pair]
        stats['region_free_energy'][pair][region2_key] = stats['region_2_free_energy'][pair]

    # Run Bayesian analysis
    try:
        print("  Running Bayesian analysis...")
        bayesian_results = bayesian_analysis_multi_regions(stats, config)

        if bayesian_results:
            # Save results
            bayesian_path = os.path.join(config.output_dir, 'bayesian_results.pkl')
            with open(bayesian_path, 'wb') as f:
                pickle.dump(bayesian_results, f)
            print(f"  Bayesian results saved: {bayesian_path}")

            # Generate plots (THIS IS THE KEY STEP)
            print("  Generating Bayesian plots...")
            plot_bayesian_results(bayesian_results, config)
            print("  Bayesian plots generated successfully!")

            # Print summary
            print("\n  Bayesian Summary:")
            for pair in bayesian_results.keys():
                if 'dimer' in bayesian_results[pair]:
                    median = bayesian_results[pair]['dimer']['median']
                    hdi_low = bayesian_results[pair]['dimer']['hdi_low']
                    hdi_high = bayesian_results[pair]['dimer']['hdi_high']
                    print(f"    {pair} dimer: {median:.2f} kJ/mol (95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}])")

                if 'region_1' in bayesian_results[pair]:
                    median = bayesian_results[pair]['region_1']['median']
                    hdi_low = bayesian_results[pair]['region_1']['hdi_low']
                    hdi_high = bayesian_results[pair]['region_1']['hdi_high']
                    print(f"    {pair} {config.binding_regions[0].name}: {median:.2f} kJ/mol (95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}])")

                if 'region_2' in bayesian_results[pair]:
                    median = bayesian_results[pair]['region_2']['median']
                    hdi_low = bayesian_results[pair]['region_2']['hdi_low']
                    hdi_high = bayesian_results[pair]['region_2']['hdi_high']
                    print(f"    {pair} {config.binding_regions[1].name}: {median:.2f} kJ/mol (95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}])")

                if 'difference' in bayesian_results[pair]:
                    median = bayesian_results[pair]['difference']['median']
                    pval = bayesian_results[pair]['difference']['pval']
                    print(f"    {pair} difference: {median:.2f} kJ/mol (p = {pval:.4f})")
        else:
            print("  Bayesian analysis returned no results")

    except Exception as e:
        print(f"  ERROR in Bayesian analysis: {e}")
        traceback.print_exc()

    # Save all results
    print("\nSaving results...")
    save_results(stats, time_series, config.output_dir, compat_config)

    # Final report
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print(f"Results saved to: {config.output_dir}/")
    print("\nGenerated files:")

    try:
        files = sorted(os.listdir(config.output_dir))
        for f in files:
            fpath = os.path.join(config.output_dir, f)
            if os.path.isfile(fpath):
                size = os.path.getsize(fpath)
                if size > 1024*1024:
                    size_str = f"{size/(1024*1024):.1f} MB"
                elif size > 1024:
                    size_str = f"{size/1024:.1f} KB"
                else:
                    size_str = f"{size} bytes"
                print(f"  - {f} ({size_str})")
    except Exception as e:
        print(f"  Could not list files: {e}")

    print("="*60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
