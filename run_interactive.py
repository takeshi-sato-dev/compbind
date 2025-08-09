#!/usr/bin/env python
"""Interactive runner for dimer mechanism analysis using modular architecture."""

import sys
import os
import time
import numpy as np
import warnings
import pickle
import traceback

# Import ALL modules
from config import AnalysisConfig, BindingRegion
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


def get_user_config():
    """Get configuration from user input."""
    
    print("\n" + "="*60)
    print("DIMER MECHANISM ANALYSIS - CONFIGURATION")
    print("="*60)
    
    config = AnalysisConfig()
    
    # File paths
    psf_file = input("\nEnter PSF file path [test_data/test_system.psf]: ").strip()
    config.psf_file = psf_file if psf_file else "test_data/test_system.psf"
    
    xtc_file = input("Enter XTC file path [test_data/test_trajectory.xtc]: ").strip()
    config.xtc_file = xtc_file if xtc_file else "test_data/test_trajectory.xtc"
    
    output_dir = input("Enter output directory [test_output]: ").strip()
    config.output_dir = output_dir if output_dir else "test_output"
    
    # Frame range
    print("\nFrame range for analysis:")
    start_frame = input("  Start frame [0]: ").strip()
    config.start_frame = int(start_frame) if start_frame else 0
    
    stop_frame = input("  Stop frame [100]: ").strip()
    config.stop_frame = int(stop_frame) if stop_frame else 100
    
    step_frame = input("  Step [5]: ").strip()
    config.step_frame = int(step_frame) if step_frame else 5
    
    # Region definitions
    print("\nDefine analysis regions (residue numbers):")
    print("Region 1 (e.g., N-terminal domain):")
    region_1_start = input("  Start residue [65]: ").strip()
    region_1_start = int(region_1_start) if region_1_start else 65
    
    region_1_end = input("  End residue [75]: ").strip()
    region_1_end = int(region_1_end) if region_1_end else 75
    
    print("Region 2 (e.g., C-terminal domain):")
    region_2_start = input("  Start residue [76]: ").strip()
    region_2_start = int(region_2_start) if region_2_start else 76
    
    region_2_end = input("  End residue [91]: ").strip()
    region_2_end = int(region_2_end) if region_2_end else 91
    
    # Add regions
    config.add_region(
        name="TM_N_term",
        resid_start=region_1_start,
        resid_end=region_1_end,
        description="N-terminal TM region",
        color='#2271B5'
    )
    
    config.add_region(
        name="TM_C_term",
        resid_start=region_2_start,
        resid_end=region_2_end,
        description="C-terminal TM region",
        color='#F58700'
    )
    
    config.figure_format = ['png', 'svg']
    config.bayesian_samples = 2000
    config.bayesian_tune = 1000
    
    os.makedirs(config.output_dir, exist_ok=True)
    
    print("\n" + "="*60)
    print("Configuration Summary:")
    print(f"  PSF: {config.psf_file}")
    print(f"  XTC: {config.xtc_file}")
    print(f"  Output: {config.output_dir}")
    print(f"  Frames: {config.start_frame} to {config.stop_frame} (step {config.step_frame})")
    for region in config.binding_regions:
        print(f"  {region.name}: residues {region.resid_start}-{region.resid_end}")
    print("="*60)
    
    return config


def convert_region_selections(proteins, region_1_selections, region_2_selections, config):
    """Convert old-style region selections to new format"""
    region_selections = {}
    
    if config.binding_regions:
        region_selections[config.binding_regions[0].name] = region_1_selections
        if len(config.binding_regions) > 1:
            region_selections[config.binding_regions[1].name] = region_2_selections
    
    return region_selections


def run_full_analysis(config):
    """Run complete analysis pipeline using all modules."""
    
    print("\n" + "="*60)
    print("RUNNING FULL ANALYSIS")
    print("="*60)
    
    # Step 1: Load trajectory
    print("\n[1/7] Loading trajectory...")
    universe = load_universe(config.psf_file, config.xtc_file)
    if not universe:
        print("  Failed to load trajectory")
        return False
    print(f"  Loaded {len(universe.trajectory)} frames")
    
    # Step 2: Select proteins and regions
    print("\n[2/7] Selecting proteins and regions...")
    
    class CompatConfig:
        def __init__(self, config):
            self.REGION_1_NAME = config.binding_regions[0].name if config.binding_regions else "Region1"
            self.REGION_2_NAME = config.binding_regions[1].name if len(config.binding_regions) > 1 else "Region2"
            self.REGION_1_RESID_START = config.binding_regions[0].resid_start if config.binding_regions else 1
            self.REGION_1_RESID_END = config.binding_regions[0].resid_end if config.binding_regions else 10
            self.REGION_2_RESID_START = config.binding_regions[1].resid_start if len(config.binding_regions) > 1 else 11
            self.REGION_2_RESID_END = config.binding_regions[1].resid_end if len(config.binding_regions) > 1 else 20
    
    compat_config = CompatConfig(config)
    proteins, region_1_selections, region_2_selections = select_proteins_and_regions(universe, compat_config)
    
    if not proteins:
        print("  No proteins found")
        return False
    print(f"  Selected {len(proteins)} proteins")
    
    region_selections = convert_region_selections(proteins, region_1_selections, region_2_selections, config)
    
    # Step 3: Analyze trajectory
    print("\n[3/7] Analyzing trajectory...")
    all_results, detailed_results = analyze_trajectory_detailed(
        universe, proteins, region_selections, config, sample_frames=50
    )
    print(f"  Analyzed {len(all_results)} frames")
    
    # Step 4: Process results
    print("\n[4/7] Processing results with statistics...")
    stats = process_results_with_statistics(all_results, config)
    time_series = collect_time_series_data(all_results)
    print(f"  Processed {len(stats['dimer_probability'])} protein pairs")
    
    # Show which pairs are dimers
    dimers = [(pair, prob) for pair, prob in stats['dimer_probability'].items() if prob > 0]
    if dimers:
        print("  Dimers found:")
        for pair, prob in dimers:
            print(f"    {pair}: {prob*100:.1f}%")
    
    # Step 5: Generate plots
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
        universe, proteins, region_selections, config, sample_frames=20
    )
    
    if residue_contact_maps:
        plot_residue_contact_heatmap(residue_contact_maps, config)
        importance_scores = calculate_residue_importance_scores(residue_contact_maps, config)
        plot_residue_importance(importance_scores, config)
    
    # Step 7: Bayesian analysis using the bayesian module
    print("\n[7/7] Bayesian analysis...")
    
    # Convert stats format for bayesian module
    stats['region_contact_frames'] = {}
    stats['region_free_energy'] = {}
    
    for pair in stats['dimer_probability']:
        stats['region_contact_frames'][pair] = {}
        stats['region_free_energy'][pair] = {}
        
        if config.binding_regions:
            region1_key = f"{config.binding_regions[0].name}-{config.binding_regions[0].name}"
            stats['region_contact_frames'][pair][region1_key] = stats['region_1_contact_frames'][pair]
            stats['region_free_energy'][pair][region1_key] = stats['region_1_free_energy'][pair]
            
            if len(config.binding_regions) > 1:
                region2_key = f"{config.binding_regions[1].name}-{config.binding_regions[1].name}"
                stats['region_contact_frames'][pair][region2_key] = stats['region_2_contact_frames'][pair]
                stats['region_free_energy'][pair][region2_key] = stats['region_2_free_energy'][pair]
    
    # Run Bayesian analysis
    try:
        print("  Running Bayesian analysis from module...")
        bayesian_results = bayesian_analysis_multi_regions(stats, config)
        
        if bayesian_results:
            # Save results
            bayesian_path = os.path.join(config.output_dir, 'bayesian_results.pkl')
            with open(bayesian_path, 'wb') as f:
                pickle.dump(bayesian_results, f)
            print(f"  Bayesian results saved: {bayesian_path}")
            
            # Generate plots
            try:
                plot_bayesian_results(bayesian_results, config)
                print("  Bayesian plots generated")
            except Exception as e:
                print(f"  Warning: Could not generate Bayesian plots: {e}")
            
            # Print summary
            print("\n  Bayesian Summary:")
            for pair in list(bayesian_results.keys())[:5]:
                if 'dimer' in bayesian_results[pair]:
                    print(f"    {pair} dimer: {bayesian_results[pair]['dimer']['median']:.2f} kJ/mol")
                for region in config.binding_regions:
                    region_key = f"{region.name}-{region.name}"
                    if region_key in bayesian_results[pair]:
                        print(f"    {pair} {region.name}: {bayesian_results[pair][region_key]['median']:.2f} kJ/mol")
        else:
            print("  Bayesian analysis returned no results")
            
    except Exception as e:
        print(f"  ERROR in Bayesian analysis: {e}")
        traceback.print_exc()
    
    # Save all results
    print("\nSaving results...")
    save_results(stats, time_series, config.output_dir, compat_config)
    
    # Final report - list actual files
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print(f"Results saved to: {config.output_dir}/")
    print("\nGenerated files:")
    
    try:
        files = sorted(os.listdir(config.output_dir))
        for f in files:
            size = os.path.getsize(os.path.join(config.output_dir, f))
            if size > 1024*1024:
                size_str = f"{size/(1024*1024):.1f} MB"
            elif size > 1024:
                size_str = f"{size/1024:.1f} KB"
            else:
                size_str = f"{size} bytes"
            print(f"  - {f} ({size_str})")
    except:
        print("  Could not list files")
    
    print("="*60)
    
    return True


def main():
    """Main interactive function."""
    
    # Check environment
    if 'VIRTUAL_ENV' in os.environ:
        venv = os.path.basename(os.environ['VIRTUAL_ENV'])
        print(f"Using virtual environment: {venv}")
        if 'bayesian' not in venv.lower():
            print("WARNING: Not in bayesian_env - Bayesian analysis may not work!")
            print("Activate with: source bayesian_env/bin/activate")
    else:
        print("WARNING: No virtual environment detected")
        print("Bayesian analysis requires: source bayesian_env/bin/activate")
    
    config = get_user_config()
    
    print("\n" + "="*60)
    print("READY FOR ANALYSIS")
    print("="*60)
    
    while True:
        print("\nOptions:")
        print("  1. Run full analysis")
        print("  2. Show configuration")
        print("  3. Save configuration")
        print("  4. Exit")
        
        choice = input("\nEnter choice [1]: ").strip()
        if not choice:
            choice = "1"
        
        if choice == "1":
            start_time = time.time()
            success = run_full_analysis(config)
            if success:
                elapsed = time.time() - start_time
                print(f"\nTotal time: {elapsed:.1f} seconds")
            break
        
        elif choice == "2":
            print("\nConfiguration:")
            print(f"  PSF: {config.psf_file}")
            print(f"  XTC: {config.xtc_file}")
            print(f"  Output: {config.output_dir}")
            print(f"  Frames: {config.start_frame}-{config.stop_frame} (step {config.step_frame})")
            for region in config.binding_regions:
                print(f"  {region.name}: residues {region.resid_start}-{region.resid_end}")
        
        elif choice == "3":
            config_file = os.path.join(config.output_dir, 'analysis_config.yaml')
            config.to_yaml(config_file)
            print(f"Configuration saved to: {config_file}")
        
        elif choice == "4":
            print("Exiting...")
            break
    
    return 0


if __name__ == "__main__":
    sys.exit(main())