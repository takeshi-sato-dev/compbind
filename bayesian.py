"""
Bayesian analysis functions for dimerization mechanism
Using the same approach as original dimer_mech.py
"""

import os
import sys
import numpy as np
from typing import Dict, Optional
from config import AnalysisConfig

# Suppress all warnings BEFORE importing PyMC3
import warnings
warnings.filterwarnings("ignore")

# Suppress logging messages
import logging
logging.getLogger("pymc3").setLevel(logging.ERROR)
logging.getLogger("pymc").setLevel(logging.ERROR)
logging.getLogger("theano").setLevel(logging.ERROR)
logging.getLogger("theano.tensor.blas").setLevel(logging.ERROR)

# Suppress stderr for library warnings
import contextlib
import io


def bayesian_analysis_multi_regions(stats: Dict, config: AnalysisConfig) -> Optional[Dict]:
    """Perform Bayesian analysis - EXACT SAME AS ORIGINAL dimer_mech.py"""
    
    # Try to import required libraries (with suppressed output)
    try:
        # Suppress import warnings
        old_stderr = sys.stderr
        sys.stderr = io.StringIO()
        
        import pymc3 as pm
        import arviz as az
        
        sys.stderr = old_stderr
        print("    Successfully imported PyMC3 and Arviz.")
    except ImportError:
        sys.stderr = old_stderr
        print("    Libraries required for Bayesian estimation are not installed.")
        print("    Install them using: pip install pymc3 arviz")
        return None
    
    # Extract only dimer pairs
    dimer_threshold = 0.2
    pairs = [pair for pair, prob in stats['dimer_probability'].items() if prob >= dimer_threshold]
    
    # Use the first dimerizing pair
    if pairs:
        target_pair = pairs[0]
        print(f"    Focusing on {target_pair} for Bayesian analysis")
    else:
        print("    No dimerizing pairs observed. Bayesian estimation will not be performed.")
        return None
    
    print(f"    Performing Bayesian estimation for {target_pair}...")
    results = {}
    
    # Prepare data
    pair_results = {
        'dimer': {},
        'region_1': {},
        'region_2': {},
        'difference': {}
    }
    
    # Get data
    dimer_count = stats['dimer_frames'][target_pair]
    total_count = stats['total_frames']
    region_1_contact = stats.get('region_1_contact_frames', {}).get(target_pair, 0)
    region_2_contact = stats.get('region_2_contact_frames', {}).get(target_pair, 0)
    dimer_frames = stats['dimer_frames'][target_pair]
    
    print(f"    Data: Total frames={total_count}, Dimer frames={dimer_count}, "
          f"Region 1 contact frames={region_1_contact}, Region 2 contact frames={region_2_contact}")
    
    # Dimer formation model (with suppressed warnings but visible progress bar)
    try:
        # Temporarily suppress stderr for warnings
        old_stderr = sys.stderr
        temp_stderr = io.StringIO()
        
        with pm.Model() as dimer_model:
            # Non-informative prior
            theta_dimer = pm.Beta('theta_dimer', alpha=1, beta=1)
            # Binomial likelihood
            y_dimer = pm.Binomial('y_dimer', n=total_count, p=theta_dimer, observed=dimer_count)
            
            # Sampling with progress bar - redirect stderr only during sampling
            print("    Sampling dimer formation model...")
            sys.stderr = temp_stderr
            trace_dimer = pm.sample(2000, tune=1000, return_inferencedata=True, 
                                   progressbar=True, cores=4)
            sys.stderr = old_stderr
        
        print(f"    Dimer formation model sampling completed")
        
        # Get posterior distribution samples
        theta_dimer_samples = trace_dimer.posterior['theta_dimer'].values.flatten()
        
        # Check mean value of samples (for debugging)
        print(f"    Posterior mean of dimer formation probability: {np.mean(theta_dimer_samples):.4f}")
        
    except Exception as e:
        sys.stderr = old_stderr
        print(f"    Error during dimer formation model estimation: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Region 1 contact model (if we have dimer frames)
    if dimer_frames > 0:
        try:
            old_stderr = sys.stderr
            temp_stderr = io.StringIO()
            
            with pm.Model() as region1_model:
                theta_r1 = pm.Beta('theta_r1', alpha=1, beta=1)
                y_r1 = pm.Binomial('y_r1', n=dimer_frames, p=theta_r1, observed=region_1_contact)
                
                print("    Sampling region 1 contact model...")
                sys.stderr = temp_stderr
                trace_r1 = pm.sample(2000, tune=1000, return_inferencedata=True, 
                                    progressbar=True, cores=4)
                sys.stderr = old_stderr
            
            print(f"    Region 1 contact model sampling completed")
            theta_r1_samples = trace_r1.posterior['theta_r1'].values.flatten()
            print(f"    Posterior mean of region 1 contact probability: {np.mean(theta_r1_samples):.4f}")
            
        except Exception as e:
            sys.stderr = old_stderr
            print(f"    Error during region 1 contact model estimation: {e}")
            theta_r1_samples = None
        
        # Region 2 contact model
        try:
            old_stderr = sys.stderr
            temp_stderr = io.StringIO()
            
            with pm.Model() as region2_model:
                theta_r2 = pm.Beta('theta_r2', alpha=1, beta=1)
                y_r2 = pm.Binomial('y_r2', n=dimer_frames, p=theta_r2, observed=region_2_contact)
                
                print("    Sampling region 2 contact model...")
                sys.stderr = temp_stderr
                trace_r2 = pm.sample(2000, tune=1000, return_inferencedata=True, 
                                    progressbar=True, cores=4)
                sys.stderr = old_stderr
            
            print(f"    Region 2 contact model sampling completed")
            theta_r2_samples = trace_r2.posterior['theta_r2'].values.flatten()
            print(f"    Posterior mean of region 2 contact probability: {np.mean(theta_r2_samples):.4f}")
            
        except Exception as e:
            sys.stderr = old_stderr
            print(f"    Error during region 2 contact model estimation: {e}")
            theta_r2_samples = None
    else:
        theta_r1_samples = None
        theta_r2_samples = None
    
    # Free energy calculation
    try:
        # Boltzmann constant and energy unit conversion factor
        RT = config.boltzmann_const * config.temperature  # RT in kJ/mol at 298K
        
        # Dimer formation free energy
        dimer_energy_samples = []
        for p in theta_dimer_samples:
            # Limit probability to 0.01-0.99 range (for numerical stability)
            p_limited = max(0.01, min(0.99, p))
            # Free energy change: ΔG = -RT ln(p/(1-p))
            energy = -RT * np.log(p_limited / (1 - p_limited))
            dimer_energy_samples.append(energy)
        
        # Debug information
        print(f"    Dimer energy mean: {np.mean(dimer_energy_samples):.4f} kJ/mol")
        
        # Region 1 contact free energy
        if theta_r1_samples is not None:
            r1_energy_samples = []
            for p in theta_r1_samples:
                p_limited = max(0.01, min(0.99, p))
                energy = -RT * np.log(p_limited / (1 - p_limited))
                r1_energy_samples.append(energy)
            print(f"    Region 1 energy mean: {np.mean(r1_energy_samples):.4f} kJ/mol")
        else:
            r1_energy_samples = []
        
        # Region 2 contact free energy
        if theta_r2_samples is not None:
            r2_energy_samples = []
            for p in theta_r2_samples:
                p_limited = max(0.01, min(0.99, p))
                energy = -RT * np.log(p_limited / (1 - p_limited))
                r2_energy_samples.append(energy)
            print(f"    Region 2 energy mean: {np.mean(r2_energy_samples):.4f} kJ/mol")
        else:
            r2_energy_samples = []
        
        # Difference in free energy
        if r1_energy_samples and r2_energy_samples:
            diff_energy_samples = [r1 - r2 for r1, r2 in zip(r1_energy_samples, r2_energy_samples)]
            print(f"    Difference energy mean: {np.mean(diff_energy_samples):.4f} kJ/mol")
        else:
            diff_energy_samples = []
        
        print(f"    Free energy calculation completed")
        
    except Exception as e:
        print(f"    Error during free energy calculation: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Calculate statistics
    try:
        # Import arviz with suppressed warnings
        old_stderr = sys.stderr
        sys.stderr = io.StringIO()
        import arviz as az
        sys.stderr = old_stderr
        
        # Dimer formation
        pair_results['dimer']['mean'] = float(np.mean(dimer_energy_samples))
        pair_results['dimer']['median'] = float(np.median(dimer_energy_samples))
        pair_results['dimer']['hdi_low'], pair_results['dimer']['hdi_high'] = \
            [float(x) for x in az.hdi(np.array(dimer_energy_samples), 0.95)]
        pair_results['dimer']['samples'] = dimer_energy_samples
        
        # Region 1
        if r1_energy_samples:
            pair_results['region_1']['mean'] = float(np.mean(r1_energy_samples))
            pair_results['region_1']['median'] = float(np.median(r1_energy_samples))
            pair_results['region_1']['hdi_low'], pair_results['region_1']['hdi_high'] = \
                [float(x) for x in az.hdi(np.array(r1_energy_samples), 0.95)]
            pair_results['region_1']['samples'] = r1_energy_samples
        
        # Region 2
        if r2_energy_samples:
            pair_results['region_2']['mean'] = float(np.mean(r2_energy_samples))
            pair_results['region_2']['median'] = float(np.median(r2_energy_samples))
            pair_results['region_2']['hdi_low'], pair_results['region_2']['hdi_high'] = \
                [float(x) for x in az.hdi(np.array(r2_energy_samples), 0.95)]
            pair_results['region_2']['samples'] = r2_energy_samples
        
        # Difference
        if diff_energy_samples:
            pair_results['difference']['mean'] = float(np.mean(diff_energy_samples))
            pair_results['difference']['median'] = float(np.median(diff_energy_samples))
            pair_results['difference']['hdi_low'], pair_results['difference']['hdi_high'] = \
                [float(x) for x in az.hdi(np.array(diff_energy_samples), 0.95)]
            pair_results['difference']['samples'] = diff_energy_samples
            
            # Calculate p-value (whether difference differs from 0)
            median_diff = pair_results['difference']['median']
            if median_diff < 0:
                pval = np.mean([1 if d >= 0 else 0 for d in diff_energy_samples])
            else:
                pval = np.mean([1 if d <= 0 else 0 for d in diff_energy_samples])
            pair_results['difference']['pval'] = float(pval)
        
        print(f"    Statistics calculation completed")
        
    except Exception as e:
        print(f"    Error during statistics calculation: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Store pair results
    results[target_pair] = pair_results
    print(f"    Bayesian estimation for pair {target_pair} completed")
    
    return results


def plot_bayesian_results(bayesian_results: Dict, config: AnalysisConfig):
    """Plot Bayesian analysis results - ALL THREE PLOTS like original"""
    import matplotlib.pyplot as plt
    
    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    if not bayesian_results:
        return
    
    # Get the target pair
    target_pair = list(bayesian_results.keys())[0]
    pair_data = bayesian_results[target_pair]
    
    # PLOT 1: Bayesian energy distribution
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        dimer_samples = pair_data['dimer']['samples']
        
        counts, bins, patches = ax.hist(dimer_samples, bins=50, density=True, 
                                       alpha=0.7, color='#4C72B0', edgecolor='black')
        
        median = pair_data['dimer']['median']
        hdi_low = pair_data['dimer']['hdi_low']
        hdi_high = pair_data['dimer']['hdi_high']
        
        ax.axvline(median, color='red', linestyle='-', linewidth=2, label=f'Median: {median:.2f} kJ/mol')
        ax.axvline(hdi_low, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(hdi_high, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        
        y_max = ax.get_ylim()[1]
        ax.fill_betweenx([0, y_max], hdi_low, hdi_high, alpha=0.2, color='red', 
                        label=f'95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}] kJ/mol')
        
        ax.set_xlabel('Dimerization Free Energy (kJ/mol)', fontsize=14)
        ax.set_ylabel('Probability Density', fontsize=14)
        ax.set_title(f'Bayesian Posterior Distribution: {target_pair}', fontsize=16)
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        for fmt in config.figure_format:
            filename = os.path.join(output_dir, f'bayesian_energy_distribution_{target_pair}.{fmt}')
            plt.savefig(filename, dpi=300)
        
        plt.close()
        print(f"    Bayesian energy distribution plot saved")
        
    except Exception as e:
        print(f"    Error creating energy distribution plot: {e}")
    
    # PLOT 2: Bayesian dimerization energy (bar plot)
    try:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        pairs = list(bayesian_results.keys())
        median_energies = [bayesian_results[p]['dimer']['median'] for p in pairs]
        
        # Sort by energy
        pair_energy = list(zip(pairs, median_energies))
        pair_energy.sort(key=lambda x: x[1])
        sorted_pairs = [item[0] for item in pair_energy]
        sorted_medians = [item[1] for item in pair_energy]
        
        # Error bars
        yerr_low = []
        yerr_high = []
        for pair in sorted_pairs:
            yerr_low.append(bayesian_results[pair]['dimer']['median'] - bayesian_results[pair]['dimer']['hdi_low'])
            yerr_high.append(bayesian_results[pair]['dimer']['hdi_high'] - bayesian_results[pair]['dimer']['median'])
        yerr = np.vstack([yerr_low, yerr_high])
        
        # Bar plot
        bars = ax.bar(range(len(sorted_pairs)), sorted_medians, color='#4C72B0', alpha=0.8)
        ax.errorbar(range(len(sorted_pairs)), sorted_medians, yerr=yerr, fmt='none', ecolor='black', capsize=5)
        
        ax.set_xlabel('Protein Pairs')
        ax.set_ylabel('Dimerization Free Energy (kJ/mol)')
        ax.set_title('Bayesian Estimate of Dimerization Free Energy')
        ax.set_xticks(range(len(sorted_pairs)))
        ax.set_xticklabels(sorted_pairs, rotation=45, ha='right')
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.axhline(y=0, color='r', linestyle='-', alpha=0.3)
        
        plt.tight_layout()
        
        for fmt in config.figure_format:
            filename = os.path.join(output_dir, f'bayesian_dimerization_energy.{fmt}')
            plt.savefig(filename, dpi=300)
        
        plt.close()
        print(f"    Bayesian dimerization energy plot saved")
        
    except Exception as e:
        print(f"    Error creating dimerization energy plot: {e}")
    
    # PLOT 3: Bayesian energy comparison (regions)
    if 'region_1' in pair_data and 'region_2' in pair_data:
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            
            # Left: Region distributions
            r1_samples = pair_data['region_1']['samples']
            r2_samples = pair_data['region_2']['samples']
            
            region1_color = config.binding_regions[0].color if config.binding_regions else '#2271B5'
            region2_color = config.binding_regions[1].color if len(config.binding_regions) > 1 else '#F58700'
            region1_name = config.binding_regions[0].name if config.binding_regions else 'Region1'
            region2_name = config.binding_regions[1].name if len(config.binding_regions) > 1 else 'Region2'
            
            ax1.hist(r1_samples, bins=30, alpha=0.6, color=region1_color, label=region1_name)
            ax1.hist(r2_samples, bins=30, alpha=0.6, color=region2_color, label=region2_name)
            
            ax1.axvline(pair_data['region_1']['median'], color=region1_color, linestyle='-', linewidth=2)
            ax1.axvline(pair_data['region_2']['median'], color=region2_color, linestyle='-', linewidth=2)
            
            ax1.set_xlabel('Free Energy (kJ/mol)')
            ax1.set_ylabel('Frequency')
            ax1.set_title(f'Region-Specific Free Energy: {target_pair}')
            ax1.legend()
            ax1.grid(True, linestyle='--', alpha=0.7)
            
            # Right: Difference distribution
            if 'difference' in pair_data:
                diff_samples = pair_data['difference']['samples']
                
                ax2.hist(diff_samples, bins=30, alpha=0.7, color='#4DBBD5')
                ax2.axvline(pair_data['difference']['median'], color='b', linestyle='-', linewidth=2, label='Median')
                ax2.axvline(0, color='r', linestyle='-', alpha=0.7, label='Zero')
                
                ax2.set_xlabel('Energy Difference (kJ/mol)')
                ax2.set_ylabel('Frequency')
                ax2.set_title(f'Difference Distribution: {target_pair}')
                ax2.legend()
                ax2.grid(True, linestyle='--', alpha=0.7)
            
            plt.tight_layout()
            
            for fmt in config.figure_format:
                filename = os.path.join(output_dir, f'bayesian_energy_comparison_{target_pair}.{fmt}')
                plt.savefig(filename, dpi=300)
            
            plt.close()
            print(f"    Bayesian energy comparison plot saved")
            
        except Exception as e:
            print(f"    Error creating energy comparison plot: {e}")