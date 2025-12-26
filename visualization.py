"""
Visualization functions for dimerization mechanism analysis
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List
from trajectory_analysis import determine_dimerizing_pairs

# Set plot parameters
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 12,
    'axes.linewidth': 1.5,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05
})


def plot_time_series(time_series: Dict, config):
    """Plot time series data"""
    os.makedirs(config.output_dir, exist_ok=True)
    
    dimerizing_pairs = determine_dimerizing_pairs(time_series)
    
    if not dimerizing_pairs:
        print("  No dimerizing pairs found for time series plot")
        return
    
    for pair in dimerizing_pairs:
        data = time_series[pair]
        if len(data['times']) == 0:
            continue
        
        times = np.array(data['times']) / 1000.0  # ps to ns
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
        
        # Use first region color or default
        dimer_color = '#1A9641'
        region1_color = config.binding_regions[0].color if config.binding_regions else '#2271B5'
        region2_color = config.binding_regions[1].color if len(config.binding_regions) > 1 else '#F58700'
        
        ax1.plot(times, data['dimer_status'], color=dimer_color, linewidth=1.5)
        ax1.set_ylabel('Dimer State')
        ax1.set_title(f'Time Evolution: {pair}')
        ax1.grid(True, linestyle='--', alpha=0.7)
        
        ax2.plot(times, data['region_1_contacts'], color=region1_color, linewidth=1.5)
        ax2.set_ylabel(f'{config.binding_regions[0].name if config.binding_regions else "Region1"}\nContacts')
        ax2.grid(True, linestyle='--', alpha=0.7)
        
        ax3.plot(times, data['region_2_contacts'], color=region2_color, linewidth=1.5)
        ax3.set_xlabel('Time (ns)')
        ax3.set_ylabel(f'{config.binding_regions[1].name if len(config.binding_regions) > 1 else "Region2"}\nContacts')
        ax3.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        
        for fmt in config.figure_format:
            plt.savefig(os.path.join(config.output_dir, f'time_series_{pair}.{fmt}'), format=fmt)
        plt.close()


def plot_free_energy_with_errors(stats: Dict, config):
    """Plot free energy with error bars"""
    os.makedirs(config.output_dir, exist_ok=True)
    
    pairs = [pair for pair, prob in stats['dimer_probability'].items() if prob >= 0.2]
    
    if not pairs:
        print("  No dimerizing pairs for free energy plot")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    x = np.arange(len(pairs))
    width = 0.35
    
    region1_color = config.binding_regions[0].color if config.binding_regions else '#2271B5'
    region2_color = config.binding_regions[1].color if len(config.binding_regions) > 1 else '#F58700'
    diff_color = '#4DBBD5'
    
    region_1_energies = [stats['region_1_free_energy'][pair] for pair in pairs]
    region_2_energies = [stats['region_2_free_energy'][pair] for pair in pairs]
    differences = [stats['free_energy_difference'][pair] for pair in pairs]
    
    # Error bars (ensure non-negative values)
    region_1_errors = [[abs(stats['region_1_free_energy'][p] - stats['region_1_ci_low'][p]) for p in pairs],
                       [abs(stats['region_1_ci_high'][p] - stats['region_1_free_energy'][p]) for p in pairs]]
    region_2_errors = [[abs(stats['region_2_free_energy'][p] - stats['region_2_ci_low'][p]) for p in pairs],
                       [abs(stats['region_2_ci_high'][p] - stats['region_2_free_energy'][p]) for p in pairs]]
    diff_errors = [[abs(stats['free_energy_difference'][p] - stats['diff_ci_low'][p]) for p in pairs],
                   [abs(stats['diff_ci_high'][p] - stats['free_energy_difference'][p]) for p in pairs]]
    
    ax1.bar(x - width/2, region_1_energies, width, color=region1_color, 
            label=config.binding_regions[0].name if config.binding_regions else 'Region1', 
            yerr=region_1_errors)
    ax1.bar(x + width/2, region_2_energies, width, color=region2_color, 
            label=config.binding_regions[1].name if len(config.binding_regions) > 1 else 'Region2', 
            yerr=region_2_errors)
    ax1.set_xlabel('Protein Pairs')
    ax1.set_ylabel('Free Energy (kJ/mol)')
    ax1.set_title('Free Energy Comparison by Region')
    ax1.set_xticks(x)
    ax1.set_xticklabels(pairs, rotation=45)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    ax2.bar(x, differences, color=diff_color, yerr=diff_errors)
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax2.set_xlabel('Protein Pairs')
    ax2.set_ylabel('Delta G Difference (kJ/mol)')
    ax2.set_title('Free Energy Difference')
    ax2.set_xticks(x)
    ax2.set_xticklabels(pairs, rotation=45)
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    for fmt in config.figure_format:
        plt.savefig(os.path.join(config.output_dir, f'free_energy_comparison.{fmt}'), format=fmt)
    plt.close()


def plot_dimerization_free_energy(stats: Dict, config):
    """Plot dimerization free energy"""
    os.makedirs(config.output_dir, exist_ok=True)
    
    pairs = [pair for pair, prob in stats['dimer_probability'].items() if prob >= 0.01]
    
    if not pairs:
        print("  No pairs for dimerization energy plot")
        return
    
    RT = config.boltzmann_const * config.temperature
    energies = []
    
    for pair in pairs:
        prob = max(0.01, min(0.99, stats['dimer_probability'][pair]))
        energy = -RT * np.log(prob / (1 - prob))
        energies.append(energy)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sorted_data = sorted(zip(pairs, energies), key=lambda x: x[1])
    sorted_pairs = [x[0] for x in sorted_data]
    sorted_energies = [x[1] for x in sorted_data]
    
    ax.bar(range(len(sorted_pairs)), sorted_energies, color='#4C72B0', alpha=0.8)
    ax.set_xlabel('Protein Pairs')
    ax.set_ylabel('Dimerization Free Energy (kJ/mol)')
    ax.set_title('Free Energy of Dimerization')
    ax.set_xticks(range(len(sorted_pairs)))
    ax.set_xticklabels(sorted_pairs, rotation=45, ha='right')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    for fmt in config.figure_format:
        plt.savefig(os.path.join(config.output_dir, f'dimerization_free_energy.{fmt}'), format=fmt)
    plt.close()