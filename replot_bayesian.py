#!/usr/bin/env python
"""
Re-plot Bayesian results with modified visual styles
Reads saved bayesian_results.pkl and regenerates plots with customizable parameters
"""

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Font settings - CUSTOMIZE HERE
matplotlib.rcParams.update({
    'font.family': 'Arial',           # Font family: 'Arial', 'Helvetica', 'Times New Roman', etc.
    'font.size': 12,                  # Base font size
    'axes.labelsize': 14,             # Axis label font size
    'axes.titlesize': 16,             # Title font size
    'xtick.labelsize': 12,            # X-axis tick label size
    'ytick.labelsize': 12,            # Y-axis tick label size
    'legend.fontsize': 12,            # Legend font size
})

# Add compbind to path
compbind_path = '/Users/takeshi/Library/CloudStorage/OneDrive-学校法人京都薬科大学/manuscript/JOSS/dimer_mechxx4done/compbind'
if os.path.exists(compbind_path):
    sys.path.insert(0, compbind_path)

from config import AnalysisConfig

def replot_bayesian_results(results_pkl_path, output_dir, plot_config):
    """
    Re-plot Bayesian results with custom visual parameters

    Parameters
    ----------
    results_pkl_path : str
        Path to bayesian_results.pkl file
    output_dir : str
        Directory to save new plots
    plot_config : dict
        Visual parameters for plots. Available options:
        - 'hist_edgecolor': color of histogram edges (default: 'none')
        - 'hist_linewidth': width of histogram edges (default: 0)
        - 'figure_formats': list of formats to save ['png', 'pdf', 'svg']
        - 'dpi': resolution (default: 300)
        - 'alpha': histogram transparency (default: 0.7)
        - 'xlabel_fontsize': x-axis label font size (default: 20)
        - 'ylabel_fontsize': y-axis label font size (default: 20)
        - 'title_fontsize': title font size (default: 20)
        - 'legend_fontsize': legend font size (default: 20)
    """

    # Load Bayesian results
    print(f"Loading Bayesian results from {results_pkl_path}...")
    with open(results_pkl_path, 'rb') as f:
        bayesian_results = pickle.load(f)

    # Load or create config
    config = AnalysisConfig()
    config.output_dir = output_dir

    # Apply plot configuration
    hist_edgecolor = plot_config.get('hist_edgecolor', 'none')
    hist_linewidth = plot_config.get('hist_linewidth', 0)
    figure_formats = plot_config.get('figure_formats', ['png', 'pdf', 'svg'])
    dpi = plot_config.get('dpi', 300)
    alpha = plot_config.get('alpha', 0.7)
    alpha_region = plot_config.get('alpha_region', 0.6)

    # Font sizes
    xlabel_fontsize = plot_config.get('xlabel_fontsize', 20)
    ylabel_fontsize = plot_config.get('ylabel_fontsize', 20)
    title_fontsize = plot_config.get('title_fontsize', 20)
    legend_fontsize = plot_config.get('legend_fontsize', 20)

    os.makedirs(output_dir, exist_ok=True)

    if not bayesian_results:
        print("ERROR: No Bayesian results found")
        return

    # Get target pair
    target_pair = list(bayesian_results.keys())[0]
    pair_data = bayesian_results[target_pair]

    print(f"\nRe-plotting Bayesian results for {target_pair}")
    print(f"Visual parameters:")
    print(f"  - Histogram edge color: {hist_edgecolor}")
    print(f"  - Histogram edge width: {hist_linewidth}")
    print(f"  - Alpha (transparency): {alpha}")
    print(f"  - DPI: {dpi}")
    print(f"  - Formats: {figure_formats}")

    # PLOT 1: Bayesian energy distribution
    print("\n[1/3] Bayesian energy distribution...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))

        dimer_samples = pair_data['dimer']['samples']

        counts, bins, patches = ax.hist(dimer_samples, bins=50, density=True,
                                       alpha=alpha, color='#4C72B0',
                                       edgecolor=hist_edgecolor, linewidth=hist_linewidth)

        median = pair_data['dimer']['median']
        hdi_low = pair_data['dimer']['hdi_low']
        hdi_high = pair_data['dimer']['hdi_high']

        ax.axvline(median, color='red', linestyle='-', linewidth=2, label=f'Median: {median:.2f} kJ/mol')
        ax.axvline(hdi_low, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(hdi_high, color='red', linestyle='--', linewidth=1.5, alpha=0.7)

        y_max = ax.get_ylim()[1]
        ax.fill_betweenx([0, y_max], hdi_low, hdi_high, alpha=0.2, color='red',
                        label=f'95% HDI: [{hdi_low:.2f}, {hdi_high:.2f}] kJ/mol')

        ax.set_xlabel('Dimerization Free Energy (kJ/mol)', fontsize=xlabel_fontsize)
        ax.set_ylabel('Probability Density', fontsize=ylabel_fontsize)
        ax.set_title(f'Bayesian Posterior Distribution: {target_pair}', fontsize=title_fontsize)
        ax.legend(fontsize=legend_fontsize)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        for fmt in figure_formats:
            filename = os.path.join(output_dir, f'bayesian_energy_distribution_{target_pair}.{fmt}')
            plt.savefig(filename, dpi=dpi)
            print(f"  Saved: {filename}")

        plt.close()

    except Exception as e:
        print(f"  ERROR: {e}")

    # PLOT 2: Bayesian dimerization energy (bar plot)
    print("\n[2/3] Bayesian dimerization energy...")
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

        for fmt in figure_formats:
            filename = os.path.join(output_dir, f'bayesian_dimerization_energy.{fmt}')
            plt.savefig(filename, dpi=dpi)
            print(f"  Saved: {filename}")

        plt.close()

    except Exception as e:
        print(f"  ERROR: {e}")

    # PLOT 3: Bayesian energy comparison (regions)
    print("\n[3/3] Bayesian energy comparison (regional)...")
    if 'region_1' in pair_data and 'region_2' in pair_data:
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

            # Left panel: Region-Specific Free Energy Distributions
            r1_samples = pair_data['region_1']['samples']
            r2_samples = pair_data['region_2']['samples']

            # Use default colors (can be customized)
            region1_color = plot_config.get('region1_color', '#2271B5')
            region2_color = plot_config.get('region2_color', '#F58700')
            region1_name = plot_config.get('region1_name', 'TM N-term')
            region2_name = plot_config.get('region2_name', 'JM')

            # Histogram for region 1
            ax1.hist(r1_samples, bins=30, alpha=alpha_region, color=region1_color,
                    label=region1_name, edgecolor=hist_edgecolor, linewidth=hist_linewidth)
            # Histogram for region 2
            ax1.hist(r2_samples, bins=30, alpha=alpha_region, color=region2_color,
                    label=region2_name, edgecolor=hist_edgecolor, linewidth=hist_linewidth)

            # Median lines
            ax1.axvline(pair_data['region_1']['median'], color=region1_color, linestyle='-', linewidth=2.5)
            ax1.axvline(pair_data['region_2']['median'], color=region2_color, linestyle='-', linewidth=2.5)

            # HDI lines (dashed)
            ax1.axvline(pair_data['region_1']['hdi_low'], color=region1_color, linestyle='--', linewidth=1.5, alpha=0.7)
            ax1.axvline(pair_data['region_1']['hdi_high'], color=region1_color, linestyle='--', linewidth=1.5, alpha=0.7)
            ax1.axvline(pair_data['region_2']['hdi_low'], color=region2_color, linestyle='--', linewidth=1.5, alpha=0.7)
            ax1.axvline(pair_data['region_2']['hdi_high'], color=region2_color, linestyle='--', linewidth=1.5, alpha=0.7)

            ax1.set_xlabel('Free Energy (kJ/mol)', fontsize=xlabel_fontsize)
            ax1.set_ylabel('Frequency', fontsize=ylabel_fontsize)
            ax1.set_title(f'Region-Specific Free Energy Distributions: {target_pair}', fontsize=title_fontsize, fontweight='bold')
            ax1.legend(fontsize=legend_fontsize)
            ax1.grid(True, linestyle='--', alpha=0.3)

            # Right panel: Difference Distribution
            if 'difference' in pair_data:
                diff_samples = pair_data['difference']['samples']
                diff_median = pair_data['difference']['median']
                diff_hdi_low = pair_data['difference']['hdi_low']
                diff_hdi_high = pair_data['difference']['hdi_high']
                pval = pair_data['difference'].get('pval', 0)

                # Significance stars
                if pval < 0.001:
                    sig_stars = "*** p < 0.001"
                elif pval < 0.01:
                    sig_stars = "** p < 0.01"
                elif pval < 0.05:
                    sig_stars = "* p < 0.05"
                else:
                    sig_stars = "ns"

                # Histogram for difference
                ax2.hist(diff_samples, bins=30, alpha=alpha, color='#4DBBD5',
                        edgecolor=hist_edgecolor, linewidth=hist_linewidth)

                # Median line (blue)
                ax2.axvline(diff_median, color='b', linestyle='-', linewidth=2.5, label='Median')
                # Zero reference line (red)
                ax2.axvline(0, color='r', linestyle='-', linewidth=2, alpha=0.7, label='Zero')
                # HDI lines (dashed black)
                ax2.axvline(diff_hdi_low, color='black', linestyle='--', linewidth=1.5, alpha=0.7, label='95% HDI')
                ax2.axvline(diff_hdi_high, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

                ax2.set_xlabel('Energy Difference (kJ/mol)', fontsize=xlabel_fontsize)
                ax2.set_ylabel('Frequency', fontsize=ylabel_fontsize)
                ax2.set_title(f'Difference Distribution (N-C): {target_pair}\n{sig_stars}',
                            fontsize=title_fontsize, fontweight='bold')
                ax2.legend(fontsize=legend_fontsize)
                ax2.grid(True, linestyle='--', alpha=0.3)

            plt.tight_layout()

            for fmt in figure_formats:
                filename = os.path.join(output_dir, f'bayesian_energy_comparison_{target_pair}.{fmt}')
                plt.savefig(filename, dpi=dpi, bbox_inches='tight')
                print(f"  Saved: {filename}")

            plt.close()

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  WARNING: Regional data not found")

    print(f"\nAll plots saved to: {output_dir}/")


def main():
    """Example usage"""

    # Example 1: Remove histogram edges
    print("="*60)
    print("Re-plotting Bayesian results with custom styles")
    print("="*60)

    # Path to saved Bayesian results - check multiple locations
    results_pkl = None
    possible_paths = [
        "bayesian_results.pkl",                           # Current directory
        "dimerization_mechanism/bayesian_results.pkl",    # Standard location
        "../dimerization_mechanism/bayesian_results.pkl"  # If running from subdirectory
    ]

    for path in possible_paths:
        if os.path.exists(path):
            results_pkl = path
            print(f"Found Bayesian results: {path}")
            break

    if results_pkl is None:
        print(f"ERROR: bayesian_results.pkl not found in any of these locations:")
        for path in possible_paths:
            print(f"  - {path}")
        print("\nPlease run the main analysis first to generate bayesian_results.pkl")
        return 1

    # Custom plot configuration
    plot_config = {
        'hist_edgecolor': 'none',      # No black edges on histograms
        'hist_linewidth': 0,            # Edge width
        'figure_formats': ['png', 'pdf', 'svg'],  # Output formats
        'dpi': 300,                     # Resolution
        'alpha': 0.7,                   # Histogram transparency (single)
        'alpha_region': 0.6,            # Histogram transparency (regional comparison)
        'region1_color': '#2271B5',     # TM N-term color
        'region2_color': '#F58700',     # JM color
        'region1_name': 'TM N-term',    # Region 1 label
        'region2_name': 'JM'            # Region 2 label
    }

    # Output directory for re-plotted figures (same directory as pkl file)
    pkl_dir = os.path.dirname(results_pkl) if os.path.dirname(results_pkl) else "."
    output_dir = os.path.join(pkl_dir, "replotted")

    # Re-plot
    replot_bayesian_results(results_pkl, output_dir, plot_config)

    print("\n" + "="*60)
    print("Re-plotting complete!")
    print("="*60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
