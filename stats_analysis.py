"""
Statistical analysis functions for dimerization mechanism
"""

import numpy as np
from typing import Dict, List, Tuple
from scipy import stats


def calculate_free_energy(contact_prob: float, temperature: float = 298.0, 
                         boltzmann_const: float = 0.001987 * 4.184) -> float:
    """Calculate free energy from contact probability"""
    contact_prob = max(0.01, min(0.99, contact_prob))
    RT = boltzmann_const * temperature
    return -RT * np.log(contact_prob / (1 - contact_prob))


def bootstrap_free_energy(contact_frames: int, total_frames: int, 
                         temperature: float = 298.0,
                         boltzmann_const: float = 0.001987 * 4.184,
                         n_bootstraps: int = 1000) -> Tuple[float, float, float]:
    """Bootstrap estimation of free energy with confidence intervals"""
    data = np.zeros(total_frames)
    data[:contact_frames] = 1
    
    energies = []
    for _ in range(n_bootstraps):
        sample = np.random.choice(data, size=total_frames, replace=True)
        prob = np.mean(sample)
        energy = calculate_free_energy(prob, temperature, boltzmann_const)
        energies.append(energy)
    
    energies = np.array(energies)
    mean_energy = np.mean(energies)
    ci_low = np.percentile(energies, 2.5)
    ci_high = np.percentile(energies, 97.5)
    
    return mean_energy, ci_low, ci_high


def calculate_p_value(region1_contacts: int, region2_contacts: int, 
                     total_frames: int) -> float:
    """Calculate p-value for difference between two regions"""
    contingency_table = np.array([
        [region1_contacts, total_frames - region1_contacts],
        [region2_contacts, total_frames - region2_contacts]
    ])
    
    _, p_value = stats.fisher_exact(contingency_table)
    return p_value


def process_results_with_statistics(all_results: List[Dict], config) -> Dict:
    """Process results with statistical analysis"""
    protein_pairs = set()
    for result in all_results:
        protein_pairs.update(result['dimers'].keys())
    
    stats_dict = {
        'total_frames': len(all_results),
        'dimer_frames': {pair: 0 for pair in protein_pairs},
        'dimer_probability': {pair: 0.0 for pair in protein_pairs},
        'region_1_contact_frames': {pair: 0 for pair in protein_pairs},
        'region_2_contact_frames': {pair: 0 for pair in protein_pairs},
        'region_1_contact_prob': {pair: 0.0 for pair in protein_pairs},
        'region_2_contact_prob': {pair: 0.0 for pair in protein_pairs},
        'region_1_free_energy': {pair: 0.0 for pair in protein_pairs},
        'region_2_free_energy': {pair: 0.0 for pair in protein_pairs},
        'free_energy_difference': {pair: 0.0 for pair in protein_pairs},
        'region_1_ci_low': {pair: 0.0 for pair in protein_pairs},
        'region_1_ci_high': {pair: 0.0 for pair in protein_pairs},
        'region_2_ci_low': {pair: 0.0 for pair in protein_pairs},
        'region_2_ci_high': {pair: 0.0 for pair in protein_pairs},
        'diff_ci_low': {pair: 0.0 for pair in protein_pairs},
        'diff_ci_high': {pair: 0.0 for pair in protein_pairs},
        'p_value': {pair: 1.0 for pair in protein_pairs}
    }
    
    # Count frames (updated for Boolean contact detection - LIPAC-compatible)
    for result in all_results:
        for pair in protein_pairs:
            if pair in result['dimers'] and result['dimers'][pair]:
                stats_dict['dimer_frames'][pair] += 1
                # Binary contact detection: True/False
                if pair in result['region_1_contacts'] and result['region_1_contacts'][pair]:
                    stats_dict['region_1_contact_frames'][pair] += 1
                if pair in result['region_2_contacts'] and result['region_2_contacts'][pair]:
                    stats_dict['region_2_contact_frames'][pair] += 1
    
    # Calculate statistics
    for pair in protein_pairs:
        stats_dict['dimer_probability'][pair] = stats_dict['dimer_frames'][pair] / stats_dict['total_frames']
        
        if stats_dict['dimer_frames'][pair] > 0:
            stats_dict['region_1_contact_prob'][pair] = stats_dict['region_1_contact_frames'][pair] / stats_dict['dimer_frames'][pair]
            stats_dict['region_2_contact_prob'][pair] = stats_dict['region_2_contact_frames'][pair] / stats_dict['dimer_frames'][pair]
            
            # Bootstrap analysis
            r1_energy, r1_ci_low, r1_ci_high = bootstrap_free_energy(
                stats_dict['region_1_contact_frames'][pair], 
                stats_dict['dimer_frames'][pair],
                config.temperature,
                config.boltzmann_const
            )
            stats_dict['region_1_free_energy'][pair] = r1_energy
            stats_dict['region_1_ci_low'][pair] = r1_ci_low
            stats_dict['region_1_ci_high'][pair] = r1_ci_high
            
            r2_energy, r2_ci_low, r2_ci_high = bootstrap_free_energy(
                stats_dict['region_2_contact_frames'][pair], 
                stats_dict['dimer_frames'][pair],
                config.temperature,
                config.boltzmann_const
            )
            stats_dict['region_2_free_energy'][pair] = r2_energy
            stats_dict['region_2_ci_low'][pair] = r2_ci_low
            stats_dict['region_2_ci_high'][pair] = r2_ci_high
            
            stats_dict['free_energy_difference'][pair] = r1_energy - r2_energy
            
            # Difference confidence intervals
            diff_samples = np.random.normal(r1_energy, (r1_ci_high - r1_ci_low) / 3.92, 1000) - \
                          np.random.normal(r2_energy, (r2_ci_high - r2_ci_low) / 3.92, 1000)
            stats_dict['diff_ci_low'][pair] = np.percentile(diff_samples, 2.5)
            stats_dict['diff_ci_high'][pair] = np.percentile(diff_samples, 97.5)
            
            # P-value
            stats_dict['p_value'][pair] = calculate_p_value(
                stats_dict['region_1_contact_frames'][pair],
                stats_dict['region_2_contact_frames'][pair],
                stats_dict['dimer_frames'][pair]
            )
    
    return stats_dict