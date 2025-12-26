"""
Input/Output utilities for dimerization mechanism analysis.
"""

import os
import pickle
import json
import traceback
import pandas as pd
import MDAnalysis as mda
import numpy as np


def load_universe(top_file, traj_file):
    """Load trajectory files into MDAnalysis Universe."""
    print(f"Loading topology: {top_file}")
    print(f"Loading trajectory: {traj_file}")
    
    try:
        universe = mda.Universe(top_file, traj_file)
        print(f"Trajectory loaded successfully with {len(universe.trajectory)} frames.")
        return universe
    except Exception as e:
        print(f"Error loading universe: {str(e)}")
        traceback.print_exc()
        return None


def save_results(stats, time_series, output_dir, config):
    """Save analysis results to CSV and pickle files."""
    REGION_1_NAME = getattr(config, 'REGION_1_NAME', 'Region1')
    REGION_2_NAME = getattr(config, 'REGION_2_NAME', 'Region2')
    
    dimer_threshold = 0.2
    pairs = [pair for pair, prob in stats['dimer_probability'].items() 
             if prob >= dimer_threshold]
    
    if not pairs:
        pairs = list(stats['dimer_probability'].keys())
    
    df = pd.DataFrame({
        'Protein_Pair': pairs,
        'Dimer_Probability': [stats['dimer_probability'][pair] for pair in pairs],
        f'{REGION_1_NAME}_Contact_Probability': [stats['region_1_contact_prob'][pair] for pair in pairs],
        f'{REGION_2_NAME}_Contact_Probability': [stats['region_2_contact_prob'][pair] for pair in pairs],
        f'{REGION_1_NAME}_Free_Energy_kJ_mol': [stats['region_1_free_energy'][pair] for pair in pairs],
        f'{REGION_2_NAME}_Free_Energy_kJ_mol': [stats['region_2_free_energy'][pair] for pair in pairs],
        'Free_Energy_Difference_kJ_mol': [stats['free_energy_difference'][pair] for pair in pairs]
    })
    
    csv_path = os.path.join(output_dir, 'contact_free_energy_analysis.csv')
    df.to_csv(csv_path, index=False)
    print(f"Results saved to {csv_path}")
    
    pickle_path = os.path.join(output_dir, 'contact_analysis_full.pkl')
    with open(pickle_path, 'wb') as f:
        pickle.dump(stats, f)
    
    time_series_path = os.path.join(output_dir, 'time_series_data.pkl')
    with open(time_series_path, 'wb') as f:
        pickle.dump(time_series, f)
    
    params = {
        'REGION_1_NAME': REGION_1_NAME,
        'REGION_1_RESID_START': getattr(config, 'REGION_1_RESID_START', 0),
        'REGION_1_RESID_END': getattr(config, 'REGION_1_RESID_END', 0),
        'REGION_2_NAME': REGION_2_NAME,
        'REGION_2_RESID_START': getattr(config, 'REGION_2_RESID_START', 0),
        'REGION_2_RESID_END': getattr(config, 'REGION_2_RESID_END', 0),
        'PROTEIN_CONTACT_CUTOFF': getattr(config, 'PROTEIN_CONTACT_CUTOFF', 6.0),
        'DIMER_CUTOFF': getattr(config, 'DIMER_CUTOFF', 20.0),
        'START_FRAME': getattr(config, 'START_FRAME', 0),
        'STOP_FRAME': getattr(config, 'STOP_FRAME', -1),
        'STEP_FRAME': getattr(config, 'STEP_FRAME', 1),
        'TEMPERATURE': getattr(config, 'TEMPERATURE', 298),
        'ENERGY_UNIT': 'kJ/mol'
    }
    
    params_path = os.path.join(output_dir, 'analysis_parameters.json')
    with open(params_path, 'w') as f:
        json.dump(params, f, indent=4)
    
    print("\nSummary of Results:")
    print(df)
    
    return df