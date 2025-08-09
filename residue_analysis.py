"""
Residue-level contact analysis functions
Provides detailed per-residue contact information
"""

import numpy as np
import MDAnalysis as mda
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from config import AnalysisConfig, BindingRegion


def calculate_residue_contacts(region1: mda.AtomGroup, region2: mda.AtomGroup,
                              box: np.ndarray, cutoff: float) -> np.ndarray:
    """Calculate residue-level contact matrix between two regions"""
    if len(region1) == 0 or len(region2) == 0:
        return np.array([])
    
    try:
        # Initialize contact matrix
        contact_matrix = np.zeros((len(region1.residues), len(region2.residues)))
        
        # Calculate minimum distance between each residue pair
        for i, res1 in enumerate(region1.residues):
            for j, res2 in enumerate(region2.residues):
                min_dist = float('inf')
                
                if len(res1.atoms) == 0 or len(res2.atoms) == 0:
                    continue
                
                # Quick screening with residue COMs
                res1_com = res1.atoms.center_of_mass()
                res2_com = res2.atoms.center_of_mass()
                
                # PBC correction
                res_diff = res1_com - res2_com
                for dim in range(3):
                    if res_diff[dim] > box[dim] * 0.5:
                        res_diff[dim] -= box[dim]
                    elif res_diff[dim] < -box[dim] * 0.5:
                        res_diff[dim] += box[dim]
                
                res_com_dist = np.sqrt(np.sum(res_diff * res_diff))
                
                # Skip if COM distance is too large
                max_atom_dist = 10.0
                if res_com_dist > (cutoff + max_atom_dist):
                    continue
                
                # Calculate actual minimum distance
                for atom1 in res1.atoms:
                    for atom2 in res2.atoms:
                        diff = atom1.position - atom2.position
                        for dim in range(3):
                            if diff[dim] > box[dim] * 0.5:
                                diff[dim] -= box[dim]
                            elif diff[dim] < -box[dim] * 0.5:
                                diff[dim] += box[dim]
                        
                        dist = np.sqrt(np.sum(diff * diff))
                        min_dist = min(min_dist, dist)
                        
                        if min_dist <= cutoff:
                            break
                    
                    if min_dist <= cutoff:
                        break
                
                # Mark as contact if within cutoff
                if min_dist <= cutoff:
                    contact_matrix[i, j] = 1
        
        return contact_matrix
    
    except Exception as e:
        print(f"Error in residue contact calculation: {str(e)}")
        return np.array([])


def analyze_residue_contacts_detailed(universe: mda.Universe, proteins: Dict,
                                     region_selections: Dict, config: AnalysisConfig,
                                     sample_frames: Optional[int] = None) -> Dict:
    """Perform detailed residue-level contact analysis"""
    
    # Determine frames to analyze
    max_frame = len(universe.trajectory)
    effective_stop = min(config.stop_frame, max_frame)
    frames = list(range(config.start_frame, effective_stop, config.step_frame))
    
    # Sample frames if requested
    if sample_frames and sample_frames < len(frames):
        sample_indices = np.linspace(0, len(frames)-1, sample_frames, dtype=int)
        frames = [frames[i] for i in sample_indices]
        print(f"Sampling {len(frames)} frames for residue-level analysis")
    
    # Initialize results
    residue_contact_maps = {}
    protein_pairs = []
    
    protein_names = list(proteins.keys())
    for i in range(len(protein_names)):
        for j in range(i + 1, len(protein_names)):
            pair_name = f"{protein_names[i]}-{protein_names[j]}"
            protein_pairs.append(pair_name)
            residue_contact_maps[pair_name] = {}
            
            # Initialize contact maps for each region
            for region in config.binding_regions:
                region_key = f"{region.name}-{region.name}"
                residue_contact_maps[pair_name][region_key] = {
                    'contact_count': None,
                    'contact_prob': None,
                    'dimer_frames': 0,
                    'resids_1': [],
                    'resids_2': []
                }
    
    # Analyze each frame
    print("Performing residue-level contact analysis...")
    for frame_idx in tqdm(frames, desc="Analyzing residue contacts"):
        try:
            universe.trajectory[frame_idx]
            box = universe.dimensions[:3]
            
            # Check each protein pair
            for i in range(len(protein_names)):
                for j in range(i + 1, len(protein_names)):
                    protein1_name = protein_names[i]
                    protein2_name = protein_names[j]
                    pair_name = f"{protein1_name}-{protein2_name}"
                    
                    # Check if dimer
                    from core import is_dimer
                    dimer_status = is_dimer(
                        proteins[protein1_name],
                        proteins[protein2_name],
                        box,
                        config.dimer_cutoff
                    )
                    
                    if dimer_status:
                        # Analyze each region
                        for region in config.binding_regions:
                            region_key = f"{region.name}-{region.name}"
                            
                            if (protein1_name in region_selections.get(region.name, {}) and
                                protein2_name in region_selections.get(region.name, {})):
                                
                                # Calculate residue contact matrix
                                contact_matrix = calculate_residue_contacts(
                                    region_selections[region.name][protein1_name],
                                    region_selections[region.name][protein2_name],
                                    box,
                                    config.protein_contact_cutoff
                                )
                                
                                if len(contact_matrix) > 0:
                                    # Initialize or accumulate
                                    if residue_contact_maps[pair_name][region_key]['contact_count'] is None:
                                        residue_contact_maps[pair_name][region_key]['contact_count'] = contact_matrix
                                        # Get residue IDs
                                        residue_contact_maps[pair_name][region_key]['resids_1'] = \
                                            region_selections[region.name][protein1_name].residues.resids
                                        residue_contact_maps[pair_name][region_key]['resids_2'] = \
                                            region_selections[region.name][protein2_name].residues.resids
                                    else:
                                        residue_contact_maps[pair_name][region_key]['contact_count'] += contact_matrix
                                    
                                    residue_contact_maps[pair_name][region_key]['dimer_frames'] += 1
        
        except Exception as e:
            print(f"Error analyzing frame {frame_idx}: {str(e)}")
            continue
    
    # Calculate contact probabilities
    for pair_name in residue_contact_maps:
        for region_key in residue_contact_maps[pair_name]:
            if residue_contact_maps[pair_name][region_key]['contact_count'] is not None:
                dimer_frames = residue_contact_maps[pair_name][region_key]['dimer_frames']
                if dimer_frames > 0:
                    residue_contact_maps[pair_name][region_key]['contact_prob'] = \
                        residue_contact_maps[pair_name][region_key]['contact_count'] / dimer_frames
    
    return residue_contact_maps


def identify_key_residue_pairs(contact_prob_matrix: np.ndarray, 
                              resids_1: np.ndarray, resids_2: np.ndarray,
                              threshold: float = 0.5) -> List[Tuple]:
    """Identify key residue pairs with high contact probability"""
    if contact_prob_matrix is None or len(contact_prob_matrix) == 0:
        return []
    
    key_pairs = []
    
    for i in range(contact_prob_matrix.shape[0]):
        for j in range(contact_prob_matrix.shape[1]):
            if contact_prob_matrix[i, j] >= threshold:
                key_pairs.append((
                    resids_1[i],
                    resids_2[j],
                    contact_prob_matrix[i, j]
                ))
    
    # Sort by probability (descending)
    key_pairs.sort(key=lambda x: x[2], reverse=True)
    
    return key_pairs


def plot_residue_contact_heatmap(residue_contact_maps: Dict, config: AnalysisConfig):
    """Plot residue-level contact heatmaps"""
    import os
    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot for each protein pair and region
    for pair_name in residue_contact_maps:
        for region in config.binding_regions:
            region_key = f"{region.name}-{region.name}"
            
            if region_key not in residue_contact_maps[pair_name]:
                continue
            
            contact_prob = residue_contact_maps[pair_name][region_key]['contact_prob']
            
            if contact_prob is None or len(contact_prob) == 0:
                continue
            
            resids_1 = residue_contact_maps[pair_name][region_key]['resids_1']
            resids_2 = residue_contact_maps[pair_name][region_key]['resids_2']
            
            # Create figure
            plt.figure(figsize=(10, 8))
            
            # Create heatmap
            sns.heatmap(contact_prob, cmap='viridis', vmin=0, vmax=1,
                       xticklabels=resids_2, yticklabels=resids_1,
                       cbar_kws={'label': 'Contact Probability'})
            
            plt.title(f'{region.name} Residue Contact Map: {pair_name}')
            plt.xlabel(f'{pair_name.split("-")[1]} Residue ID')
            plt.ylabel(f'{pair_name.split("-")[0]} Residue ID')
            
            plt.tight_layout()
            
            # Save plot
            for fmt in config.figure_format:
                filename = os.path.join(output_dir, 
                                      f'residue_heatmap_{pair_name}_{region.name}.{fmt}')
                plt.savefig(filename, format=fmt)
                print(f"Saved residue heatmap: {filename}")
            
            plt.close()
            
            # Identify and print key residue pairs
            key_pairs = identify_key_residue_pairs(contact_prob, resids_1, resids_2)
            if key_pairs:
                print(f"\nKey residue pairs for {pair_name} - {region.name}:")
                for res1, res2, prob in key_pairs[:10]:  # Top 10
                    print(f"  Residue {res1} - Residue {res2}: {prob:.2f}")


def calculate_residue_importance_scores(residue_contact_maps: Dict, 
                                       config: AnalysisConfig) -> Dict:
    """Calculate importance scores for each residue based on contact patterns"""
    importance_scores = {}
    
    for region in config.binding_regions:
        region_key = f"{region.name}-{region.name}"
        scores = {}
        
        # Aggregate across all protein pairs
        for pair_name in residue_contact_maps:
            if region_key not in residue_contact_maps[pair_name]:
                continue
            
            contact_prob = residue_contact_maps[pair_name][region_key]['contact_prob']
            if contact_prob is None:
                continue
            
            resids_1 = residue_contact_maps[pair_name][region_key]['resids_1']
            resids_2 = residue_contact_maps[pair_name][region_key]['resids_2']
            
            # Score for protein 1 residues
            for i, resid in enumerate(resids_1):
                if resid not in scores:
                    scores[resid] = []
                scores[resid].append(np.mean(contact_prob[i, :]))
            
            # Score for protein 2 residues
            for j, resid in enumerate(resids_2):
                if resid not in scores:
                    scores[resid] = []
                scores[resid].append(np.mean(contact_prob[:, j]))
        
        # Calculate average scores
        importance_scores[region.name] = {}
        for resid, score_list in scores.items():
            importance_scores[region.name][resid] = np.mean(score_list)
    
    return importance_scores


def plot_residue_importance(importance_scores: Dict, config: AnalysisConfig):
    """Plot residue importance scores"""
    import os
    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    n_regions = len(config.binding_regions)
    if n_regions == 0:
        return
    
    fig, axes = plt.subplots(n_regions, 1, figsize=(12, 4 * n_regions))
    if n_regions == 1:
        axes = [axes]
    
    for i, region in enumerate(config.binding_regions):
        ax = axes[i]
        
        if region.name not in importance_scores or not importance_scores[region.name]:
            ax.text(0.5, 0.5, f"No data for {region.name}", 
                   transform=ax.transAxes, ha='center')
            continue
        
        # Sort residues by ID
        resids = sorted(importance_scores[region.name].keys())
        scores = [importance_scores[region.name][resid] for resid in resids]
        
        # Create bar plot
        bars = ax.bar(range(len(resids)), scores, color=region.color, alpha=0.7)
        
        # Highlight important residues
        threshold = np.mean(scores) + np.std(scores)
        for j, (resid, score) in enumerate(zip(resids, scores)):
            if score > threshold:
                bars[j].set_color(region.color)
                bars[j].set_alpha(1.0)
                ax.text(j, score, str(resid), ha='center', va='bottom', fontsize=8)
        
        ax.set_xlabel('Residue Index')
        ax.set_ylabel('Importance Score')
        ax.set_title(f'{region.name} - Residue Importance for Dimerization')
        ax.set_xticks(range(0, len(resids), max(1, len(resids)//20)))
        ax.set_xticklabels([resids[i] for i in range(0, len(resids), max(1, len(resids)//20))],
                          rotation=45)
        ax.axhline(y=threshold, color='r', linestyle='--', alpha=0.5, 
                  label=f'Threshold (mean + std)')
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    for fmt in config.figure_format:
        filename = os.path.join(output_dir, f'residue_importance.{fmt}')
        plt.savefig(filename, format=fmt)
        print(f"Saved residue importance plot: {filename}")
    
    plt.close()


def export_residue_data_for_visualization(residue_contact_maps: Dict,
                                         importance_scores: Dict,
                                         config: AnalysisConfig):
    """Export residue data for external visualization tools"""
    import os
    import pandas as pd
    
    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Export importance scores as B-factors for PyMOL
    for region in config.binding_regions:
        if region.name not in importance_scores:
            continue
        
        # Create PDB B-factor file
        bfactor_lines = []
        bfactor_lines.append(f"# B-factor values for {region.name}")
        bfactor_lines.append("# Load in PyMOL with: alter resi X, b=VALUE")
        
        for resid, score in importance_scores[region.name].items():
            # Scale score to B-factor range (0-100)
            bfactor = min(100, score * 100)
            bfactor_lines.append(f"alter resi {resid}, b={bfactor:.2f}")
        
        bfactor_file = os.path.join(output_dir, f'bfactors_{region.name}.pml')
        with open(bfactor_file, 'w') as f:
            f.write('\n'.join(bfactor_lines))
        print(f"Exported B-factors for {region.name}: {bfactor_file}")
    
    # Export contact matrices as CSV
    for pair_name in residue_contact_maps:
        for region in config.binding_regions:
            region_key = f"{region.name}-{region.name}"
            
            if region_key not in residue_contact_maps[pair_name]:
                continue
            
            contact_prob = residue_contact_maps[pair_name][region_key]['contact_prob']
            if contact_prob is None:
                continue
            
            resids_1 = residue_contact_maps[pair_name][region_key]['resids_1']
            resids_2 = residue_contact_maps[pair_name][region_key]['resids_2']
            
            # Create DataFrame
            df = pd.DataFrame(contact_prob, 
                            index=resids_1,
                            columns=resids_2)
            
            csv_file = os.path.join(output_dir, 
                                   f'contact_matrix_{pair_name}_{region.name}.csv')
            df.to_csv(csv_file)
            print(f"Exported contact matrix: {csv_file}")