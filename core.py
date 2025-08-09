"""
Core analysis functions for dimerization mechanism analysis.
"""

import numpy as np
import MDAnalysis as mda
import traceback


def select_proteins_and_regions(universe, config):
    """Select protein chains and define regions for analysis."""
    
    # Select all proteins
    proteins = {}
    
    # Get unique segments/chains
    unique_segments = np.unique(universe.atoms.segids)
    
    print(f"  Segments in system: {unique_segments}")
    
    # Filter for protein segments - USE ACTUAL SEGMENT NAMES
    for segid in unique_segments:
        seg = universe.select_atoms(f"segid {segid}")
        # Check if this is a protein (has protein atoms)
        protein_atoms = seg.select_atoms("protein")
        if len(protein_atoms) > 0:
            # Use the actual segment ID (PROA, PROB, etc.)
            proteins[segid] = seg
            print(f"    {segid}: {len(seg.residues)} residues, {len(protein_atoms)} atoms (protein)")
        else:
            non_protein_atoms = len(seg.atoms)
            print(f"    {segid}: {non_protein_atoms} atoms (non-protein)")
    
    if not proteins:
        print("  No proteins found in the system!")
        return {}, {}, {}
    
    print(f"  Found {len(proteins)} protein chains: {sorted(proteins.keys())}")
    
    # Select regions for each protein
    region_1_selections = {}
    region_2_selections = {}
    
    for prot_name, protein in proteins.items():
        # Region 1
        try:
            region_1 = protein.select_atoms(
                f"resid {config.REGION_1_RESID_START}:{config.REGION_1_RESID_END}"
            )
            if len(region_1) > 0:
                region_1_selections[prot_name] = region_1
                print(f"    {prot_name} {config.REGION_1_NAME}: {len(region_1)} atoms")
        except:
            pass
        
        # Region 2
        try:
            region_2 = protein.select_atoms(
                f"resid {config.REGION_2_RESID_START}:{config.REGION_2_RESID_END}"
            )
            if len(region_2) > 0:
                region_2_selections[prot_name] = region_2
                print(f"    {prot_name} {config.REGION_2_NAME}: {len(region_2)} atoms")
        except:
            pass
    
    return proteins, region_1_selections, region_2_selections


def is_dimer(protein1, protein2, box, cutoff=20.0):
    """Check if two proteins form a dimer based on COM distance."""
    if len(protein1) == 0 or len(protein2) == 0:
        return False
    
    try:
        com1 = protein1.center_of_mass()
        com2 = protein2.center_of_mass()
        
        diff = com1 - com2
        for dim in range(3):
            if diff[dim] > box[dim] * 0.5:
                diff[dim] -= box[dim]
            elif diff[dim] < -box[dim] * 0.5:
                diff[dim] += box[dim]
        
        dist = np.sqrt(np.sum(diff * diff))
        return dist <= cutoff
    except Exception as e:
        print(f"Error during dimer check: {str(e)}")
        return False


def calculate_region_contacts(region1, region2, box, cutoff=6.0):
    """Calculate number of contacts between two regions."""
    if len(region1) == 0 or len(region2) == 0:
        return 0
    
    try:
        r1_com = region1.center_of_mass()
        r2_com = region2.center_of_mass()
        
        diff = r1_com - r2_com
        for dim in range(3):
            if diff[dim] > box[dim] * 0.5:
                diff[dim] -= box[dim]
            elif diff[dim] < -box[dim] * 0.5:
                diff[dim] += box[dim]
        
        com_dist = np.sqrt(np.sum(diff * diff))
        
        if com_dist > 30.0:
            return 0
        
        contact_count = 0
        
        for res1 in region1.residues:
            for res2 in region2.residues:
                min_dist = float('inf')
                
                if len(res1.atoms) == 0 or len(res2.atoms) == 0:
                    continue
                
                for atom1 in res1.atoms[:5]:  # Limit for speed
                    for atom2 in res2.atoms[:5]:
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
                
                if min_dist <= cutoff:
                    contact_count += 1
        
        return contact_count
    
    except Exception as e:
        print(f"Error during contact calculation: {str(e)}")
        return 0


def analyze_frame_with_time(frame_idx, universe, proteins, region_1_selections, region_2_selections):
    """Analyze a single frame for dimers and region contacts."""
    try:
        universe.trajectory[frame_idx]
        time = universe.trajectory[frame_idx].time
        box = universe.dimensions[:3]
        
        results = {
            'frame': frame_idx,
            'time': time,
            'dimers': {},
            'region_1_contacts': {},
            'region_2_contacts': {}
        }
        
        protein_names = list(proteins.keys())
        
        for i in range(len(protein_names)):
            for j in range(i + 1, len(protein_names)):
                protein1_name = protein_names[i]
                protein2_name = protein_names[j]
                
                pair_name = f"{protein1_name}-{protein2_name}"
                
                dimer_status = is_dimer(proteins[protein1_name], proteins[protein2_name], box)
                results['dimers'][pair_name] = dimer_status
                
                region_1_contacts = 0
                region_2_contacts = 0
                
                if dimer_status:
                    if protein1_name in region_1_selections and protein2_name in region_1_selections:
                        region_1_contacts = calculate_region_contacts(
                            region_1_selections[protein1_name], 
                            region_1_selections[protein2_name], 
                            box
                        )
                    
                    if protein1_name in region_2_selections and protein2_name in region_2_selections:
                        region_2_contacts = calculate_region_contacts(
                            region_2_selections[protein1_name], 
                            region_2_selections[protein2_name], 
                            box
                        )
                
                results['region_1_contacts'][pair_name] = region_1_contacts
                results['region_2_contacts'][pair_name] = region_2_contacts
        
        return results
    
    except Exception as e:
        print(f"Error analyzing frame {frame_idx}: {str(e)}")
        return {'frame': frame_idx, 'time': 0, 'dimers': {}, 
                'region_1_contacts': {}, 'region_2_contacts': {}}