"""
Trajectory analysis functions for dimerization mechanism
"""

import numpy as np
from typing import Dict, List, Tuple
from tqdm import tqdm
from core import analyze_frame_with_time


def analyze_trajectory_detailed(universe, proteins: Dict, region_selections: Dict, 
                               config, sample_frames: int = 50) -> Tuple[List, List]:
    """Analyze trajectory with detailed sampling"""
    
    max_frame = len(universe.trajectory)
    effective_stop = min(config.stop_frame, max_frame)
    frames = list(range(config.start_frame, effective_stop, config.step_frame))
    
    if sample_frames and sample_frames < len(frames):
        sample_indices = np.linspace(0, len(frames)-1, sample_frames, dtype=int)
        sample_frames_list = [frames[i] for i in sample_indices]
    else:
        sample_frames_list = frames
    
    print(f"  Total frames to analyze: {len(frames)}")
    print(f"  Detailed analysis frames: {len(sample_frames_list)}")
    
    all_results = []
    detailed_results = []
    
    # Get region selections
    region_1_selections = {}
    region_2_selections = {}
    
    if config.binding_regions:
        if config.binding_regions[0].name in region_selections:
            region_1_selections = region_selections[config.binding_regions[0].name]
        if len(config.binding_regions) > 1 and config.binding_regions[1].name in region_selections:
            region_2_selections = region_selections[config.binding_regions[1].name]
    
    for frame_idx in tqdm(frames, desc="  Analyzing frames"):
        frame_results = analyze_frame_with_time(frame_idx, universe, proteins, 
                                               region_1_selections, region_2_selections)
        all_results.append(frame_results)
        
        if frame_idx in sample_frames_list:
            detailed_results.append(frame_results)
    
    return all_results, detailed_results


def collect_time_series_data(all_results: List[Dict]) -> Dict:
    """Collect time series data from results"""
    protein_pairs = set()
    for result in all_results:
        protein_pairs.update(result['dimers'].keys())
    
    time_series = {pair: {
        'times': [],
        'dimer_status': [],
        'region_1_contacts': [],
        'region_2_contacts': []
    } for pair in protein_pairs}
    
    for result in all_results:
        for pair in protein_pairs:
            if pair in result['dimers']:
                time_series[pair]['times'].append(result['time'])
                time_series[pair]['dimer_status'].append(1 if result['dimers'][pair] else 0)
                time_series[pair]['region_1_contacts'].append(
                    result.get('region_1_contacts', {}).get(pair, 0))
                time_series[pair]['region_2_contacts'].append(
                    result.get('region_2_contacts', {}).get(pair, 0))
    
    return time_series


def determine_dimerizing_pairs(time_series: Dict, threshold: float = 0.2) -> List[str]:
    """Identify pairs that form dimers"""
    dimerizing_pairs = []
    
    for pair, data in time_series.items():
        if len(data['dimer_status']) > 0:
            dimer_probability = sum(data['dimer_status']) / len(data['dimer_status'])
            if dimer_probability >= threshold:
                dimerizing_pairs.append(pair)
    
    return dimerizing_pairs