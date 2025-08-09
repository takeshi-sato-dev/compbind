"""
Configuration management for dimerization mechanism analysis
Supports multiple binding regions/motifs
"""

import json
import yaml
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import numpy as np


@dataclass
class BindingRegion:
    """Definition of a binding region/motif"""
    name: str
    resid_start: int
    resid_end: int
    description: Optional[str] = None
    color: Optional[str] = None  # For visualization
    
    def __post_init__(self):
        """Assign default color if not provided"""
        if self.color is None:
            # Default color palette
            default_colors = ['#2271B5', '#F58700', '#1A9641', '#E31A1C', 
                            '#6A3D9A', '#FF7F00', '#B2DF8A', '#FB9A99']
            # Use hash of name to get consistent color
            color_idx = hash(self.name) % len(default_colors)
            self.color = default_colors[color_idx]


@dataclass
class AnalysisConfig:
    """Configuration for dimerization analysis"""
    # Input files
    psf_file: str = 'step5_assembly.psf'
    xtc_file: str = 'md_wrapped.xtc'
    
    # Output settings
    output_dir: str = 'dimerization_mechanism'
    
    # Frame selection
    start_frame: int = 20000
    stop_frame: int = 80000
    step_frame: int = 5
    
    # Analysis parameters
    protein_contact_cutoff: float = 6.0  # Å
    dimer_cutoff: float = 20.0  # Å
    
    # Thermodynamic parameters
    temperature: float = 298.0  # K
    boltzmann_const: float = 0.001987 * 4.184  # kJ/(mol⋅K)
    
    # Binding regions (now supports multiple)
    binding_regions: List[BindingRegion] = field(default_factory=list)
    
    # Protein selection
    n_proteins: int = 4
    segids: Optional[List[str]] = None
    
    # Statistical parameters
    bootstrap_samples: int = 1000
    bayesian_samples: int = 2000
    bayesian_tune: int = 1000
    
    # Visualization settings
    dpi: int = 300
    figure_format: List[str] = field(default_factory=lambda: ['png', 'pdf'])
    
    @classmethod
    def from_yaml(cls, yaml_file: str) -> 'AnalysisConfig':
        """Load configuration from YAML file"""
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
        
        # Parse binding regions
        regions = []
        if 'binding_regions' in data:
            for region_data in data['binding_regions']:
                regions.append(BindingRegion(**region_data))
        
        # Remove binding_regions from data to avoid duplicate
        config_data = {k: v for k, v in data.items() if k != 'binding_regions'}
        
        return cls(binding_regions=regions, **config_data)
    
    @classmethod
    def from_json(cls, json_file: str) -> 'AnalysisConfig':
        """Load configuration from JSON file"""
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Parse binding regions
        regions = []
        if 'binding_regions' in data:
            for region_data in data['binding_regions']:
                regions.append(BindingRegion(**region_data))
        
        # Remove binding_regions from data to avoid duplicate
        config_data = {k: v for k, v in data.items() if k != 'binding_regions'}
        
        return cls(binding_regions=regions, **config_data)
    
    def to_yaml(self, yaml_file: str):
        """Save configuration to YAML file"""
        data = {
            'psf_file': self.psf_file,
            'xtc_file': self.xtc_file,
            'output_dir': self.output_dir,
            'start_frame': self.start_frame,
            'stop_frame': self.stop_frame,
            'step_frame': self.step_frame,
            'protein_contact_cutoff': self.protein_contact_cutoff,
            'dimer_cutoff': self.dimer_cutoff,
            'temperature': self.temperature,
            'n_proteins': self.n_proteins,
            'bootstrap_samples': self.bootstrap_samples,
            'bayesian_samples': self.bayesian_samples,
            'bayesian_tune': self.bayesian_tune,
            'binding_regions': [
                {
                    'name': r.name,
                    'resid_start': r.resid_start,
                    'resid_end': r.resid_end,
                    'description': r.description,
                    'color': r.color
                }
                for r in self.binding_regions
            ]
        }
        
        if self.segids:
            data['segids'] = self.segids
        
        with open(yaml_file, 'w') as f:
            yaml.dump(data, f, default_flow_style=False)
    
    def add_region(self, name: str, resid_start: int, resid_end: int, 
                   description: Optional[str] = None, color: Optional[str] = None):
        """Add a binding region"""
        region = BindingRegion(name, resid_start, resid_end, description, color)
        self.binding_regions.append(region)
    
    def get_region_pairs(self) -> List[Tuple[BindingRegion, BindingRegion]]:
        """Get all unique pairs of binding regions for comparison"""
        pairs = []
        n_regions = len(self.binding_regions)
        for i in range(n_regions):
            for j in range(i + 1, n_regions):
                pairs.append((self.binding_regions[i], self.binding_regions[j]))
        return pairs
    
    def validate(self) -> bool:
        """Validate configuration"""
        errors = []
        
        # Check files exist (optional, can be skipped if files not yet available)
        # import os
        # if not os.path.exists(self.psf_file):
        #     errors.append(f"PSF file not found: {self.psf_file}")
        # if not os.path.exists(self.xtc_file):
        #     errors.append(f"XTC file not found: {self.xtc_file}")
        
        # Check frame range
        if self.start_frame < 0:
            errors.append("start_frame must be >= 0")
        if self.stop_frame <= self.start_frame:
            errors.append("stop_frame must be > start_frame")
        if self.step_frame <= 0:
            errors.append("step_frame must be > 0")
        
        # Check cutoffs
        if self.protein_contact_cutoff <= 0:
            errors.append("protein_contact_cutoff must be > 0")
        if self.dimer_cutoff <= 0:
            errors.append("dimer_cutoff must be > 0")
        
        # Check regions
        if len(self.binding_regions) < 1:
            errors.append("At least one binding region must be defined")
        
        for region in self.binding_regions:
            if region.resid_start < 1:
                errors.append(f"Region {region.name}: resid_start must be >= 1")
            if region.resid_end < region.resid_start:
                errors.append(f"Region {region.name}: resid_end must be >= resid_start")
        
        # Check for overlapping regions (warning, not error)
        for i, r1 in enumerate(self.binding_regions):
            for j, r2 in enumerate(self.binding_regions[i+1:], i+1):
                if (r1.resid_start <= r2.resid_end and r2.resid_start <= r1.resid_end):
                    print(f"Warning: Regions {r1.name} and {r2.name} overlap")
        
        if errors:
            print("Configuration errors:")
            for error in errors:
                print(f"  - {error}")
            return False
        
        return True


def create_example_config() -> AnalysisConfig:
    """Create an example configuration with multiple motifs"""
    config = AnalysisConfig()
    
    # Add example binding regions
    config.add_region(
        name="GxxxG_motif",
        resid_start=65,
        resid_end=70,
        description="GxxxG dimerization motif",
        color="#2271B5"
    )
    
    config.add_region(
        name="Leu_zipper_N",
        resid_start=71,
        resid_end=77,
        description="N-terminal part of leucine zipper",
        color="#F58700"
    )
    
    config.add_region(
        name="Leu_zipper_C",
        resid_start=78,
        resid_end=84,
        description="C-terminal part of leucine zipper",
        color="#1A9641"
    )
    
    config.add_region(
        name="TM_C_term",
        resid_start=85,
        resid_end=91,
        description="C-terminal TM region",
        color="#E31A1C"
    )
    
    return config


def create_example_yaml(filename: str = "example_config.yaml"):
    """Create an example YAML configuration file"""
    config = create_example_config()
    config.to_yaml(filename)
    print(f"Example configuration saved to {filename}")
    return config