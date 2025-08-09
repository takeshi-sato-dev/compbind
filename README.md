# CompBind - Competitive Binding Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![Platform](https://img.shields.io/badge/platform-macOS%20(Apple%20Silicon)-lightgrey)](https://www.apple.com/)

A Python package for analyzing competitive binding between different protein regions from molecular dynamics trajectories, optimized for Apple Silicon Macs.

## ⚠️ Platform Requirements

**This package is currently optimized for Apple Silicon (M1/M2/M3) Macs running macOS 14+**

The Bayesian analysis module requires specific versions of PyMC3 and Theano that are compatible with ARM64 architecture. Users on Intel Macs or Linux systems may need to adapt the installation process.

## Features

- 🎯 **Competitive binding analysis** between multiple protein regions
- 🔍 **Automatic protein detection** from PSF topology files
- ⚖️ **Differential binding assessment** to quantify competitive advantages
- 📊 **Comprehensive statistics** including bootstrap and Bayesian inference
- 📈 **Publication-ready visualizations** for competitive binding landscapes
- ⚡ **Optimized for Apple Silicon** with native ARM64 support
- 📝 **Flexible configuration** via YAML/JSON or interactive interface

## Installation

### Prerequisites

- macOS 14+ on Apple Silicon (M1/M2/M3)
- Python 3.10.14 (tested version)
- Xcode Command Line Tools

### Step 1: Create Python Environment

Using pyenv (recommended):

    pyenv install 3.10.14
    pyenv virtualenv 3.10.14 bayesian_env
    pyenv activate bayesian_env

Or using venv:

    python3.10 -m venv bayesian_env
    source bayesian_env/bin/activate

### Step 2: Install Dependencies

    # Upgrade pip
    pip install --upgrade pip

    # Install core packages
    pip install numpy==1.22.4
    pip install scipy==1.7.3
    pip install MDAnalysis==2.7.0
    pip install matplotlib==3.8.4
    pip install seaborn pandas pyyaml tqdm

    # Install Bayesian packages (specific versions for Apple Silicon)
    pip install pymc3==3.11.6
    pip install theano-pymc==1.1.2
    pip install arviz==0.12.1

### Step 3: Install CompBind

    git clone https://github.com/yourusername/compbind.git
    cd compbind
    pip install -e .

### Step 4: Verify Installation

    python test_environment.py

You should see all packages marked with ✓. The BLAS warning about architecture mismatch is expected and doesn't affect functionality.

## Known Issues

### BLAS Architecture Warning

You may see:

    ld: warning: ignoring file '/opt/local/lib/libopenblas-r1.dylib': found architecture 'x86_64', required architecture 'arm64'

This warning is harmless and the analysis will run correctly.

### PyMC3 Version Warning

The package uses PyMC3 (v3.11.6) instead of the newer PyMC (v5) for stability on Apple Silicon. The warning about outdated version can be ignored.

## Quick Start

### Interactive Mode

Always activate the environment first:

    source bayesian_env/bin/activate  # or pyenv activate bayesian_env
    python run_interactive.py

Follow the prompts to:
1. Specify your PSF and trajectory files
2. Define competing binding regions
3. Run the competitive binding analysis pipeline

### Configuration File Mode

Create a configuration file (config.yaml):

    psf_file: path/to/your/system.psf
    xtc_file: path/to/your/trajectory.xtc
    output_dir: results
    start_frame: 0
    stop_frame: 10000
    step_frame: 10

    # Define competing binding regions
    binding_regions:
      - name: TM_N_term
        resid_start: 65
        resid_end: 75
        description: N-terminal TM region
        color: '#2271B5'
      - name: TM_C_term
        resid_start: 76
        resid_end: 91
        description: C-terminal TM region
        color: '#F58700'

Run analysis:

    from config import AnalysisConfig
    from run_interactive import run_full_analysis

    config = AnalysisConfig.from_yaml('config.yaml')
    run_full_analysis(config)

## Output Files

### Data Files

- `contact_free_energy_analysis.csv` - Competitive binding statistics
- `bayesian_results.pkl` - Bayesian analysis of binding preferences
- `time_series_data.pkl` - Time-dependent competitive binding data
- `analysis_parameters.json` - Analysis configuration

### Visualizations

- `time_series_*.png/svg` - Competitive binding dynamics over time
- `free_energy_comparison.png/svg` - Energy comparisons between competing regions
- `bayesian_energy_distribution_*.png/svg` - Posterior distributions of binding preferences
- `bayesian_dimerization_energy.png/svg` - Overall binding energy summary
- `bayesian_energy_comparison_*.png/svg` - Competitive region comparison
- `residue_heatmap_*.png/svg` - Residue-level contact maps for each region
- `residue_importance.png/svg` - Key residues in competitive binding

## Module Structure

    compbind/
    ├── config.py              # Configuration management
    ├── core.py                # Core analysis functions
    ├── trajectory_analysis.py # Trajectory processing
    ├── stats_analysis.py      # Statistical calculations
    ├── bayesian.py           # Bayesian inference
    ├── visualization.py      # Plotting functions
    ├── residue_analysis.py   # Residue-level analysis
    ├── io_utils.py          # Input/output utilities
    └── run_interactive.py   # Interactive interface

## Example: Analyzing specific regions

    from config import AnalysisConfig

    config = AnalysisConfig()
    config.add_region("GxxxG_motif", 65, 70, "GxxxG dimerization motif")
    config.add_region("Leucine_zipper", 71, 84, "Leucine zipper region")
    config.add_region("C_terminal", 85, 91, "C-terminal region")

## Example: Custom analysis pipeline

    from io_utils import load_universe
    from trajectory_analysis import analyze_trajectory_detailed
    from stats_analysis import process_results_with_statistics

    universe = load_universe(psf_file, xtc_file)
    results = analyze_trajectory_detailed(universe, proteins, regions, config)
    stats = process_results_with_statistics(results, config)

## Troubleshooting

### ImportError with NumPy attributes

If you see errors about `np.bool`, `np.int`, etc., make sure you're using the correct environment:

    source bayesian_env/bin/activate

### Multiprocessing hangs

The Bayesian sampling uses multicore processing. If it hangs, you may need to:

    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1

### For Intel Mac or Linux Users

You'll need to create a different environment. Try:

    # Create environment with conda
    conda create -n compbind python=3.10
    conda activate compbind
    conda install -c conda-forge pymc arviz mdanalysis numpy scipy matplotlib seaborn pandas

## Citation

If you use CompBind in your research, please cite:

    @article{compbind2024,
      title={CompBind: A Python package for competitive binding analysis of protein interactions from molecular dynamics trajectories},
      author={Takeshi Sato},
      journal={Journal of Open Source Software},
      year={2024},
      doi={10.21105/joss.xxxxx}
    }

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Contact

- Issues: [GitHub Issues](https://github.com/takeshi-sato-dev/compbind/issues)
- Email: your.email@institution.edu

## Acknowledgments

- MDAnalysis community for trajectory analysis tools
- PyMC developers for Bayesian inference framework
- Apple Silicon community for ARM64 compatibility solutions