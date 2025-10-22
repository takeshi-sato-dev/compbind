# Contributing to Protein Dimerization Analysis

Thank you for your interest in contributing to this project! We welcome contributions from the community.

## How to Contribute

### Reporting Issues

If you find a bug or have a suggestion for improvement:

1. Check if the issue already exists in [GitHub Issues](https://github.com/yourusername/protein-dimer-analysis/issues)
2. If not, create a new issue with:
   - Clear description of the problem
   - Steps to reproduce (if it's a bug)
   - Expected vs actual behavior
   - Your environment (Python version, OS, etc.)

### Submitting Changes

1. **Fork the repository** and create your branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Add or update functionality
   - Follow the existing code style
   - Add docstrings to new functions
   - Update README.md if needed

3. **Test your changes**:
   ```bash
   python run_interactive.py
   # Ensure it runs without errors
   ```

4. **Commit your changes**:
   ```bash
   git commit -m "Add: brief description of changes"
   ```
   
   Commit message prefixes:
   - `Add:` for new features
   - `Fix:` for bug fixes
   - `Update:` for updates to existing features
   - `Docs:` for documentation only

5. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Submit a Pull Request** through GitHub

## Code Style Guidelines

- Follow PEP 8 Python style guide
- Use meaningful variable and function names
- Add type hints where appropriate
- Include docstrings for all functions:
  ```python
  def calculate_contacts(region1, region2, cutoff=6.0):
      """
      Calculate contacts between two regions.
      
      Parameters
      ----------
      region1 : MDAnalysis.AtomGroup
          First region selection
      region2 : MDAnalysis.AtomGroup
          Second region selection
      cutoff : float
          Distance cutoff in Angstroms
      
      Returns
      -------
      int
          Number of contacts
      """
  ```

## Adding New Features

When adding new analysis methods:

1. Add the function to the appropriate module:
   - `core.py` for basic analysis functions
   - `stats_analysis.py` for statistical methods
   - `visualization.py` for new plot types
   - `bayesian.py` for Bayesian methods

2. Update `run_interactive.py` if it's a major feature

3. Document the feature in README.md

## Questions?

Feel free to open an issue for any questions about contributing.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.