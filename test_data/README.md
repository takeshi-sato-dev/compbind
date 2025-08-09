# Test Data

This directory contains small test datasets for unit testing.

## Files

- `test_system.psf`: Topology file (6.62 MB)
- `test_trajectory.xtc`: Trajectory with 50 frames (9.63 MB)

## Generation

These files were generated from the full trajectory using `create_test_data.py`.
- Original PSF: step5_assembly.psf
- Original XTC: step7_production.xtc
- Original trajectory: 100301 frames
- Test trajectory: frames 80000 to 84900 (step 100)

## Usage

```python
import MDAnalysis as mda
u = mda.Universe('test_data/test_system.psf', 'test_data/test_trajectory.xtc')
print(f'Loaded {len(u.atoms)} atoms, {len(u.trajectory)} frames')
```
