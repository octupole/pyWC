# Pytim Scripts Installation Guide

## Summary

Three Willard-Chandler analysis scripts have been successfully installed as command-line tools in the pytim conda environment.

## Installed Scripts

### 1. `pytim-compare-wc-backends`
**Purpose:** Compare CPU vs GPU performance for Willard-Chandler density evaluation

**Source:** `scripts/compare_cpu_gpu.py`

**Usage:**
```bash
pytim-compare-wc-backends --help
pytim-compare-wc-backends --particles 10000 --grid 50000
pytim-compare-wc-backends --particles 5000 --chunk 4096 --max-mem 268435456
```

**Key Features:**
- Benchmarks CPU vs GPU (CuPy) backends
- Configurable particle count and grid size
- Adjustable GPU chunk size and memory limits
- Reports speedup and accuracy metrics

---

### 2. `pytim-wc-area`
**Purpose:** Compute membrane-solvent interface area and thickness using Willard-Chandler surfaces

**Source:** `scripts/compute_willard_chandler_area.py`

**Usage:**
```bash
pytim-wc-area --help
pytim-wc-area -s topology.tpr -x trajectory.xtc --selection "resname POPC"
pytim-wc-area -s system.gro -x traj.xtc --alpha 3.0 --mesh 2.0 --step 10
pytim-wc-area -s input.psf -x output.dcd --start-time 1000 --end-time 5000
```

**Key Features:**
- Analyzes trajectories frame-by-frame
- Computes triangulated surface area
- Generates thickness maps on XY grid
- Flexible atom selection (direct, file, or environment variable)
- Time window and frame stride controls
- Outputs:
  - Per-frame area statistics (CSV)
  - Averaged thickness map (CSV)
  - Surface roughness metrics (RMS, peak-to-valley)

**Parameters:**
- `--alpha`: Gaussian kernel width (default: 3.0 Å)
- `--mesh`: Grid spacing (default: 2.0 Å)
- `--density-level`: Relative density threshold 0-2 (default: 1.0)
- `--grid-size`: XY bins for thickness (default: 20)
- `--step`: Analyze every Nth frame (default: 1)

---

### 3. `pytim-bending-rigidity`
**Purpose:** Estimate membrane bending rigidity from Willard-Chandler surfaces via Helfrich fitting

**Source:** `scripts/compute_bending_rigidity.py`

**Usage:**
```bash
pytim-bending-rigidity --help
pytim-bending-rigidity -s topology.tpr -x trajectory.xtc --bilayer-selection "resname POPC"
pytim-bending-rigidity -s system.gro -x traj.xtc --leaflet-a "resid 1:128" --leaflet-b "resid 129:256"
pytim-bending-rigidity -s input.tpr -x output.xtc --temperature 310 --qmin 0.1 --qmax 0.5
```

**Key Features:**
- Builds separate WC surfaces for each leaflet
- Computes mid-surface height fluctuations
- Fourier analysis of surface undulations
- Helfrich fit: `kBT / <|h(q)|²> = σq² + κq⁴`
- Extracts bending rigidity (κ) and surface tension (σ)

**Workflow:**
1. Leaflet assignment (manual or automatic via LeafletFinder)
2. Per-frame WC surface construction
3. Height field rasterization on XY grid
4. 2D FFT to compute undulation spectrum
5. Radial averaging in q-space
6. Linear regression for κ and σ

**Parameters:**
- `--alpha`: Gaussian width (default: 3.0 Å)
- `--mesh`: WC grid spacing (default: 2.0 Å)
- `--grid-size`: XY bins for FFT (default: 32)
- `--temperature`: Simulation temperature (default: 310 K)
- `--qmin`, `--qmax`: Wave-vector range for fitting
- `--detrend-plane`: Remove best-fit plane before FFT

**Outputs:**
- Radially averaged spectrum: `wc_undulation_spectrum.csv`
- Mean thickness map: `wc_thickness_mean.csv`
- Optional 2D spectrum: `wc_modes_raw.npz`
- Fitted parameters: κ (bending rigidity), σ (surface tension)

---

## Installation Details

### Changes Made

1. **Created `scripts/__init__.py`**
   - Makes `scripts/` a Python package
   - Exports all three script modules

2. **Updated `pyproject.toml`**
   - Added `[project.scripts]` section with three entry points:
     ```toml
     [project.scripts]
     pytim-compare-wc-backends = "scripts.compare_cpu_gpu:main"
     pytim-wc-area = "scripts.compute_willard_chandler_area:cli"
     pytim-bending-rigidity = "scripts.compute_bending_rigidity:cli"
     ```
   - Modified `[tool.setuptools.packages.find]` to include `scripts*`

3. **Added CLI wrappers**
   - Added `cli()` function to `compute_willard_chandler_area.py`
   - Added `cli()` function to `compute_bending_rigidity.py`
   - These wrap the `main(args)` functions for proper entry point handling

### Installation Command

```bash
# Using the pytim-git conda environment
/opt/miniforge3/envs/pytim-git/bin/python -m pip install -e .
```

### Verification

All scripts successfully installed to:
```
/opt/miniforge3/envs/pytim-git/bin/pytim-compare-wc-backends
/opt/miniforge3/envs/pytim-git/bin/pytim-wc-area
/opt/miniforge3/envs/pytim-git/bin/pytim-bending-rigidity
```

---

## Script Comparison

| Feature | compare-wc-backends | wc-area | bending-rigidity |
|---------|---------------------|---------|------------------|
| **Purpose** | Performance testing | Surface area analysis | Mechanical properties |
| **Input** | Random particles | MD trajectory | MD trajectory |
| **Output** | Timing & speedup | Area & thickness CSV | κ, σ, spectrum CSV |
| **GPU Support** | Yes (CuPy) | No | No |
| **Trajectory** | Not required | Required | Required |
| **Complexity** | Simple benchmark | Moderate | Advanced |
| **Use Case** | Testing/optimization | Interface characterization | Membrane mechanics |

---

## Example Workflows

### Quick Performance Test
```bash
# Test CPU vs GPU with default settings
pytim-compare-wc-backends

# Larger system for better statistics
pytim-compare-wc-backends --particles 20000 --grid 100000
```

### Membrane Interface Analysis
```bash
# Analyze lipid bilayer interface
pytim-wc-area \
  -s membrane.tpr \
  -x trajectory.xtc \
  --selection "resname POPC POPE" \
  --alpha 3.0 \
  --mesh 1.5 \
  --step 5 \
  --grid-size 30 \
  -o area_results.csv \
  --thickness-output thickness_map.csv
```

### Bending Rigidity Calculation
```bash
# Calculate bending modulus for DPPC bilayer at 323K
pytim-bending-rigidity \
  -s dppc_bilayer.tpr \
  -x production.xtc \
  --bilayer-selection "resname DPPC" \
  --alpha 3.0 \
  --mesh 2.0 \
  --grid-size 32 \
  --temperature 323.0 \
  --qmin 0.05 \
  --qmax 0.4 \
  --stride 10 \
  --spectrum-output dppc_spectrum.csv
```

---

## Technical Notes

### Dependencies
All scripts require:
- MDAnalysis
- NumPy
- SciPy
- scikit-image
- pytim (with Willard-Chandler support)

GPU backend additionally requires:
- CuPy (with appropriate CUDA version)

### Performance Considerations

**pytim-compare-wc-backends:**
- GPU speedup typically 5-20x for large systems
- Memory usage scales with grid size
- Use `--chunk` and `--max-mem` to control GPU memory

**pytim-wc-area:**
- Time complexity: O(N_frames × N_atoms × grid³)
- Memory: density field + triangulated surface
- Typical runtime: ~1-10 s/frame depending on system size

**pytim-bending-rigidity:**
- Most computationally intensive script
- FFT complexity: O(N_frames × grid² log grid)
- Memory: accumulates full spectrum across trajectory
- Recommended: grid-size = 32-64 for good q-space resolution

### Troubleshooting

**Script not found:**
```bash
# Ensure conda environment is activated
conda activate pytim-git

# Or use full path
/opt/miniforge3/envs/pytim-git/bin/pytim-wc-area --help
```

**Import errors:**
```bash
# Reinstall in editable mode
cd /home/marchi/git/pytim
/opt/miniforge3/envs/pytim-git/bin/python -m pip install -e .
```

**GPU backend unavailable:**
```bash
# Check CuPy installation
python -c "import cupy; print(cupy.__version__)"

# Install if needed
conda install -c conda-forge cupy
```

---

## Files Modified

1. `pyproject.toml` - Added script entry points and package configuration
2. `scripts/__init__.py` - Created package initialization
3. `scripts/compute_willard_chandler_area.py` - Added `cli()` wrapper
4. `scripts/compute_bending_rigidity.py` - Added `cli()` wrapper

## Files Added

- `SCRIPTS_INSTALLATION.md` - This documentation file

---

## Future Enhancements

Potential improvements:
- Add progress bars to long-running analyses
- Parallel frame processing with multiprocessing/dask
- Interactive plotting options
- Configuration file support for complex workflows
- Integration with Jupyter notebooks via magic commands
- Automatic parameter optimization based on system properties
