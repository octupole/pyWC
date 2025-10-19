# Command-Line Tools

pyWC provides production-ready command-line tools for membrane analysis. These tools require no Python programming and can be integrated into analysis pipelines or batch processing workflows.

## Overview

| Tool | Purpose | Output |
|------|---------|--------|
| `pywc-wc-area` | Compute surface area and thickness profiles | CSV files with per-frame statistics |
| `pywc-bending-rigidity` | Estimate bending modulus from fluctuations | Bending rigidity value |
| `pywc-compare-wc-backends` | Benchmark CPU vs GPU backends | Performance comparison |

## pywc-wc-area

The main tool for membrane surface area and thickness analysis across MD trajectories.

### Quick Start

```bash
pywc-wc-area -s membrane.tpr -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 --mesh 2.5 \
    --backend cpp \
    -b 0 -e 10000 --step 10
```

### Key Features

- **Dual outputs**: Surface area per frame + averaged thickness map
- **Multiple backends**: Python, C++ (CPU), or GPU (CuPy)
- **Time/frame selection**: Analyze specific portions of trajectories
- **Automatic centering**: Handle periodic boundaries correctly
- **Performance timing**: Built-in benchmarking

### Output Files

1. **Area file** (`willard_chandler_area.csv`): Per-frame surface area statistics
2. **Thickness file** (`willard_chandler_thickness.csv`): Spatially-resolved thickness map

### Command-Line Options

#### Required Arguments

```bash
-s, --topology FILE        # Topology file (TPR, PSF, GRO, etc.)
-x, --trajectory FILE      # Trajectory file (XTC, TRR, DCD, etc.)
```

#### Atom Selection

```bash
--selection "EXPR"         # MDAnalysis selection string
--selection-file FILE      # Read selection from text file (recommended)
--selection-env VAR        # Read selection from environment variable
```

**Example selection file** (`atoms.txt`):
```
resname DPPC DOPC and name PO4 C31 C32
```

#### Willard-Chandler Parameters

```bash
--alpha FLOAT              # Gaussian kernel width (Å) [default: 3.0]
--mesh FLOAT               # Grid spacing (Å) [default: 2.0]
--density-cutoff FLOAT     # Absolute density threshold (overrides --density-level)
--density-level FLOAT      # Relative density level [0-2] [default: 1.0]
                          #   0 = minimum density
                          #   1 = midpoint (typical)
                          #   2 = maximum density
```

!!! tip "Parameter Selection"
    - **`alpha`**: Use 2.5-3.5 Å for lipid membranes, 3.0-4.0 Å for protein surfaces
    - **`mesh`**: Use 2.0-2.5 Å for most applications (balance between accuracy and speed)
    - **`density-level`**: Use 1.0 for typical liquid/vapor interfaces

#### Backend Selection

```bash
--backend {cpu,cpp,cupy,gpu,python}
                          # Surface computation backend
                          #   cpp/cpu: C++ with OpenMP (recommended)
                          #   cupy/gpu: GPU acceleration (requires CuPy)
                          #   python: Pure Python (testing only)
```

#### Trajectory Control

```bash
--start INDEX             # First frame index (inclusive)
--stop INDEX              # Last frame index (exclusive)
-b, --start-time TIME     # Start time in picoseconds
-e, --end-time TIME       # End time in picoseconds
--step N                  # Analyze every Nth frame [default: 1]
```

#### Thickness Analysis

```bash
--grid-size N             # XY grid bins for thickness map [default: 20]
                          # Set to 0 to disable thickness calculation
--center                  # Center analysis group before evaluation
                          # RECOMMENDED for membrane systems
```

#### Output Control

```bash
-o, --output FILE         # CSV file for area data
                          # [default: willard_chandler_area.csv]
--thickness-output FILE   # CSV file for thickness map
                          # [default: willard_chandler_thickness.csv]
```

#### Performance Options

```bash
--enable-timing           # Enable detailed timing
--print-timing            # Print timing breakdown at end
--skip-frames N           # Skip first N frames in timing stats [default: 2]
--no-autoassign           # Disable automatic surface recomputation
```

### Example 1: Basic Membrane Analysis

Analyze a DPPC membrane over the entire trajectory:

```bash
# Create selection file
echo "resname DPPC and name PO4" > dppc_headgroups.txt

# Run analysis with C++ backend
pywc-wc-area \
    -s dppc_membrane.tpr \
    -x production.xtc \
    --selection-file dppc_headgroups.txt \
    --alpha 3.0 \
    --mesh 2.5 \
    --backend cpp \
    --center \
    -o dppc_area.csv \
    --thickness-output dppc_thickness.csv
```

**Output:**
```
frame      1 time       10.0 ps → area    15234.45 Ų (plane    16000.00 Ų)
frame      2 time       20.0 ps → area    15198.32 Ų (plane    16000.00 Ų)
...
Total time: 25.34 ms (average per frame, excluding first 2 frame(s))

Completed. Area statistics (Ų): mean=15216.78, std=45.23, min=15098.12, max=15334.56
Per-frame data saved to dppc_area.csv
Thickness grid written to dppc_thickness.csv (mean thickness 38.4 Å)
Upper surface roughness: mean RMS 2.15 Å, mean P-V 8.73 Å
Lower surface roughness: mean RMS 2.08 Å, mean P-V 8.21 Å
```

### Example 2: Time Window Selection

Analyze only the equilibrated portion (50-100 ns):

```bash
pywc-wc-area \
    -s system.tpr \
    -x trajectory.xtc \
    --selection-file selection.txt \
    --alpha 3.0 \
    --mesh 2.0 \
    -b 50000 \
    -e 100000 \
    --step 10 \
    --backend cpp
```

### Example 3: GPU Acceleration for Large Systems

Use GPU backend for a system with >100k atoms:

```bash
pywc-wc-area \
    -s large_membrane.tpr \
    -x trajectory.xtc \
    --selection-file lipids.txt \
    --alpha 3.0 \
    --mesh 2.5 \
    --backend cupy \
    --center \
    --enable-timing \
    --print-timing
```

### Example 4: Custom Density Threshold

Use a specific density cutoff instead of relative level:

```bash
pywc-wc-area \
    -s system.tpr \
    -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 \
    --mesh 2.0 \
    --density-cutoff 0.015 \
    --backend cpp
```

### Example 5: High-Resolution Thickness Map

Generate a detailed thickness map with fine spatial resolution:

```bash
pywc-wc-area \
    -s membrane.tpr \
    -x trajectory.xtc \
    --selection-file headgroups.txt \
    --alpha 3.0 \
    --mesh 2.0 \
    --grid-size 50 \
    --center \
    --backend cpp \
    --thickness-output thickness_highres.csv
```

### Understanding the Outputs

#### Area Output CSV

The area CSV file contains per-frame statistics:

| Column | Description |
|--------|-------------|
| `frame` | Frame number |
| `time_ps` | Simulation time (ps) |
| `area_angstrom2` | Willard-Chandler surface area (Ų) |
| `box_x`, `box_y`, `box_z` | Simulation box dimensions (Å) |
| `plane_area` | Projected XY plane area (Ų) |
| `density_cutoff` | Density threshold used for isosurface |
| `alpha` | Gaussian width parameter |
| `mesh` | Grid spacing used |
| `upper_rms`, `lower_rms` | RMS roughness of upper/lower surfaces |
| `upper_peak_to_valley`, `lower_peak_to_valley` | Peak-to-valley roughness |

**Example:**
```csv
frame,time_ps,area_angstrom2,box_x,box_y,box_z,plane_area,density_cutoff,alpha,mesh,upper_rms,lower_rms,upper_peak_to_valley,lower_peak_to_valley
1,10.0,15234.45,80.0,80.0,100.0,6400.0,0.0145,3.0,2.5,2.15,2.08,8.73,8.21
2,20.0,15198.32,80.1,79.9,100.1,6399.99,0.0145,3.0,2.5,2.18,2.11,8.91,8.45
```

#### Thickness Output CSV

The thickness CSV file contains a spatially-resolved thickness map:

| Column | Description |
|--------|-------------|
| `x` | X coordinate of grid cell center (Å) |
| `y` | Y coordinate of grid cell center (Å) |
| `thickness` | Local membrane thickness (Å) |

**Example:**
```csv
x,y,thickness
2.0,2.0,38.42
2.0,6.0,38.67
2.0,10.0,37.89
...
```

This can be visualized as a heatmap or contour plot:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load thickness map
df = pd.read_csv('willard_chandler_thickness.csv')

# Pivot for heatmap
grid = df.pivot(index='x', columns='y', values='thickness')

# Plot
plt.figure(figsize=(10, 8))
plt.imshow(grid, cmap='viridis', origin='lower', aspect='auto')
plt.colorbar(label='Thickness (Å)')
plt.xlabel('Y (Å)')
plt.ylabel('X (Å)')
plt.title('Membrane Thickness Map')
plt.savefig('thickness_map.png', dpi=300)
```

### Performance Considerations

#### Backend Selection

For a 105,000 atom GOLO/GOLH membrane system:

| Backend | Time per frame | Speedup |
|---------|----------------|---------|
| Python | 427 ms | 1.0x |
| C++ (CPU) | 25 ms | **17.1x** |
| GPU (CuPy) | 20 ms | **21.4x** |

**What is being measured:**

The timing captures the complete Willard-Chandler surface computation per frame, including:

1. **`center()`**: System centering in the unit cell
2. **`prepare_box()`**: Density field grid preparation
3. **`define_cluster_group()`**: KDE computation - evaluating Gaussian kernels at all grid points
4. **`compute_surface()`**: Marching cubes isosurface extraction and area calculation

**Performance bottleneck varies by backend:**
- **Python backend**: `compute_surface` dominates (~97% of time) - the marching cubes and KDE steps are slow
- **C++ and GPU backends**: `center` becomes the bottleneck (~38-46% of time) - surface computation is heavily accelerated, so centering dominates

These timings are automatically collected by the script when using the built-in timing system:

```python
# From compute_willard_chandler_area.py, lines 354-364
interface = pywc.WillardChandler(
    universe,
    group=analysis_group,
    alpha=args.alpha,
    mesh=args.mesh,
    centered=args.center,
    density_cutoff=args.density_cutoff,
    autoassign=args.no_autoassign,
    surface_backend=args.backend,
    centering_backend='cpu',
    enable_timing=True,  # Enables automatic timing
)
```

The reported time is the **mean per-frame computation time** (excluding the first 2 frames to remove initialization overhead):

```python
# From compute_willard_chandler_area.py, lines 518-526
timings = interface.get_detailed_timings(skip_frames=args.skip_frames)
if timings and 'total' in timings:
    stats = timings['total']
    if args.skip_frames > 0:
        print(f"Total time: {stats['mean']*1000:.2f} ms "
              f"(average per frame, excluding first {args.skip_frames} frame(s))")
```

**Recommendations:**
- Use `--backend cpp` for most applications
- Use `--backend cupy` for systems >100k atoms with NVIDIA GPU
- Use `--backend python` only for testing or debugging

#### Frame Stepping

For long trajectories, use `--step` to sample every Nth frame:

```bash
# Analyze every 10th frame (saves 90% computation time)
pywc-wc-area ... --step 10
```

This is often sufficient for converged statistics while dramatically reducing runtime.

#### Grid Size vs Performance

The `--grid-size` parameter affects thickness calculation:

| Grid Size | Spatial Resolution | Performance Impact |
|-----------|-------------------|-------------------|
| 10 | ~8 Å bins | Fast |
| 20 | ~4 Å bins | Standard (recommended) |
| 50 | ~1.6 Å bins | Slower, high detail |

### The `--center` Option

!!! warning "Important for Membranes"
    Always use `--center` for membrane systems to avoid artifacts from periodic boundaries.

The `--center` flag centers the analysis group in the unit cell before computing the density field. This prevents the membrane from being split across periodic boundaries, which would create artificial discontinuities in the surface.

**Without `--center`:**
```
Membrane wraps around box edge → Density field split → Incorrect surface
```

**With `--center`:**
```
Membrane centered in box → Continuous density field → Correct surface
```

### Troubleshooting

#### Empty surface / No triangles

**Problem:** `area_angstrom2` is 0.0 or NaN

**Solutions:**
- Check that your selection contains atoms: `--selection "resname DPPC" | wc -l`
- Adjust `--density-level` (try 0.5 or 1.5)
- Verify trajectory is not corrupted
- Check that `--alpha` and `--mesh` are reasonable

#### Memory errors with GPU backend

**Problem:** `cupy.cuda.memory.OutOfMemoryError`

**Solutions:**
- Use smaller `--mesh` value (e.g., 3.0 instead of 2.0)
- Use `--backend cpp` instead
- Reduce selection to fewer atoms

#### Surface area much larger than box area

**Problem:** `area_angstrom2 >> plane_area`

**Solutions:**
- Use `--center` flag
- Check that selection doesn't include bulk solvent
- Verify `--density-level` is appropriate (try 1.0)

---

## pywc-bending-rigidity

Estimate membrane bending modulus from surface fluctuations.

### Basic Usage

```bash
pywc-bending-rigidity \
    -s membrane.tpr \
    -x trajectory.xtc \
    --selection-file lipids.txt \
    --alpha 3.0 \
    --mesh 2.5
```

!!! note "Coming Soon"
    Detailed documentation for this tool is under development. See `pywc-bending-rigidity --help` for current options.

---

## pywc-compare-wc-backends

Benchmark computational backends for your system.

### Basic Usage

```bash
pywc-compare-wc-backends --help
```

!!! note "Coming Soon"
    Detailed documentation for this tool is under development. See `pywc-compare-wc-backends --help` for current options.

---

## Integration with Analysis Pipelines

### Bash Script Example

```bash
#!/bin/bash
# Batch process multiple membrane simulations

for system in DPPC DOPC POPC; do
    echo "Processing ${system}..."

    pywc-wc-area \
        -s ${system}/system.tpr \
        -x ${system}/production.xtc \
        --selection-file ${system}/selection.txt \
        --alpha 3.0 \
        --mesh 2.5 \
        --backend cpp \
        --center \
        -b 50000 -e 100000 \
        -o ${system}_area.csv \
        --thickness-output ${system}_thickness.csv
done
```

### Python Wrapper Example

```python
import subprocess
import pandas as pd

systems = ['DPPC', 'DOPC', 'POPC']
results = {}

for system in systems:
    # Run analysis
    subprocess.run([
        'pywc-wc-area',
        '-s', f'{system}/system.tpr',
        '-x', f'{system}/production.xtc',
        '--selection-file', f'{system}/selection.txt',
        '--alpha', '3.0',
        '--mesh', '2.5',
        '--backend', 'cpp',
        '-o', f'{system}_area.csv'
    ])

    # Load results
    df = pd.read_csv(f'{system}_area.csv')
    results[system] = {
        'mean_area': df['area_angstrom2'].mean(),
        'std_area': df['area_angstrom2'].std()
    }

# Compare
for system, stats in results.items():
    print(f"{system}: {stats['mean_area']:.1f} ± {stats['std_area']:.1f} Ų")
```

---

## Next Steps

- Learn about the [Python API](api/willard_chandler.md) for custom analysis
- Explore [Examples](examples.md) for common use cases
- Understand [Backend Performance](api/backends.md) tradeoffs
