# compute_willard_chandler_area.py

The `compute_willard_chandler_area.py` script (installed as `pywc-wc-area`) is the main command-line tool for analyzing membrane surface area and thickness profiles across MD trajectories.

## Overview

This script computes the Willard-Chandler dividing surface for each frame of a trajectory and outputs:

1. **Surface area per frame** - CSV file with area, roughness, and system parameters
2. **Thickness map** - Spatially-resolved membrane thickness averaged over the trajectory

## Basic Usage

```bash
pywc-wc-area -s topology.tpr -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 \
    --mesh 2.5 \
    --backend cpp
```

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-s, --topology FILE` | Topology file (TPR, PSF, GRO, PDB, etc.) |
| `-x, --trajectory FILE` | Trajectory file (XTC, TRR, DCD, etc.) |

## Atom Selection

You have three ways to specify which atoms to use for surface calculation:

### Method 1: Selection File (Recommended)

```bash
pywc-wc-area -s system.tpr -x trajectory.xtc \
    --selection-file selection.txt
```

**selection.txt:**
```
resname DPPC DOPC and name PO4 C31 C32
```

### Method 2: Command Line

```bash
pywc-wc-area -s system.tpr -x trajectory.xtc \
    --selection "resname DPPC and name PO4"
```

### Method 3: Environment Variable

```bash
export WC_SELECTION="resname DPPC and name PO4"
pywc-wc-area -s system.tpr -x trajectory.xtc \
    --selection-env WC_SELECTION
```

## Willard-Chandler Parameters

### Alpha (Gaussian Width)

```bash
--alpha 3.0  # Width of Gaussian kernel in Ångströms
```

**Guidelines:**
- Lipid membranes: 2.5-3.5 Å
- Protein surfaces: 3.0-4.0 Å
- Smaller values → sharper, noisier surfaces
- Larger values → smoother surfaces

### Mesh (Grid Spacing)

```bash
--mesh 2.5  # Grid spacing in Ångströms
```

**Guidelines:**
- Standard: 2.0-2.5 Å (good balance)
- High detail: 1.5-2.0 Å (slower)
- Fast preview: 3.0-4.0 Å

### Density Cutoff

**Option 1: Relative level (default)**
```bash
--density-level 1.0  # Midpoint between min and max density
```
- `0.0` = minimum density
- `1.0` = midpoint (typical for interfaces)
- `2.0` = maximum density

**Option 2: Absolute cutoff**
```bash
--density-cutoff 0.015  # Specific density value
```

## Backend Selection

```bash
--backend {cpp,cpu,cupy,gpu,python}
```

| Backend | When to Use | Performance |
|---------|-------------|-------------|
| `cpp` or `cpu` | **Default** - most systems | 17x faster than Python |
| `cupy` or `gpu` | Large systems (>100k atoms) with NVIDIA GPU | 21x faster than Python |
| `python` | Testing, debugging, or no C++ compiler | Baseline |

## Trajectory Control

### Frame Range

```bash
--start 100       # Start at frame 100
--stop 500        # Stop before frame 500
--step 10         # Analyze every 10th frame
```

### Time Window

```bash
-b 50000          # Start time in ps (50 ns)
-e 100000         # End time in ps (100 ns)
--step 10         # Analyze every 10th frame
```

## Centering (Critical for Membranes!)

```bash
--center          # Center the analysis group before each frame
```

**Always use `--center` for membrane systems!**

Without centering, the membrane can wrap across periodic boundaries, creating artifacts in the density field and incorrect surfaces.

## Thickness Analysis

### Grid Size

```bash
--grid-size 20    # 20x20 grid for thickness map (default)
```

- `10` → ~8 Å spatial resolution (fast)
- `20` → ~4 Å spatial resolution (recommended)
- `50` → ~1.6 Å spatial resolution (detailed, slower)
- `0` → Disable thickness calculation

## Output Files

### Area Output

```bash
-o willard_chandler_area.csv  # Per-frame statistics
```

**Columns:**
- `frame` - Frame number
- `time_ps` - Simulation time (ps)
- `area_angstrom2` - Surface area (Å^2)
- `box_x, box_y, box_z` - Box dimensions (Å)
- `plane_area` - Projected XY area (Å^2)
- `density_cutoff` - Density threshold used
- `alpha, mesh` - Parameters used
- `upper_rms, lower_rms` - Surface roughness
- `upper_peak_to_valley, lower_peak_to_valley` - Peak-to-valley roughness

### Thickness Output

```bash
--thickness-output willard_chandler_thickness.csv  # Thickness map
```

**Columns:**
- `x, y` - Grid cell coordinates (Å)
- `thickness` - Local membrane thickness (Å)

## Performance Options

```bash
--enable-timing   # Enable detailed timing
--print-timing    # Print full timing breakdown at end
--skip-frames 2   # Skip first N frames in timing stats (default: 2)
```

## Complete Examples

### Example 1: Basic DPPC Membrane Analysis

```bash
# Create selection file
cat > dppc_headgroups.txt << EOF
resname DPPC and name PO4
EOF

# Run analysis
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

**Expected output:**
```
frame      1 time       10.0 ps → area    15234.45 Å^2 (plane    16000.00 Å^2)
frame      2 time       20.0 ps → area    15198.32 Å^2 (plane    16000.00 Å^2)
...
Total time: 25.34 ms (average per frame, excluding first 2 frame(s))

Completed. Area statistics (Å^2): mean=15216.78, std=45.23, min=15098.12, max=15334.56
Per-frame data saved to dppc_area.csv
Thickness grid written to dppc_thickness.csv (mean thickness 38.4 Å)
Upper surface roughness: mean RMS 2.15 Å, mean P-V 8.73 Å
Lower surface roughness: mean RMS 2.08 Å, mean P-V 8.21 Å
```

### Example 2: Time Window Analysis (50-100 ns)

```bash
pywc-wc-area \
    -s system.tpr \
    -x trajectory.xtc \
    --selection-file selection.txt \
    --alpha 3.0 \
    --mesh 2.5 \
    -b 50000 \
    -e 100000 \
    --step 10 \
    --backend cpp \
    --center
```

### Example 3: GPU Acceleration for Large System

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

### Example 4: High-Resolution Thickness Map

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

### Example 5: Custom Density Threshold

```bash
pywc-wc-area \
    -s system.tpr \
    -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 \
    --mesh 2.0 \
    --density-cutoff 0.015 \
    --backend cpp \
    --center
```

## Analyzing Results

### Plotting Surface Area vs Time

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('willard_chandler_area.csv')

# Plot
plt.figure(figsize=(10, 6))
plt.plot(df['time_ps'] / 1000, df['area_angstrom2'])
plt.xlabel('Time (ns)')
plt.ylabel('Surface Area (Å^2)')
plt.title('Membrane Surface Area Evolution')
plt.grid(True, alpha=0.3)
plt.savefig('area_vs_time.png', dpi=300)
```

### Visualizing Thickness Map

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load thickness data
df = pd.read_csv('willard_chandler_thickness.csv')

# Pivot for heatmap
grid = df.pivot(index='x', columns='y', values='thickness')

# Plot
plt.figure(figsize=(10, 8))
im = plt.imshow(grid.T, cmap='viridis', origin='lower',
                aspect='auto', interpolation='bilinear')
plt.colorbar(im, label='Thickness (Å)')
plt.xlabel('X (Å)')
plt.ylabel('Y (Å)')
plt.title('Membrane Thickness Map')
plt.savefig('thickness_map.png', dpi=300)

# Statistics
print(f"Mean thickness: {grid.values.flatten().mean():.2f} Å")
print(f"Std deviation: {grid.values.flatten().std():.2f} Å")
```

### Calculating Area Per Lipid (APL)

```python
import pandas as pd

# Load data
df = pd.read_csv('willard_chandler_area.csv')

# Number of lipids in upper leaflet
n_lipids_upper = 512

# Calculate APL (divide by 2 for one leaflet)
df['apl'] = df['area_angstrom2'] / 2 / n_lipids_upper

print(f"Mean APL: {df['apl'].mean():.2f} ± {df['apl'].std():.2f} Å^2")
```

## Performance Benchmarks

Real-world performance on a 105,000 atom GOLO/GOLH interdigitated membrane (α=3.0 Å, mesh=2.5 Å):

| Backend | Time per frame | Speedup |
|---------|----------------|---------|
| Python | 427 ms | 1.0x |
| C++ (CPU) | 25 ms | **17.1x** |
| GPU (CuPy) | 20 ms | **21.4x** |

**Bottleneck analysis:**
- **Python**: `compute_surface` dominates (~97% of time)
- **C++ & GPU**: `center` becomes bottleneck (~38-46%) since surface computation is accelerated

## Troubleshooting

### Problem: No surface / area = 0

**Possible causes:**
1. Selection returned no atoms
2. Density cutoff too high or too low
3. Corrupted trajectory frame

**Solutions:**
```bash
# Check selection
echo "resname DPPC and name PO4" | \
    gmx select -s system.tpr -select -

# Try different density level
--density-level 0.5  # or 1.5

# Try absolute cutoff
--density-cutoff 0.01
```

### Problem: Surface area >> box area

**Cause:** Membrane wrapping across periodic boundaries

**Solution:** Use `--center` flag
```bash
--center
```

### Problem: GPU out of memory

**Solutions:**
```bash
# Increase mesh spacing (reduces grid size)
--mesh 3.0  # instead of 2.0

# Use CPU backend
--backend cpp

# Reduce selection to fewer atoms
```

### Problem: Thickness map has gaps

**Cause:** Low grid resolution or small trajectory

**Solutions:**
```bash
# Increase grid size
--grid-size 50

# Analyze more frames
--step 1  # instead of 10
```

## Integration with Workflows

### Batch Processing Multiple Systems

```bash
#!/bin/bash
for lipid in DPPC DOPC POPC; do
    echo "Processing ${lipid}..."
    pywc-wc-area \
        -s ${lipid}/system.tpr \
        -x ${lipid}/production.xtc \
        --selection-file ${lipid}/selection.txt \
        --alpha 3.0 \
        --mesh 2.5 \
        --backend cpp \
        --center \
        -o ${lipid}_area.csv \
        --thickness-output ${lipid}_thickness.csv
done
```

### Using with Snakemake

```python
# Snakefile
rule compute_area:
    input:
        tpr = "{system}/system.tpr",
        xtc = "{system}/production.xtc",
        sel = "{system}/selection.txt"
    output:
        area = "{system}_area.csv",
        thickness = "{system}_thickness.csv"
    shell:
        """
        pywc-wc-area \
            -s {input.tpr} \
            -x {input.xtc} \
            --selection-file {input.sel} \
            --alpha 3.0 --mesh 2.5 \
            --backend cpp --center \
            -o {output.area} \
            --thickness-output {output.thickness}
        """
```

## Tips and Best Practices

1. **Always use `--center` for membranes** - Prevents periodic boundary artifacts
2. **Start with `--step 10`** - Analyze every 10th frame for faster testing
3. **Use `--backend cpp` by default** - Best performance/compatibility balance
4. **Save selection to file** - Easier to document and reproduce
5. **Check first frame output** - Run on 1 frame first to verify parameters
6. **Monitor timing** - Use `--enable-timing` to track performance

## Command-Line Help

For the complete list of options:

```bash
pywc-wc-area --help
```

## See Also

- [CLI Tools Guide](cli-tools.md) - Comprehensive command-line tools documentation
- [API Reference](api/willard_chandler.md) - Python API for WillardChandler class
- [Examples](examples.md) - More usage examples
- [Backends](api/backends.md) - Backend performance comparison
