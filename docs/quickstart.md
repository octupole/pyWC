# Quick Start Guide

This guide will help you get started with pyWC in just a few minutes.

## Your First Surface Calculation

### Step 1: Prepare Your Data

pyWC works with any trajectory format supported by MDAnalysis (GROMACS, CHARMM, LAMMPS, AMBER, etc.).

For this tutorial, we'll use the test trajectory included with pyWC:

```python
import MDAnalysis as mda
from pywc import WillardChandler
from pywc.datafiles import NPT_RUN_TPR, TRAJ_TEST_XTC, SELECTION_TXT

# Load the test trajectory
u = mda.Universe(NPT_RUN_TPR, TRAJ_TEST_XTC)
```

### Step 2: Select Atoms

Select the group of atoms you want to analyze. For this example, we'll load a selection from file:

```python
# Load selection from file
with open(SELECTION_TXT) as f:
    selection_text = f.read().strip()

# Apply selection
group = u.select_atoms(selection_text)

print(f"Selected {len(group)} atoms")
```

### Step 3: Create WillardChandler Object

```python
wc = WillardChandler(
    u,
    group=group,
    alpha=3.0,           # Gaussian width in Ångströms
    mesh=2.5,            # Grid spacing in Ångströms
    centered=True,       # Center the system
    surface_backend='cpp' # Use C++ backend
)
```

### Step 4: Compute and Analyze Surface

```python
# Get triangulated surface
vertices, faces, normals = wc.triangulated_surface

# Get surface area
print(f"Surface area: {wc.surface_area:.2f} Å^2")

# Export to VTK for visualization
wc.writevtk.surface("water_surface.vtk")
```

That's it! You've computed your first Willard-Chandler surface.

## Understanding the Parameters

### `alpha` - Gaussian Width

Controls the smoothness of the density field:

- **Smaller values (1.5-2.5 Å)**: Sharper, more detailed surfaces
- **Larger values (3.0-4.0 Å)**: Smoother surfaces, less noise

**Rule of thumb**: Start with 2-3 Å for water, 3-4 Å for lipids.

```python
# Sharp surface
wc_sharp = WillardChandler(u, group=water, alpha=2.0, mesh=2.0)

# Smooth surface
wc_smooth = WillardChandler(u, group=water, alpha=4.0, mesh=2.0)
```

### `mesh` - Grid Spacing

Controls the resolution of the density grid:

- **Smaller values (1.0-1.5 Å)**: Higher resolution, slower computation
- **Larger values (2.0-3.0 Å)**: Lower resolution, faster computation

**Rule of thumb**: Use 2.0 Å for most applications.

```python
# High resolution (slower)
wc_highres = WillardChandler(u, group=water, alpha=2.4, mesh=1.0)

# Standard resolution
wc_standard = WillardChandler(u, group=water, alpha=2.4, mesh=2.0)
```

### `centered` - System Centering

Whether to center the selected group in the simulation box:

- `True`: Recommended for most cases, avoids periodic boundary artifacts
- `False`: Use if system is already properly positioned

```python
wc = WillardChandler(u, group=water, alpha=2.4, mesh=2.0, centered=True)
```

### `surface_backend` - Computational Backend

Choose the backend based on your system size and available hardware:

| Backend | When to Use | Performance |
|---------|-------------|-------------|
| `'cpp'` | Default, best for most systems | ~35x faster than Python |
| `'cupy'` | Large systems (>10k atoms) with NVIDIA GPU | ~5-6x faster than C++ |
| `'python'` | Testing, debugging, or when C++ unavailable | Baseline |

```python
# CPU backend (default, recommended)
wc_cpu = WillardChandler(u, group=water, alpha=2.4, mesh=2.0,
                         surface_backend='cpp')

# GPU backend (requires CuPy)
wc_gpu = WillardChandler(u, group=water, alpha=2.4, mesh=2.0,
                         surface_backend='cupy')
```

## Common Use Cases

### Lipid Bilayer Surfaces

```python
import MDAnalysis as mda
from pywc import WillardChandler

# Load lipid bilayer simulation
u = mda.Universe("bilayer.gro", "bilayer.xtc")

# Select phosphate groups (headgroups)
headgroups = u.select_atoms("name P or name PO4")

# Create surface
wc = WillardChandler(
    u,
    group=headgroups,
    alpha=3.0,        # Larger alpha for lipids
    mesh=2.0,
    centered=True
)

# Analyze
print(f"Membrane surface area: {wc.surface_area:.2f} Å^2")
wc.writevtk.surface("membrane.vtk")
```

### Protein Surface

```python
# Select protein
protein = u.select_atoms("protein")

# Create surface
wc = WillardChandler(
    u,
    group=protein,
    alpha=3.5,        # Larger for proteins
    mesh=2.5,
    centered=True
)

print(f"Protein surface area: {wc.surface_area:.2f} Å^2")
```

### Micelle Analysis

```python
# Select surfactant headgroups
surfactant = u.select_atoms("resname DPC and name N")

wc = WillardChandler(
    u,
    group=surfactant,
    alpha=3.0,
    mesh=2.0,
    centered=True
)

# Export for visualization
wc.writeobj("micelle.obj")
```

## Trajectory Analysis

Process multiple frames to analyze surface evolution:

```python
import numpy as np

# Track surface area over time
areas = []
times = []

for ts in u.trajectory[::10]:  # Every 10th frame
    # Surface recomputes automatically per frame
    areas.append(wc.surface_area)
    times.append(ts.time)
    print(f"Frame {ts.frame}: Area = {wc.surface_area:.2f} Å^2")

# Calculate statistics
mean_area = np.mean(areas)
std_area = np.std(areas)
print(f"\nMean area: {mean_area:.2f} ± {std_area:.2f} Å^2")

# Plot
import matplotlib.pyplot as plt
plt.plot(times, areas)
plt.xlabel("Time (ps)")
plt.ylabel("Surface Area (Å^2)")
plt.savefig("area_vs_time.png")
```

## Exporting Results

### VTK Format (ParaView, VMD)

```python
# Export surface
wc.writevtk.surface("surface.vtk")

# Export density field
wc.writevtk.density("density.vtk")

# Export particles
wc.writevtk.particles("atoms.vtk")
```

**Visualize in ParaView:**
1. Open `surface.vtk`
2. Apply "Extract Surface" filter
3. Color by normals or other properties

### Wavefront OBJ (Blender, 3D Software)

```python
wc.writeobj("surface.obj")
```

**Open in Blender:**
1. File → Import → Wavefront (.obj)
2. Select `surface.obj`

### Gaussian Cube (Quantum Chemistry)

```python
wc.writecube("density.cube")
```

Compatible with Gaussian, ORCA, VMD, etc.

### PDB (Molecular Viewers)

```python
wc.writepdb("surface.pdb")
```

## Performance Tips

### Enable Timing

Track performance of different components:

```python
wc = WillardChandler(
    u, group=water,
    alpha=2.4, mesh=2.0,
    enable_timing=True
)

# After computation
print(f"Total time: {wc.get_timing():.4f} s")
```

### Choose Appropriate Backend

```python
# For small systems (<1000 atoms)
wc = WillardChandler(u, group=small_group, surface_backend='cpp')

# For large systems (>10000 atoms) with GPU
wc = WillardChandler(u, group=large_group, surface_backend='cupy')
```

### Optimize Parameters

Balance accuracy and performance:

```python
# Fast computation (lower accuracy)
wc_fast = WillardChandler(u, group=water, alpha=3.0, mesh=3.0)

# High accuracy (slower)
wc_accurate = WillardChandler(u, group=water, alpha=2.0, mesh=1.0)
```

### Skip Warmup Frames

For trajectory analysis, skip initial equilibration:

```python
for ts in u.trajectory[100:]:  # Skip first 100 frames
    # Surface recomputes automatically per frame
    # Analysis...
```

## Next Steps

- [API Reference](api/willard_chandler.md) - Complete API documentation
- [Examples](examples.md) - More detailed examples
- [Backends](api/backends.md) - Deep dive into computational backends
- [Citation](citation.md) - How to cite pyWC in publications
