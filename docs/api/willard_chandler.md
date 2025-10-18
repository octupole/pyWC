# WillardChandler API Reference

The `WillardChandler` class is the main interface for computing intrinsic density surfaces.

## Class Definition

```python
from pywc import WillardChandler

wc = WillardChandler(
    universe,
    group=None,
    alpha=2.4,
    mesh=2.0,
    centered=False,
    surface_backend='cpp',
    enable_timing=False,
    **kwargs
)
```

## Parameters

### `universe` (MDAnalysis.Universe)
**Required.** The MDAnalysis Universe object containing the trajectory.

```python
import MDAnalysis as mda
u = mda.Universe("topology.pdb", "trajectory.xtc")
```

### `group` (MDAnalysis.AtomGroup, optional)
The group of atoms to analyze. If `None`, uses all atoms in the universe.

```python
# Select specific atoms
group = u.select_atoms("resname DPPC and name P")

# Or use all atoms
group = None  # Will use u.atoms
```

### `alpha` (float, default=2.4)
Gaussian width parameter in Ångströms. Controls the smoothness of the density field.

- **Range**: Typically 1.5 - 4.0 Å
- **Smaller values**: Sharper, more detailed surfaces
- **Larger values**: Smoother surfaces, less noise

**Recommendations:**
- Water: 2.0 - 2.5 Å
- Lipids: 3.0 - 3.5 Å
- Proteins: 3.0 - 4.0 Å

```python
wc = WillardChandler(u, group=water, alpha=2.4)
```

### `mesh` (float, default=2.0)
Grid spacing in Ångströms. Controls the resolution of the density grid.

- **Range**: Typically 1.0 - 3.0 Å
- **Smaller values**: Higher resolution, slower
- **Larger values**: Lower resolution, faster

```python
wc = WillardChandler(u, group=water, mesh=2.0)
```

### `centered` (bool, default=False)
Whether to center the selected group in the simulation box.

- `True`: Recommended for most cases
- `False`: Use if system is already positioned correctly

```python
wc = WillardChandler(u, group=water, centered=True)
```

### `surface_backend` (str, default='cpp')
Backend for density computation.

**Options:**
- `'cpp'`: C++ with OpenMP (default, recommended)
- `'cupy'`: GPU via CuPy/CUDA
- `'python'`: Pure Python (slow, for testing)

```python
# CPU backend
wc_cpu = WillardChandler(u, group=water, surface_backend='cpp')

# GPU backend
wc_gpu = WillardChandler(u, group=water, surface_backend='cupy')
```

### `enable_timing` (bool, default=False)
Enable performance timing diagnostics.

```python
wc = WillardChandler(u, group=water, enable_timing=True)
# ... after computation ...
print(f"Time: {wc.get_timing():.4f} s")
```

## Attributes

### `triangulated_surface`
Returns the triangulated surface as a tuple of (vertices, faces, normals).

**Returns:**
- `vertices` (ndarray, shape (N, 3)): Surface vertex coordinates in Ångströms
- `faces` (ndarray, shape (M, 3)): Triangle indices
- `normals` (ndarray, shape (N, 3)): Vertex normal vectors

```python
vertices, faces, normals = wc.triangulated_surface

print(f"Number of vertices: {len(vertices)}")
print(f"Number of faces: {len(faces)}")
```

### `surface_area`
Total surface area in Ų.

```python
area = wc.surface_area
print(f"Surface area: {area:.2f} Ų")
```

### `density_field`
The computed 3D density field on the grid.

**Returns:** ndarray, shape (nx, ny, nz)

```python
density = wc.density_field
print(f"Density shape: {density.shape}")
print(f"Min/Max density: {density.min():.3f} / {density.max():.3f}")
```

### `grid`
Grid points where density is evaluated.

**Returns:** ndarray, shape (nx*ny*nz, 3)

```python
grid_points = wc.grid
```

## Methods

### `assign_surface()`
Compute the surface for the current frame.

```python
wc.assign_surface()
```

**Use case:** Call this when iterating over trajectory frames:

```python
for ts in u.trajectory:
    wc.assign_surface()
    area = wc.surface_area
    print(f"Frame {ts.frame}: {area:.2f} Ų")
```

### `get_timing()`
Get the total computation time (if `enable_timing=True`).

**Returns:** float - Time in seconds

```python
wc = WillardChandler(u, group=water, enable_timing=True)
wc.assign_surface()
time = wc.get_timing()
print(f"Computation time: {time:.4f} s")
```

### `writevtk.surface(filename)`
Export surface to VTK format.

**Parameters:**
- `filename` (str): Output file path

```python
wc.writevtk.surface("surface.vtk")
```

### `writevtk.density(filename)`
Export density field to VTK format.

```python
wc.writevtk.density("density.vtk")
```

### `writevtk.atoms(filename)`
Export atom positions to VTK format.

```python
wc.writevtk.atoms("atoms.vtk")
```

### `writeobj(filename)`
Export surface to Wavefront OBJ format.

```python
wc.writeobj("surface.obj")
```

### `writecube(filename)`
Export density to Gaussian Cube format.

```python
wc.writecube("density.cube")
```

### `writepdb.surface(filename)`
Export surface vertices as PDB.

```python
wc.writepdb.surface("surface.pdb")
```

## Complete Example

```python
import MDAnalysis as mda
from pywc import WillardChandler
import numpy as np

# Load trajectory
u = mda.Universe("system.pdb", "traj.xtc")
lipids = u.select_atoms("resname DPPC and name P")

# Create WillardChandler object
wc = WillardChandler(
    u,
    group=lipids,
    alpha=3.0,
    mesh=2.0,
    centered=True,
    surface_backend='cpp',
    enable_timing=True
)

# Analyze trajectory
areas = []
for ts in u.trajectory[::10]:  # Every 10th frame
    wc.assign_surface()
    areas.append(wc.surface_area)

    # Export every 100th frame
    if ts.frame % 100 == 0:
        wc.writevtk.surface(f"surface_frame_{ts.frame}.vtk")

# Statistics
mean_area = np.mean(areas)
std_area = np.std(areas)
print(f"Mean area: {mean_area:.2f} ± {std_area:.2f} Ų")
print(f"Average time per frame: {wc.get_timing():.4f} s")

# Get final surface
vertices, faces, normals = wc.triangulated_surface
print(f"Final surface: {len(vertices)} vertices, {len(faces)} faces")
```

## See Also

- [Backends](backends.md) - Backend selection and performance
- [Quick Start](../quickstart.md) - Basic usage examples
- [Examples](../examples.md) - Real-world applications
