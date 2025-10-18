# Examples

This page provides real-world examples of using pyWC for different molecular systems.

## Water/Vapor Interface

Analyze a liquid-vapor interface:

```python
import MDAnalysis as mda
from pywc import WillardChandler
import numpy as np

# Load water system
u = mda.Universe("water.gro", "water.xtc")
water = u.select_atoms("resname SOL")

# Create WillardChandler object
wc = WillardChandler(
    u,
    group=water,
    alpha=2.4,        # Standard for water
    mesh=2.0,
    centered=True,
    surface_backend='cpp'
)

# Analyze trajectory
areas = []
for ts in u.trajectory:
    wc.assign_surface()
    areas.append(wc.surface_area)

# Statistics
print(f"Mean area: {np.mean(areas):.2f} ± {np.std(areas):.2f} Ų")

# Export final frame
wc.writevtk.surface("water_surface.vtk")
```

## Lipid Bilayer

Characterize membrane surfaces:

```python
import MDAnalysis as mda
from pywc import WillardChandler

# Load DPPC bilayer
u = mda.Universe("dppc.gro", "dppc.xtc")

# Select phosphate groups (headgroups)
upper_leaflet = u.select_atoms("resname DPPC and name P and prop z > 0")
lower_leaflet = u.select_atoms("resname DPPC and name P and prop z < 0")

# Analyze upper leaflet
wc_upper = WillardChandler(
    u,
    group=upper_leaflet,
    alpha=3.0,        # Larger for lipids
    mesh=2.0,
    centered=False    # Leaflets already positioned
)

wc_upper.assign_surface()
print(f"Upper leaflet area: {wc_upper.surface_area:.2f} Ų")

# Analyze lower leaflet
wc_lower = WillardChandler(
    u,
    group=lower_leaflet,
    alpha=3.0,
    mesh=2.0
)

wc_lower.assign_surface()
print(f"Lower leaflet area: {wc_lower.surface_area:.2f} Ų")

# Export both surfaces
wc_upper.writevtk.surface("upper_leaflet.vtk")
wc_lower.writevtk.surface("lower_leaflet.vtk")
```

## Micelle Analysis

Analyze surfactant micelle:

```python
import MDAnalysis as mda
from pywc import WillardChandler
import matplotlib.pyplot as plt

# Load DPC micelle
u = mda.Universe("micelle.pdb")
headgroups = u.select_atoms("resname DPC and name N")

# Create surface
wc = WillardChandler(
    u,
    group=headgroups,
    alpha=3.0,
    mesh=2.0,
    centered=True
)

# Get surface
vertices, faces, normals = wc.triangulated_surface

# Calculate radius distribution
center = vertices.mean(axis=0)
radii = np.linalg.norm(vertices - center, axis=1)

print(f"Mean radius: {radii.mean():.2f} ± {radii.std():.2f} Å")

# Plot radius distribution
plt.hist(radii, bins=30)
plt.xlabel("Radius (Å)")
plt.ylabel("Count")
plt.savefig("micelle_radius_dist.png")

# Export
wc.writeobj("micelle.obj")
```

## Protein Surface

Compute protein surface:

```python
import MDAnalysis as mda
from pywc import WillardChandler

# Load protein system
u = mda.Universe("protein.pdb", "protein.xtc")
protein = u.select_atoms("protein")

# Create surface
wc = WillardChandler(
    u,
    group=protein,
    alpha=3.5,        # Larger for proteins
    mesh=2.5,
    centered=True,
    surface_backend='cpp'
)

# Process trajectory
for ts in u.trajectory[::10]:  # Every 10th frame
    wc.assign_surface()
    area = wc.surface_area
    print(f"Frame {ts.frame:4d}: Area = {area:.2f} Ų")

    # Save snapshots
    if ts.frame % 100 == 0:
        wc.writevtk.surface(f"protein_{ts.frame:04d}.vtk")
```

## GPU Backend for Large Systems

Use GPU acceleration for large systems:

```python
import MDAnalysis as mda
from pywc import WillardChandler
import time

# Load large system
u = mda.Universe("large_system.gro", "large_system.xtc")
atoms = u.select_atoms("all")

print(f"System size: {len(atoms)} atoms")

# Compare backends
for backend in ['cpp', 'cupy']:
    try:
        wc = WillardChandler(
            u,
            group=atoms,
            alpha=2.4,
            mesh=2.0,
            surface_backend=backend,
            enable_timing=True
        )

        start = time.time()
        wc.assign_surface()
        elapsed = time.time() - start

        print(f"{backend:10s}: {elapsed:.4f} s, Area = {wc.surface_area:.2f} Ų")

    except Exception as e:
        print(f"{backend:10s}: Not available ({e})")
```

## Trajectory Analysis with Plotting

Analyze surface evolution over time:

```python
import MDAnalysis as mda
from pywc import WillardChandler
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory
u = mda.Universe("system.gro", "trajectory.xtc")
group = u.select_atoms("resname DPPC and name P")

# Create WillardChandler
wc = WillardChandler(
    u, group=group,
    alpha=3.0, mesh=2.0,
    centered=True
)

# Collect data
times = []
areas = []
n_vertices = []

for ts in u.trajectory[::5]:  # Every 5th frame
    wc.assign_surface()

    times.append(ts.time)
    areas.append(wc.surface_area)

    vertices, _, _ = wc.triangulated_surface
    n_vertices.append(len(vertices))

# Convert to arrays
times = np.array(times)
areas = np.array(areas)
n_vertices = np.array(n_vertices)

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Surface area
ax1.plot(times, areas, 'b-', linewidth=0.5, alpha=0.7)
ax1.axhline(areas.mean(), color='r', linestyle='--',
            label=f'Mean: {areas.mean():.1f} Ų')
ax1.fill_between(times,
                 areas.mean()-areas.std(),
                 areas.mean()+areas.std(),
                 alpha=0.2, color='r')
ax1.set_ylabel('Surface Area (Ų)')
ax1.legend()
ax1.grid(alpha=0.3)

# Number of vertices
ax2.plot(times, n_vertices, 'g-', linewidth=0.5, alpha=0.7)
ax2.set_xlabel('Time (ps)')
ax2.set_ylabel('Number of Vertices')
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('surface_analysis.png', dpi=300)

# Statistics
print(f"Area: {areas.mean():.2f} ± {areas.std():.2f} Ų")
print(f"Vertices: {n_vertices.mean():.0f} ± {n_vertices.std():.0f}")
```

## Computing Bending Rigidity

Use the built-in script:

```bash
pywc-bending-rigidity \
    --topology system.gro \
    --trajectory traj.xtc \
    --selection "resname DPPC and name P" \
    --alpha 3.0 \
    --mesh 2.0 \
    --output bending_rigidity.csv
```

Or from Python:

```python
from scripts.compute_bending_rigidity import analyze_bending_rigidity

results = analyze_bending_rigidity(
    topology="system.gro",
    trajectory="traj.xtc",
    selection="resname DPPC and name P",
    alpha=3.0,
    mesh=2.0
)

print(f"Bending modulus: {results['kappa']:.3f} kT")
```

## Thickness Mapping

Compute membrane thickness:

```python
from scripts.compute_willard_chandler_area import compute_thickness_map

thickness_map = compute_thickness_map(
    topology="bilayer.gro",
    trajectory="bilayer.xtc",
    selection_upper="resname DPPC and name P and prop z > 0",
    selection_lower="resname DPPC and name P and prop z < 0",
    grid_spacing=2.0
)

# Visualize
import matplotlib.pyplot as plt

plt.imshow(thickness_map, cmap='viridis', origin='lower')
plt.colorbar(label='Thickness (Å)')
plt.xlabel('X (grid points)')
plt.ylabel('Y (grid points)')
plt.title('Membrane Thickness Map')
plt.savefig('thickness_map.png', dpi=300)
```

## Export to Multiple Formats

```python
from pywc import WillardChandler

# ... create and compute surface ...

# VTK for ParaView
wc.writevtk.surface("surface.vtk")
wc.writevtk.density("density.vtk")
wc.writevtk.atoms("atoms.vtk")

# Wavefront OBJ for Blender
wc.writeobj("surface.obj")

# Gaussian Cube for quantum chemistry
wc.writecube("density.cube")

# PDB for molecular viewers
wc.writepdb.surface("surface.pdb")
```

## Custom Analysis

Extract surface properties for custom analysis:

```python
from pywc import WillardChandler
import numpy as np

# Compute surface
wc = WillardChandler(u, group=atoms, alpha=2.4, mesh=2.0)
vertices, faces, normals = wc.triangulated_surface

# Calculate local curvature
def compute_curvature(vertices, faces):
    # Custom curvature calculation
    # (simplified example)
    curvatures = []
    for i, vertex in enumerate(vertices):
        # Find neighboring faces
        neighbor_faces = faces[np.any(faces == i, axis=1)]
        # Compute local curvature
        # ... your calculation ...
        curvatures.append(curvature_value)
    return np.array(curvatures)

curvatures = compute_curvature(vertices, faces)

# Visualize
import matplotlib.pyplot as plt
plt.hist(curvatures, bins=50)
plt.xlabel('Curvature (1/Å)')
plt.ylabel('Count')
plt.savefig('curvature_distribution.png')
```

## See Also

- [API Reference](api/willard_chandler.md) - Complete API documentation
- [Quick Start](quickstart.md) - Basic usage
- [Backends](api/backends.md) - Backend selection guide
