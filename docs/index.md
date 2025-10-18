# pyWC: Willard-Chandler Surface Analysis Toolkit

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![GitHub](https://img.shields.io/badge/github-pyWC-blue.svg)](https://github.com/octupole/pyWC)

</div>

**pyWC** is a high-performance Python toolkit for computing and analyzing **Willard-Chandler intrinsic density surfaces** from molecular dynamics simulations.

## About This Project

This package is a **specialized fork** of [pytim](https://github.com/Marcello-Sega/pytim) by Marcello Sega and collaborators, focusing exclusively on the Willard-Chandler method with significant performance improvements:

- **~35x faster** KDE evaluation through optimized C++ implementation
- **GPU acceleration** via CuPy for large systems (>10k atoms)
- **10-20x faster** system centering with OpenMP parallelization
- **Streamlined codebase** by removing unused modules

!!! info "Upstream Project"
    - **Original:** [Marcello-Sega/pytim](https://github.com/Marcello-Sega/pytim)
    - **Paper:** Sega, M., Fabian, B., & Jedlovszky, P. (2018). *Pytim: A python package for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry, 39(25), 2118-2125. [DOI: 10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384)

## Scientific Background

The Willard-Chandler method identifies intrinsic molecular surfaces by:

1. Computing a smooth density field via Gaussian kernel density estimation
2. Extracting the isosurface at critical density using marching cubes
3. Providing a triangulated surface representation for analysis

!!! quote "Reference"
    Willard, A. P., & Chandler, D. (2010). *Instantaneous liquid interfaces.* The Journal of Physical Chemistry B, 114(5), 1954-1958. [DOI: 10.1021/jp909219k](https://doi.org/10.1021/jp909219k)

## Key Features

### Multi-Backend Performance
- **CPU (C++)**: OpenMP-parallelized pybind11 extensions (~35x speedup)
- **GPU (CUDA)**: CuPy implementation with custom CUDA kernels
- **Python**: Pure NumPy/SciPy fallback for testing

### Core Capabilities
- ✓ Compute intrinsic density fields on structured 3D grids
- ✓ Extract triangulated surfaces via marching cubes
- ✓ Calculate surface area, bending rigidity, and thickness profiles
- ✓ Process entire MD trajectories with frame-to-frame continuity
- ✓ Accurate periodic boundary condition handling

### Output Formats
- **VTK**: For visualization in ParaView, VMD
- **Wavefront OBJ**: For 3D graphics software (Blender, etc.)
- **Gaussian Cube**: For quantum chemistry packages
- **PDB**: For molecular visualization

## Quick Example

```python
import MDAnalysis as mda
from pywc import WillardChandler

# Load trajectory
u = mda.Universe("system.pdb", "trajectory.xtc")
group = u.select_atoms("resname DPPC")

# Create interface with CPU backend
wc = WillardChandler(
    u,
    group=group,
    alpha=3.0,           # Gaussian width (Å)
    mesh=2.0,            # Grid spacing (Å)
    centered=True,       # Center the group
    surface_backend='cpp'
)

# Get surface
vertices, faces, normals = wc.triangulated_surface
print(f"Surface area: {wc.surface_area:.2f} Ų")

# Export
wc.writevtk.surface("surface.vtk")
```

## Performance Benchmarks

| System Size | Backend | Time (s) | Speedup |
|-------------|---------|----------|---------|
| 1,000 atoms | Python  | 2.45     | 1x      |
| 1,000 atoms | C++     | 0.12     | 20x     |
| 10,000 atoms| C++     | 1.85     | -       |
| 10,000 atoms| GPU     | 0.32     | 5.8x    |

*Benchmarks on AMD Ryzen 5950X (C++) and NVIDIA RTX 3090 (GPU)*

## Use Cases

- **Membrane biophysics**: Lipid bilayer surface characterization
- **Liquid interfaces**: Water/vapor, oil/water surfaces
- **Protein-membrane interactions**: Binding site identification
- **Nanoparticle analysis**: Micelle and vesicle shapes
- **Material science**: Polymer surfaces, nanostructures

## Getting Started

Ready to start? Check out the [Installation Guide](installation.md) and [Quick Start Tutorial](quickstart.md).

## Community

- **GitHub Repository**: [octupole/pyWC](https://github.com/octupole/pyWC)
- **Report Issues**: [GitHub Issues](https://github.com/octupole/pyWC/issues)
- **Contribute**: See our [Contributing Guide](contributing.md)
- **Contact**: Massimo Marchi (massimo@octupole.org)

## License

This project is licensed under the **GNU General Public License v3.0** (GPL-3.0) - the same license as the upstream pytim project.
