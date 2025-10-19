# pyWC: Membrane Interface Analysis Toolkit

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![GitHub](https://img.shields.io/badge/github-pyWC-blue.svg)](https://github.com/octupole/pyWC)

</div>

**pyWC** is a high-performance Python toolkit for **analyzing membrane interfaces** in molecular dynamics simulations using the Willard-Chandler method. It computes intrinsic density surfaces, enabling quantitative analysis of **membrane surfaces**, **thickness profiles**, and **bending rigidity**.

## Main Contributions

### 1. Command-Line Tools
Ready-to-use scripts in `./scripts/` for membrane analysis:

- **`pywc-wc-area`**: Compute surface area over trajectories
- **`pywc-wc-thickness`**: Calculate membrane thickness profiles
- **`pywc-bending-rigidity`**: Estimate bending modulus from fluctuations
- **`pywc-compare-wc-backends`**: Benchmark computational backends

### 2. Performance Optimizations
Significant improvements over the original pytim implementation:

- **~35x faster** KDE evaluation through optimized C++ with OpenMP
- **GPU acceleration** via CuPy for large membrane systems (>10k atoms)
- **10-20x faster** system centering with parallel algorithms
- **Efficient memory management** for trajectory processing

!!! info "Upstream Project"
    - **Original:** [Marcello-Sega/pytim](https://github.com/Marcello-Sega/pytim)
    - **Paper:** Sega, M., Fabian, B., & Jedlovszky, P. (2018). *Pytim: A python package for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry, 39(25), 2118-2125. [DOI: 10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384)

## Scientific Background

The Willard-Chandler method identifies **intrinsic membrane surfaces** by:

1. Computing a smooth density field from membrane atoms via Gaussian kernel density estimation
2. Extracting the isosurface at critical density using marching cubes
3. Providing a triangulated surface representation for geometric and mechanical analysis

This is particularly powerful for **biological membranes** where thermal fluctuations create complex surfaces and local thickness variations encode mechanical properties.

!!! quote "Reference"
    Willard, A. P., & Chandler, D. (2010). *Instantaneous liquid interfaces.* The Journal of Physical Chemistry B, 114(5), 1954-1958. [DOI: 10.1021/jp909219k](https://doi.org/10.1021/jp909219k)

## Key Features

### Multi-Backend Performance
- **CPU (C++)**: OpenMP-parallelized pybind11 extensions
- **GPU (CUDA)**: CuPy implementation for large membranes
- **Python**: NumPy/SciPy fallback

### Membrane Analysis Capabilities
- ✓ **Surface extraction**: Intrinsic membrane surfaces from lipid density
- ✓ **Thickness profiles**: Local membrane thickness with spatial resolution
- ✓ **Surface area**: Track area changes during transitions or pore formation
- ✓ **Bending rigidity**: Estimate mechanical properties from undulation spectra
- ✓ **Trajectory processing**: Analyze MD simulations with continuity
- ✓ **Periodic boundaries**: Accurate handling of membrane systems

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
