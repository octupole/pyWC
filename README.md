# pyWC: Willard-Chandler Surface Analysis Toolkit

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**pyWC** is a high-performance Python toolkit for computing and analyzing **Willard-Chandler intrinsic density surfaces** from molecular dynamics simulations. This package provides a focused, optimized implementation with multiple computational backends (CPU/GPU) for maximum performance on large molecular systems.

## About This Fork

This project is a **specialized fork** of [pytim](https://github.com/Marcello-Sega/pytim) by Marcello Sega and collaborators. While pytim provides a comprehensive suite of interfacial analysis methods (ITIM, GITIM, SASA, Willard-Chandler, and DBSCAN), pyWC focuses exclusively on the Willard-Chandler method with significant performance improvements:

- **~35x faster** KDE evaluation through optimized C++ implementation
- **GPU acceleration** via CuPy for large systems (>10k atoms)
- **10-20x faster** system centering with OpenMP parallelization
- **Streamlined codebase** by removing unused ITIM/GITIM modules

**Upstream Project:** [Marcello-Sega/pytim](https://github.com/Marcello-Sega/pytim)
**Original Publication:** Sega, M., Fabian, B., & Jedlovszky, P. (2018). *Pytim: A python package for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry, 39(25), 2118-2125. [DOI: 10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384)

## Scientific Background

The Willard-Chandler method identifies intrinsic molecular surfaces by:
1. Computing a smooth density field via Gaussian kernel density estimation
2. Extracting the isosurface at critical density using marching cubes
3. Providing a triangulated surface representation for analysis

**Reference:** Willard, A. P., & Chandler, D. (2010). *Instantaneous liquid interfaces.* The Journal of Physical Chemistry B, 114(5), 1954-1958. [DOI: 10.1021/jp909219k](https://doi.org/10.1021/jp909219k)

## Key Features

### Multi-Backend Performance
- **CPU (C++)**: OpenMP-parallelized pybind11 extensions (~35x speedup)
- **GPU (CUDA)**: CuPy implementation with custom CUDA kernels
- **Python**: Pure NumPy/SciPy fallback for testing

### Core Capabilities
- Compute intrinsic density fields on structured 3D grids
- Extract triangulated surfaces via marching cubes
- Calculate surface area, bending rigidity, and thickness profiles
- Process entire MD trajectories with frame-to-frame continuity
- Accurate periodic boundary condition handling

### Output Formats
- **VTK**: For visualization in ParaView, VMD
- **Wavefront OBJ**: For 3D graphics software (Blender, etc.)
- **Gaussian Cube**: For quantum chemistry packages
- **PDB**: For molecular visualization

## Installation

### Basic Installation (CPU only)

```bash
pip install .
```

This installs the C++ extensions with OpenMP parallelization.

### GPU Acceleration (Optional)

For GPU support, install CuPy matching your CUDA toolkit:

```bash
# For CUDA 12.x
pip install cupy-cuda12x

# Or install with GPU extras
pip install .[gpu]
```

### Dependencies

**Required:**
- Python ≥ 3.10
- NumPy ≥ 2.1.3
- SciPy ≥ 1.11.3
- scikit-image ≥ 0.24.0
- MDAnalysis ≥ 2.8.0
- pybind11 ≥ 2.11 (build-time)

**Optional:**
- CuPy ≥ 12.0 (for GPU acceleration)

## Quick Start

### Basic Usage

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
wc.writeobj("surface.obj")
wc.writecube("density.cube")
```

### GPU Acceleration

```python
# Use GPU backend for large systems
wc = WillardChandler(
    u,
    group=group,
    alpha=3.0,
    mesh=2.0,
    surface_backend='cupy',  # GPU backend
    enable_timing=True
)

# Check performance
print(f"Computation time: {wc.get_timing():.4f} s")
```

### Trajectory Analysis

```python
# Process entire trajectory
areas = []
for ts in u.trajectory:
    wc.assign_surface()
    areas.append(wc.surface_area)

# Analyze results
import numpy as np
print(f"Mean area: {np.mean(areas):.2f} ± {np.std(areas):.2f} Ų")
```

## Command-Line Tools

pyWC includes several analysis scripts:

```bash
# Compare CPU vs GPU backends
pywc-compare-wc-backends --help

# Compute surface area over trajectory
pywc-wc-area -h

# Analyze bending rigidity
pywc-bending-rigidity -h
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

## Documentation

### Parameter Selection

- **`alpha`** (Gaussian width): Typically 2-4 Å. Smaller values = sharper surfaces, larger values = smoother surfaces
- **`mesh`** (grid spacing): Typically 1-3 Å. Smaller values = more detail but slower computation
- **`centered`**: Set `True` to avoid periodic boundary artifacts

### Backend Selection

| Backend | When to Use |
|---------|-------------|
| `'cpp'` | Default for most applications (best CPU performance) |
| `'cupy'`| Large systems (>10k atoms) with NVIDIA GPU |
| `'python'` | Testing, debugging, or when C++ compilation fails |

## Project Status

**Current Version:** 1.0.4+wc

**Maintainer:** Massimo Marchi (massimo@octupole.org)

**Original pytim Authors:** Marcello Sega, Balazs Fabian, Gyorgy Hantal, Pal Jedlovszky

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

Report issues at: https://github.com/octupole/pyWC/issues

## Citation

If you use pyWC in your research, please cite both this fork and the original pytim paper:

```bibtex
@software{pywc2025,
  author = {Marchi, Massimo},
  title = {pyWC: Willard-Chandler Surface Analysis Toolkit},
  year = {2025},
  url = {https://github.com/octupole/pyWC},
  note = {Fork of pytim focused on Willard-Chandler method with performance enhancements}
}

@article{pytim2018,
  author = {Sega, Marcello and Fabian, Balazs and Jedlovszky, Pál},
  title = {Pytim: A python package for the interfacial analysis of molecular simulations},
  journal = {Journal of Computational Chemistry},
  volume = {39},
  number = {25},
  pages = {2118-2125},
  year = {2018},
  doi = {10.1002/jcc.25384}
}
```

For the Willard-Chandler method itself:

```bibtex
@article{willard2010,
  author = {Willard, Adam P. and Chandler, David},
  title = {Instantaneous liquid interfaces},
  journal = {The Journal of Physical Chemistry B},
  volume = {114},
  number = {5},
  pages = {1954-1958},
  year = {2010},
  doi = {10.1021/jp909219k}
}
```

## License

This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.

This is the same license as the upstream pytim project, ensuring compatibility and compliance with derivative work requirements.

## Acknowledgments

This work builds upon the excellent [pytim](https://github.com/Marcello-Sega/pytim) package. We gratefully acknowledge:

- **Marcello Sega** (original pytim author and maintainer)
- **Balazs Fabian**, **Gyorgy Hantal**, **Pal Jedlovszky** (pytim co-authors)
- The pytim community for the foundational implementation

Special thanks to the MDAnalysis team for providing the molecular analysis framework that both pytim and pyWC build upon.

## Related Projects

- **[pytim](https://github.com/Marcello-Sega/pytim)**: Comprehensive interfacial analysis toolkit (ITIM, GITIM, SASA, Willard-Chandler)
- **[MDAnalysis](https://www.mdanalysis.org/)**: Python library for analyzing molecular dynamics trajectories
- **[gitesei/willard-chandler](https://github.com/gitesei/willard-chandler)**: Alternative Willard-Chandler implementation with orientation analysis

---

**Questions?** Contact Massimo Marchi at massimo@octupole.org or open an issue on GitHub.
