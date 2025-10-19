# pyWC: Membrane Interface Analysis with Willard-Chandler Surfaces

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**pyWC** is a high-performance Python toolkit for **analyzing membrane interfaces** in molecular dynamics simulations. It computes intrinsic density surfaces using the Willard-Chandler method, enabling quantitative analysis of **membrane surfaces**, **thickness profiles**, and **bending rigidity**. The package is optimized for large biological membrane systems with multiple computational backends (CPU/GPU).

## Main Contributions

### 1. Command-Line Tools for Membrane Analysis
**pyWC provides ready-to-use scripts** in `./scripts/` for common membrane analysis tasks:

- **`pywc-wc-area`**: Compute surface area and thickness profiles over trajectories
- **`pywc-bending-rigidity`**: Estimate bending modulus from surface fluctuations
- **`pywc-compare-wc-backends`**: Benchmark different computational backends

These tools enable immediate analysis without writing custom Python code.

### 2. Optimized Willard-Chandler Implementation
Significant performance improvements over the original pytim implementation:

- **~35x faster** KDE evaluation through optimized C++ implementation with OpenMP
- **GPU acceleration** via CuPy for large membrane systems (>10k atoms)
- **10-20x faster** system centering with parallel algorithms
- **Efficient memory management** for trajectory processing

### 3. Membrane-Specific Features
- **Thickness profiles**: Compute local membrane thickness along interface normals
- **Bending rigidity**: Extract mechanical properties from surface fluctuations
- **Surface curvature**: Analyze membrane deformations and undulations
- **Interface tracking**: Maintain surface continuity across MD frames

## About This Fork

This project is a **specialized fork** of [pytim](https://github.com/Marcello-Sega/pytim) by Marcello Sega and collaborators. While pytim provides a comprehensive suite of interfacial analysis methods (ITIM, GITIM, SASA, Willard-Chandler, and DBSCAN), pyWC focuses exclusively on the Willard-Chandler method optimized for membrane analysis.

**Upstream Project:** [Marcello-Sega/pytim](https://github.com/Marcello-Sega/pytim)
**Original Publication:** Sega, M., Fabian, B., & Jedlovszky, P. (2018). *Pytim: A python package for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry, 39(25), 2118-2125. [DOI: 10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384)

## Scientific Background

The Willard-Chandler method identifies **intrinsic membrane surfaces** by:
1. Computing a smooth density field from membrane atoms via Gaussian kernel density estimation (KDE)
2. Extracting the isosurface at critical density using marching cubes algorithm
3. Providing a triangulated surface representation for geometric and mechanical analysis

This approach is particularly powerful for **biological membranes** where:
- Thermal fluctuations create complex, non-planar surfaces
- Local thickness variations reflect membrane composition and interactions
- Surface curvature and undulations encode mechanical properties

**Reference:** Willard, A. P., & Chandler, D. (2010). *Instantaneous liquid interfaces.* The Journal of Physical Chemistry B, 114(5), 1954-1958. [DOI: 10.1021/jp909219k](https://doi.org/10.1021/jp909219k)

## Key Features

### Multi-Backend Performance
- **CPU (C++)**: OpenMP-parallelized pybind11 extensions (~35x speedup)
- **GPU (CUDA)**: CuPy implementation with custom CUDA kernels for large membranes
- **Python**: Pure NumPy/SciPy fallback for testing

### Membrane Analysis Capabilities
- **Surface extraction**: Compute intrinsic membrane surfaces from lipid density fields
- **Thickness profiles**: Calculate local membrane thickness with spatial resolution
- **Surface area**: Track membrane area changes during phase transitions or pore formation
- **Bending rigidity**: Estimate mechanical properties from surface undulation spectra
- **Curvature analysis**: Quantify local and global membrane deformations
- **Trajectory processing**: Analyze entire MD simulations with frame-to-frame continuity

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

### Command-Line Analysis

For quick membrane analysis without coding:

```bash
# Compute membrane surface area and thickness over trajectory
pywc-wc-area -s membrane.tpr -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 --mesh 2.5 \
    -b 0 -e 1000 --step 10

# Estimate bending rigidity
pywc-bending-rigidity -s membrane.tpr -x trajectory.xtc \
    --selection-file atoms.txt \
    --alpha 3.0 --mesh 2.5
```

### Python API for Membrane Analysis

```python
import MDAnalysis as mda
from pywc import WillardChandler

# Load membrane simulation
u = mda.Universe("membrane.tpr", "trajectory.xtc")

# Select lipid headgroups for surface calculation
lipids = u.select_atoms("resname DPPC DOPC and name PO4")

# Compute Willard-Chandler surface
wc = WillardChandler(
    u,
    group=lipids,
    alpha=3.0,           # Gaussian width (Å)
    mesh=2.5,            # Grid spacing (Å)
    centered=True,       # Center membrane in box
    surface_backend='cpp'  # Use optimized C++ backend
)

# Analyze membrane properties
print(f"Membrane surface area: {wc.surface_area:.2f} Ų")
print(f"Average thickness: {wc.thickness_profile.mean():.2f} Å")

# Export for visualization
wc.writevtk.surface("membrane_surface.vtk")  # ParaView
wc.writeobj("membrane_surface.obj")          # Blender
wc.writecube("density.cube")                 # VMD
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

Real-world membrane analysis using `compute_willard_chandler_area.py`:

| System | Atoms | Backend | Time (ms) | Speedup vs Python |
|--------|-------|---------|-----------|-------------------|
| GOLO/GOLH interdigitated membrane | 105,000 | Python | 427 | 1.0x |
| GOLO/GOLH interdigitated membrane | 105,000 | C++ (cpu) | 25 | **17.1x** |
| GOLO/GOLH interdigitated membrane | 105,000 | GPU (cupy) | 20 | **21.4x** |

*Benchmarks measured with the `pywc-wc-area` command-line tool on a large biological membrane system (α=3.0 Å, mesh=2.5 Å)*

**What's measured:** Complete Willard-Chandler surface computation per frame, including:

1. System centering and grid preparation (`prepare_box`)
2. KDE density field evaluation (`define_cluster_group`) - **most compute-intensive step**
3. Marching cubes isosurface extraction (`compute_surface`)

The timing is built into the script via `enable_timing=True` (see [compute_willard_chandler_area.py:364](scripts/compute_willard_chandler_area.py#L364)). Time reported is mean per-frame (excluding first 2 frames to remove initialization overhead).

**Key insight**: The C++ backend provides ~17x speedup for production membrane analysis, while GPU provides ~21x speedup for large systems.

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
