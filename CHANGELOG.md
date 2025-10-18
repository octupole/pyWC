# Changelog

All notable changes to pyWC will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.4+wc] - 2025-10-18

### Added - Initial Fork Release

This is the initial release of **pyWC**, a specialized fork of [pytim](https://github.com/Marcello-Sega/pytim) focusing exclusively on the Willard-Chandler method with significant performance improvements.

#### New Features
- **GPU Acceleration**: Full CuPy/CUDA backend implementation
  - Custom CUDA kernels for Gaussian accumulation
  - Cell-list neighbor search on GPU
  - Chunked processing for memory management
  - ~5-6x speedup over C++ for large systems (>10k atoms)

- **Enhanced C++ Backend**: Optimized KDE evaluation
  - `evaluate_pbc_fast_auto()`: Self-contained cell-list generation (~35x faster than Python)
  - OpenMP parallelization with atomic updates
  - Efficient periodic boundary condition handling

- **Fast Centering Algorithms**: Multiple implementation levels
  - Pure Python fallback
  - Partially optimized C++ (histogram acceleration)
  - Fully optimized C++ (bypasses Python overhead)
  - GPU-based centering (~10-20x faster than Python)

- **Comprehensive Documentation**:
  - Complete README.md with usage examples
  - CITATION.cff for academic citations
  - AUTHORS.md documenting contributors
  - Performance benchmarks

- **Command-Line Tools**:
  - `pywc-compare-wc-backends`: Benchmark CPU vs GPU
  - `pywc-wc-area`: Compute surface area over trajectories
  - `pywc-bending-rigidity`: Analyze bending rigidity properties

#### Performance Improvements
- **~35x faster** KDE evaluation (C++ vs Python)
- **~10-20x faster** system centering (C++ vs Python)
- **GPU support** for systems with >10k atoms
- OpenMP parallelization throughout

#### Removed from pytim
- ITIM (Identification of Truly Interfacial Molecules) module
- GITIM (Generalized ITIM) module
- SASA (Solvent Accessible Surface Area) module
- Legacy observables not relevant to Willard-Chandler analysis
- Unused utility functions and dependencies
- ~60% code reduction for focused, maintainable codebase

#### Changed from pytim
- Renamed package from `pytim` to `pywc`
- Streamlined API focusing on `WillardChandler` class
- Updated dependencies to modern versions:
  - NumPy ≥ 2.1.3
  - SciPy ≥ 1.11.3
  - scikit-image ≥ 0.24.0
  - MDAnalysis ≥ 2.8.0
  - Python ≥ 3.10

- Enhanced backend selection system:
  - `surface_backend='cpp'` (default, optimized C++)
  - `surface_backend='cupy'` (GPU acceleration)
  - `surface_backend='python'` (fallback)

- Improved output format support:
  - VTK files for ParaView/VMD
  - Wavefront OBJ for 3D graphics
  - Gaussian Cube for quantum chemistry
  - PDB for molecular visualization

#### Fixed
- Removed debug print statement in centering implementation
- Improved error handling in GPU backend
- Better memory management for large systems
- Fixed periodic boundary condition edge cases

### Technical Details

#### C++ Extensions (pybind11)
- `_wc_kde.cpp`: KDE evaluation with cell-list neighbor search
- `center_fast.cpp`: Parallelized histogram-based centering
- `center_fast_full.cpp`: Fully optimized centering
- All use OpenMP for CPU parallelization

#### GPU Implementation (CuPy)
- `center_gpu.py`: GPU-accelerated centering
- `wc_core/gpu.py`: CUDA kernels for density evaluation
- Automatic fallback to CPU if GPU unavailable

#### Build System
- Modern `pyproject.toml` configuration
- Automatic CUDA detection
- Optional CuPy installation
- Cross-platform support (Linux, macOS, Windows)

### Fork Rationale

pyWC was created to:

1. **Focus**: Provide a specialized, high-performance toolkit for Willard-Chandler analysis
2. **Performance**: Optimize for large molecular systems with GPU acceleration
3. **Maintainability**: Reduce codebase complexity by removing unused modules
4. **Modularity**: Clean separation of CPU/GPU backends

### Acknowledgments

This fork builds upon the excellent work of:
- **Marcello Sega** (original pytim author)
- **Balazs Fabian**, **Gyorgy Hantal**, **Pál Jedlovszky** (pytim co-authors)

See [AUTHORS.md](AUTHORS.md) for complete attribution.

### Migration from pytim

If you're migrating from pytim's Willard-Chandler implementation:

```python
# pytim (old)
from pytim import WillardChandler
wc = WillardChandler(u, group=group, alpha=3.0, mesh=2.0)

# pyWC (new) - same interface!
from pywc import WillardChandler
wc = WillardChandler(u, group=group, alpha=3.0, mesh=2.0,
                     surface_backend='cpp')  # or 'cupy' for GPU
```

The API is backward-compatible for Willard-Chandler functionality.

---

## pytim History (Pre-Fork)

This project was forked from pytim at commit `a78b2e1` (first commit to pyWC repository).

For the complete history of pytim development, see:
- [pytim repository](https://github.com/Marcello-Sega/pytim)
- [pytim releases](https://github.com/Marcello-Sega/pytim/releases)
- [pytim changelog](https://github.com/Marcello-Sega/pytim/blob/master/CHANGELOG.md)

### pytim Reference

Sega, M., Fabian, B., & Jedlovszky, P. (2018). Pytim: A python package for the interfacial analysis of molecular simulations. *Journal of Computational Chemistry*, 39(25), 2118-2125. [DOI: 10.1002/jcc.25384](https://doi.org/10.1002/jcc.25384)

---

[Unreleased]: https://github.com/octupole/pyWC/compare/v1.0.4+wc...HEAD
[1.0.4+wc]: https://github.com/octupole/pyWC/releases/tag/v1.0.4+wc
