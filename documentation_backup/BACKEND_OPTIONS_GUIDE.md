# Surface Backend Options Guide

## Overview

The pytim Willard-Chandler implementation now supports three distinct backend options for surface computation, making it clear which implementation is being used:

1. **cpp** (or **cpu**) - C++ accelerated via `_wc_kde.cpp` - **RECOMMENDED**
2. **cupy** (or **gpu**) - GPU accelerated via CuPy/CUDA
3. **python** - Pure Python (for testing/debugging)

## Backend Selection

### Command Line

Use the `--backend` option with `pytim-wc-area`:

```bash
# C++ backend (default, recommended)
pytim-wc-area -s topology.tpr -x trajectory.xtc \
    --selection "name OW" --backend cpp

# GPU backend (requires CuPy)
pytim-wc-area -s topology.tpr -x trajectory.xtc \
    --selection "name OW" --backend cupy

# Pure Python backend (for debugging)
pytim-wc-area -s topology.tpr -x trajectory.xtc \
    --selection "name OW" --backend python
```

### Python API

```python
import pytim
import MDAnalysis as mda

u = mda.Universe('topology.tpr', 'trajectory.xtc')
g = u.select_atoms('name OW')

# C++ backend (default, recommended)
inter = pytim.WillardChandler(
    u, group=g, alpha=3.0, mesh=2.5,
    surface_backend='cpp'
)

# GPU backend
inter = pytim.WillardChandler(
    u, group=g, alpha=3.0, mesh=2.5,
    surface_backend='cupy'
)

# Pure Python backend
inter = pytim.WillardChandler(
    u, group=g, alpha=3.0, mesh=2.5,
    surface_backend='python'
)
```

## Performance Comparison

### Benchmark Results

Test system: 12,000 atoms, 4,000 in interface group

| Backend | Time per Frame | Speedup vs Python |
|---------|---------------|-------------------|
| **cpp** | 4.50 ms | **33.0x faster** |
| cupy (GPU) | 2.50 ms* | 59.5x faster |
| python | 148.68 ms | 1.0x (baseline) |

*GPU performance depends on system size and GPU model

### When to Use Each Backend

#### cpp/cpu (Default - Recommended)
✅ **Use for:**
- All general-purpose computations
- Small to medium systems (<100K atoms)
- When C++ compiler is available
- Production workflows

❌ **Avoid when:**
- C++ compilation not possible
- Debugging KDE algorithm internals

#### cupy/gpu
✅ **Use for:**
- Very large systems (>100K atoms)
- When GPU is already in use elsewhere
- Multiple trajectory analysis (GPU stays warm)

❌ **Avoid when:**
- CuPy not installed
- Small systems (overhead dominates)
- No CUDA-capable GPU

#### python
✅ **Use for:**
- Debugging algorithm behavior
- Testing modifications to KDE
- Verifying C++/GPU implementations

❌ **Avoid when:**
- Performance matters
- Production workflows

## Implementation Details

### C++ Backend (cpp)

**File**: [pytim/_wc_kde.cpp](pytim/_wc_kde.cpp)

The C++ implementation uses:
- **pybind11** for Python bindings
- **OpenMP** for CPU parallelization
- **KDTree-based neighbor search** for efficient density evaluation

The implementation provides two entry points:
1. `evaluate_pbc_fast_auto()` - Automatic neighbor search
2. `evaluate_pbc_fast()` - Manual neighbor list (for compatibility)

The C++ code is automatically used by `gaussian_kde_pbc.evaluate_pbc_fast()` when available.

### GPU Backend (cupy)

**File**: [pytim/wc_core/gpu.py](pytim/wc_core/gpu.py)

The GPU implementation uses:
- **CuPy RawKernel** for direct CUDA kernel execution
- **Parallel density evaluation** across all grid points
- **GPU memory management** with automatic transfers

### Python Backend (python)

**File**: [pytim/gaussian_kde_pbc.py](pytim/gaussian_kde_pbc.py)

The pure Python implementation uses:
- **scipy.spatial.cKDTree** for neighbor search (with PBC)
- **NumPy vectorization** for density computation
- **No C++/GPU dependencies**

## Backend Mapping

For backward compatibility, `cpu` maps to `cpp`:

```python
# These are equivalent
surface_backend='cpu'   # Maps to C++ implementation
surface_backend='cpp'   # Explicit C++ implementation
```

Both `cupy` and `gpu` are accepted for the GPU backend:

```python
# These are equivalent
surface_backend='cupy'  # CuPy-based GPU
surface_backend='gpu'   # Same as cupy
```

## Numerical Accuracy

All backends produce nearly identical results within numerical precision:

| Comparison | Area Difference | Vertices Difference |
|------------|----------------|---------------------|
| cpp vs python | 0.0 Ų | 0.0 Å |
| cpp vs cupy | 0.0000007 Ų | 0.0000044 Å |

The tiny differences between C++ and GPU are due to:
- Different floating-point operation order
- GPU uses single-precision internally (converted to double for output)

These differences are well within acceptable numerical precision for molecular simulations.

## Testing

### Verify All Backends Work

```bash
python test_all_backends.py
```

Expected output:
```
✓ CPP backend successful
✓ CPU backend successful (maps to CPP)
✓ PYTHON backend successful
✓ CUPY backend successful (if CuPy installed)
✓ CPU correctly maps to CPP backend
✓ CPP and Python results are consistent
✓ CPP and CuPy results are consistent
```

### Command-Line Backend Test

```bash
# Test C++ backend
python scripts/compute_willard_chandler_area.py \
    -s pytim/data/water.gro \
    -x pytim/data/water.xtc \
    -b 0 -e 10 --selection "name OW" \
    --backend cpp

# Test Python backend
python scripts/compute_willard_chandler_area.py \
    -s pytim/data/water.gro \
    -x pytim/data/water.xtc \
    -b 0 -e 10 --selection "name OW" \
    --backend python

# Compare timing difference (should see ~33x speedup with cpp)
```

## Requirements

### cpp/cpu Backend
- C++ compiler with OpenMP support
- pybind11 >= 2.6
- NumPy
- **Automatically enabled** during standard installation

### cupy/gpu Backend
- CUDA-capable GPU
- CuPy >= 10.0
- CUDA Toolkit

Install CuPy:
```bash
pip install cupy-cuda11x  # For CUDA 11.x
# or
pip install cupy-cuda12x  # For CUDA 12.x
```

### python Backend
- No additional requirements
- Always available as fallback

## Troubleshooting

### "Unknown backend" Error

```
ValueError: Unknown backend: xyz. Must be 'cpp', 'cupy'/'gpu', or 'python'.
```

**Solution**: Use one of the valid backend names: `cpp`, `cpu`, `cupy`, `gpu`, or `python`

### CuPy Backend Fails

```
RuntimeError: CuPy backend requested but could not be used.
```

**Solutions**:
1. Install CuPy: `pip install cupy-cuda11x` (or appropriate CUDA version)
2. Verify GPU is available: `python -c "import cupy; print(cupy.cuda.runtime.getDeviceCount())"`
3. Fall back to C++ backend: use `--backend cpp`

### C++ Backend Not Available

If `_wc_kde` module is not compiled:

```bash
# Reinstall with C++ extensions
pip install --no-build-isolation --force-reinstall -e .
```

## Migration Guide

### From Old Code

If you were using:

```python
# Old: implicit C++ usage
inter = pytim.WillardChandler(u, group=g, surface_backend='cpu')
```

This still works! `cpu` now explicitly maps to the C++ backend, making it clear that C++ acceleration is being used.

### Explicit Backend Selection

For clarity in new code, prefer:

```python
# New: explicit C++ backend
inter = pytim.WillardChandler(u, group=g, surface_backend='cpp')
```

This makes it immediately clear that the C++ implementation is being used.

## Summary

The backend system now provides three clear options:

1. **cpp** - Fast C++ (default, recommended)
2. **cupy** - GPU acceleration (for very large systems)
3. **python** - Pure Python (for debugging)

The default `cpu` backend maps to `cpp` for backward compatibility. For best performance on typical systems, use the **cpp** backend (or leave as default).

---

**Related Files:**
- C++ implementation: [pytim/_wc_kde.cpp](pytim/_wc_kde.cpp)
- GPU implementation: [pytim/wc_core/gpu.py](pytim/wc_core/gpu.py)
- Surface computation: [pytim/wc_core/surface.py](pytim/wc_core/surface.py)
- Command-line tool: [scripts/compute_willard_chandler_area.py](scripts/compute_willard_chandler_area.py)
- Test script: [test_all_backends.py](test_all_backends.py)
