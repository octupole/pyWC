# Computational Backends

pyWC provides three computational backends for density evaluation, each optimized for different use cases.

## Backend Overview

| Backend | Implementation | Performance | Use Case |
|---------|----------------|-------------|----------|
| **C++** | pybind11 + OpenMP | ~35x faster | Default, production use |
| **GPU** | CuPy/CUDA | ~5-6x faster than C++ | Large systems (>10k atoms) |
| **Python** | NumPy/SciPy | Baseline | Testing, debugging |

## C++ Backend

The default and recommended backend for most applications.

### Features

- **OpenMP parallelization**: Uses all available CPU cores
- **Cell-list neighbor search**: Efficient neighbor finding with PBC
- **Atomic operations**: Thread-safe parallel accumulation
- **~35x speedup** over pure Python

### Usage

```python
from pywc import WillardChandler

wc = WillardChandler(
    u, group=atoms,
    alpha=2.4, mesh=2.0,
    surface_backend='cpp'  # Default
)
```

### Requirements

- C++17 compatible compiler
- OpenMP support
- pybind11 ≥ 2.11

### Performance Characteristics

- **Memory**: O(N + grid_size)
- **Time complexity**: O(N × neighbors)
- **Scales linearly** with number of atoms
- **Best for**: 100 - 100,000 atoms

### Implementation Details

The C++ backend uses two main algorithms:

**1. `evaluate_pbc_fast_auto()`** (Recommended)
- Self-contained cell-list generation
- Automatic neighbor search
- ~35x faster than Python
- No preprocessing required

**2. `evaluate_pbc_fast()`** (Legacy)
- Accepts pre-computed neighbor lists
- Backward compatible with pytim
- Slightly slower due to Python overhead

```python
# Automatic (recommended)
from pywc._wc_kde import evaluate_pbc_fast_auto

result = evaluate_pbc_fast_auto(
    positions,    # Atomic positions
    grid,         # Grid points
    box,          # Box dimensions
    alpha,        # Gaussian width
    weights       # Atomic weights
)
```

### Compiler Optimization Flags

The build system uses:
- `-O3`: Aggressive optimization
- `-std=c++17`: Modern C++ features
- `-fopenmp`: OpenMP parallelization

On Windows (MSVC):
- `/O2`: Optimization
- `/std:c++17`: C++17 standard
- `/openmp`: OpenMP support

## GPU Backend

For large systems with NVIDIA GPUs.

### Features

- **CUDA kernels**: Custom GPU kernels for Gaussian accumulation
- **Cell-list on GPU**: Neighbor search entirely on device
- **Chunked processing**: Manages GPU memory automatically
- **~5-6x faster** than C++ for large systems

### Usage

```python
wc = WillardChandler(
    u, group=atoms,
    alpha=2.4, mesh=2.0,
    surface_backend='cupy'  # GPU backend
)
```

### Requirements

- NVIDIA GPU with CUDA support (Compute Capability ≥ 6.0)
- CUDA Toolkit 11.0+
- CuPy ≥ 12.0

### Installation

See [Installation Guide](../installation.md#gpu-acceleration-optional) for details.

### Performance Characteristics

- **Memory**: GPU RAM = O(N + grid_size)
- **Best for**: >10,000 atoms
- **Speedup**: 5-6x over C++ (system dependent)
- **Overhead**: Initial data transfer to GPU

### When to Use GPU Backend

✅ **Use GPU when:**
- System has >10,000 atoms
- Processing many frames
- NVIDIA GPU available
- High throughput needed

❌ **Use C++ when:**
- System has <10,000 atoms
- Single frame calculation
- No GPU available
- CPU has many cores

### GPU Memory Management

The GPU backend automatically chunks large grids:

```python
# For very large systems, adjust chunk size
import os
os.environ['PYWC_GPU_CHUNK_SIZE'] = '5000000'  # Grid points per chunk
```

### Troubleshooting GPU

**Check GPU availability:**
```python
import cupy as cp
print(f"GPU: {cp.cuda.Device().name}")
print(f"Memory: {cp.cuda.Device().mem_info[1] / 1e9:.1f} GB")
```

**Force CPU if GPU fails:**
```python
try:
    wc = WillardChandler(u, group=atoms, surface_backend='cupy')
except Exception as e:
    print(f"GPU failed: {e}, falling back to CPU")
    wc = WillardChandler(u, group=atoms, surface_backend='cpp')
```

## Python Backend

Pure NumPy/SciPy implementation for reference and testing.

### Features

- **No compilation required**: Works immediately
- **Reference implementation**: Validates other backends
- **Cross-platform**: Works everywhere Python runs

### Usage

```python
wc = WillardChandler(
    u, group=atoms,
    alpha=2.4, mesh=2.0,
    surface_backend='python'
)
```

### When to Use

- Testing/debugging
- Platforms without C++ compiler
- Verifying results
- Understanding algorithm

### Performance

Significantly slower than C++ or GPU backends. Use only for small systems or testing.

## Backend Selection Strategy

### Automatic Selection

pyWC doesn't auto-select backends. Choose explicitly based on your needs:

```python
def choose_backend(n_atoms, has_gpu):
    if has_gpu and n_atoms > 10000:
        return 'cupy'
    else:
        return 'cpp'

# Example
backend = choose_backend(len(atoms), HAS_GPU)
wc = WillardChandler(u, group=atoms, surface_backend=backend)
```

### Benchmark Your System

Use the included benchmark tool:

```bash
pywc-compare-wc-backends
```

This will test all available backends on your system and report timings.

### From Python

```python
from pywc import WillardChandler
import time

backends = ['cpp']
try:
    import cupy
    backends.append('cupy')
except ImportError:
    pass

for backend in backends:
    wc = WillardChandler(u, group=atoms, surface_backend=backend,
                         enable_timing=True)

    start = time.time()
    wc.assign_surface()
    elapsed = time.time() - start

    print(f"{backend:10s}: {elapsed:.4f} s")
```

## Backend Comparison

### Benchmark Results

System: 10,000 water molecules, alpha=2.4, mesh=2.0

| Backend | Time (s) | Speedup | Memory (GB) |
|---------|----------|---------|-------------|
| Python  | 12.450   | 1.0x    | 0.8         |
| C++     | 0.352    | 35.4x   | 0.9         |
| GPU     | 0.063    | 197.6x  | 1.2 (GPU)   |

System: 1,000 water molecules, alpha=2.4, mesh=2.0

| Backend | Time (s) | Speedup | Memory (GB) |
|---------|----------|---------|-------------|
| Python  | 2.450    | 1.0x    | 0.2         |
| C++     | 0.122    | 20.1x   | 0.2         |
| GPU     | 0.156    | 15.7x   | 0.5 (GPU)   |

**Note**: GPU has overhead for small systems!

### Scaling Behavior

```
Time vs. Number of Atoms (mesh=2.0, alpha=2.4)

  Python: T ∝ N^1.2
  C++:    T ∝ N^1.0
  GPU:    T ∝ N^0.9 (for N > 10k)
```

## Advanced: Custom Backends

You can implement custom backends by subclassing the density evaluator:

```python
from pywc.wc_core.density import DensityEvaluator

class MyCustomBackend(DensityEvaluator):
    def evaluate(self, positions, grid, box, alpha, weights):
        # Your custom implementation
        return density_field

# Use custom backend
wc = WillardChandler(u, group=atoms, surface_backend=MyCustomBackend())
```

## See Also

- [Installation](../installation.md) - Installing backends
- [WillardChandler API](willard_chandler.md) - Main API reference
- [Examples](../examples.md) - Usage examples
