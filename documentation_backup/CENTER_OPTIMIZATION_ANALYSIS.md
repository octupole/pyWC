# Center() Method Optimization Analysis

## Current Performance Issue

Based on timing analysis, the `center()` method in `interface.py` (lines 402-409) is a bottleneck. This method calls `Interface._center()` which performs iterative centering of the liquid slab.

## Current Implementation Analysis

### Call Chain
```
center() → Interface.center_system() → Interface._center()
```

### Bottleneck Location: `Interface._center()` (lines 301-373)

**Key operations:**
1. **Multiple histogram calculations** (lines 344-345, 361-362)
2. **Iterative shifting** in while loop (lines 352-362)
3. **Multiple array operations** (lines 336-338, 369-373)
4. **Repeated `_attempt_shift()` calls** (line 358-359)

### Performance Hotspots

```python
# Called in EVERY iteration of while loop:
while (histo[0] > delta or histo[-1] > delta):
    total_shift += shift
    _pos_group += shift  # Array operation

    Interface._attempt_shift(...)  # Calls centerbox()

    histo, _ = np.histogram(      # EXPENSIVE
        _pos_group, bins=10, range=_range, density=True)
```

**Problems:**
1. **Histogram calculation in tight loop** - O(N) per iteration
2. **`centerbox()` calls** - Modifies full universe positions
3. **Multiple array copies** - `get_x()`, `get_y()`, `get_z()` create copies
4. **Column stacking** - Rebuilds position array (line 373)

### Typical Iteration Count
- Usually converges in 5-20 iterations
- Each iteration: histogram (O(N)) + array ops (O(N))
- Total: O(N × iterations) ≈ O(10N) to O(20N)

## Optimization Strategies

### Strategy 1: **Pure NumPy Optimization** (Easiest, 2-3x speedup)

**Changes:**
- Pre-allocate arrays
- Avoid repeated copies
- Cache box dimensions
- Vectorize shift operations

**Implementation:**
```python
@staticmethod
def _center_optimized(group, direction, halfbox_shift=False):
    """Optimized pure NumPy version"""
    dim = group.universe.coord.dimensions
    _dir_idx = {'x': 0, 'y': 1, 'z': 2}[direction]

    # Pre-allocate and cache
    positions = group.universe.atoms.positions
    n_atoms = len(group)

    if halfbox_shift:
        _range = (-dim[_dir_idx] / 2., dim[_dir_idx] / 2.)
    else:
        _range = (0., dim[_dir_idx])

    shift = dim[_dir_idx] / 100.
    total_shift = 0.0

    # Get positions once
    _pos_group = group.positions[:, _dir_idx].copy()

    # Pre-compute histogram edges
    bins = 10
    bin_edges = np.linspace(_range[0], _range[1], bins + 1)

    # Initial histogram
    histo = np.histogram(_pos_group, bins=bin_edges, density=True)[0]
    max_val, min_val = histo.max(), histo.min()
    delta = min_val + (max_val - min_val) / 3.

    # Iterative centering
    max_iter = 100
    for _ in range(max_iter):
        if not (histo[0] > delta or histo[-1] > delta):
            break

        total_shift += shift
        if total_shift >= dim[_dir_idx]:
            raise ValueError("Centering failure")

        # Shift in-place
        _pos_group += shift

        # Recompute histogram (still expensive but unavoidable)
        histo = np.histogram(_pos_group, bins=bin_edges, density=True)[0]

    # Final centering
    _center_ = _pos_group.mean()
    box_half = 0. if halfbox_shift else dim[_dir_idx] / 2.
    final_shift = total_shift - _center_ + box_half

    # Apply shift to full positions array
    positions[:, _dir_idx] += final_shift
    group.universe.atoms.positions = positions
```

**Speedup:** ~2-3x (avoids copies, streamlines logic)

---

### Strategy 2: **Numba JIT Compilation** (Medium effort, 5-10x speedup)

Use Numba to JIT-compile the histogram and shift logic.

**Implementation:**
```python
import numba

@numba.jit(nopython=True)
def _compute_histogram_1d(data, bins, range_min, range_max):
    """Fast 1D histogram computation"""
    counts = np.zeros(bins, dtype=np.int64)
    bin_width = (range_max - range_min) / bins

    for val in data:
        if range_min <= val < range_max:
            bin_idx = int((val - range_min) / bin_width)
            if 0 <= bin_idx < bins:
                counts[bin_idx] += 1

    # Normalize to density
    total = counts.sum()
    if total > 0:
        density = counts.astype(np.float64) / (total * bin_width)
    else:
        density = np.zeros(bins, dtype=np.float64)

    return density

@numba.jit(nopython=True)
def _center_loop_numba(pos_group, shift, range_min, range_max, delta, max_shift):
    """Numba-accelerated centering loop"""
    total_shift = 0.0
    bins = 10

    for iteration in range(100):
        histo = _compute_histogram_1d(pos_group, bins, range_min, range_max)

        if not (histo[0] > delta or histo[-1] > delta):
            break

        total_shift += shift
        if total_shift >= max_shift:
            return pos_group, -1.0  # Signal failure

        # Shift positions
        for i in range(len(pos_group)):
            pos_group[i] += shift

    return pos_group, total_shift
```

**Speedup:** ~5-10x (JIT compilation eliminates Python overhead)

---

### Strategy 3: **C++/Cython Extension** (Most effort, 10-20x speedup)

Create a compiled extension for the entire centering operation.

**File: `pytim/center_fast.pyx`**
```cython
# cython: language_level=3, boundscheck=False, wraparound=False
import numpy as np
cimport numpy as cnp
cimport cython
from libc.math cimport fabs

@cython.boundscheck(False)
@cython.wraparound(False)
def center_fast(cnp.ndarray[cnp.float64_t, ndim=1] pos_group,
                double range_min,
                double range_max,
                double delta,
                double shift,
                double max_shift,
                int bins=10):
    """
    Fast centering using Cython.

    Parameters
    ----------
    pos_group : ndarray (N,)
        Positions along centering direction
    range_min, range_max : float
        Histogram range
    delta : float
        Density threshold
    shift : float
        Shift increment
    max_shift : float
        Maximum allowed shift
    bins : int
        Number of histogram bins

    Returns
    -------
    pos_group : ndarray
        Shifted positions
    total_shift : float
        Total applied shift
    """
    cdef:
        int i, j, bin_idx
        int n = pos_group.shape[0]
        double total_shift = 0.0
        double bin_width = (range_max - range_min) / bins
        double val
        cnp.ndarray[cnp.int64_t, ndim=1] counts = np.zeros(bins, dtype=np.int64)
        cnp.ndarray[cnp.float64_t, ndim=1] density = np.zeros(bins, dtype=np.float64)
        int total_count
        int max_iter = 100

    for iteration in range(max_iter):
        # Reset counts
        for i in range(bins):
            counts[i] = 0

        # Compute histogram
        for i in range(n):
            val = pos_group[i]
            if range_min <= val < range_max:
                bin_idx = <int>((val - range_min) / bin_width)
                if 0 <= bin_idx < bins:
                    counts[bin_idx] += 1

        # Normalize to density
        total_count = 0
        for i in range(bins):
            total_count += counts[i]

        if total_count > 0:
            for i in range(bins):
                density[i] = <double>counts[i] / (total_count * bin_width)
        else:
            for i in range(bins):
                density[i] = 0.0

        # Check convergence
        if not (density[0] > delta or density[bins-1] > delta):
            break

        # Apply shift
        total_shift += shift
        if total_shift >= max_shift:
            raise ValueError("Centering failure: maximum shift exceeded")

        for i in range(n):
            pos_group[i] += shift

    return pos_group, total_shift
```

**Build configuration in `setup.py`:**
```python
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "pytim.center_fast",
        ["pytim/center_fast.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=['-O3', '-march=native']
    )
]

setup(
    ext_modules=cythonize(ext_modules, compiler_directives={
        'language_level': "3",
        'boundscheck': False,
        'wraparound': False
    })
)
```

**Speedup:** ~10-20x (compiled code, no Python overhead, cache-friendly)

---

### Strategy 4: **Hybrid OpenMP Parallelization** (Expert, 20-50x speedup)

Add OpenMP parallelization to Cython implementation.

**File: `pytim/center_fast_omp.pyx`**
```cython
# cython: language_level=3
from cython.parallel import prange
cimport openmp

@cython.boundscheck(False)
@cython.wraparound(False)
def center_fast_omp(cnp.ndarray[cnp.float64_t, ndim=1] pos_group,
                    double range_min,
                    double range_max,
                    double delta,
                    double shift,
                    double max_shift,
                    int bins=10,
                    int num_threads=4):
    """Parallel version using OpenMP"""
    cdef:
        int i, bin_idx, thread_id
        int n = pos_group.shape[0]
        double bin_width = (range_max - range_min) / bins
        cnp.ndarray[cnp.int64_t, ndim=2] thread_counts = np.zeros((num_threads, bins), dtype=np.int64)
        cnp.ndarray[cnp.int64_t, ndim=1] counts = np.zeros(bins, dtype=np.int64)

    openmp.omp_set_num_threads(num_threads)

    # Parallel histogram computation
    with nogil:
        for i in prange(n, schedule='static'):
            thread_id = openmp.omp_get_thread_num()
            val = pos_group[i]
            if range_min <= val < range_max:
                bin_idx = <int>((val - range_min) / bin_width)
                if 0 <= bin_idx < bins:
                    thread_counts[thread_id, bin_idx] += 1

    # Reduce thread counts
    for thread_id in range(num_threads):
        for i in range(bins):
            counts[i] += thread_counts[thread_id, i]

    # ... rest of algorithm
```

**Speedup:** ~20-50x on multi-core systems

---

## Recommendation

### **Recommended Approach: Strategy 3 (Cython)**

**Reasons:**
1. **Best effort/benefit ratio**
   - 10-20x speedup achievable
   - Fits existing codebase (already uses Cython for other modules)
   - No new dependencies (Numba would add dependency)

2. **Maintainability**
   - Cython code is readable
   - Easy to debug
   - Type hints help catch errors

3. **Compatibility**
   - Works on all platforms
   - No runtime JIT compilation
   - Predictable performance

### **Implementation Plan**

1. **Create `pytim/center_fast.pyx`** with Cython implementation
2. **Modify `interface.py`**:
   ```python
   try:
       from .center_fast import center_fast
       HAS_CENTER_FAST = True
   except ImportError:
       HAS_CENTER_FAST = False

   @staticmethod
   def _center(group, direction, halfbox_shift=False):
       if HAS_CENTER_FAST:
           return Interface._center_fast(group, direction, halfbox_shift)
       else:
           return Interface._center_python(group, direction, halfbox_shift)
   ```

3. **Update `setup.py`** to compile extension
4. **Add fallback** for systems without Cython compiler

### **Expected Results**

| Implementation | Speedup | Effort | Dependencies |
|---------------|---------|--------|--------------|
| Current | 1x | - | - |
| NumPy Optimized | 2-3x | Low | None |
| Numba JIT | 5-10x | Medium | Numba |
| **Cython** | **10-20x** | **Medium** | **Cython (build-time)** |
| Cython+OpenMP | 20-50x | High | Cython, OpenMP |

### **Timeline**

- **NumPy optimization**: 2-4 hours
- **Cython implementation**: 1-2 days
- **Testing & integration**: 1 day
- **Total**: ~3-4 days for production-ready solution

---

## Alternative: Algorithm Improvement

Instead of optimizing the implementation, **reconsider if centering is always necessary**:

### **Option A: Cache-aware centering**
Only re-center when box dimensions change significantly:
```python
def center(self, force=False):
    if force or self._box_changed():
        Interface.center_system(...)
        self._last_box = self.universe.dimensions.copy()
```

### **Option B: User-controlled centering**
Make centering optional with a warning:
```python
WillardChandler(u, group=g, alpha=3.0, centered=False)  # Skip expensive centering
```

### **Option C: Approximate centering**
Use faster, approximate method:
- Center based on COM instead of density profile
- Use coarser histogram (5 bins instead of 10)
- Reduce convergence threshold

---

## Benchmarking Plan

To validate optimizations:

```python
import time
import numpy as np

# Generate test case
n_atoms = 10000
positions = np.random.rand(n_atoms, 3) * 50

# Benchmark original
start = time.perf_counter()
for _ in range(100):
    Interface._center(group, 'z', halfbox_shift=True)
time_original = time.perf_counter() - start

# Benchmark optimized
start = time.perf_counter()
for _ in range(100):
    Interface._center_optimized(group, 'z', halfbox_shift=True)
time_optimized = time.perf_counter() - start

speedup = time_original / time_optimized
print(f"Speedup: {speedup:.1f}x")
```

---

## Conclusion

The `center()` method can be significantly optimized:

1. **Quick win**: NumPy optimization (2-3x) - implement today
2. **Best overall**: Cython implementation (10-20x) - implement this week
3. **Future**: OpenMP parallelization (20-50x) - if still needed

The Cython approach is **strongly recommended** as it:
- Provides excellent speedup (10-20x)
- Integrates well with existing codebase
- Has acceptable implementation effort
- Maintains code readability

Would you like me to implement the Cython optimization?
