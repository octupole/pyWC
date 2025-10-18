# WillardChandler Backend Timing and Comparison Features

## Summary of Changes

This document describes the modifications made to `pytim/willard_chandler.py` to add backend comparison capabilities with timing measurements.

## New Features

### 1. **Timing Support**

The WillardChandler class now tracks computation time for surface calculations.

#### New Parameters:
- `enable_timing` (bool): Enable timing measurements (default: False)

#### Usage:
```python
import pytim
import MDAnalysis as mda

u = mda.Universe(MICELLE_PDB)
g = u.select_atoms('resname DPC')

inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)
print(f"Computation time: {inter.get_timing():.4f} s")
```

### 2. **Backend Selection**

You can now explicitly choose between CPU and GPU backends.

#### Existing Parameters (now documented):
- `surface_backend` (str): 'cpu' or 'cupy' (default: 'cpu')
- `surface_backend_options` (dict): Backend-specific options

#### Usage:
```python
# CPU backend
inter_cpu = pytim.WillardChandler(u, group=g, alpha=3.0, surface_backend='cpu')

# GPU backend with options
inter_gpu = pytim.WillardChandler(
    u, group=g, alpha=3.0,
    surface_backend='cupy',
    surface_backend_options={
        'chunk_size': 4096,
        'max_chunk_mem': 256 * 1024 * 1024
    }
)
```

### 3. **Backend Comparison**

New method to compare performance of different backends.

#### New Method: `compare_backends()`

```python
def compare_backends(self, backends=['cpu', 'cupy'], backend_options=None, verbose=True)
```

**Parameters:**
- `backends` (list): List of backend names to compare
- `backend_options` (dict): Dictionary mapping backend names to their options
- `verbose` (bool): Print comparison results

**Returns:**
- Dictionary with timing results and status for each backend

#### Usage:
```python
inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0)

# Compare backends with default options
results = inter.compare_backends(['cpu', 'cupy'])

# Output:
# ============================================================
# Backend Performance Comparison
# ============================================================
# cpu         : 0.1234 s
# cupy        : 0.0456 s (2.71x)
# ============================================================

# Compare with custom options for each backend
results = inter.compare_backends(
    backends=['cpu', 'cupy'],
    backend_options={
        'cupy': {'chunk_size': 8192}
    }
)

# Access detailed results
for backend, result in results.items():
    if result['success']:
        print(f"{backend}: {result['time']:.4f} s")
    else:
        print(f"{backend} failed: {result['error']}")
```

### 4. **Get Timing Method**

New method to retrieve the last computation time.

#### New Method: `get_timing()`

```python
def get_timing(self)
```

**Returns:**
- Float: Time in seconds, or None if timing was not enabled

## Implementation Details

### Internal Attributes Added:
- `_surface_computation_time`: Stores the last computation time
- `_backend_timings`: Dictionary storing times for different backends
- `_enable_timing`: Flag to enable/disable timing

### Modified Methods:
- `__init__()`: Added `enable_timing` support
- `_assign_layers()`: Added timing measurements using `time.perf_counter()`

### New Imports:
- `import time`: For high-resolution timing

## Code Changes

### File: `pytim/willard_chandler.py`

1. **Import statement** (line 10):
   ```python
   import time
   ```

2. **Class attributes** (lines 86-88):
   ```python
   _surface = None
   _surface_computation_time = None
   _backend_timings = {}
   ```

3. **Updated `__init__` parameters** (lines 140-141):
   ```python
   self.surface_backend_options = surface_backend_options or {}
   self._enable_timing = kargs.get('enable_timing', False)
   ```

4. **Timing in `_assign_layers()`** (lines 250-263):
   ```python
   start_time = time.perf_counter()
   result = compute_surface(...)
   elapsed_time = time.perf_counter() - start_time

   if self._enable_timing:
       self._surface_computation_time = elapsed_time
       self._backend_timings[self.surface_backend] = elapsed_time
   ```

5. **New methods** (lines 225-333):
   - `get_timing()`
   - `compare_backends()`

6. **Updated docstring** (lines 63-65):
   - Documented `surface_backend`, `surface_backend_options`, `enable_timing`

## Testing

### Syntax Validation:
```bash
python -m py_compile pytim/willard_chandler.py
```
✓ No syntax errors

### Example Usage:
See `example_backend_usage.py` for complete examples.

## Compatibility

- **Backward Compatible**: All changes are additive; existing code will work unchanged
- **Optional Features**: Timing and comparison features are opt-in
- **Graceful Degradation**: GPU backend fails gracefully if CuPy is not installed

## Performance Considerations

- Timing uses `time.perf_counter()` for high-resolution measurements
- Minimal overhead when `enable_timing=False` (default)
- `compare_backends()` runs each backend independently for fair comparison
- Original backend configuration is restored after comparison

## Future Enhancements

Potential future additions:
- Warmup runs for more accurate GPU timing
- Memory usage tracking
- Support for additional backends
- Automated performance reporting
