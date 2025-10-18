# Willard-Chandler Timing - Quick Reference

## Enable Timing

```python
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)
```

## Get Last Computation Time

```python
time_sec = inter.get_timing()
print(f"Last computation: {time_sec:.4f} s")
```

## Print Detailed Breakdown

```python
# Process multiple frames
for ts in u.trajectory[:10]:
    inter._assign_layers()

# Show detailed breakdown
inter.print_timing_breakdown()

# Or compact version
inter.print_timing_breakdown(verbose=False)
```

## Get Raw Statistics

```python
timings = inter.get_detailed_timings()

# Access specific component
compute_stats = timings['compute_surface']
print(f"Mean: {compute_stats['mean']:.4f} s")
print(f"Std:  {compute_stats['std']:.4f} s")
print(f"Min:  {compute_stats['min']:.4f} s")
print(f"Max:  {compute_stats['max']:.4f} s")

# Get raw data for plotting
all_times = compute_stats['times']  # List of floats (seconds)
```

## Components Timed

- `prepare_box` - Box preparation
- `define_cluster_group` - Cluster group definition
- `center` - Centering (if enabled)
- `get_positions` - Position extraction
- `compute_surface` - **Main computation** (KDE + marching cubes)
- `total` - Sum of all components

## Compare Backends

```python
# Method 1: Using compare_backends()
inter = pytim.WillardChandler(u, group=g, alpha=3.0)
results = inter.compare_backends(['cpu', 'cupy'])

# Method 2: Manual detailed comparison
cpu_inter = pytim.WillardChandler(u, group=g, surface_backend='cpu', enable_timing=True)
gpu_inter = pytim.WillardChandler(u, group=g, surface_backend='cupy', enable_timing=True)

for ts in u.trajectory[:10]:
    cpu_inter._assign_layers()
    gpu_inter._assign_layers()

cpu_inter.print_timing_breakdown()
gpu_inter.print_timing_breakdown()
```

## Typical Output

```
======================================================================
Willard-Chandler Timing Breakdown
======================================================================

compute_surface
  Mean time:  145.234 ms  ( 95.2%)
  Std dev:      2.145 ms
  Min/Max:    142.567 / 148.901 ms
  Count:             10 calls

define_cluster_group
  Mean time:    4.123 ms  (  2.7%)
  Std dev:      0.234 ms
  Min/Max:      3.890 / 4.456 ms
  Count:             10 calls

get_positions
  Mean time:    1.234 ms  (  0.8%)
  Std dev:      0.089 ms
  Min/Max:      1.123 / 1.389 ms
  Count:             10 calls

prepare_box
  Mean time:    0.234 ms  (  0.2%)
  Std dev:      0.012 ms
  Min/Max:      0.221 / 0.256 ms
  Count:             10 calls

TOTAL
  Mean time:  152.567 ms
  Total time:   1.526 s
  Frames:            10
======================================================================
```

## Expected Distribution

- `compute_surface`: 90-95% (normal)
- `define_cluster_group`: 2-5%
- `get_positions`: 1-2%
- `prepare_box`: <1%

## Optimization Tips

**If `compute_surface` > 95%:**
- ✅ Normal - this is the main computation
- Consider GPU backend for speedup
- Reduce mesh density if acceptable

**If `define_cluster_group` > 10%:**
- Check `cluster_cut` parameter
- Consider disabling cluster analysis

**If `get_positions` > 5%:**
- Large atom group - consider reducing selection

## Examples

See:
- `example_detailed_timing.py` - Complete usage example
- `DETAILED_TIMING_GUIDE.md` - Full documentation with plots
