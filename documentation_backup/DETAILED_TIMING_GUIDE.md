# Detailed Timing Guide for Willard-Chandler

## Overview

The WillardChandler class now includes detailed timing capabilities for debugging and performance analysis. When `enable_timing=True`, the code tracks the time spent in each component of the surface calculation.

## Components Timed

When timing is enabled, the following components are tracked individually:

1. **prepare_box** - Box preparation and validation
2. **define_cluster_group** - Cluster group definition
3. **center** - Centering of atom group (if enabled)
4. **get_positions** - Position extraction and filtering
5. **compute_surface** - Main surface computation (KDE + marching cubes)
6. **total** - Sum of all components

## Usage

### Basic Timing

```python
import MDAnalysis as mda
import pytim

u = mda.Universe('topology.tpr', 'trajectory.xtc')
g = u.select_atoms('resname POPC')

# Enable timing
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)

# Get last computation time
print(f"Computation time: {inter.get_timing():.4f} s")
```

### Detailed Breakdown

```python
# Process multiple frames
for ts in u.trajectory[:10]:
    inter._assign_layers()

# Print formatted breakdown
inter.print_timing_breakdown()
```

**Output:**
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

### Programmatic Access

```python
# Get raw timing statistics
timings = inter.get_detailed_timings()

# Access specific component
compute_surface_stats = timings['compute_surface']
print(f"Mean: {compute_surface_stats['mean']:.4f} s")
print(f"Std:  {compute_surface_stats['std']:.4f} s")
print(f"Min:  {compute_surface_stats['min']:.4f} s")
print(f"Max:  {compute_surface_stats['max']:.4f} s")

# Get raw data for plotting
times = compute_surface_stats['times']  # List of all measurements
```

## Advanced Usage

### Backend Comparison with Details

```python
# Compare backends
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)

# First, do backend comparison
results = inter.compare_backends(['cpu', 'cupy'])

# Then get detailed breakdown for each
for backend in ['cpu', 'cupy']:
    inter_temp = pytim.WillardChandler(
        u, group=g, alpha=3.0,
        surface_backend=backend,
        enable_timing=True
    )

    # Process frames
    for ts in u.trajectory[:10]:
        inter_temp._assign_layers()

    print(f"\n{backend.upper()} Backend:")
    inter_temp.print_timing_breakdown(verbose=False)
```

### Plotting Timing Data

```python
import matplotlib.pyplot as plt
import numpy as np

# Collect timing data
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)
for ts in u.trajectory[:50]:
    inter._assign_layers()

timings = inter.get_detailed_timings()

# Plot time series
fig, ax = plt.subplots(figsize=(10, 6))

for component in ['prepare_box', 'define_cluster_group', 'compute_surface']:
    times = np.array(timings[component]['times']) * 1000  # Convert to ms
    ax.plot(times, label=component, marker='o', markersize=3)

ax.set_xlabel('Frame')
ax.set_ylabel('Time (ms)')
ax.set_title('Willard-Chandler Component Timing')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('wc_timing.png', dpi=150)
```

### Stacked Bar Chart

```python
import matplotlib.pyplot as plt
import numpy as np

timings = inter.get_detailed_timings()

# Prepare data
components = ['prepare_box', 'define_cluster_group', 'get_positions', 'compute_surface']
means = [timings[c]['mean'] * 1000 for c in components]  # ms
stds = [timings[c]['std'] * 1000 for c in components]

# Create stacked bar
fig, ax = plt.subplots(figsize=(8, 6))
colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(components)))

bottom = 0
for i, (component, mean, std) in enumerate(zip(components, means, stds)):
    ax.bar(0, mean, bottom=bottom, color=colors[i], label=component,
           width=0.5, yerr=std, capsize=5)
    bottom += mean

ax.set_ylabel('Time (ms)')
ax.set_title('Willard-Chandler Timing Breakdown')
ax.set_xticks([])
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig('wc_timing_stacked.png', dpi=150)
```

### Identify Bottlenecks

```python
timings = inter.get_detailed_timings()

# Calculate percentages
total = sum(t['mean'] for c, t in timings.items() if c != 'total')
percentages = {c: (t['mean']/total*100) for c, t in timings.items() if c != 'total'}

# Sort by percentage
sorted_components = sorted(percentages.items(), key=lambda x: x[1], reverse=True)

print("\nBottleneck Analysis:")
for component, pct in sorted_components:
    mean_ms = timings[component]['mean'] * 1000
    print(f"{component:25s}: {pct:5.1f}% ({mean_ms:7.2f} ms)")

# Identify if optimization is worthwhile
if percentages['compute_surface'] > 90:
    print("\n⚠ compute_surface dominates (>90%) - consider GPU backend")
elif percentages['define_cluster_group'] > 10:
    print("\n⚠ define_cluster_group is significant (>10%) - check cluster settings")
```

## Interpreting Results

### Typical Time Distribution

For a well-optimized calculation:
- **compute_surface**: 90-95% (expected bottleneck)
- **define_cluster_group**: 2-5%
- **get_positions**: 1-2%
- **prepare_box**: <1%

### Optimization Strategies

**If compute_surface dominates (>95%):**
- ✅ Normal - this is the main computation
- Consider GPU backend for speedup
- Reduce mesh density if acceptable
- Increase alpha (smooths surface, faster)

**If define_cluster_group is high (>10%):**
- Check `cluster_cut` parameter
- Consider disabling cluster analysis if not needed
- Reduce `extra_cluster_groups` if applicable

**If get_positions is high (>5%):**
- May indicate large atom groups
- Consider reducing selection size
- Check if `include_zero_radius` is needed

**If prepare_box is high (>5%):**
- Unusual - may indicate trajectory I/O issues
- Check if box dimensions are changing

## Performance Tips

### 1. Minimize Overhead
```python
# Good: Enable timing only when needed
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)

# Bad: Timing adds ~1-2% overhead
# Don't enable for production runs unless debugging
```

### 2. Sample Representative Frames
```python
# Don't need to time every frame
for i, ts in enumerate(u.trajectory):
    if i % 10 == 0:  # Sample every 10th frame
        inter._assign_layers()
```

### 3. Compare Backend Performance
```python
# Quantify GPU speedup with breakdown
cpu_inter = pytim.WillardChandler(u, group=g, surface_backend='cpu', enable_timing=True)
gpu_inter = pytim.WillardChandler(u, group=g, surface_backend='cupy', enable_timing=True)

for ts in u.trajectory[:10]:
    cpu_inter._assign_layers()
    gpu_inter._assign_layers()

cpu_timings = cpu_inter.get_detailed_timings()
gpu_timings = gpu_inter.get_detailed_timings()

speedup = cpu_timings['compute_surface']['mean'] / gpu_timings['compute_surface']['mean']
print(f"GPU speedup: {speedup:.2f}x")
```

## API Reference

### Methods

#### `get_timing()`
Returns the computation time of the last surface calculation.

**Returns:** `float` or `None` - Time in seconds

#### `get_detailed_timings()`
Returns detailed statistics for all components.

**Returns:** `dict` with structure:
```python
{
    'component_name': {
        'mean': float,    # Average time (seconds)
        'std': float,     # Standard deviation (seconds)
        'min': float,     # Minimum time (seconds)
        'max': float,     # Maximum time (seconds)
        'total': float,   # Total accumulated time (seconds)
        'count': int,     # Number of measurements
        'times': list     # Raw measurements (seconds)
    }
}
```

#### `print_timing_breakdown(verbose=True)`
Prints formatted timing breakdown.

**Parameters:**
- `verbose` (bool): Include detailed statistics (default: True)

**Returns:** None (prints to stdout)

## Examples

### Example 1: Quick Performance Check
```python
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)
for ts in u.trajectory[:5]:
    inter._assign_layers()
inter.print_timing_breakdown(verbose=False)
```

### Example 2: Export Timing Data
```python
import json

inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)
for ts in u.trajectory[:100]:
    inter._assign_layers()

timings = inter.get_detailed_timings()

# Remove raw times for cleaner JSON
export_data = {
    comp: {k: v for k, v in stats.items() if k != 'times'}
    for comp, stats in timings.items()
}

with open('wc_timings.json', 'w') as f:
    json.dump(export_data, f, indent=2)
```

### Example 3: Trajectory Analysis
```python
# Analyze timing across trajectory
inter = pytim.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)

frame_times = []
for ts in u.trajectory:
    inter._assign_layers()
    frame_times.append(inter.get_timing())

# Identify slow frames
import numpy as np
mean_time = np.mean(frame_times)
std_time = np.std(frame_times)
slow_frames = [i for i, t in enumerate(frame_times) if t > mean_time + 2*std_time]

print(f"Slow frames (>2σ): {slow_frames}")
```

## Troubleshooting

### No timing data
```python
timings = inter.get_detailed_timings()
if timings is None:
    print("Enable timing with enable_timing=True")
```

### High variability
- Check if trajectory has varying box sizes
- Ensure consistent frame processing
- May indicate system memory pressure

### Unexpected component times
- Use `verbose=True` to see min/max ranges
- Check for outliers in raw data
- Consider warmup frames (JIT compilation)

## Summary

The detailed timing feature helps you:
- ✅ Identify performance bottlenecks
- ✅ Compare CPU vs GPU effectiveness
- ✅ Optimize parameter choices
- ✅ Debug slow calculations
- ✅ Monitor performance across trajectories

For production runs, disable timing (`enable_timing=False`, default) to minimize overhead.
