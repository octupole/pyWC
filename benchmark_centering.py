#!/usr/bin/env python3
"""Benchmark centering implementations."""

import time
import numpy as np
import MDAnalysis as mda
import pytim
from pytim.datafiles import MICELLE_PDB
from pytim._center_impl import HAS_CENTER_FAST, HAS_CENTER_FULL

print("=" * 70)
print("Centering Performance Benchmark")
print("=" * 70)

print(f"\nOptimization status:")
print(f"  HAS_CENTER_FAST: {HAS_CENTER_FAST}")
print(f"  HAS_CENTER_FULL: {HAS_CENTER_FULL}")

# Load test system
u = mda.Universe(MICELLE_PDB)
g = u.select_atoms('resname DPC')

print(f"\nTest system: {len(g)} atoms in centering group")
print(f"Total atoms: {len(u.atoms)}")

# Warm-up
inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0, centered=True)

# Test with C++ optimization (current default)
print("\n" + "-" * 70)
print("Testing OPTIMIZED C++ implementation (default)")
print("-" * 70)

cpp_times = []
for i in range(10):
    u.trajectory[0]  # Reset to first frame
    start = time.perf_counter()
    inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                                   centered=True, enable_timing=True)
    elapsed = time.perf_counter() - start
    cpp_times.append(elapsed * 1000)

    # Get center timing
    timings = inter.get_detailed_timings()
    if i == 0 and timings and 'center' in timings:
        print(f"First run - Center component: {timings['center']['mean']*1000:.3f} ms")

print(f"\nC++ Results (10 runs):")
print(f"  Mean total time: {np.mean(cpp_times):.2f} ms")
print(f"  Std dev:         {np.std(cpp_times):.2f} ms")
print(f"  Min:             {np.min(cpp_times):.2f} ms")
print(f"  Max:             {np.max(cpp_times):.2f} ms")

# Get detailed breakdown
timings = inter.get_detailed_timings()
if timings:
    print(f"\n  Component breakdown:")
    for component in ['prepare_box', 'define_cluster_group', 'center',
                     'get_positions', 'compute_surface']:
        if component in timings:
            stats = timings[component]
            print(f"    {component:25s}: {stats['mean']*1000:8.3f} ms")

# Test with Python fallback
print("\n" + "-" * 70)
print("Testing PYTHON fallback implementation")
print("-" * 70)

# Temporarily disable C++ by using force_python in _center
from pytim import interface
original_center = interface.Interface._center

def python_center_wrapper(group, direction, halfbox_shift=False):
    from pytim._center_impl import _center_python
    return _center_python(group, direction, halfbox_shift)

interface.Interface._center = staticmethod(python_center_wrapper)

python_times = []
for i in range(10):
    u.trajectory[0]  # Reset to first frame
    start = time.perf_counter()
    inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                                   centered=True, enable_timing=True)
    elapsed = time.perf_counter() - start
    python_times.append(elapsed * 1000)

print(f"\nPython Results (10 runs):")
print(f"  Mean total time: {np.mean(python_times):.2f} ms")
print(f"  Std dev:         {np.std(python_times):.2f} ms")
print(f"  Min:             {np.min(python_times):.2f} ms")
print(f"  Max:             {np.max(python_times):.2f} ms")

# Restore original
interface.Interface._center = original_center

# Comparison
print("\n" + "=" * 70)
print("COMPARISON")
print("=" * 70)

speedup = np.mean(python_times) / np.mean(cpp_times)
time_saved = np.mean(python_times) - np.mean(cpp_times)

print(f"\nOverall WillardChandler initialization:")
print(f"  Python mean:     {np.mean(python_times):.2f} ms")
print(f"  C++ mean:        {np.mean(cpp_times):.2f} ms")
print(f"  Speedup:         {speedup:.2f}x")
print(f"  Time saved:      {time_saved:.2f} ms per call")

print("\n" + "=" * 70)
