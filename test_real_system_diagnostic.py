#!/usr/bin/env python3
"""Diagnostic for real system centering."""

import MDAnalysis as mda
import pywc
import time
import sys

# Check optimization status
from pywc._center_impl import HAS_CENTER_FAST, HAS_CENTER_FULL

print("=" * 70)
print("Real System Centering Diagnostic")
print("=" * 70)
print(f"\nOptimization status:")
print(f"  HAS_CENTER_FAST: {HAS_CENTER_FAST}")
print(f"  HAS_CENTER_FULL: {HAS_CENTER_FULL}")

# Load real system
print("\nLoading system...")
u = mda.Universe('/home/marchi/pyTimTest/npt_run.tpr',
                 '/home/marchi/pyTimTest/traj_npt.xtc')

# Read selection
with open('/home/marchi/pyTimTest/selection.txt', 'r') as f:
    selection = f.read().strip()

print(f"Selection: {selection[:50]}...")
g = u.select_atoms(selection)
print(f"Selected {len(g)} atoms (total system: {len(u.atoms)} atoms)")

# Test with C++ (current default)
print("\n" + "-" * 70)
print("Testing with C++ optimization (current default)")
print("-" * 70)

u.trajectory[0]  # Go to first frame
start = time.perf_counter()
inter = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.5,
                              centered=True, enable_timing=True)
elapsed = time.perf_counter() - start

print(f"Total time: {elapsed*1000:.2f} ms")

timings = inter.get_detailed_timings()
if timings:
    print("\nComponent breakdown:")
    for component in ['prepare_box', 'define_cluster_group', 'center',
                     'get_positions', 'compute_surface']:
        if component in timings:
            stats = timings[component]
            print(f"  {component:25s}: {stats['mean']*1000:8.3f} ms ({stats['count']} calls)")

# Now test with forced Python version
print("\n" + "-" * 70)
print("Testing with PYTHON fallback (forced)")
print("-" * 70)

from pywc import interface
original_center = interface.Interface._center

def python_center_wrapper(group, direction, halfbox_shift=False):
    from pywc._center_impl import _center_python
    return _center_python(group, direction, halfbox_shift)

interface.Interface._center = staticmethod(python_center_wrapper)

u.trajectory[0]  # Reset
start = time.perf_counter()
inter_py = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.5,
                                 centered=True, enable_timing=True)
elapsed_py = time.perf_counter() - start

print(f"Total time: {elapsed_py*1000:.2f} ms")

timings_py = inter_py.get_detailed_timings()
if timings_py:
    print("\nComponent breakdown:")
    for component in ['prepare_box', 'define_cluster_group', 'center',
                     'get_positions', 'compute_surface']:
        if component in timings_py:
            stats = timings_py[component]
            print(f"  {component:25s}: {stats['mean']*1000:8.3f} ms ({stats['count']} calls)")

# Restore
interface.Interface._center = original_center

# Compare
print("\n" + "=" * 70)
print("COMPARISON")
print("=" * 70)

if timings and timings_py and 'center' in timings and 'center' in timings_py:
    cpp_center = timings['center']['mean'] * 1000
    py_center = timings_py['center']['mean'] * 1000
    speedup = py_center / cpp_center

    print(f"\nCenter component:")
    print(f"  Python:    {py_center:.2f} ms")
    print(f"  C++:       {cpp_center:.2f} ms")
    print(f"  Speedup:   {speedup:.2f}x")
    print(f"  Saved:     {py_center - cpp_center:.2f} ms")
else:
    print("\nCould not compare - timing data incomplete")

print(f"\nOverall:")
print(f"  Python:    {elapsed_py*1000:.2f} ms")
print(f"  C++:       {elapsed*1000:.2f} ms")
print(f"  Speedup:   {elapsed_py/elapsed:.2f}x")
print(f"  Saved:     {(elapsed_py - elapsed)*1000:.2f} ms")

print("\n" + "=" * 70)
