#!/usr/bin/env python3
"""Test if there's a warm-up effect."""

import MDAnalysis as mda
import pytim
import time
import numpy as np

print("=" * 70)
print("Warm-up Effect Test")
print("=" * 70)

# Load system
u = mda.Universe('/home/marchi/pyTimTest/npt_run.tpr',
                 '/home/marchi/pyTimTest/traj_npt.xtc')

with open('/home/marchi/pyTimTest/selection.txt', 'r') as f:
    selection = f.read().strip()

g = u.select_atoms(selection)
print(f"\nSystem: {len(g)} atoms in group, {len(u.atoms)} total")

# Test multiple runs
print("\nRunning 10 iterations with C++ optimization...")
cpp_times = []
for i in range(10):
    u.trajectory[0]
    start = time.perf_counter()
    inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.5,
                                  centered=True, enable_timing=True)
    elapsed = time.perf_counter() - start
    cpp_times.append(elapsed * 1000)
    print(f"  Run {i+1:2d}: {elapsed*1000:7.2f} ms", end="")

    # Get center time
    timings = inter.get_detailed_timings()
    if timings and 'center' in timings:
        center_time = timings['center']['mean'] * 1000
        print(f"  (center: {center_time:6.2f} ms)")
    else:
        print()

print(f"\nC++ Results:")
print(f"  First run:  {cpp_times[0]:.2f} ms")
print(f"  Mean:       {np.mean(cpp_times):.2f} ms")
print(f"  Std:        {np.std(cpp_times):.2f} ms")
print(f"  Mean (2-10):{np.mean(cpp_times[1:]):.2f} ms")

# Now test Python
print("\n" + "-" * 70)
print("Running 10 iterations with Python fallback...")

from pytim import interface
original_center = interface.Interface._center

def python_center_wrapper(group, direction, halfbox_shift=False):
    from pytim._center_impl import _center_python
    return _center_python(group, direction, halfbox_shift)

interface.Interface._center = staticmethod(python_center_wrapper)

py_times = []
for i in range(10):
    u.trajectory[0]
    start = time.perf_counter()
    inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.5,
                                  centered=True, enable_timing=True)
    elapsed = time.perf_counter() - start
    py_times.append(elapsed * 1000)
    print(f"  Run {i+1:2d}: {elapsed*1000:7.2f} ms", end="")

    timings = inter.get_detailed_timings()
    if timings and 'center' in timings:
        center_time = timings['center']['mean'] * 1000
        print(f"  (center: {center_time:6.2f} ms)")
    else:
        print()

interface.Interface._center = original_center

print(f"\nPython Results:")
print(f"  First run:  {py_times[0]:.2f} ms")
print(f"  Mean:       {np.mean(py_times):.2f} ms")
print(f"  Std:        {np.std(py_times):.2f} ms")
print(f"  Mean (2-10):{np.mean(py_times[1:]):.2f} ms")

# Compare
print("\n" + "=" * 70)
print("COMPARISON (excluding first run)")
print("=" * 70)
print(f"  Python mean:  {np.mean(py_times[1:]):.2f} ms")
print(f"  C++ mean:     {np.mean(cpp_times[1:]):.2f} ms")
print(f"  Speedup:      {np.mean(py_times[1:])/np.mean(cpp_times[1:]):.2f}x")
print("=" * 70)
