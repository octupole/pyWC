#!/usr/bin/env python3
"""Diagnostic script to check centering optimization status."""

import time
import numpy as np

print("="*70)
print("Centering Optimization Diagnostic")
print("="*70)

# Check if modules are available
print("\n1. Checking C++ module availability:")
print("-" * 70)

try:
    from pywc._center_impl import HAS_CENTER_FAST, HAS_CENTER_FULL
    print(f"HAS_CENTER_FAST:  {HAS_CENTER_FAST}")
    print(f"HAS_CENTER_FULL:  {HAS_CENTER_FULL}")
except ImportError as e:
    print(f"ERROR importing _center_impl: {e}")
    HAS_CENTER_FAST = False
    HAS_CENTER_FULL = False

if HAS_CENTER_FAST:
    try:
        from pywc.center_fast import center_fast
        print("✓ center_fast module loaded successfully")
    except ImportError as e:
        print(f"✗ Failed to import center_fast: {e}")

if HAS_CENTER_FULL:
    try:
        from pywc.center_fast_full import center_full_optimized
        print("✓ center_fast_full module loaded successfully")
    except ImportError as e:
        print(f"✗ Failed to import center_fast_full: {e}")

# Check which implementation Interface._center will use
print("\n2. Checking Interface._center implementation:")
print("-" * 70)

try:
    from pywc.interface import Interface, HAS_CENTER_FAST as IF_HAS_FAST
    from pywc.interface import center_wrapper

    print(f"Interface sees HAS_CENTER_FAST: {IF_HAS_FAST}")
    print(f"Interface has center_wrapper: {center_wrapper is not None}")

except ImportError as e:
    print(f"ERROR: {e}")

# Test with actual data
print("\n3. Testing with sample data:")
print("-" * 70)

try:
    import MDAnalysis as mda
    import pywc
    from pywc.datafiles import MICELLE_PDB

    u = mda.Universe(MICELLE_PDB)
    g = u.select_atoms('resname DPC')

    print(f"Loaded system with {len(g)} atoms")

    # Test centering with timing
    print("\nTiming centering operation...")

    # Create interface with centering enabled
    start = time.perf_counter()
    inter = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                                  centered=True, enable_timing=True)
    elapsed = time.perf_counter() - start

    print(f"Total WillardChandler init time: {elapsed*1000:.2f} ms")

    # Get detailed timings
    timings = inter.get_detailed_timings()
    if timings:
        print("\nDetailed timing breakdown:")
        for component in ['prepare_box', 'define_cluster_group', 'center',
                         'get_positions', 'compute_surface']:
            if component in timings:
                stats = timings[component]
                print(f"  {component:25s}: {stats['mean']*1000:8.3f} ms")

    # Test multiple frames to see if optimization helps
    print("\n4. Testing across multiple frames:")
    print("-" * 70)

    frame_times = []
    for i in range(5):
        start = time.perf_counter()
        inter._assign_layers()
        frame_times.append((time.perf_counter() - start) * 1000)

    print(f"Frame times (ms): {frame_times}")
    print(f"Average: {np.mean(frame_times):.2f} ms")

    # Get breakdown again
    timings = inter.get_detailed_timings()
    if timings and 'center' in timings:
        center_stats = timings['center']
        print(f"\nCenter component:")
        print(f"  Mean:   {center_stats['mean']*1000:.3f} ms")
        print(f"  Std:    {center_stats['std']*1000:.3f} ms")
        print(f"  Calls:  {center_stats['count']}")

except Exception as e:
    print(f"ERROR during testing: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("Diagnostic complete")
print("="*70)
