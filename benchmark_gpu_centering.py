#!/usr/bin/env python
"""
Benchmark script comparing CPU vs GPU centering performance.

This script compares the performance of:
1. C++ CPU centering (center_fast_full)
2. GPU centering (CuPy CUDA kernels)
"""

import time
import numpy as np
import MDAnalysis as mda
import pytim
from pytim._center_impl import HAS_CENTER_GPU, HAS_CENTER_FULL

def benchmark_centering(u, g, backend, n_frames=100):
    """Benchmark centering for a given backend."""
    print(f"\n{'='*70}")
    print(f"Benchmarking {backend.upper()} backend...")
    print(f"{'='*70}")

    # Warm-up (3 frames)
    for _ in range(3):
        u.trajectory[0]
        inter = pytim.WillardChandler(
            u,
            group=g,
            alpha=3.0,
            mesh=2.5,
            centered=True,
            centering_backend=backend,
            enable_timing=False
        )

    # Actual benchmark
    times = []
    for i in range(n_frames):
        u.trajectory[0]

        start = time.perf_counter()
        inter = pytim.WillardChandler(
            u,
            group=g,
            alpha=3.0,
            mesh=2.5,
            centered=True,
            centering_backend=backend,
            enable_timing=False
        )
        elapsed = time.perf_counter() - start
        times.append(elapsed)

        if (i + 1) % 20 == 0:
            print(f"  Progress: {i+1}/{n_frames} frames...")

    times = np.array(times) * 1000  # Convert to ms

    print(f"\n  Results (excluding warm-up):")
    print(f"    Mean:   {np.mean(times):.3f} ms")
    print(f"    Median: {np.median(times):.3f} ms")
    print(f"    Std:    {np.std(times):.3f} ms")
    print(f"    Min:    {np.min(times):.3f} ms")
    print(f"    Max:    {np.max(times):.3f} ms")

    return times

def main():
    print("=" * 70)
    print("GPU vs CPU Centering Benchmark")
    print("=" * 70)

    # Check backend availability
    print(f"\nBackend availability:")
    print(f"  C++ centering (HAS_CENTER_FULL): {HAS_CENTER_FULL}")
    print(f"  GPU centering (HAS_CENTER_GPU):  {HAS_CENTER_GPU}")

    if not HAS_CENTER_GPU:
        print("\nWARNING: GPU centering not available (CuPy not installed)")
        print("Cannot run benchmark...")
        return

    if not HAS_CENTER_FULL:
        print("\nWARNING: C++ centering not available")
        print("Cannot run benchmark...")
        return

    # Load test data - use larger system for meaningful benchmark
    import os
    data_dir = os.path.join(os.path.dirname(pytim.__file__), 'data')

    # Try to use a larger system
    test_files = [
        ('water.gro', 'name OW'),
        ('water-small.gro', 'name OW')
    ]

    u = None
    for filename, selection in test_files:
        try:
            filepath = f"{data_dir}/{filename}"
            u = mda.Universe(filepath, filepath)
            g = u.select_atoms(selection)
            print(f"\nTest system: {filename}")
            print(f"  Total atoms: {len(u.atoms)}")
            print(f"  Group atoms: {len(g)}")
            print(f"  Box: {u.dimensions[:3]}")
            break
        except:
            continue

    if u is None:
        print("ERROR: Could not load test data")
        return

    n_frames = 50

    # Benchmark CPU
    times_cpu = benchmark_centering(u, g, 'cpu', n_frames)

    # Benchmark GPU
    times_gpu = benchmark_centering(u, g, 'gpu', n_frames)

    # Compare results
    print(f"\n{'='*70}")
    print("Performance Comparison")
    print(f"{'='*70}")

    mean_cpu = np.mean(times_cpu)
    mean_gpu = np.mean(times_gpu)
    speedup = mean_cpu / mean_gpu

    print(f"\n  CPU (C++) mean time: {mean_cpu:.3f} ms")
    print(f"  GPU (CuPy) mean time: {mean_gpu:.3f} ms")
    print(f"  Speedup: {speedup:.2f}x", end="")

    if speedup > 1.0:
        print(" (GPU is faster)")
    elif speedup < 1.0:
        print(f" (CPU is {1/speedup:.2f}x faster)")
    else:
        print(" (equal)")

    print(f"\n{'='*70}")
    print("Benchmark complete!")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
