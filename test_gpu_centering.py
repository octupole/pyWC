#!/usr/bin/env python
"""
Test script to verify GPU centering integration.

This script tests:
1. That the GPU centering backend can be selected
2. That it runs without errors
3. Comparison between CPU and GPU backends (if both available)
"""

import sys
import numpy as np
import MDAnalysis as mda
import pytim
from pytim._center_impl import HAS_CENTER_GPU, HAS_CENTER_FULL

def main():
    print("=" * 70)
    print("GPU Centering Integration Test")
    print("=" * 70)

    # Check backend availability
    print(f"\nBackend availability:")
    print(f"  C++ centering (HAS_CENTER_FULL): {HAS_CENTER_FULL}")
    print(f"  GPU centering (HAS_CENTER_GPU):  {HAS_CENTER_GPU}")

    if not HAS_CENTER_GPU:
        print("\nWARNING: GPU centering not available (CuPy not installed)")
        print("Skipping GPU tests...")
        return

    # Load test data
    import os
    data_dir = os.path.join(os.path.dirname(pytim.__file__), 'data')
    u = mda.Universe(f"{data_dir}/water-small.gro",
                     f"{data_dir}/water-small.gro")

    g = u.select_atoms('name OW')
    print(f"\nTest system:")
    print(f"  Total atoms: {len(u.atoms)}")
    print(f"  Group atoms: {len(g)}")
    print(f"  Box: {u.dimensions[:3]}")

    # Test CPU backend (C++)
    print(f"\n{'='*70}")
    print("Testing CPU (C++) centering backend...")
    print(f"{'='*70}")
    try:
        inter_cpu = pytim.WillardChandler(
            u,
            group=g,
            alpha=3.0,
            mesh=2.5,
            centered=True,
            centering_backend='cpu',
            enable_timing=True
        )
        print("✓ CPU backend successful")

        # Get positions after CPU centering
        pos_cpu = np.copy(u.atoms.positions)

        # Reset positions
        u.trajectory[0]

    except Exception as e:
        print(f"✗ CPU backend failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # Test GPU backend
    print(f"\n{'='*70}")
    print("Testing GPU (CuPy) centering backend...")
    print(f"{'='*70}")
    try:
        inter_gpu = pytim.WillardChandler(
            u,
            group=g,
            alpha=3.0,
            mesh=2.5,
            centered=True,
            centering_backend='gpu',
            enable_timing=True
        )
        print("✓ GPU backend successful")

        # Get positions after GPU centering
        pos_gpu = np.copy(u.atoms.positions)

    except Exception as e:
        print(f"✗ GPU backend failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # Compare results
    print(f"\n{'='*70}")
    print("Comparing CPU vs GPU results...")
    print(f"{'='*70}")

    pos_diff = np.abs(pos_cpu - pos_gpu)
    max_diff = np.max(pos_diff)
    mean_diff = np.mean(pos_diff)

    print(f"  Maximum position difference: {max_diff:.10f} Å")
    print(f"  Mean position difference:    {mean_diff:.10f} Å")

    # Check if results are close (within numerical precision)
    if max_diff < 1e-6:
        print("\n✓ CPU and GPU results match within numerical precision")
    else:
        print(f"\n✗ WARNING: Large difference detected (> 1e-6 Å)")

    print(f"\n{'='*70}")
    print("Test complete!")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
