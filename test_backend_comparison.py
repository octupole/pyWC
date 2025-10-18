#!/usr/bin/env python3
"""Test script to verify backend comparison functionality in WillardChandler."""

import numpy as np
import MDAnalysis as mda
import pytim
from pytim.datafiles import MICELLE_PDB

print("="*70)
print("Testing WillardChandler Backend Comparison")
print("="*70)

# Load test system
print("\n1. Loading test system...")
u = mda.Universe(MICELLE_PDB)
g = u.select_atoms('resname DPC')
print(f"   Loaded {len(g)} atoms")

# Test 1: Create interface with timing enabled
print("\n2. Creating WillardChandler interface with CPU backend and timing...")
inter = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                              surface_backend='cpu', enable_timing=True)
print(f"   Surface area: {inter.surface_area:.2f} Å²")
timing = inter.get_timing()
if timing is not None:
    print(f"   Computation time: {timing:.4f} s")
else:
    print("   Warning: Timing not recorded")

# Test 2: Compare backends
print("\n3. Comparing CPU and CuPy backends...")
try:
    results = inter.compare_backends(backends=['cpu', 'cupy'])

    print("\n4. Results summary:")
    for backend, result in results.items():
        if result['success']:
            print(f"   {backend}: {result['time']:.4f} s - SUCCESS")
        else:
            print(f"   {backend}: FAILED - {result['error']}")

except Exception as e:
    print(f"   Backend comparison failed: {e}")
    print("   (This is expected if CuPy is not installed)")

# Test 3: Create interface with GPU backend (if available)
print("\n5. Testing direct GPU backend creation...")
try:
    inter_gpu = pytim.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                                     surface_backend='cupy',
                                     surface_backend_options={'chunk_size': 4096},
                                     enable_timing=True)
    print(f"   GPU Surface area: {inter_gpu.surface_area:.2f} Å²")
    gpu_timing = inter_gpu.get_timing()
    if gpu_timing is not None:
        print(f"   GPU Computation time: {gpu_timing:.4f} s")
except Exception as e:
    print(f"   GPU backend not available: {e}")

print("\n" + "="*70)
print("Testing complete!")
print("="*70)
