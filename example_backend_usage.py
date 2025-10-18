#!/usr/bin/env python3
"""
Example usage of the new backend comparison features in WillardChandler.

This script demonstrates:
1. How to enable timing for surface computation
2. How to use different backends (CPU vs GPU)
3. How to compare backend performance
"""

# Example 1: Enable timing for a single backend
# ------------------------------------------------
import MDAnalysis as mda
import pywc
from pywc.datafiles import MICELLE_PDB

u = mda.Universe(MICELLE_PDB)
g = u.select_atoms('resname DPC')

# Create interface with timing enabled
inter = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0,
                              surface_backend='cpu', enable_timing=True)

# Get the computation time
timing = inter.get_timing()
print(f"CPU computation time: {timing:.4f} s")


# Example 2: Use GPU backend with custom options
# ------------------------------------------------
try:
    inter_gpu = pywc.WillardChandler(
        u, group=g,
        alpha=3.0,
        mesh=2.0,
        surface_backend='cupy',
        surface_backend_options={
            'chunk_size': 4096,
            'max_chunk_mem': 256 * 1024 * 1024
        },
        enable_timing=True
    )
    print(f"GPU computation time: {inter_gpu.get_timing():.4f} s")
except Exception as e:
    print(f"GPU backend not available: {e}")


# Example 3: Compare multiple backends
# ------------------------------------------------
inter = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0)

# Compare CPU vs GPU with custom options for each
results = inter.compare_backends(
    backends=['cpu', 'cupy'],
    backend_options={
        'cupy': {'chunk_size': 4096, 'max_chunk_mem': 256 * 1024 * 1024}
    },
    verbose=True
)

# Output will show:
# ============================================================
# Backend Performance Comparison
# ============================================================
# cpu         : 0.1234 s
# cupy        : 0.0456 s (2.71x)
# ============================================================


# Example 4: Access detailed timing results
# ------------------------------------------------
for backend, result in results.items():
    if result['success']:
        print(f"{backend}: {result['time']:.4f} s")
    else:
        print(f"{backend} failed: {result['error']}")
