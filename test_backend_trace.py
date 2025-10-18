#!/usr/bin/env python
"""
Test to trace which backend is actually being used.
"""

import numpy as np
import MDAnalysis as mda
import pywc
import os

# Patch to trace which method is being called
original_evaluate_pbc_fast = pywc.gaussian_kde_pbc.gaussian_kde_pbc.evaluate_pbc_fast

def traced_evaluate_pbc_fast(self, points):
    print("\n>>> evaluate_pbc_fast called")
    from pywc.gaussian_kde_pbc import _wc_kde
    print(f">>> _wc_kde available: {_wc_kde is not None}")
    if _wc_kde is not None:
        print(f">>> Has evaluate_pbc_fast_auto: {hasattr(_wc_kde, 'evaluate_pbc_fast_auto')}")
    result = original_evaluate_pbc_fast(self, points)
    print(f">>> Result shape: {result.shape}")
    return result

pywc.gaussian_kde_pbc.gaussian_kde_pbc.evaluate_pbc_fast = traced_evaluate_pbc_fast

# Now run a test
data_dir = os.path.join(os.path.dirname(pywc.__file__), 'data')
u = mda.Universe(f"{data_dir}/water-small.gro", f"{data_dir}/water-small.gro")
g = u.select_atoms('name OW')

print("=" * 70)
print("Testing CPU backend")
print("=" * 70)

inter = pywc.WillardChandler(
    u,
    group=g,
    alpha=3.0,
    mesh=2.5,
    surface_backend='cpu',
    centered=False
)

print("\n" + "=" * 70)
print("Test complete")
print("=" * 70)
