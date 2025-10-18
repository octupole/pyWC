#!/usr/bin/env python
"""
Test script to verify all three surface backends work correctly.

Tests:
1. cpp (C++ via _wc_kde.cpp) - default
2. cupy (GPU via CuPy)
3. python (pure Python, no C++ acceleration)
"""

import sys
import numpy as np
import MDAnalysis as mda
import pytim
import os

def test_backend(backend_name, u, g):
    """Test a specific backend."""
    print(f"\n{'='*70}")
    print(f"Testing {backend_name.upper()} backend")
    print(f"{'='*70}")

    try:
        u.trajectory[0]  # Reset to first frame
        inter = pytim.WillardChandler(
            u,
            group=g,
            alpha=3.0,
            mesh=2.5,
            surface_backend=backend_name,
            centered=False
        )
        print(f"✓ {backend_name.upper()} backend successful")
        print(f"  Surface area: {inter.surface_area:.2f} Ų")
        return inter.surface_area, inter.triangulated_surface[0]  # Return area and vertices
    except Exception as e:
        print(f"✗ {backend_name.upper()} backend failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def main():
    print("=" * 70)
    print("Surface Backend Comparison Test")
    print("=" * 70)

    # Check availability
    from pytim.gaussian_kde_pbc import _wc_kde
    print(f"\nBackend availability:")
    print(f"  C++ (_wc_kde):  {_wc_kde is not None}")

    try:
        import cupy
        print(f"  CuPy (GPU):     True")
        cupy_available = True
    except ImportError:
        print(f"  CuPy (GPU):     False")
        cupy_available = False

    # Load test data
    data_dir = os.path.join(os.path.dirname(pytim.__file__), 'data')
    u = mda.Universe(f"{data_dir}/water-small.gro", f"{data_dir}/water-small.gro")
    g = u.select_atoms('name OW')

    print(f"\nTest system: water-small.gro")
    print(f"  Total atoms: {len(u.atoms)}")
    print(f"  Group atoms: {len(g)}")
    print(f"  Box: {u.dimensions[:3]}")

    # Test all backends
    results = {}

    # Test C++ backend
    area_cpp, verts_cpp = test_backend('cpp', u, g)
    if area_cpp is not None:
        results['cpp'] = (area_cpp, verts_cpp)

    # Test CPU backend (should map to cpp)
    area_cpu, verts_cpu = test_backend('cpu', u, g)
    if area_cpu is not None:
        results['cpu'] = (area_cpu, verts_cpu)

    # Test Python backend
    area_py, verts_py = test_backend('python', u, g)
    if area_py is not None:
        results['python'] = (area_py, verts_py)

    # Test CuPy backend if available
    if cupy_available:
        area_cupy, verts_cupy = test_backend('cupy', u, g)
        if area_cupy is not None:
            results['cupy'] = (area_cupy, verts_cupy)

    # Compare results
    print(f"\n{'='*70}")
    print("Results Comparison")
    print(f"{'='*70}\n")

    if len(results) < 2:
        print("Not enough successful backends to compare")
        return

    print("Surface areas:")
    for backend, (area, _) in results.items():
        print(f"  {backend:10s}: {area:.4f} Ų")

    # Check if cpu and cpp give identical results
    if 'cpu' in results and 'cpp' in results:
        area_diff = abs(results['cpu'][0] - results['cpp'][0])
        print(f"\nCPU vs CPP difference: {area_diff:.10f} Ų")
        if area_diff < 1e-10:
            print("✓ CPU correctly maps to CPP backend")
        else:
            print("✗ WARNING: CPU and CPP differ!")

    # Compare cpp vs python
    if 'cpp' in results and 'python' in results:
        area_diff = abs(results['cpp'][0] - results['python'][0])
        vert_diff = np.max(np.abs(results['cpp'][1] - results['python'][1]))
        print(f"\nCPP vs Python:")
        print(f"  Area difference:     {area_diff:.10f} Ų")
        print(f"  Vertices max diff:   {vert_diff:.10f} Å")
        if area_diff < 0.01:  # Allow small numerical differences
            print("✓ CPP and Python results are consistent")
        else:
            print("✗ WARNING: Large difference between CPP and Python")

    # Compare cpp vs cupy if available
    if 'cpp' in results and 'cupy' in results:
        area_diff = abs(results['cpp'][0] - results['cupy'][0])
        vert_diff = np.max(np.abs(results['cpp'][1] - results['cupy'][1]))
        print(f"\nCPP vs CuPy:")
        print(f"  Area difference:     {area_diff:.10f} Ų")
        print(f"  Vertices max diff:   {vert_diff:.10f} Å")
        if area_diff < 0.01:
            print("✓ CPP and CuPy results are consistent")
        else:
            print("✗ WARNING: Large difference between CPP and CuPy")

    print(f"\n{'='*70}")
    print("Test complete!")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
