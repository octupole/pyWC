#!/usr/bin/env python3
"""Example demonstrating detailed timing breakdown for Willard-Chandler calculations.

This script shows how to use the debugging timing features to understand
where time is spent in the Willard-Chandler surface calculation.
"""

import MDAnalysis as mda
import pytim
from pytim.datafiles import MICELLE_PDB

print("="*70)
print("Willard-Chandler Detailed Timing Example")
print("="*70)

# Load test system
u = mda.Universe(MICELLE_PDB)
g = u.select_atoms('resname DPC')
print(f"\nLoaded system with {len(g)} atoms")

# Create interface with timing enabled
print("\n1. Creating WillardChandler with CPU backend and timing...")
inter = pytim.WillardChandler(
    u,
    group=g,
    alpha=3.0,
    mesh=2.0,
    surface_backend='cpu',
    enable_timing=True
)

print(f"   Surface area: {inter.surface_area:.2f} Å²")
print(f"   Initial computation time: {inter.get_timing():.4f} s")

# Process multiple frames to collect statistics
print("\n2. Processing multiple frames to collect timing statistics...")
n_frames = min(5, len(u.trajectory))
for i, ts in enumerate(u.trajectory[:n_frames]):
    inter._assign_layers()
    print(f"   Frame {ts.frame}: {inter.get_timing()*1000:.2f} ms")

# Get detailed timing breakdown
print("\n3. Detailed timing breakdown:")
inter.print_timing_breakdown(verbose=True)

# Get raw timing data for custom analysis
print("\n4. Accessing raw timing data:")
timings = inter.get_detailed_timings()
if timings:
    print("\nAvailable timing components:")
    for component in sorted(timings.keys()):
        if component != 'total':
            stats = timings[component]
            print(f"  {component:25s}: {stats['mean']*1000:7.3f} ms (avg)")

# Example: Compare CPU vs GPU backends
print("\n5. Comparing backends (if GPU available)...")
try:
    inter_comparison = pytim.WillardChandler(
        u, group=g, alpha=3.0, mesh=2.0, enable_timing=True
    )

    backend_results = inter_comparison.compare_backends(
        backends=['cpu', 'cupy'],
        verbose=True
    )

    # Show detailed breakdown for each backend
    for backend in ['cpu', 'cupy']:
        if backend_results[backend]['success']:
            print(f"\n{backend.upper()} Backend Details:")
            inter_temp = pytim.WillardChandler(
                u, group=g, alpha=3.0, mesh=2.0,
                surface_backend=backend,
                enable_timing=True
            )
            # Process a few frames
            for ts in u.trajectory[:3]:
                inter_temp._assign_layers()
            inter_temp.print_timing_breakdown(verbose=False)

except Exception as e:
    print(f"   Backend comparison skipped: {e}")

print("\n" + "="*70)
print("Example complete!")
print("="*70)
