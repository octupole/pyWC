#!/usr/bin/env python3
"""Compare CPU vs GPU density evaluation for the Willard–Chandler kernel."""

from __future__ import annotations

import argparse
import time

import numpy as np

from pywc.wc_core.surface import compute_surface


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=10000,
                        help="Number of random particles (default: 10000)")
    parser.add_argument("--grid", type=int, default=50000,
                        help="Number of random grid points (default: 50000)")
    parser.add_argument("--box", type=float, nargs=3, default=(60.0, 60.0, 60.0),
                        metavar=("LX", "LY", "LZ"),
                        help="Box lengths (default: 60 60 60)")
    parser.add_argument("--sigma", type=float, default=2.0,
                        help="Gaussian width (alpha)")
    parser.add_argument("--mesh", type=float, default=2.5,
                        help="Target mesh spacing")
    parser.add_argument("--seed", type=int, default=0,
                        help="Random seed")
    parser.add_argument("--chunk", type=int, default=4096,
                        help="Max particles per GPU chunk")
    parser.add_argument("--max-mem", type=int, default=256 * 1024 * 1024,
                        help="Max chunk memory (bytes) on GPU")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    box = np.asarray(args.box, dtype=float)
    pos = rng.random((args.particles, 3)) * box
    grid = rng.random((args.grid, 3)) * box

    print(f"Particles: {args.particles}, grid points: {args.grid}")
    print(f"Box: {box}, sigma: {args.sigma}, mesh: {args.mesh}\n")

    start = time.perf_counter()
    cpu_res = compute_surface(pos, box, mesh=args.mesh, sigma=args.sigma,
                              density_cutoff=None, backend='cpu')
    cpu_time = time.perf_counter() - start

    start = time.perf_counter()
    gpu_res = compute_surface(pos, box, mesh=args.mesh, sigma=args.sigma,
                              density_cutoff=None, backend='cupy',
                              backend_options={'chunk_size': args.chunk,
                                               'max_chunk_mem': args.max_mem})
    gpu_time = time.perf_counter() - start

    diff = gpu_res.density_field - cpu_res.density_field
    max_abs = np.max(np.abs(diff))
    rms = np.sqrt(np.mean(diff ** 2))

    print("CPU time : %.3f s" % cpu_time)
    print("GPU time : %.3f s" % gpu_time)
    print("Speedup  : %.2fx" % (cpu_time / gpu_time if gpu_time > 0 else float('inf')))
    print("Max abs diff: %.3e" % max_abs)
    print("RMS diff   : %.3e" % rms)


if __name__ == "__main__":
    main()
