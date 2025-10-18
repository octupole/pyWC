"""Quick benchmark comparing Python vs pybind11 KDE accumulation.

Usage example:

    python benchmarks/willard_chandler_kde.py --points 4000 --grid 64000 --repeat 5

The script relies only on NumPy/Scipy and does not need trajectory files.
"""
from __future__ import annotations

import argparse
import contextlib
import time
from typing import Callable, Tuple

import numpy as np
from pytim.gaussian_kde_pbc import gaussian_kde_pbc as _gaussian_kde_pbc
from pytim.wc_core.surface import compute_surface


@contextlib.contextmanager
def force_python_impl():
    from pytim import gaussian_kde_pbc as module

    saved = module._wc_kde  # type: ignore[attr-defined]
    module._wc_kde = None
    try:
        yield
    finally:
        module._wc_kde = saved


def _create_inputs(n_particles: int, n_grid: int, box: Tuple[float, float, float], sigma: float):
    rng = np.random.default_rng(0)
    box_arr = np.asarray(box, dtype=float)
    positions = rng.random((n_particles, 3)) * box_arr

    grid_x = np.linspace(0.0, box[0], int(round(n_grid ** (1.0 / 3))), endpoint=False)
    grid_points = np.stack(np.meshgrid(grid_x, grid_x, grid_x, indexing="ij"), axis=-1).reshape(-1, 3)
    grid_points = grid_points[:n_grid]

    data = positions.T
    kde = _gaussian_kde_pbc(data, bw_method=sigma / data.std(ddof=1))
    kde.pos = positions.copy()
    kde.box = np.asarray(box, dtype=float)
    kde.sigma = sigma

    grid_eval = grid_points.T.copy()

    mesh = box_arr[0] / int(round(n_grid ** (1.0 / 3)))
    return kde, grid_eval, positions, box_arr, mesh


def _run(kde,
         grid,
         use_accel: bool,
         backend: str,
         backend_options: dict,
         positions: np.ndarray,
         box: np.ndarray,
         mesh: float,
         sigma: float) -> float:
    if use_accel:
        if backend == "cupy":
            result = compute_surface(positions,
                                      box,
                                      mesh=mesh,
                                      sigma=sigma,
                                      density_cutoff=None,
                                      backend=backend,
                                      backend_options=backend_options)
            return result.density_field
        return kde.evaluate_pbc_fast(grid)
    with force_python_impl():
        return kde.evaluate_pbc_fast(grid)


def benchmark(step_fn: Callable[[], float], repeat: int) -> float:
    timings = []
    for _ in range(repeat):
        start = time.perf_counter()
        step_fn()
        timings.append(time.perf_counter() - start)
    return float(np.median(timings))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=4000, help="Number of particles")
    parser.add_argument("--grid", type=int, default=64000, help="Number of grid points")
    parser.add_argument("--sigma", type=float, default=2.0, help="Gaussian width")
    parser.add_argument("--repeat", type=int, default=5, help="Number of repetitions")
    parser.add_argument("--backend",
                        choices=["cpu", "cupy"],
                        default="cpu",
                        help="Accelerated backend to benchmark (default: cpu)")
    parser.add_argument("--backend-option",
                        action="append",
                        default=[],
                        metavar="KEY=VALUE",
                        help="Backend-specific option (may be given multiple times)")
    parser.add_argument(
        "--min-speedup",
        type=float,
        default=None,
        help="If set, fail (exit 1) when accelerated speedup falls below this threshold.")
    args = parser.parse_args()

    box = (60.0, 60.0, 60.0)
    kde, grid, positions, box_arr, mesh = _create_inputs(args.points, args.grid, box, args.sigma)

    backend_options = {}
    for entry in args.backend_option:
        if "=" not in entry:
            raise SystemExit(f"Invalid backend option '{entry}', expected KEY=VALUE")
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        for caster in (int, float):
            try:
                value_casted = caster(value)
                break
            except ValueError:
                continue
        else:
            value_casted = value
        backend_options[key] = value_casted

    accel = benchmark(lambda: _run(kde,
                                   grid,
                                   True,
                                   args.backend,
                                   backend_options,
                                   positions,
                                   box_arr,
                                   mesh,
                                   args.sigma), args.repeat)
    python = benchmark(lambda: _run(kde,
                                    grid,
                                    False,
                                    args.backend,
                                    backend_options,
                                    positions,
                                    box_arr,
                                    mesh,
                                    args.sigma), args.repeat)

    speedup = python / accel if accel > 0 else float("inf")
    print(f"Accelerated median: {accel:.4f} s")
    print(f"Python median    : {python:.4f} s")
    print(f"Speedup          : {speedup:.2f}x")

    if args.min_speedup is not None and speedup < args.min_speedup:
        raise SystemExit(
            f"Speedup {speedup:.2f}x is below required minimum ({args.min_speedup:.2f}x)")


if __name__ == "__main__":
    main()
