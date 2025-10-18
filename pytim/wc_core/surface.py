"""High-level Willard–Chandler surface construction routines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Optional

import numpy as np
from skimage import measure

from .grid import build_grid
from .density import density_map


@dataclass
class SurfaceResult:
    """Container for the data produced by :func:`compute_surface`."""

    density_field: np.ndarray
    grid: np.ndarray
    spacing: np.ndarray
    ngrid: np.ndarray
    triangulated_surface: Tuple[np.ndarray, np.ndarray, np.ndarray]
    surface_area: float


def compute_surface(pos: np.ndarray,
                    box: np.ndarray,
                    mesh: float,
                    sigma: float,
                    density_cutoff: float | None,
                    order: str = "xyz",
                    backend: str = "cpu",
                    backend_options: Optional[Dict[str, object]] = None) -> SurfaceResult:
    """Compute the Willard–Chandler density field and triangulated surface.

    Parameters
    ----------
    pos : ndarray, shape (N, 3)
        Coordinates of the particles used to construct the intrinsic density.
    box : ndarray, shape (3,)
        Simulation box lengths.
    mesh : float
        Target mesh spacing for the grid used in the KDE evaluation.
    sigma : float
        Gaussian width (interface ``alpha``).
    density_cutoff : float or None
        Isosurface value forwarded to :func:`skimage.measure.marching_cubes`.
    order : str, optional
        Axis ordering used when building the grid (default 'xyz').
    backend : str, optional
        Backend for density computation (default 'cpu'):
        - 'cpp' or 'cpu': C++ accelerated (via _wc_kde.cpp) - recommended
        - 'cupy' or 'gpu': GPU accelerated (requires CuPy)
        - 'python': Pure Python (for testing/debugging)
    backend_options : dict, optional
        Additional options passed to the backend.

    Returns
    -------
    SurfaceResult
        Data structure containing the density field and triangulated surface.
    """
    ngrid, spacing, grid = build_grid(box, mesh, order=order)
    backend = backend.lower()
    backend_options = backend_options or {}

    # Map 'cpu' to 'cpp' for backward compatibility
    if backend == "cpu":
        backend = "cpp"

    if backend == "cupy" or backend == "gpu":
        try:
            from .gpu import evaluate_density_gpu
            density_field = evaluate_density_gpu(pos, grid.T, sigma, box,
                                                 **backend_options)
        except Exception as exc:  # pragma: no cover - optional and user controlled
            raise RuntimeError(
                "CuPy backend requested but could not be used."
            ) from exc
    elif backend == "cpp":
        # Use C++ accelerated KDE (via _wc_kde.cpp)
        kernel, _ = density_map(pos, grid, sigma, box)
        kernel.pos = pos.copy()
        density_field = kernel.evaluate_pbc_fast(grid)
    elif backend == "python":
        # Force pure Python implementation
        from ..gaussian_kde_pbc import _wc_kde
        original_wc_kde = _wc_kde
        # Temporarily disable C++ extension
        import pytim.gaussian_kde_pbc
        pytim.gaussian_kde_pbc._wc_kde = None
        try:
            kernel, _ = density_map(pos, grid, sigma, box)
            kernel.pos = pos.copy()
            density_field = kernel.evaluate_pbc_fast(grid)
        finally:
            # Restore C++ extension
            pytim.gaussian_kde_pbc._wc_kde = original_wc_kde
    else:
        raise ValueError(f"Unknown backend: {backend}. Must be 'cpp', 'cupy'/'gpu', or 'python'.")

    volume = density_field.reshape(tuple(np.array(ngrid[::-1]).astype(int)))
    verts, faces, normals, _ = measure.marching_cubes(
        volume, density_cutoff, spacing=tuple(spacing))
    verts += spacing[::-1] / 2.
    triangulated_surface = (np.fliplr(verts), faces, np.fliplr(normals))
    surface_area = measure.mesh_surface_area(verts, faces)
    return SurfaceResult(
        density_field=density_field,
        grid=grid,
        spacing=spacing,
        ngrid=ngrid,
        triangulated_surface=triangulated_surface,
        surface_area=surface_area,
    )
