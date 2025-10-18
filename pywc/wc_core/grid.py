"""Grid helpers for Willard–Chandler calculations."""

from __future__ import annotations

from typing import Tuple

import numpy as np

from ..utilities_mesh import compute_compatible_mesh_params, generate_grid_in_box


def build_grid(box: np.ndarray, mesh: float, order: str = "xyz") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return grid metadata for the supplied box and mesh size.

    Parameters
    ----------
    box : array-like, shape (3,)
        Box lengths along x, y, z.
    mesh : float
        Target mesh spacing used to discretise the box.
    order : str, optional
        Axis ordering passed to :func:`generate_grid_in_box` (default 'xyz').

    Returns
    -------
    ngrid : ndarray, shape (3,)
        Number of grid points along each axis.
    spacing : ndarray, shape (3,)
        Actual spacing realised along each axis.
    grid : ndarray, shape (3, N)
        3×N transposed grid coordinates suitable for :mod:`gaussian_kde_pbc`.
    """
    ngrid, spacing = compute_compatible_mesh_params(mesh, box)
    grid = generate_grid_in_box(box, ngrid, order=order)
    return ngrid, spacing, grid
