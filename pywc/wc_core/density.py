"""Density helpers for Willard–Chandler calculations."""

from __future__ import annotations

import numpy as np

from ..gaussian_kde_pbc import gaussian_kde_pbc


def density_map(pos: np.ndarray, grid: np.ndarray, sigma: float, box: np.ndarray):
    """Return a periodic Gaussian KDE object evaluated on ``grid``.

    This is a light wrapper around :class:`pywc.gaussian_kde_pbc.gaussian_kde_pbc`
    mirroring the historical behaviour of :func:`pywc.utilities.density_map`.
    The function is kept here so that both the legacy utilities module and the
    refactored Willard–Chandler core can share the same implementation.

    Parameters
    ----------
    pos : ndarray
        Particle coordinates with shape (N, 3).
    grid : ndarray
        Grid coordinates as returned by :func:`wc_core.grid.build_grid`.
    sigma : float
        Gaussian width (the ``alpha`` parameter in the interface definition).
    box : ndarray
        Simulation box lengths, shape (3,).

    Returns
    -------
    kernel : :class:`gaussian_kde_pbc`
        KDE object with ``box`` and ``sigma`` attributes already configured.
    std : float
        Standard deviation of the KDE input, matching the historic return
        signature.
    """
    values = np.vstack([pos[:, 0], pos[:, 1], pos[:, 2]])
    kernel = gaussian_kde_pbc(values, bw_method=sigma / values.std(ddof=1))
    kernel.box = box
    kernel.sigma = sigma
    return kernel, values.std(ddof=1)
