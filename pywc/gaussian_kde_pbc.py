# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
""" Module: pywc.gaussian_kde_pbc
    ==============================
"""
from __future__ import print_function
import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial import cKDTree

try:
    from . import _wc_kde  # type: ignore[attr-defined]
except ImportError:  # pragma: no cover - optional acceleration
    _wc_kde = None


class gaussian_kde_pbc(gaussian_kde):
    # note that here "points" are those on the grid

    def evaluate_pbc_fast(self, points):
        grid = points
        pos = np.ascontiguousarray(self.pos, dtype=float)
        box = np.ascontiguousarray(np.asarray(self.box, dtype=float)[:3], dtype=float)
        gridT = np.ascontiguousarray(grid[::-1].T, dtype=float)

        if _wc_kde is not None:
            try:
                return _wc_kde.evaluate_pbc_fast_auto(
                    gridT,
                    pos,
                    box,
                    float(self.sigma))
            except AttributeError:
                # older extension without the auto neighbor search
                pass
            except Exception:
                # fall back to the python implementation if the accelerator fails
                pass

        d = self.sigma * 2.5
        tree = cKDTree(gridT, boxsize=box)
        # the indices of grid elements within distance d from each of the pos
        indlist = tree.query_ball_point(pos, d)

        if _wc_kde is not None:
            try:
                return _wc_kde.evaluate_pbc_fast(
                    gridT,
                    pos,
                    box,
                    float(self.sigma),
                    indlist)
            except Exception:
                # fall back to the python implementation if the accelerator fails
                pass

        results = np.zeros(grid.shape[1], dtype=float)
        scale = 2. * self.sigma**2
        for n, ind in enumerate(indlist):
            dr = gridT[ind, :] - pos[n]
            cond = np.where(dr > box / 2.)
            dr[cond] -= box[cond[1]]
            cond = np.where(dr < -box / 2.)
            dr[cond] += box[cond[1]]
            dens = np.exp(-np.sum(dr * dr, axis=1) / scale)
            results[ind] += dens

        return results
