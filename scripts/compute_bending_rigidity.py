#!/usr/bin/env python3
r"""Estimate membrane bending rigidity from Willard–Chandler surfaces.

This helper script follows the workflow sketched in ``IMG_20250219_145156.jpg``:

1. Build two Willard–Chandler (WC) surfaces, one per leaflet.
2. Convert the triangulated surfaces into gridded height fields and average
   them to obtain a mid-surface height map as well as the local thickness.
3. Fourier transform the mid-surface to obtain the undulation spectrum
   :math:`\langle |h(\mathbf{q})|^2 \rangle`.
4. Perform a Helfrich fit ``k_BT / <|h(q)|^2> = sigma * q^2 + kappa * q^4`` to
   extract the surface tension ``sigma`` and bending rigidity ``kappa``.

Outputs are written as CSV files for further inspection and plotting.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import MDAnalysis as mda
import numpy as np
import pywc
from skimage import measure


def split_leaflets_by_height(bilayer_group: mda.AtomGroup) -> Tuple[mda.AtomGroup, mda.AtomGroup]:
    residues = bilayer_group.residues
    if len(residues) < 2:
        raise ValueError("Bilayer selection must contain at least two residues for leaflet partitioning")

    heights = np.array([res.atoms.center_of_geometry(wrap=True)[2] for res in residues])
    median_z = np.median(heights)
    top_mask = heights >= median_z
    bottom_mask = ~top_mask

    top_residues = residues[top_mask]
    bottom_residues = residues[bottom_mask]

    if len(top_residues) == 0 or len(bottom_residues) == 0:
        raise ValueError("Automatic leaflet split failed; check the bilayer selection")

    return top_residues.atoms, bottom_residues.atoms


KBOLTZ_KJMOL_PER_K = 0.00831446261815324  # kJ mol^-1 K^-1


def positive_float(value: str) -> float:
    number = float(value)
    if number <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate bending rigidity via WC surfaces",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-s", "--topology", required=True, help="Topology file (TPR/GRO/PSF…)")
    parser.add_argument("-x", "--trajectory", required=True, help="Trajectory file (XTC/TRR/DCD…)")
    parser.add_argument(
        "--leaflet-a",
        default=None,
        help="Selection string for the first leaflet (if omitted, use --bilayer-selection)",
    )
    parser.add_argument(
        "--leaflet-b",
        default=None,
        help="Selection string for the second leaflet (if omitted, use --bilayer-selection)",
    )
    parser.add_argument(
        "--bilayer-selection",
        default=None,
        help="Single selection covering both leaflets; the script splits it using MDAnalysis.LeafletFinder",
    )
    parser.add_argument(
        "--alpha", type=positive_float, default=3.0, help="Gaussian width α (Å) for WC surfaces"
    )
    parser.add_argument(
        "--mesh", type=positive_float, default=2.0, help="WC grid spacing (Å)"
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=32,
        help="Number of bins along x/y when rasterising the surfaces",
    )
    parser.add_argument("--stride", type=int, default=1, help="Frame stride during analysis")
    parser.add_argument("--start", type=int, default=None, help="First frame index (inclusive)")
    parser.add_argument("--stop", type=int, default=None, help="Last frame index (exclusive)")
    parser.add_argument(
        "--temperature",
        type=positive_float,
        default=310.0,
        help="Simulation temperature in Kelvin (used for k_BT)",
    )
    parser.add_argument(
        "--qmin",
        type=float,
        default=0.0,
        help="Minimum wave-vector magnitude q to include in Helfrich fit (Å⁻¹)",
    )
    parser.add_argument(
        "--qmax",
        type=float,
        default=None,
        help="Maximum q to include in Helfrich fit (Å⁻¹). None → all available modes",
    )
    parser.add_argument(
        "--min-modes",
        type=int,
        default=4,
        help="Minimum number of discrete modes per q-bin required for spectrum averaging",
    )
    parser.add_argument(
        "--q-bins",
        type=int,
        default=25,
        help="Number of radial bins in reciprocal-space for spectrum averaging",
    )
    parser.add_argument(
        "--force-density",
        type=float,
        default=None,
        help="Explicit density cutoff for WC surfaces (overrides automatic level selection)",
    )
    parser.add_argument(
        "--density-level",
        type=float,
        default=1.0,
        help="Relative density level in [0, 2] for the WC isosurface when --force-density is omitted",
    )
    parser.add_argument(
        "--detrend-plane",
        action="store_true",
        help="Remove the best-fit plane from each mid-surface before the FFT",
    )
    parser.add_argument(
        "--thickness-output",
        type=Path,
        default=Path("wc_thickness_mean.csv"),
        help="CSV file storing the time-averaged thickness map",
    )
    parser.add_argument(
        "--spectrum-output",
        type=Path,
        default=Path("wc_undulation_spectrum.csv"),
        help="CSV file where averaged |h(q)|^2 per radial bin is written",
    )
    parser.add_argument(
        "--mode-output",
        type=Path,
        default=Path("wc_modes_raw.npz"),
        help="Optional NPZ file with the full 2D average spectrum (for debugging)",
    )
    return parser


def resolve_selection(expr: str) -> str:
    """Allow referencing environment variables via ${VAR}."""

    if expr.startswith("$"):
        key = expr[1:]
        value = os.environ.get(key)
        if not value:
            raise ValueError(f"Environment variable {key!r} referenced by selection is not set")
        return value
    return expr


def reassign_surface(interface: pywc.WillardChandler, force: bool) -> None:
    if force or not getattr(interface, "autoassign", True):
        interface._assign_layers()


def build_surface_with_cutoff(interface: pywc.WillardChandler, cutoff: float) -> Optional[float]:
    density_field = getattr(interface, "density_field", None)
    spacing = getattr(interface, "spacing", None)
    grid = getattr(interface, "ngrid", None)
    if density_field is None or spacing is None or grid is None:
        interface.density_cutoff = cutoff
        return None

    volume = density_field.reshape(tuple(np.array(grid[::-1]).astype(int)))
    verts, faces, normals, _ = measure.marching_cubes(volume, cutoff, spacing=tuple(spacing))
    interface.triangulated_surface = [np.fliplr(verts), faces, np.fliplr(normals)]
    interface.surface_area = measure.mesh_surface_area(verts, faces)
    interface.density_cutoff = cutoff
    return cutoff


def apply_density_level(interface: pywc.WillardChandler, level: float) -> Optional[float]:
    density_field = getattr(interface, "density_field", None)
    spacing = getattr(interface, "spacing", None)
    grid = getattr(interface, "ngrid", None)
    if density_field is None or spacing is None or grid is None:
        return None

    rho_min = float(np.nanmin(density_field))
    rho_max = float(np.nanmax(density_field))
    if np.isclose(rho_min, rho_max):
        interface.density_cutoff = rho_min
        return rho_min

    cutoff = rho_min + 0.5 * level * (rho_max - rho_min)
    span = rho_max - rho_min
    if span > 0.0:
        eps = span * 1e-6
        cutoff = min(max(cutoff, rho_min + eps), rho_max - eps)

    return build_surface_with_cutoff(interface, cutoff)


def ensure_surface(
    interface: pywc.WillardChandler,
    *,
    force_reassign: bool,
    absolute_level: Optional[float],
    relative_level: float,
) -> float:
    reassign_surface(interface, force_reassign)
    if absolute_level is not None:
        cutoff = build_surface_with_cutoff(interface, absolute_level)
        return float(absolute_level if cutoff is None else cutoff)

    cutoff = apply_density_level(interface, relative_level)
    if cutoff is None:
        return float(interface.density_cutoff)
    return float(cutoff)


def rasterise_surface(
    vertices: np.ndarray,
    lx: float,
    ly: float,
    grid_size: int,
    *,
    reduce: str = "mean",
) -> np.ndarray:
    """Project triangulated surface vertices onto an XY grid and aggregate z values."""

    if grid_size <= 0:
        raise ValueError("grid_size must be positive")

    heights = np.full((grid_size, grid_size), np.nan, dtype=float)
    counts = np.zeros((grid_size, grid_size), dtype=int)

    if vertices.size == 0:
        return heights

    x = np.mod(vertices[:, 0], lx)
    y = np.mod(vertices[:, 1], ly)
    z = vertices[:, 2]

    ix = np.floor(x / lx * grid_size).astype(int)
    iy = np.floor(y / ly * grid_size).astype(int)
    ix = np.clip(ix, 0, grid_size - 1)
    iy = np.clip(iy, 0, grid_size - 1)

    if reduce == "mean":
        accum = np.zeros((grid_size, grid_size), dtype=float)
        for xi, yi, zi in zip(ix, iy, z):
            accum[xi, yi] += zi
            counts[xi, yi] += 1
        mask = counts > 0
        heights[mask] = accum[mask] / counts[mask]
    elif reduce == "max":
        for xi, yi, zi in zip(ix, iy, z):
            counts[xi, yi] += 1
            current = heights[xi, yi]
            if np.isnan(current) or zi > current:
                heights[xi, yi] = zi
    elif reduce == "min":
        for xi, yi, zi in zip(ix, iy, z):
            counts[xi, yi] += 1
            current = heights[xi, yi]
            if np.isnan(current) or zi < current:
                heights[xi, yi] = zi
    else:
        raise ValueError(f"Unknown reduction '{reduce}'")

    return heights


def fill_nan_nearest(field: np.ndarray) -> np.ndarray:
    """Replace NaNs with the nearest non-NaN neighbour (Euclidean distance)."""

    output = field.copy()
    nan_mask = np.isnan(output)
    if not np.any(nan_mask):
        return output
    valid_coords = np.argwhere(~nan_mask)
    valid_values = output[~nan_mask]
    if valid_values.size == 0:
        raise ValueError("Surface rasterisation produced only NaN values; check selections")

    for (ix, iy) in np.argwhere(nan_mask):
        distances = np.sum((valid_coords - np.array([ix, iy])) ** 2, axis=1)
        nearest = np.argmin(distances)
        output[ix, iy] = valid_values[nearest]
    return output


def detrend_plane(field: np.ndarray, lx: float, ly: float) -> np.ndarray:
    """Subtract best-fit plane ax + by + c from field (in-place safe)."""

    nx, ny = field.shape
    xs = (np.arange(nx) + 0.5) * (lx / nx)
    ys = (np.arange(ny) + 0.5) * (ly / ny)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    A = np.column_stack((X.ravel(), Y.ravel(), np.ones(nx * ny)))
    b = field.ravel()
    coeffs, *_ = np.linalg.lstsq(A, b, rcond=None)
    plane = (coeffs[0] * X + coeffs[1] * Y + coeffs[2]).reshape(field.shape)
    return field - plane


def radial_average(
    qx: np.ndarray,
    qy: np.ndarray,
    power: np.ndarray,
    *,
    bins: int,
    min_modes: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (q_centres, power_mean, counts) for radial bins."""

    q_mod = np.sqrt(qx ** 2 + qy ** 2)
    mask = q_mod > 0.0

    q_vals = q_mod[mask]
    power_vals = power[mask]

    q_max = q_vals.max()
    edges = np.linspace(0.0, q_max, bins + 1)
    bin_indices = np.digitize(q_vals, edges) - 1

    q_mean = np.zeros(bins)
    p_mean = np.zeros(bins)
    counts = np.zeros(bins, dtype=int)

    for idx in range(bins):
        mask_bin = bin_indices == idx
        n = np.count_nonzero(mask_bin)
        if n == 0:
            continue
        counts[idx] = n
        q_mean[idx] = q_vals[mask_bin].mean()
        p_mean[idx] = power_vals[mask_bin].mean()

    valid = counts >= min_modes
    return q_mean[valid], p_mean[valid], counts[valid]


@dataclass
class SpectralResult:
    q: np.ndarray
    power: np.ndarray
    counts: np.ndarray
    sigma: float
    sigma_err: float
    kappa: float
    kappa_err: float


def helfrich_fit(
    q: np.ndarray,
    power: np.ndarray,
    temperature: float,
) -> Tuple[float, float, float, float]:
    """Fit kBT / power = sigma * q^2 + kappa * q^4 via linear least squares."""

    if q.size < 2:
        raise ValueError("Not enough q-points for a Helfrich fit (need ≥2)")

    kbt = KBOLTZ_KJMOL_PER_K * temperature
    y = kbt / power
    A = np.column_stack((q ** 2, q ** 4))
    coeffs, residuals, rank, svals = np.linalg.lstsq(A, y, rcond=None)
    sigma, kappa = coeffs

    dof = max(0, q.size - 2)
    if residuals.size == 0 or dof == 0:
        sigma_err = float("nan")
        kappa_err = float("nan")
    else:
        variance = residuals[0] / dof
        cov = variance * np.linalg.inv(A.T @ A)
        sigma_err = math.sqrt(cov[0, 0])
        kappa_err = math.sqrt(cov[1, 1])

    return sigma, sigma_err, kappa, kappa_err


def assign_leaflets(
    universe: mda.Universe,
    *,
    sel_a: Optional[str],
    sel_b: Optional[str],
    bilayer_selection: Optional[str],
) -> Tuple[mda.AtomGroup, mda.AtomGroup]:
    if sel_a and sel_b:
        group_a = universe.select_atoms(resolve_selection(sel_a))
        group_b = universe.select_atoms(resolve_selection(sel_b))
        return group_a, group_b

    if bilayer_selection is None:
        raise ValueError(
            "Either --leaflet-a/--leaflet-b must be provided or --bilayer-selection must be specified"
        )

    try:
        from MDAnalysis.analysis.leaflet import LeafletFinder
    except ImportError as exc:
        raise ImportError(
            "MDAnalysis.analysis.leaflet.LeafletFinder is required for --bilayer-selection"
        ) from exc

    selection = resolve_selection(bilayer_selection)
    bilayer_group = universe.select_atoms(selection)
    if len(bilayer_group) == 0:
        raise ValueError("--bilayer-selection did not match any atoms")

    leaflet = LeafletFinder(universe, selection, cutoff=10.0)
    groups = sorted(leaflet.groups(), key=lambda g: g.n_atoms, reverse=True)
    if len(groups) >= 2:
        return groups[0], groups[1]

    # Fallback: monotonic split by height
    return split_leaflets_by_height(bilayer_group)


def main(args: argparse.Namespace) -> None:
    if args.grid_size <= 0:
        raise ValueError("--grid-size must be positive")
    if args.stride <= 0:
        raise ValueError("--stride must be positive")
    if args.density_level < 0.0 or args.density_level > 2.0:
        raise ValueError("--density-level must lie in [0, 2]")

    u = mda.Universe(args.topology, args.trajectory)
    group_a, group_b = assign_leaflets(
        u,
        sel_a=args.leaflet_a,
        sel_b=args.leaflet_b,
        bilayer_selection=args.bilayer_selection,
    )
    if len(group_a) == 0 or len(group_b) == 0:
        raise ValueError("Leaflet assignment produced empty atom groups")

    wc_kwargs = dict(alpha=args.alpha, mesh=args.mesh, density_cutoff=args.force_density)
    wc_a = pywc.WillardChandler(u, group=group_a, centered=True, autoassign=False, **wc_kwargs)
    wc_b = pywc.WillardChandler(u, group=group_b, centered=True, autoassign=False, **wc_kwargs)

    grid_n = args.grid_size
    thickness_accum = np.zeros((grid_n, grid_n), dtype=float)
    thickness_counts = np.zeros_like(thickness_accum, dtype=int)
    box_x_sum = 0.0
    box_y_sum = 0.0

    spectra_accum = None
    frames_used = 0

    for ts in u.trajectory[args.start:args.stop:args.stride]:
        lx, ly, lz = map(float, ts.dimensions[:3])
        cutoff_a = ensure_surface(
            wc_a,
            force_reassign=True,
            absolute_level=args.force_density,
            relative_level=args.density_level,
        )
        cutoff_b = ensure_surface(
            wc_b,
            force_reassign=True,
            absolute_level=args.force_density,
            relative_level=args.density_level,
        )

        verts_a = wc_a.triangulated_surface[0]
        verts_b = wc_b.triangulated_surface[0]
        top_map = fill_nan_nearest(rasterise_surface(verts_a, lx, ly, grid_n, reduce="max"))
        bottom_map = fill_nan_nearest(rasterise_surface(verts_b, lx, ly, grid_n, reduce="min"))

        mid_map = 0.5 * (top_map + bottom_map)
        thickness_map = top_map - bottom_map
        thickness_accum += thickness_map
        thickness_counts += 1
        box_x_sum += lx
        box_y_sum += ly

        if args.detrend_plane:
            mid_map = detrend_plane(mid_map, lx, ly)

        mid_map -= np.mean(mid_map)

        fft = np.fft.fft2(mid_map) / (grid_n * grid_n)
        power = np.abs(fft) ** 2

        if spectra_accum is None:
            spectra_accum = power
            qx = 2.0 * math.pi * np.fft.fftfreq(grid_n, d=lx / grid_n)
            qy = 2.0 * math.pi * np.fft.fftfreq(grid_n, d=ly / grid_n)
            QX, QY = np.meshgrid(qy, qx)  # note ordering to match FFT layout
        else:
            spectra_accum += power

        frames_used += 1

    if frames_used == 0 or spectra_accum is None:
        raise RuntimeError("No frames were processed – adjust --start/--stop/--stride")

    spectra_mean = spectra_accum / frames_used

    if args.mode_output:
        np.savez(args.mode_output, power=spectra_mean, qx=QX, qy=QY)

    q_vals, power_vals, counts = radial_average(
        QX, QY, spectra_mean, bins=args.q_bins, min_modes=args.min_modes
    )

    if args.qmax is not None:
        mask = q_vals <= args.qmax
        q_vals = q_vals[mask]
        power_vals = power_vals[mask]
        counts = counts[mask]
    if args.qmin > 0.0:
        mask = q_vals >= args.qmin
        q_vals = q_vals[mask]
        power_vals = power_vals[mask]
        counts = counts[mask]

    sigma, sigma_err, kappa, kappa_err = helfrich_fit(q_vals, power_vals, args.temperature)

    args.spectrum_output.parent.mkdir(parents=True, exist_ok=True)
    with args.spectrum_output.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["q_A^-1", "power_A2", "mode_count"])
        for q, p, c in zip(q_vals, power_vals, counts):
            writer.writerow([q, p, int(c)])

    print(
        f"Processed {frames_used} frames. Bending rigidity κ = {kappa:.4f} ± {kappa_err:.4f} kJ/mol, "
        f"tension σ = {sigma:.4f} ± {sigma_err:.4f} kJ/mol/Å²"
    )
    print(f"Radially averaged spectrum saved to {args.spectrum_output}")

    if thickness_counts.any():
        thickness_mean = thickness_accum / thickness_counts
        lx_mean = box_x_sum / frames_used
        ly_mean = box_y_sum / frames_used
        dx = lx_mean / grid_n
        dy = ly_mean / grid_n
        args.thickness_output.parent.mkdir(parents=True, exist_ok=True)
        with args.thickness_output.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["x_A", "y_A", "thickness_A"])
            for ix in range(grid_n):
                for iy in range(grid_n):
                    x = (ix + 0.5) * dx
                    y = (iy + 0.5) * dy
                    writer.writerow([x, y, thickness_mean[ix, iy]])
        print(f"Mean thickness grid written to {args.thickness_output}")


def cli():
    """Command-line interface entry point."""
    main(build_parser().parse_args())


if __name__ == "__main__":
    cli()
