#!/usr/bin/env python3
"""Compute membrane–solvent interface area and thickness using pytim's Willard–Chandler surface.

This script iterates over a trajectory, builds the Willard–Chandler dividing
surface for a user-provided atom selection, and records the total triangulated
surface area per frame. The `alpha` Gaussian width, mesh spacing, frame/time
window, and selections are configurable from the command line. Additionally, a
planar grid is used to average the upper and lower interfacial surfaces and
report a membrane thickness map.
"""

import argparse
import csv
import os
import time
from pathlib import Path
from typing import Optional

import MDAnalysis as mda
import numpy as np
import pytim
from skimage import measure

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None


def positive_float(value: str) -> float:
    val = float(value)
    if val <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return val


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute Willard–Chandler interface area along a trajectory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--topology",
        required=True,
        help="Topology file understood by MDAnalysis (e.g. TPR, PSF, GRO)",
    )
    parser.add_argument(
        "-x",
        "--trajectory",
        required=True,
        help="Trajectory file (e.g. XTC, TRR, DCD).",
    )
    parser.add_argument(
        "--selection",
        default="all",
        help=(
            "MDAnalysis selection string that defines the atoms used to build "
            "the interface. Combine membrane and solvent selections here if "
            "needed, e.g. 'resname C18 and not name H* or resname SOL'."
        ),
    )
    parser.add_argument(
        "--backend",
        default="cpu",
        choices=["cpu", "cpp", "cupy", "gpu", "python"],
        help=(
            "Willard-Chandler surface computation backend:\n"
            "  cpu/cpp: C++ accelerated (via _wc_kde.cpp) - recommended\n"
            "  cupy/gpu: GPU accelerated (requires CuPy)\n"
            "  python: Pure Python (for testing/debugging)"
        ),
    )
    parser.add_argument(
        "--selection-file",
        type=Path,
        default=None,
        help="Path to a text file containing the MDAnalysis selection string (overrides --selection and --selection-env).",
    )
    parser.add_argument(
        "--selection-env",
        default=None,
        help="Name of an environment variable holding the MDAnalysis selection string (overrides --selection).",
    )
    parser.add_argument(
        "--alpha",
        type=positive_float,
        default=3.0,
        help="Gaussian kernel width (Willard–Chandler α parameter, in Å).",
    )
    parser.add_argument(
        "--mesh",
        type=positive_float,
        default=2.0,
        help="Grid spacing (Å) for the density field discretization.",
    )
    parser.add_argument(
        "--density-cutoff",
        type=float,
        default=None,
        help="Optional absolute density threshold for the isosurface (overrides --density-level).",
    )
    parser.add_argument(
        "--density-level",
        type=float,
        default=1.0,
        help=(
            "Relative density level in [0, 2]: 0 = min density, 1 = midpoint, 2 = max. "
            "Ignored when --density-cutoff is provided."
        ),
    )
    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help="First frame index to analyze (inclusive).",
    )
    parser.add_argument(
        "--stop",
        type=int,
        default=None,
        help="Stop frame index (exclusive). None → full trajectory.",
    )
    parser.add_argument(
        "-b",
        "--start-time",
        type=float,
        default=None,
        help=(
            "Lower bound for the simulation time window in picoseconds. "
            "Frames with a smaller time stamp are skipped."
        ),
    )
    parser.add_argument(
        "-e",
        "--end-time",
        type=float,
        default=None,
        help=(
            "Upper bound for the simulation time window in picoseconds. "
            "Frames with a larger time stamp stop the iteration."
        ),
    )
    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Analyze every Nth frame.",
    )
    parser.add_argument(
        "--center",
        action="store_true",
        help="Center the analysis group in the unit cell before each evaluation.",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=20,
        help="Number of bins along x and y for the thickness grid (set to 0 to disable).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("willard_chandler_area.csv"),
        help="CSV file where per-frame areas will be written.",
    )
    parser.add_argument(
        "--thickness-output",
        type=Path,
        default=Path("willard_chandler_thickness.csv"),
        help="CSV file containing the averaged thickness map (x, y, thickness).",
    )
    parser.add_argument(
        "--no-autoassign",
        action="store_true",
        help=(
            "Disable automatic recomputation upon frame change. When set, the "
            "script does not explicitly reassigns the interface each frame."
        ),
    )
    parser.add_argument(
        "--enable-timing",
        action="store_true",
        help="Enable detailed timing of Willard-Chandler components (prepare_box, define_cluster_group, compute_surface).",
    )
    parser.add_argument(
        "--print-timing",
        action="store_true",
        help="Print detailed timing breakdown at the end of the run. Implies --enable-timing.",
    )
    parser.add_argument(
        "--skip-frames",
        type=int,
        default=2,
        help="Number of initial frames to skip when computing timing statistics (default: 2). "
             "This excludes warm-up overhead from timing measurements.",
    )
    return parser


def format_time(ts) -> Optional[float]:
    try:
        return float(ts.time)
    except AttributeError:
        return None


def ensure_assign(interface: pytim.WillardChandler, force: bool = False) -> None:
    """Recompute the interface if `force` is True or autoassign is disabled."""
    if force or not getattr(interface, "autoassign", True):
        interface._assign_layers()


def resolve_selection(args: argparse.Namespace) -> str:
    selection = args.selection
    if args.selection_env:
        value = os.environ.get(args.selection_env)
        if value is None:
            raise ValueError(
                f"Environment variable {args.selection_env!r} is not set but was requested via --selection-env"
            )
        selection = value.strip()

    if args.selection_file:
        try:
            selection = args.selection_file.read_text(encoding="utf-8").strip()
        except OSError as exc:
            raise ValueError(
                f"Failed to read selection from file {args.selection_file}: {exc}"
            ) from exc

    if not selection or selection.isspace():
        raise ValueError(
            "The resolved selection string is empty; please provide a valid selection expression"
        )

    return selection


def compute_frame_extrema(vertices: np.ndarray, lx: float, ly: float, grid_size: int):
    """Return per-cell extrema of the triangulated surface for the XY grid."""

    top = np.full((grid_size, grid_size), np.nan)
    bottom = np.full((grid_size, grid_size), np.nan)

    if vertices.size == 0 or lx <= 0.0 or ly <= 0.0:
        return top, bottom

    x_coords = np.mod(vertices[:, 0], lx)
    y_coords = np.mod(vertices[:, 1], ly)
    z_coords = vertices[:, 2]

    ix = np.floor(x_coords / lx * grid_size).astype(int)
    iy = np.floor(y_coords / ly * grid_size).astype(int)
    np.clip(ix, 0, grid_size - 1, out=ix)
    np.clip(iy, 0, grid_size - 1, out=iy)

    for xi, yi, zi in zip(ix, iy, z_coords):
        current_top = top[xi, yi]
        if np.isnan(current_top) or zi > current_top:
            top[xi, yi] = zi
        current_bottom = bottom[xi, yi]
        if np.isnan(current_bottom) or zi < current_bottom:
            bottom[xi, yi] = zi

    return top, bottom


def apply_density_level(
    interface: pytim.WillardChandler, level: float
) -> Optional[float]:
    """Reconstruct the triangulated surface for a given relative density level."""

    density_field = getattr(interface, "density_field", None)
    spacing = getattr(interface, "spacing", None)
    ngrid = getattr(interface, "ngrid", None)
    if density_field is None or spacing is None or ngrid is None:
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

    volume = density_field.reshape(tuple(np.array(ngrid[::-1]).astype(int)))
    verts, faces, normals, _ = measure.marching_cubes(
        volume, cutoff, spacing=tuple(spacing)
    )
    interface.triangulated_surface = [
        np.fliplr(verts), faces, np.fliplr(normals)]
    interface.surface_area = measure.mesh_surface_area(verts, faces)
    interface.density_cutoff = cutoff
    return cutoff


def fill_missing_with_nearest(data: np.ndarray) -> np.ndarray:
    """Fill NaN cells with nearest valid neighbor (Manhattan distance)."""

    if np.all(np.isnan(data)):
        return data
    mask = ~np.isnan(data)
    if np.all(mask):
        return data

    filled = data.copy()
    valid_coords = np.array(np.nonzero(mask)).T
    valid_values = data[mask]

    for idx in zip(*np.where(~mask)):
        distances = np.sum((valid_coords - idx) ** 2, axis=1)
        nearest = valid_values[np.argmin(distances)]
        filled[idx] = nearest
    return filled


def enforce_time_metadata(
    ts, *, start_time: Optional[float], end_time: Optional[float]
):
    if (start_time is not None or end_time is not None) and format_time(ts) is None:
        raise ValueError(
            "Trajectory does not supply time information but --start-time/--end-time were requested."
        )


def main(args: argparse.Namespace) -> None:
    universe = mda.Universe(args.topology, args.trajectory)
    selection_expr = resolve_selection(args)
    analysis_group = universe.select_atoms(selection_expr)
    if len(analysis_group) == 0:
        raise ValueError(
            "The provided selection did not match any atoms. Please adjust --selection."
        )

    if args.step <= 0:
        raise ValueError("--step must be a positive integer")
    if args.grid_size < 0:
        raise ValueError("--grid-size must be non-negative")
    if args.density_level < 0.0 or args.density_level > 2.0:
        raise ValueError("--density-level must lie within [0, 2]")

    use_density_level = args.density_cutoff is None

    # Always enable timing to provide performance feedback
    interface = pytim.WillardChandler(
        universe,
        group=analysis_group,
        alpha=args.alpha,
        mesh=args.mesh,
        centered=args.center,
        density_cutoff=args.density_cutoff,
        autoassign=args.no_autoassign,
        surface_backend=args.backend,
        centering_backend='cpu',
        enable_timing=True,
    )

    ensure_assign(interface, force=True)
    density_cutoff_value: Optional[float] = None
    if use_density_level:
        density_cutoff_value = apply_density_level(
            interface, args.density_level)

    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "frame",
        "time_ps",
        "area_angstrom2",
        "box_x",
        "box_y",
        "box_z",
        "plane_area",
        "density_cutoff",
        "alpha",
        "mesh",
        "upper_rms",
        "lower_rms",
        "upper_peak_to_valley",
        "lower_peak_to_valley",
    ]

    records = []
    frame_counter = -1

    grid_size = args.grid_size
    accumulate_thickness = grid_size > 0

    if accumulate_thickness:
        grid_shape = (grid_size, grid_size)
        top_sum = np.zeros(grid_shape)
        bottom_sum = np.zeros(grid_shape)
        thickness_counts = np.zeros(grid_shape, dtype=int)
        x_sum = np.zeros(grid_shape)
        y_sum = np.zeros(grid_shape)
        box_x_sum = 0.0
        box_y_sum = 0.0
        box_samples = 0
        upper_rms_values = []
        lower_rms_values = []
        upper_pv_values = []
        lower_pv_values = []
    else:
        upper_rms_values = []
        lower_rms_values = []
        upper_pv_values = []
        lower_pv_values = []
    tt0 = universe.trajectory[1]
    my_step = int(args.start_time / float(tt0.time) - 1)

    for ts in universe.trajectory[my_step:]:
        if args.start is not None and ts.frame < args.start:
            continue
        if args.stop is not None and ts.frame >= args.stop:
            break

        enforce_time_metadata(
            ts, start_time=args.start_time, end_time=args.end_time)
        time_ps = format_time(ts)
        if args.start_time is not None and time_ps is not None:
            if time_ps < args.start_time:
                continue
        if args.end_time is not None and time_ps is not None:
            if time_ps > args.end_time:
                break

        frame_counter += 1
        if args.step > 1 and frame_counter % args.step != 0:
            continue
        ensure_assign(interface, force=args.no_autoassign)
        if use_density_level:
            level_value = apply_density_level(interface, args.density_level)
            if level_value is not None:
                density_cutoff_value = level_value
        area = float(interface.surface_area)
        upper_rms = None
        lower_rms = None
        upper_pv = None
        lower_pv = None

        if accumulate_thickness:
            lx, ly = float(ts.dimensions[0]), float(ts.dimensions[1])
            top_frame, bottom_frame = compute_frame_extrema(
                interface.triangulated_surface[0], lx, ly, grid_size
            )
            valid = ~np.isnan(top_frame) & ~np.isnan(bottom_frame)
            if np.any(valid):
                top_sum[valid] += top_frame[valid]
                bottom_sum[valid] += bottom_frame[valid]
                thickness_counts[valid] += 1

                x_centers = (np.arange(grid_size) + 0.5) * (lx / grid_size)
                y_centers = (np.arange(grid_size) + 0.5) * (ly / grid_size)
                x_center_grid = np.broadcast_to(x_centers[:, None], grid_shape)
                y_center_grid = np.broadcast_to(y_centers[None, :], grid_shape)
                x_sum[valid] += x_center_grid[valid]
                y_sum[valid] += y_center_grid[valid]

                box_x_sum += lx
                box_y_sum += ly
                box_samples += 1

                top_values = top_frame[valid]
                bottom_values = bottom_frame[valid]
                if top_values.size > 0:
                    top_mean = float(np.mean(top_values))
                    upper_rms = float(
                        np.sqrt(np.mean((top_values - top_mean) ** 2)))
                    upper_pv = float(np.max(top_values) - np.min(top_values))
                    upper_rms_values.append(upper_rms)
                    upper_pv_values.append(upper_pv)
                if bottom_values.size > 0:
                    bottom_mean = float(np.mean(bottom_values))
                    lower_rms = float(
                        np.sqrt(np.mean((bottom_values - bottom_mean) ** 2))
                    )
                    lower_pv = float(np.max(bottom_values) -
                                     np.min(bottom_values))
                    lower_rms_values.append(lower_rms)
                    lower_pv_values.append(lower_pv)

        plane_area = float(ts.dimensions[0]) * float(ts.dimensions[1])
        record = {
            "frame": ts.frame,
            "time_ps": time_ps,
            "area_angstrom2": area,
            "box_x": float(ts.dimensions[0]),
            "box_y": float(ts.dimensions[1]),
            "box_z": float(ts.dimensions[2]),
            "plane_area": plane_area,
            "density_cutoff": density_cutoff_value
            if density_cutoff_value is not None
            else interface.density_cutoff,
            "alpha": args.alpha,
            "mesh": args.mesh,
            "upper_rms": upper_rms,
            "lower_rms": lower_rms,
            "upper_peak_to_valley": upper_pv,
            "lower_peak_to_valley": lower_pv,
        }
        records.append(record)
        print(
            f"frame {ts.frame:6d} time {record['time_ps']!s:>10} ps → "
            f"area {area:12.3f} Å² (plane {2 * plane_area:12.3f} Å²)"
        )

    # Print total time from WillardChandler's timing (skip initial frames)
    timings = interface.get_detailed_timings(skip_frames=args.skip_frames)
    if timings and 'total' in timings:
        stats = timings['total']
        if args.skip_frames > 0:
            print(
                f"Total time: {stats['mean']*1000:.2f} ms (average per frame, excluding first {args.skip_frames} frame(s))")
        else:
            print(
                f"Total time: {stats['mean']*1000:.2f} ms (average per frame)")
    with output_path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    if records:
        areas = np.array([row["area_angstrom2"] for row in records])
        plane_areas = np.array([row["plane_area"] for row in records])
        print(
            "\nCompleted. Area statistics (Å²): "
            f"mean={areas.mean():.3f}, std={areas.std(ddof=1) if len(areas) > 1 else 0.0:.3f}, "
            f"min={areas.min():.3f}, max={areas.max():.3f}"
        )
        print(f"Per-frame data saved to {output_path}")
        print(
            "Plane area statistics (Å²): "
            f"mean={plane_areas.mean():.3f}, std={plane_areas.std(ddof=1) if len(plane_areas) > 1 else 0.0:.3f}, "
            f"min={plane_areas.min():.3f}, max={plane_areas.max():.3f}"
        )
    else:
        print("No frames processed; check start/stop/step parameters.")

    if accumulate_thickness and box_samples > 0:
        avg_top = np.full_like(top_sum, np.nan, dtype=float)
        avg_bottom = np.full_like(bottom_sum, np.nan, dtype=float)
        avg_x = np.full_like(top_sum, np.nan, dtype=float)
        avg_y = np.full_like(bottom_sum, np.nan, dtype=float)
        valid = thickness_counts > 0
        avg_top[valid] = top_sum[valid] / thickness_counts[valid]
        avg_bottom[valid] = bottom_sum[valid] / thickness_counts[valid]
        avg_x[valid] = x_sum[valid] / thickness_counts[valid]
        avg_y[valid] = y_sum[valid] / thickness_counts[valid]
        thickness_map = avg_top - avg_bottom
        mean_thickness = np.nanmean(thickness_map)
        thickness_map_filled = fill_missing_with_nearest(thickness_map)

        thickness_output = args.thickness_output
        thickness_output.parent.mkdir(parents=True, exist_ok=True)

        lx_mean = box_x_sum / box_samples
        ly_mean = box_y_sum / box_samples
        dx = lx_mean / grid_size
        dy = ly_mean / grid_size

        with thickness_output.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["x", "y", "thickness"])
            for ix in range(grid_size):
                for iy in range(grid_size):
                    x = avg_x[ix, iy]
                    y = avg_y[ix, iy]
                    if np.isnan(x):
                        x = (ix + 0.5) * dx
                    if np.isnan(y):
                        y = (iy + 0.5) * dy
                    thickness = thickness_map_filled[ix, iy]
                    writer.writerow([x, y, thickness])

        print(
            "Thickness grid written to {} (mean thickness {:.3f} Å)".format(
                thickness_output, mean_thickness
            )
        )

        if upper_rms_values:
            print(
                "Upper surface roughness: mean RMS {:.3f} Å, mean P-V {:.3f} Å".format(
                    float(np.mean(upper_rms_values)), float(
                        np.mean(upper_pv_values))
                )
            )
        if lower_rms_values:
            print(
                "Lower surface roughness: mean RMS {:.3f} Å, mean P-V {:.3f} Å".format(
                    float(np.mean(lower_rms_values)), float(
                        np.mean(lower_pv_values))
                )
            )
        if not upper_rms_values and not lower_rms_values:
            print("No roughness statistics accumulated; check grid coverage.")
    elif accumulate_thickness:
        print("No thickness data accumulated; check grid size or trajectory content.")

    # Print timing breakdown if requested
    if args.print_timing:
        interface.print_timing_breakdown(skip_frames=args.skip_frames)


def cli():
    """Command-line interface entry point."""
    main(build_parser().parse_args())


if __name__ == "__main__":
    cli()
