#!/usr/bin/env python3
"""Plot the Willard–Chandler triangulated surface and save to an image file.

This script loads a topology and trajectory, moves to a user-specified frame,
computes the Willard–Chandler dividing surface, and renders it as a 3D plot
using matplotlib. The result is saved to a PNG (or other format) file.
"""

import argparse
import os
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import pywc


def positive_float(value: str) -> float:
    val = float(value)
    if val <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return val


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot the Willard–Chandler surface for a given frame and save to an image",
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
        help="MDAnalysis selection string for the atoms used to build the interface.",
    )
    parser.add_argument(
        "--selection-file",
        type=Path,
        default=None,
        help="Path to a text file containing the MDAnalysis selection string (overrides --selection).",
    )
    parser.add_argument(
        "--selection-env",
        default=None,
        help="Name of an environment variable holding the MDAnalysis selection string (overrides --selection).",
    )
    parser.add_argument(
        "--frame",
        type=int,
        default=0,
        help="Trajectory frame index to visualize.",
    )
    parser.add_argument(
        "--alpha",
        type=positive_float,
        default=3.0,
        help="Gaussian kernel width (Willard–Chandler alpha parameter, in Angstrom).",
    )
    parser.add_argument(
        "--mesh",
        type=positive_float,
        default=2.0,
        help="Grid spacing (Angstrom) for the density field discretization.",
    )
    parser.add_argument(
        "--density-cutoff",
        type=float,
        default=None,
        help="Absolute density threshold for the isosurface.",
    )
    parser.add_argument(
        "--backend",
        default="cpu",
        choices=["cpu", "cpp", "cupy", "gpu", "python"],
        help="Willard-Chandler surface computation backend.",
    )
    parser.add_argument(
        "--center",
        action="store_true",
        help="Center the analysis group in the unit cell before evaluation.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("wc_surface.png"),
        help="Output image file (format inferred from extension: png, pdf, svg, eps, ...).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolution of the output image in dots per inch.",
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=[10, 8],
        metavar=("WIDTH", "HEIGHT"),
        help="Figure size in inches.",
    )
    parser.add_argument(
        "--elev",
        type=float,
        default=25,
        help="Elevation angle (degrees) for the 3D view.",
    )
    parser.add_argument(
        "--azim",
        type=float,
        default=-60,
        help="Azimuth angle (degrees) for the 3D view.",
    )
    parser.add_argument(
        "--colormap",
        default="coolwarm",
        help="Matplotlib colormap for coloring the surface by z-coordinate.",
    )
    parser.add_argument(
        "--surface-alpha",
        type=float,
        default=0.8,
        help="Opacity of the surface (0 = transparent, 1 = opaque).",
    )
    parser.add_argument(
        "--show-edges",
        action="store_true",
        help="Draw triangle edges on the surface.",
    )
    parser.add_argument(
        "--show-box",
        action="store_true",
        help="Draw the simulation box outline.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom plot title. If not set, an automatic title with frame and area info is used.",
    )
    return parser


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


def draw_box_edges(ax, box):
    """Draw the edges of the orthorhombic simulation box."""
    lx, ly, lz = box[0], box[1], box[2]
    corners = np.array([
        [0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0],
        [0, 0, lz], [lx, 0, lz], [lx, ly, lz], [0, ly, lz],
    ])
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for i, j in edges:
        ax.plot3D(
            *zip(corners[i], corners[j]),
            color="gray", linewidth=0.5, linestyle="--",
        )


def main(args: argparse.Namespace) -> None:
    universe = mda.Universe(args.topology, args.trajectory)
    selection_expr = resolve_selection(args)
    analysis_group = universe.select_atoms(selection_expr)
    if len(analysis_group) == 0:
        raise ValueError(
            "The provided selection did not match any atoms. Please adjust --selection."
        )

    n_frames = len(universe.trajectory)
    if args.frame < 0 or args.frame >= n_frames:
        raise ValueError(
            f"Frame {args.frame} is out of range. Trajectory has {n_frames} frames (0-{n_frames - 1})."
        )

    universe.trajectory[args.frame]
    time_ps = universe.trajectory.ts.time
    print(f"Frame {args.frame} loaded (time = {time_ps:.1f} ps)")

    interface = pywc.WillardChandler(
        universe,
        group=analysis_group,
        alpha=args.alpha,
        mesh=args.mesh,
        centered=args.center,
        density_cutoff=args.density_cutoff,
        autoassign=True,
        surface_backend=args.backend,
    )

    vertices = interface.triangulated_surface[0]
    faces = interface.triangulated_surface[1]
    area = interface.surface_area
    box = universe.dimensions[:3]

    print(
        f"Surface computed: {len(vertices)} vertices, {len(faces)} triangles, "
        f"area = {area:.2f} Ang^2"
    )

    # Build the polygon collection for plotting
    tri_verts = vertices[faces]
    z_centers = tri_verts[:, :, 2].mean(axis=1)

    cmap = plt.get_cmap(args.colormap)
    z_min, z_max = z_centers.min(), z_centers.max()
    if np.isclose(z_min, z_max):
        z_norm = np.zeros_like(z_centers)
    else:
        z_norm = (z_centers - z_min) / (z_max - z_min)
    face_colors = cmap(z_norm)
    face_colors[:, 3] = args.surface_alpha

    fig = plt.figure(figsize=tuple(args.figsize))
    ax = fig.add_subplot(111, projection="3d")

    poly = Poly3DCollection(
        tri_verts,
        facecolors=face_colors,
        edgecolors="k" if args.show_edges else "none",
        linewidths=0.1 if args.show_edges else 0,
    )
    ax.add_collection3d(poly)

    if args.show_box:
        draw_box_edges(ax, box)

    # Set axis limits
    margin = 2.0
    ax.set_xlim(vertices[:, 0].min() - margin, vertices[:, 0].max() + margin)
    ax.set_ylim(vertices[:, 1].min() - margin, vertices[:, 1].max() + margin)
    ax.set_zlim(vertices[:, 2].min() - margin, vertices[:, 2].max() + margin)

    ax.set_xlabel("x (Ang)")
    ax.set_ylabel("y (Ang)")
    ax.set_zlabel("z (Ang)")
    ax.view_init(elev=args.elev, azim=args.azim)

    if args.title is not None:
        ax.set_title(args.title)
    else:
        ax.set_title(
            f"Willard-Chandler surface — frame {args.frame} "
            f"(t = {time_ps:.1f} ps, A = {area:.1f} Ang²)"
        )

    # Add a colorbar for z-coordinate
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=z_min, vmax=z_max)
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("z (Ang)")

    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Image saved to {output_path}")


def cli():
    """Command-line interface entry point."""
    main(build_parser().parse_args())


if __name__ == "__main__":
    cli()
