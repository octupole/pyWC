#!/usr/bin/env python3
"""Write the Willard–Chandler triangulated surface to a VTK file for a given trajectory frame.

This script loads a topology and trajectory, moves to a user-specified frame,
computes the Willard–Chandler dividing surface, and writes the result as a
VTK unstructured grid that can be visualized with ParaView, VisIt, or similar tools.

Optionally, the density field and/or atom positions can also be written to
separate VTK files.
"""

import argparse
import os
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pywc


def positive_float(value: str) -> float:
    val = float(value)
    if val <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return val


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write the Willard–Chandler surface to VTK for a given frame",
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
        default=Path("wc_surface.vtk"),
        help="Output VTK file for the triangulated surface.",
    )
    parser.add_argument(
        "--density-output",
        type=Path,
        default=None,
        help="If provided, also write the density field to this VTK file.",
    )
    parser.add_argument(
        "--particles-output",
        type=Path,
        default=None,
        help="If provided, also write the atom positions to this VTK file.",
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

    # Move to the requested frame
    universe.trajectory[args.frame]
    print(f"Frame {args.frame} loaded (time = {universe.trajectory.ts.time:.1f} ps)")

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

    # Write the triangulated surface
    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    interface.writevtk.surface(str(output_path))
    vertices = interface.triangulated_surface[0]
    faces = interface.triangulated_surface[1]
    print(
        f"Surface written to {output_path} "
        f"({len(vertices)} vertices, {len(faces)} triangles, "
        f"area = {interface.surface_area:.2f} Ang^2)"
    )

    # Optionally write the density field
    if args.density_output is not None:
        args.density_output.parent.mkdir(parents=True, exist_ok=True)
        interface.writevtk.density(str(args.density_output))
        print(f"Density field written to {args.density_output}")

    # Optionally write atom positions
    if args.particles_output is not None:
        args.particles_output.parent.mkdir(parents=True, exist_ok=True)
        interface.writevtk.particles(str(args.particles_output), group=analysis_group)
        print(f"Particle positions written to {args.particles_output}")


def cli():
    """Command-line interface entry point."""
    main(build_parser().parse_args())


if __name__ == "__main__":
    cli()
