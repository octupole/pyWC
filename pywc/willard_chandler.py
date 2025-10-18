#!/usr/bin/python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
""" Module: willard-chandler
    ========================
"""

from __future__ import print_function
import numpy as np
import time

from . import messages
from . import utilities, cube, wavefront_obj
from .wc_core.surface import compute_surface
from .sanity_check import SanityCheck
from .vtk import Writevtk

from .interface import Interface
from .patches import patchTrajectory, patchOpenMM, patchMDTRAJ

np.set_printoptions(legacy=False)  # fixes problem with skimage


class WillardChandler(Interface):
    """ Identifies the dividing surface using the Willard-Chandler method
        NOTE that this method does *not* identify surface atoms

        *(Willard, A. P.; Chandler, D. J. Phys. Chem. B 2010, 114, 1954–1958)*

        :param Object universe:   The MDAnalysis_ Universe, MDTraj_ trajectory
                                  or OpenMM_ Simulation objects.
        :param Object group:      An AtomGroup, or an array-like object with
                                  the indices of the atoms in the group.
                                  Will identify the interfacial molecules from
                                  this group
        :param float alpha:       The width of the Gaussian kernel
        :param float mesh:        The grid spacing for the density calculation
        :param float density_cutoff: The density value used to define the
                                  isosurface. `None` (default) uses the average
                                  of the minimum and maximum density.
        :param AtomGroup group:   Compute the density using this group
        :param dict radii_dict:   Dictionary with the atomic radii of
                                  the elements in the group.
                                  If None is supplied, the default one
                                  (from GROMOS 43a1) will be used.
        :param float cluster_cut: Cutoff used for neighbors or density-based
                                  cluster search (default: None disables the
                                  cluster analysis)
        :param float cluster_threshold_density: Number density threshold for
                                  the density-based cluster search. 'auto'
                                  determines the threshold automatically.
                                  Default: None uses simple neighbors cluster
                                  search, if cluster_cut is not None
        :param Object extra_cluster_groups: Additional groups, to allow for
                                  mixed interfaces

        :param bool include_zero_radius: if false (default) exclude atoms with zero radius
                                  from the surface analysis (they are always included
                                  in the cluster search, if present in the relevant
                                  group) to avoid some artefacts.
        :param bool centered:     Center the  :py:obj:`group`
        :param bool warnings:     Print warnings
        :param str surface_backend: Backend for density computation: 'cpu' (default) or 'cupy'
        :param dict surface_backend_options: Backend-specific options (e.g., chunk_size, max_chunk_mem for cupy)
        :param bool enable_timing: Enable timing measurements for surface computation (default: False)

        Example:

        >>> import MDAnalysis as mda
        >>> import pywc
        >>> from pywc.datafiles import *
        >>>
        >>> u = mda.Universe(MICELLE_PDB)
        >>> g = u.select_atoms('resname DPC')
        >>>
        >>> radii = pywc_data.vdwradii(G43A1_TOP)
        >>>
        >>> inter= pywc.WillardChandler(u, group=g, alpha=3.0, fast=True)
        >>> R, _, _, _ = pywc.utilities.fit_sphere(inter.triangulated_surface[0])
        >>> print ("Radius={:.3f}".format(R))
        Radius=19.970


        .. _MDAnalysis: http://www.mdanalysis.org/
        .. _MDTraj: http://www.mdtraj.org/
        .. _OpenMM: http://www.openmm.org/
    """

    _surface = None
    _surface_computation_time = None
    _backend_timings = {}
    _detailed_timings = {}

    @property
    def layers(self):
        """ The method does not identify layers.

        Example:

        >>> import MDAnalysis as mda
        >>> import pywc
        >>> from pywc.datafiles import *
        >>>
        >>> u = mda.Universe(MICELLE_PDB)
        >>> g = u.select_atoms('resname DPC')
        >>> inter= pywc.WillardChandler(u, group=g, alpha=3.0, fast=True)
        >>> inter.layers
        <AtomGroup with 0 atoms>
        """
        return self._layers

    def _sanity_checks(self):
        """ Basic checks to be performed after the initialization.

            >>> import pytest
            >>> with pytest.raises(Exception):
            ...     pywc.WillardChandler(u,mesh=-1)

        """

    def __init__(self,
                 universe,
                 group=None,
                 alpha=2.0,
                 radii_dict=None,
                 mesh=2.0,
                 symmetry='spherical',
                 cluster_cut=None,
                 include_zero_radius=False,
                 cluster_threshold_density=None,
                 extra_cluster_groups=None,
                 centered=False,
                 warnings=False,
                 autoassign=True,
                 density_cutoff=None,
                 surface_backend='cpu',
                 surface_backend_options=None,
                 centering_backend='cpu',
                 **kargs):

        self.autoassign, self.do_center = autoassign, centered
        self.include_zero_radius = include_zero_radius
        self.density_cutoff = density_cutoff
        self.surface_backend = surface_backend
        self.surface_backend_options = surface_backend_options or {}
        self.centering_backend = centering_backend
        self._enable_timing = kargs.get('enable_timing', False)
        if kargs:
            for key in list(kargs):
                setattr(self, key, kargs.pop(key))
        sanity = SanityCheck(self, warnings=warnings)
        sanity.assign_universe(universe, group)
        sanity.assign_alpha(alpha)

        if mesh <= 0:
            raise ValueError(messages.MESH_NEGATIVE)
        self.mesh, self.spacing, self.ngrid, self.PDB = mesh, None, None, {}

        sanity.assign_radii(radii_dict=radii_dict)

        sanity.assign_cluster_params(cluster_cut,
                                     cluster_threshold_density, extra_cluster_groups)

        self._assign_symmetry(symmetry)

        patchTrajectory(self.universe.trajectory, self)
        self._assign_layers()
        self._atoms = self._layers[:]  # this is an empty AtomGroup
        self.writevtk = Writevtk(self)

    def writecube(self, filename="pywc.cube", group=None, sequence=False, order='zyx', shift=[0,0,0], normalize=True):
        """ Write to cube files (sequences) the volumentric density and the
            atomic positions.

            :param str filename  : the file name
            :param bool sequence : if true writes a sequence of files adding
                                   the frame to the filename
            :param string order  : 'xyz' or 'zyx', to adapt to the software used
                                   for visualization. Default: 'zyx'
            :param array shift   : add shift in grid units. Default: [0,0,0]
            :param bool normalize: if true normalizes the density field to the range [0, 1]

            >>> import MDAnalysis as mda
            >>> import pywc
            >>> from pywc.datafiles import MICELLE_PDB
            >>> u = mda.Universe(MICELLE_PDB)
            >>> g = u.select_atoms('resname DPC')
            >>> inter= pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0)
            >>> inter.writecube('dens.cube') # writes on dens.cube
            >>> inter.writecube('dens.cube',group=g) # writes also  particles
            >>> inter.writecube('dens.cube',sequence=True) # dens.<frame>.cube
            >>> inter.writecube('dens.cube', normalize=False) # writes density values without normalization
        """
        if sequence is True:
            filename = cube.consecutive_filename(self.universe, filename)
        # TODO handle optional atomic_numbers
        cube.write_file(
            filename,
            group,
            self.ngrid,
            self.spacing,
            self.density_field,
            atomic_numbers=None,
            order=order,
            normalize=normalize)

    def writeobj(self, filename="pywc.obj", sequence=False):
        """ Write to wavefront obj files (sequences) the triangulated surface

            :param str filename:  the file name
            :param bool sequence: if true writes a sequence of files adding
                                  the frame to the filename

            >>> import MDAnalysis as mda
            >>> import pywc
            >>> from pywc.datafiles import MICELLE_PDB
            >>> u = mda.Universe(MICELLE_PDB)
            >>> g = u.select_atoms('resname DPC')
            >>> inter= pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0)
            >>> inter.writeobj('surf.obj') # writes on surf.obj
            >>> inter.writeobj('surf.obj',sequence=True) # surf.<frame>.obj
        """

        if sequence is True:
            filename = wavefront_obj.consecutive_filename(
                self.universe, filename)

        vert, surf = list(self.triangulated_surface[0:2])
        wavefront_obj.write_file(filename, vert, surf)

    def get_timing(self):
        """ Get the computation time for the last surface calculation.

            :return: Time in seconds, or None if timing was not enabled

            Example:

            >>> import MDAnalysis as mda  # doctest: +SKIP
            >>> import pywc  # doctest: +SKIP
            >>> from pywc.datafiles import MICELLE_PDB  # doctest: +SKIP
            >>> u = mda.Universe(MICELLE_PDB)  # doctest: +SKIP
            >>> g = u.select_atoms('resname DPC')  # doctest: +SKIP
            >>> inter = pywc.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)  # doctest: +SKIP
            >>> print(f"Computation time: {inter.get_timing():.3f} s")  # doctest: +SKIP
        """
        return self._surface_computation_time

    def get_detailed_timings(self, skip_frames=0):
        """ Get detailed timing breakdown for all components.

            :param int skip_frames: Number of initial frames to skip in statistics (default: 0).
                                   Useful to exclude warm-up overhead.
            :return: Dictionary with timing statistics for each component

            Example:

            >>> import MDAnalysis as mda  # doctest: +SKIP
            >>> import pywc  # doctest: +SKIP
            >>> from pywc.datafiles import MICELLE_PDB  # doctest: +SKIP
            >>> u = mda.Universe(MICELLE_PDB)  # doctest: +SKIP
            >>> g = u.select_atoms('resname DPC')  # doctest: +SKIP
            >>> inter = pywc.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)  # doctest: +SKIP
            >>> # Process multiple frames
            >>> for ts in u.trajectory[:10]:  # doctest: +SKIP
            ...     inter._assign_layers()  # doctest: +SKIP
            >>> timings = inter.get_detailed_timings(skip_frames=2)  # Skip first 2 frames  # doctest: +SKIP
            >>> print(f"Average compute_surface time: {timings['compute_surface']['mean']:.4f} s")  # doctest: +SKIP
        """
        if not self._enable_timing or not self._detailed_timings:
            return None

        import numpy as np
        stats = {}
        for component, times in self._detailed_timings.items():
            if times:
                # Skip initial frames if requested
                times_to_use = times[skip_frames:] if skip_frames > 0 else times
                if not times_to_use:
                    continue  # Skip if no data remains after skipping

                times_array = np.array(times_to_use)
                stats[component] = {
                    'mean': float(np.mean(times_array)),
                    'std': float(np.std(times_array)),
                    'min': float(np.min(times_array)),
                    'max': float(np.max(times_array)),
                    'total': float(np.sum(times_array)),
                    'count': len(times_array),
                    'skipped': skip_frames,
                    'times': times_to_use  # Keep raw data for further analysis
                }
        return stats

    def print_timing_breakdown(self, verbose=True, skip_frames=0):
        """ Print a formatted breakdown of timing statistics.

            :param bool verbose: Include detailed statistics (default: True)
            :param int skip_frames: Number of initial frames to skip in statistics (default: 0).
                                   Useful to exclude warm-up overhead.

            Example:

            >>> import MDAnalysis as mda  # doctest: +SKIP
            >>> import pywc  # doctest: +SKIP
            >>> from pywc.datafiles import MICELLE_PDB  # doctest: +SKIP
            >>> u = mda.Universe(MICELLE_PDB)  # doctest: +SKIP
            >>> g = u.select_atoms('resname DPC')  # doctest: +SKIP
            >>> inter = pywc.WillardChandler(u, group=g, alpha=3.0, enable_timing=True)  # doctest: +SKIP
            >>> for ts in u.trajectory[:10]:  # doctest: +SKIP
            ...     inter._assign_layers()  # doctest: +SKIP
            >>> inter.print_timing_breakdown(skip_frames=2)  # Skip first 2 frames  # doctest: +SKIP
        """
        timings = self.get_detailed_timings(skip_frames=skip_frames)
        if not timings:
            print("No timing data available. Enable timing with enable_timing=True")
            return

        print("\n" + "="*70)
        print("Willard-Chandler Timing Breakdown")
        if skip_frames > 0:
            # Get total frame count from any component
            total_frames = next(iter(timings.values()))['count'] + skip_frames
            print(f"(excluding first {skip_frames} frame(s) for warm-up, showing {total_frames - skip_frames}/{total_frames} frames)")
        print("="*70)

        # Sort by mean time (descending)
        sorted_components = sorted(timings.items(),
                                   key=lambda x: x[1]['mean'],
                                   reverse=True)

        # Calculate total for percentage
        total_mean = sum(t['mean'] for c, t in sorted_components if c != 'total')

        for component, stats in sorted_components:
            if component == 'total':
                continue  # Skip total, we'll show it separately

            mean_ms = stats['mean'] * 1000
            percentage = (stats['mean'] / total_mean * 100) if total_mean > 0 else 0

            print(f"\n{component:25s}")
            print(f"  Mean time: {mean_ms:8.3f} ms  ({percentage:5.1f}%)")

            if verbose:
                print(f"  Std dev:   {stats['std']*1000:8.3f} ms")
                print(f"  Min/Max:   {stats['min']*1000:8.3f} / {stats['max']*1000:8.3f} ms")
                print(f"  Count:     {stats['count']:8d} calls")

        # Show total
        if 'total' in timings:
            total_stats = timings['total']
            print(f"\n{'TOTAL':25s}")
            print(f"  Mean time: {total_stats['mean']*1000:8.3f} ms")
            if verbose:
                print(f"  Total time: {total_stats['total']:7.3f} s")
                print(f"  Frames:     {total_stats['count']:7d}")

        print("="*70 + "\n")

    def compare_backends(self, backends=['cpu', 'cupy'], backend_options=None, verbose=True):
        """ Compare performance of different backends for surface computation.

            :param list backends: List of backend names to compare (default: ['cpu', 'cupy'])
            :param dict backend_options: Dictionary mapping backend names to their options
            :param bool verbose: Print comparison results (default: True)
            :return: Dictionary with backend names as keys and timing results as values

            Example:

            >>> import MDAnalysis as mda  # doctest: +SKIP
            >>> import pywc  # doctest: +SKIP
            >>> from pywc.datafiles import MICELLE_PDB  # doctest: +SKIP
            >>> u = mda.Universe(MICELLE_PDB)  # doctest: +SKIP
            >>> g = u.select_atoms('resname DPC')  # doctest: +SKIP
            >>> inter = pywc.WillardChandler(u, group=g, alpha=3.0, mesh=2.0)  # doctest: +SKIP
            >>> timings = inter.compare_backends(['cpu', 'cupy'])  # doctest: +SKIP
            >>> # Backend comparison:
            >>> # cpu: 0.123 s
            >>> # cupy: 0.045 s (2.73x speedup)
        """
        if backend_options is None:
            backend_options = {}

        results = {}
        original_backend = self.surface_backend
        original_options = self.surface_backend_options

        # Get positions and box for comparison
        if self.include_zero_radius:
            pos = self.cluster_group.positions
        else:
            pos = self.cluster_group.positions[self.cluster_group.radii > 0.0]
        box = self.universe.dimensions[:3]

        for backend in backends:
            opts = backend_options.get(backend, {})
            try:
                start_time = time.perf_counter()
                result = compute_surface(pos,
                                       box,
                                       mesh=self.mesh,
                                       sigma=self.alpha,
                                       density_cutoff=self.density_cutoff,
                                       order='xyz',
                                       backend=backend,
                                       backend_options=opts)
                elapsed_time = time.perf_counter() - start_time
                results[backend] = {
                    'time': elapsed_time,
                    'success': True,
                    'error': None
                }
                self._backend_timings[backend] = elapsed_time
            except Exception as e:
                results[backend] = {
                    'time': None,
                    'success': False,
                    'error': str(e)
                }

        # Restore original backend
        self.surface_backend = original_backend
        self.surface_backend_options = original_options

        if verbose:
            print("\n" + "="*60)
            print("Backend Performance Comparison")
            print("="*60)

            # Find baseline (first successful backend)
            baseline_time = None
            for backend in backends:
                if results[backend]['success']:
                    baseline_time = results[backend]['time']
                    break

            for backend in backends:
                result = results[backend]
                if result['success']:
                    elapsed = result['time']
                    if baseline_time and elapsed > 0:
                        speedup = baseline_time / elapsed
                        speedup_str = f" ({speedup:.2f}x)" if speedup != 1.0 else ""
                    else:
                        speedup_str = ""
                    print(f"{backend:12s}: {elapsed:.4f} s{speedup_str}")
                else:
                    print(f"{backend:12s}: FAILED - {result['error']}")
            print("="*60 + "\n")

        return results

    def _assign_layers(self):
        """ There are no layers in the Willard-Chandler method.

            This function identifies the dividing surface and stores the
            triangulated isosurface, the density and the particles.

        """
        # Initialize detailed timing dictionary for this frame
        if self._enable_timing:
            frame_timings = {}

        self.reset_labels()
        # we assign an empty group for consistency
        self._layers, self.normal = self.universe.atoms[:0], None

        # Time prepare_box
        t0 = time.perf_counter()
        self.prepare_box()
        t1 = time.perf_counter()
        if self._enable_timing:
            frame_timings['prepare_box'] = t1 - t0

        # Time _define_cluster_group
        t0 = time.perf_counter()
        self._define_cluster_group()
        t1 = time.perf_counter()
        if self._enable_timing:
            frame_timings['define_cluster_group'] = t1 - t0

        # Time centering (if enabled)
        self.centered_positions = None
        if self.do_center is True:
            t0 = time.perf_counter()
            self.center()
            t1 = time.perf_counter()
            if self._enable_timing:
                frame_timings['center'] = t1 - t0

        # Time position extraction
        t0 = time.perf_counter()
        if self.include_zero_radius:
            pos = self.cluster_group.positions
        else:
            pos = self.cluster_group.positions[self.cluster_group.radii > 0.0]
        box = self.universe.dimensions[:3]
        t1 = time.perf_counter()
        if self._enable_timing:
            frame_timings['get_positions'] = t1 - t0

        # Time compute_surface
        start_time = time.perf_counter()
        result = compute_surface(pos,
                                 box,
                                 mesh=self.mesh,
                                 sigma=self.alpha,
                                 density_cutoff=self.density_cutoff,
                                 order='xyz',
                                 backend=self.surface_backend,
                                 backend_options=self.surface_backend_options)
        elapsed_time = time.perf_counter() - start_time

        if self._enable_timing:
            frame_timings['compute_surface'] = elapsed_time
            frame_timings['total'] = sum(frame_timings.values())

            # Accumulate timings
            if not hasattr(self, '_detailed_timings') or not self._detailed_timings:
                self._detailed_timings = {key: [] for key in frame_timings.keys()}

            for key, value in frame_timings.items():
                if key not in self._detailed_timings:
                    self._detailed_timings[key] = []
                self._detailed_timings[key].append(value)

            self._surface_computation_time = elapsed_time
            self._backend_timings[self.surface_backend] = elapsed_time

        self.spacing = result.spacing
        self.ngrid = result.ngrid
        self.density_field = result.density_field
        self.triangulated_surface = list(result.triangulated_surface)
        self.surface_area = result.surface_area
