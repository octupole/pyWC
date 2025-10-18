Willard–Chandler Surface Toolkit
================================

This fork of ``pytim`` provides a focused toolkit for building and analysing
Willard–Chandler intrinsic density surfaces from molecular simulations.  The
code keeps the original high-performance kernels (OpenMP-enabled pybind11,
CuPy GPU backend, cell-list neighbour search) without the legacy ITIM/GITIM
interfaces or unrelated observables.

Main features
-------------

* Compute intrinsic density fields on structured grids with either the C++
  (OpenMP) or CuPy backend.
* Extract triangulated surfaces via ``skimage.measure.marching_cubes`` and
  export them as VTK, Wavefront OBJ, or Gaussian cube files.
* Optional timing diagnostics and backend benchmarking scripts for CPU/GPU
  comparison, surface area evaluation, and bending-rigidity analysis.

Installation
------------

.. code-block:: bash

   pip install .

The CPU path requires NumPy, SciPy, scikit-image, MDAnalysis, and packaging.
For GPU acceleration install a CuPy wheel matching your CUDA toolkit, or use:

.. code-block:: bash

   pip install pytim[gpu]

Quick start
-----------

.. code-block:: python

   import MDAnalysis as mda
   from pytim import WillardChandler

   u = mda.Universe("system.pdb", "trajectory.xtc")
   group = u.select_atoms("resname DPPC")
   wc = WillardChandler(u, group=group, alpha=3.0, mesh=2.0, centered=True)
   wc.surface_backend = "cpp"        # or "cupy" if CuPy is available
   wc._assign_layers()               # triggers the surface calculation

   verts, faces, normals = wc.triangulated_surface
   print("Surface area:", wc.surface_area)
   wc.writevtk.surface("surface.vtk")

Project status
--------------

Maintained by Massimo (massimo@octupole.org).  Original pytim authors
Marcello Sega, Balazs Fabian, Gyorgy Hantal, and Pal Jedlovszky are gratefully
acknowledged for creating the upstream project that this fork builds upon.

