pyWC: Willard-Chandler Surface Analysis Toolkit
================================================

**Note:** For complete documentation, please see `README.md <README.md>`_

pyWC is a high-performance Python toolkit for computing Willard-Chandler intrinsic
density surfaces from molecular dynamics simulations.

**This is a specialized fork of** `pytim <https://github.com/Marcello-Sega/pytim>`_
**with ~35x performance improvements and GPU acceleration.**

Upstream Project
----------------

* Original Project: `Marcello-Sega/pytim <https://github.com/Marcello-Sega/pytim>`_
* Original Paper: Sega, M., Fabian, B., & Jedlovszky, P. (2018). *Pytim: A python package
  for the interfacial analysis of molecular simulations.* Journal of Computational Chemistry,
  39(25), 2118-2125. DOI: `10.1002/jcc.25384 <https://doi.org/10.1002/jcc.25384>`_

Quick Start
-----------

Installation::

    pip install .

Basic usage::

    import MDAnalysis as mda
    from pywc import WillardChandler

    u = mda.Universe("system.pdb", "trajectory.xtc")
    group = u.select_atoms("resname DPPC")

    wc = WillardChandler(u, group=group, alpha=3.0, mesh=2.0,
                         surface_backend='cpp')  # or 'cupy' for GPU

    vertices, faces, normals = wc.triangulated_surface
    print(f"Surface area: {wc.surface_area:.2f} Ų")
    wc.writevtk.surface("surface.vtk")

Performance Improvements
------------------------

* **~35x faster** KDE evaluation (C++ vs Python)
* **~10-20x faster** system centering
* **GPU acceleration** for large systems (>10k atoms)
* OpenMP parallelization throughout

Key Features
------------

* Multiple computational backends (CPU C++, GPU CUDA, Python fallback)
* Export to VTK, Wavefront OBJ, Gaussian Cube, PDB formats
* Trajectory analysis with frame-to-frame continuity
* Command-line tools for batch processing

Documentation
-------------

* `README.md <README.md>`_ - Complete documentation
* `CITATION.cff <CITATION.cff>`_ - Citation information
* `AUTHORS.md <AUTHORS.md>`_ - Contributors
* `CHANGELOG.md <CHANGELOG.md>`_ - Version history
* `CONTRIBUTING.md <CONTRIBUTING.md>`_ - Contribution guidelines

Maintainer
----------

Massimo Marchi (massimo@octupole.org)

Original pytim authors: Marcello Sega, Balazs Fabian, Gyorgy Hantal, and
Pal Jedlovszky are gratefully acknowledged for creating the upstream project.

License
-------

GNU General Public License v3.0 (GPL-3.0-only)

Repository
----------

https://github.com/octupole/pyWC
