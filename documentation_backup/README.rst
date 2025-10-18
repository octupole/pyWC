A Python Tool for Interfacial Molecules Analysis
================================================

Pytim is a cross-platform python implementation of several methods
for the detection of fluid interfaces in molecular simulations: it
allows to identify interfacial molecules from the trajectories of
major molecular dynamics simulation packages, and run some analyses
specifically conceived for interfacial molecules, such as intrinsic
profiles.

About this Fork
---------------

This fork is maintained by **Massimo** (massimo@octupole.org) and builds upon
the original ``pytim`` project authored by Marcello Sega and collaborators.
The upstream developers are gratefully acknowledged for their foundational
work; this fork focuses on documenting and extending the Willard–Chandler
surface analysis pipeline while preserving compatibility with the original
code base.

So far the following methods have been implemented:

* ITIM
* GITIM 
* SASA
* Willard Chandler — evaluates a Gaussian KDE under periodic boundary
  conditions to construct a smooth intrinsic surface, with optional GPU
  acceleration.
* DBSCAN filtering

Willard–Chandler Density Surfaces
---------------------------------

The Willard–Chandler functionality follows the approach of Willard and
Chandler (J. Phys. Chem. B 2010, 114, 1954–1958). Particle coordinates are
mapped onto a regular mesh, a Gaussian kernel density estimate is computed in
periodic space, and marching cubes extracts the dividing surface from the
resulting intrinsic density field. ``pytim`` accelerates the calculation with
pybind11/OpenMP extensions and, when available, CuPy-powered GPU kernels. This
fork provides expanded documentation of those stages while keeping behaviour
aligned with the upstream implementation.

----

Pytim relies on the MDAnalysis package for reading/writing trajectories,
and work therefore seamlessly for a number of popular trajectory
formats, including:

* GROMACS
* CHARMM/NAMD
* LAMMPS
* AMBER
* DL_Poly

