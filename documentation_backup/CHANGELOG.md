# Changelog
Changelog for the pytim project

## [Unreleased]
### Changed

- Added DistributionFunction observable in distributionfunction.py. 
  This is aims at being a class from which one can derive both density 
  profiles/maps and correlation functions.

- Added ReferenceFrame observable in basic_observables.py

- Willard-Chandler interface determination now uses a pybind11/OpenMP
  implementation of the Gaussian KDE with an internal neighbour search,
  yielding substantial speed-ups while keeping a transparent Python fallback.
- ``pip`` installation requires a C++17 compiler with OpenMP and the pybind11
  headers; documentation explains the new requirements and how to run the
  benchmark harness.
- Added ``benchmarks/willard_chandler_kde.py`` and wired it into CI to ensure
  the accelerated path maintains a minimum performance gain over the legacy
  implementation.
- Extracted the Willard–Chandler building blocks into ``pytim.wc_core`` so the
  KDE/grid/surface routines can be reused (and later extended with GPU
  back-ends) without touching the high-level interface class.
- Added an experimental CuPy-backed density evaluator (selectable via
  ``surface_backend='cupy'``) to demonstrate future GPU integration paths.
