"""
Implementation of centering algorithms for pywc.

Provides GPU (CuPy), C++ (pybind11), and pure Python implementations.
"""

import numpy as np
from . import utilities
from . import messages

# Try to import the fast C++ implementations
try:
    from .center_fast import center_fast, compute_mean, apply_shift
    HAS_CENTER_FAST = True
except ImportError:
    HAS_CENTER_FAST = False
    center_fast = None
    compute_mean = None
    apply_shift = None

try:
    from .center_fast_full import center_full_optimized
    HAS_CENTER_FULL = True
except ImportError:
    HAS_CENTER_FULL = False
    center_full_optimized = None
    print(' None 1')

# Try to import GPU implementation
try:
    from .center_gpu import center_gpu, GPUBackendUnavailable
    HAS_CENTER_GPU = True
except (ImportError, Exception):
    HAS_CENTER_GPU = False
    center_gpu = None
    GPUBackendUnavailable = None


def _center_python(group, direction, halfbox_shift=False):
    """
    Pure Python implementation of centering (original algorithm).

    This is the fallback implementation used when the C++ extension
    is not available.
    """
    dim = group.universe.coord.dimensions
    total_shift = 0

    if not (direction in {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}):
        raise ValueError(messages.WRONG_DIRECTION)

    _dir_map = {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}
    _dir = _dir_map[direction]

    _xyz_funcs = {
        0: (0, utilities.get_x),
        1: (1, utilities.get_y),
        2: (2, utilities.get_z)
    }

    direction_idx = _xyz_funcs[_dir][0]
    _pos_group = (_xyz_funcs[_dir][1])(group)

    shift = dim[direction_idx] / 100.

    _x = utilities.get_x(group.universe.atoms)
    _y = utilities.get_y(group.universe.atoms)
    _z = utilities.get_z(group.universe.atoms)

    _range = (0., dim[direction_idx])
    if halfbox_shift is True:
        _range = (-dim[direction_idx] / 2., dim[direction_idx] / 2.)

    histo, _ = np.histogram(_pos_group, bins=10, range=_range, density=True)

    max_val, min_val = np.amax(histo), np.amin(histo)
    delta = min_val + (max_val - min_val) / 3.

    # Iterative centering loop
    while (histo[0] > delta or histo[-1] > delta):
        total_shift += shift
        if total_shift >= dim[direction_idx]:
            raise ValueError(messages.CENTERING_FAILURE)
        _pos_group += shift

        from .interface import Interface
        Interface._attempt_shift(group, _pos_group, direction_idx,
                                 halfbox_shift, _dir)

        histo, _ = np.histogram(_pos_group, bins=10,
                                range=_range, density=True)

    _center_ = np.average(_pos_group)
    if halfbox_shift is False:
        box_half = dim[direction_idx] / 2.
    else:
        box_half = 0.

    _pos = {0: _x, 1: _y, 2: _z}
    if _dir in _pos.keys():
        _pos[_dir] += total_shift - _center_ + box_half

    group.universe.atoms.positions = np.column_stack((_x, _y, _z))


def _center_fast_cpp(group, direction, halfbox_shift=False):
    """
    Fast C++ implementation of centering using pybind11.

    This version is ~10-20x faster than the pure Python implementation.
    """
    dim = group.universe.coord.dimensions

    if not (direction in {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}):
        raise ValueError(messages.WRONG_DIRECTION)

    _dir_map = {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}
    _dir = _dir_map[direction]

    _xyz_funcs = {
        0: (0, utilities.get_x),
        1: (1, utilities.get_y),
        2: (2, utilities.get_z)
    }

    direction_idx = _xyz_funcs[_dir][0]

    # Get positions along centering direction (makes a copy)
    _pos_group = (_xyz_funcs[_dir][1])(group).copy()

    shift = dim[direction_idx] / 100.

    # Get full position arrays
    _x = utilities.get_x(group.universe.atoms)
    _y = utilities.get_y(group.universe.atoms)
    _z = utilities.get_z(group.universe.atoms)

    _range_min = -dim[direction_idx] / 2. if halfbox_shift else 0.
    _range_max = dim[direction_idx] / \
        2. if halfbox_shift else dim[direction_idx]

    # Compute initial histogram to get delta threshold
    histo, _ = np.histogram(_pos_group, bins=10, range=(
        _range_min, _range_max), density=True)
    max_val, min_val = np.amax(histo), np.amin(histo)
    delta = min_val + (max_val - min_val) / 3.

    # Call fast C++ implementation for the centering loop
    try:
        total_shift = center_fast(
            _pos_group,
            _range_min,
            _range_max,
            delta,
            shift,
            dim[direction_idx],
            bins=10,
            max_iter=100
        )
    except RuntimeError as e:
        raise ValueError(messages.CENTERING_FAILURE) from e

    # Compute final center and apply to all atoms
    _center_ = compute_mean(_pos_group)
    box_half = 0. if halfbox_shift else dim[direction_idx] / 2.
    final_shift = total_shift - _center_ + box_half

    # Apply final shift to the appropriate dimension
    _pos = {0: _x, 1: _y, 2: _z}
    if _dir in _pos.keys():
        apply_shift(_pos[_dir], final_shift)

    # Update universe positions
    group.universe.atoms.positions = np.column_stack((_x, _y, _z))


def _center_full_cpp(group, direction, halfbox_shift=False):
    """
    Fully optimized C++ centering that bypasses all MDAnalysis overhead.

    This version operates directly on position arrays without any Python
    callbacks during the centering loop.
    """
    if not (direction in {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}):
        raise ValueError(messages.WRONG_DIRECTION)

    _dir_map = {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}
    direction_idx = _dir_map[direction]

    # Get all atom positions (will be modified in-place)
    all_positions = group.universe.atoms.positions

    # Get group atom indices
    group_indices = group.indices.astype(np.int64)

    # Get box dimensions
    box = group.universe.dimensions[:3].astype(np.float64)

    # Call fully optimized C++ implementation
    try:
        total_shift = center_full_optimized(
            all_positions,
            group_indices,
            direction_idx,
            box,
            halfbox_shift
        )
    except RuntimeError as e:
        raise ValueError(messages.CENTERING_FAILURE) from e

    # Positions are already modified in-place, no need to copy back


def _center_gpu(group, direction, halfbox_shift=False):
    """
    GPU-accelerated centering using CuPy/CUDA.

    This version uses CUDA kernels for parallel computation of histograms
    and position updates. Requires CuPy to be installed.
    """
    if not HAS_CENTER_GPU:
        raise RuntimeError(
            "GPU centering not available. Install CuPy or use CPU backend.")

    if not (direction in {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}):
        raise ValueError(messages.WRONG_DIRECTION)

    _dir_map = {'x': 0, 'y': 1, 'z': 2, 0: 0, 1: 1, 2: 2}
    direction_idx = _dir_map[direction]

    # Get all atom positions (will be modified in-place on GPU)
    all_positions = group.universe.atoms.positions

    # Get group atom indices
    group_indices = group.indices.astype(np.int64)

    # Get box dimensions
    box = group.universe.dimensions[:3].astype(np.float64)

    # Call GPU implementation
    try:
        total_shift = center_gpu(
            all_positions,
            group_indices,
            direction_idx,
            box,
            halfbox_shift
        )
    except Exception as e:
        raise ValueError(f"GPU centering failed: {e}")


def center_wrapper(group, direction, halfbox_shift=False, force_python=False, use_gpu=False):
    """
    Wrapper that selects the appropriate centering implementation.

    Parameters
    ----------
    group : AtomGroup
        The group to center
    direction : int or str
        Centering direction (0/'x', 1/'y', 2/'z')
    halfbox_shift : bool
        If True, center to origin; if False, center to box middle
    force_python : bool
        If True, use Python implementation even if C++/GPU is available
    use_gpu : bool
        If True, use GPU implementation (requires CuPy)

    Notes
    -----
    Automatically selects the best available implementation:
    1. GPU (center_gpu) - if use_gpu=True and CuPy available
    2. Fully optimized C++ (center_fast_full) - bypasses all Python overhead
    3. Partially optimized C++ (center_fast) - only optimizes histogram
    4. Pure Python fallback
    """
    if force_python:
        return _center_python(group, direction, halfbox_shift)

    # Use GPU if requested and available
    if use_gpu:
        if HAS_CENTER_GPU:
            return _center_gpu(group, direction, halfbox_shift)
        else:
            raise RuntimeError(
                "GPU centering requested but CuPy is not available. Install CuPy or use CPU backend.")

    # Prefer the fully optimized C++ version
    if HAS_CENTER_FULL:
        return _center_full_cpp(group, direction, halfbox_shift)
    elif HAS_CENTER_FAST:
        return _center_fast_cpp(group, direction, halfbox_shift)
    else:
        return _center_python(group, direction, halfbox_shift)
