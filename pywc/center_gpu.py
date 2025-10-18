"""CuPy-backed GPU-accelerated centering algorithm for Willard-Chandler."""

from __future__ import annotations

import numpy as np

try:
    import cupy as cp
except Exception:
    cp = None


class GPUBackendUnavailable(RuntimeError):
    """Raised when the CuPy backend cannot be used."""


def _ensure_cupy() -> None:
    if cp is None:
        raise GPUBackendUnavailable(
            "CuPy is not available. Install cupy-cudaXX or use CPU centering.")


# CUDA kernel for parallel histogram computation
_HISTOGRAM_KERNEL = r"""
extern "C" __global__
void compute_histogram(
    const double* __restrict__ data,
    const int n_data,
    const double range_min,
    const double range_max,
    const int bins,
    int* __restrict__ counts)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= n_data) return;

    double val = data[idx];
    if (val >= range_min && val < range_max) {
        double bin_width = (range_max - range_min) / bins;
        int bin_idx = (int)((val - range_min) / bin_width);
        if (bin_idx >= 0 && bin_idx < bins) {
            atomicAdd(&counts[bin_idx], 1);
        }
    }
}
"""

# CUDA kernel for applying shifts with PBC
_SHIFT_PBC_KERNEL = r"""
extern "C" __global__
void apply_shift_pbc(
    double* __restrict__ positions,
    const int n_atoms,
    const double shift,
    const double box_length,
    const double pbc_shift)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= n_atoms) return;

    double val = positions[idx] + shift;
    val -= pbc_shift;
    val -= box_length * floor(val / box_length);
    positions[idx] = val + pbc_shift;
}
"""

# CUDA kernel for final centering
_FINAL_CENTER_KERNEL = r"""
extern "C" __global__
void apply_final_center(
    double* __restrict__ positions,
    const int n_atoms,
    const double total_shift,
    const double group_mean,
    const double box_half,
    const double box_length,
    const double pbc_shift)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= n_atoms) return;

    double val = positions[idx] + total_shift;
    // Apply PBC
    val -= pbc_shift;
    val -= box_length * floor(val / box_length);
    val += pbc_shift;
    // Apply final centering
    positions[idx] = val - group_mean + box_half;
}
"""

_HISTOGRAM_KERNEL_OBJ = None
_SHIFT_PBC_KERNEL_OBJ = None
_FINAL_CENTER_KERNEL_OBJ = None


def _get_histogram_kernel() -> "cp.RawKernel":
    global _HISTOGRAM_KERNEL_OBJ
    if _HISTOGRAM_KERNEL_OBJ is None:
        _HISTOGRAM_KERNEL_OBJ = cp.RawKernel(_HISTOGRAM_KERNEL, "compute_histogram")
    return _HISTOGRAM_KERNEL_OBJ


def _get_shift_pbc_kernel() -> "cp.RawKernel":
    global _SHIFT_PBC_KERNEL_OBJ
    if _SHIFT_PBC_KERNEL_OBJ is None:
        _SHIFT_PBC_KERNEL_OBJ = cp.RawKernel(_SHIFT_PBC_KERNEL, "apply_shift_pbc")
    return _SHIFT_PBC_KERNEL_OBJ


def _get_final_center_kernel() -> "cp.RawKernel":
    global _FINAL_CENTER_KERNEL_OBJ
    if _FINAL_CENTER_KERNEL_OBJ is None:
        _FINAL_CENTER_KERNEL_OBJ = cp.RawKernel(_FINAL_CENTER_KERNEL, "apply_final_center")
    return _FINAL_CENTER_KERNEL_OBJ


def center_gpu(all_positions: np.ndarray,
               group_indices: np.ndarray,
               direction: int,
               box: np.ndarray,
               halfbox_shift: bool = False,
               max_iter: int = 100) -> float:
    """
    GPU-accelerated centering algorithm using CuPy.

    This is the GPU equivalent of center_fast_full.cpp, using CuPy raw CUDA
    kernels for maximum performance.

    Parameters
    ----------
    all_positions : ndarray (N_atoms, 3)
        All atom positions (modified in-place)
    group_indices : ndarray (N_group,) of int64
        Indices of atoms in the group to center
    direction : int
        Centering dimension (0=x, 1=y, 2=z)
    box : ndarray (3,)
        Box dimensions [Lx, Ly, Lz]
    halfbox_shift : bool
        If True, center to origin; if False, to box middle
    max_iter : int
        Maximum iterations for centering loop

    Returns
    -------
    float
        Total shift applied

    Notes
    -----
    This function uses CuPy RawKernel for direct CUDA kernel execution,
    following the same pattern as gpu.py's evaluate_density_gpu().

    The algorithm:
    1. Extract group positions along centering direction
    2. Compute histogram to determine density threshold
    3. Iteratively shift group positions until centered (GPU)
    4. Apply final shift to all atoms (GPU)

    Performance: Expected 5-10x faster than CPU for large systems (>50K atoms)
    """
    _ensure_cupy()

    # Validate inputs
    if all_positions.ndim != 2 or all_positions.shape[1] != 3:
        raise ValueError("all_positions must have shape (N, 3)")
    if group_indices.ndim != 1:
        raise ValueError("group_indices must be 1D")
    if box.ndim != 1 or box.shape[0] != 3:
        raise ValueError("box must have length 3")
    if direction < 0 or direction > 2:
        raise ValueError("direction must be 0, 1, or 2")

    n_atoms = all_positions.shape[0]
    n_group = group_indices.shape[0]

    box_length = float(box[direction])
    shift_increment = box_length / 100.0

    range_min = -box_length / 2.0 if halfbox_shift else 0.0
    range_max = box_length / 2.0 if halfbox_shift else box_length
    pbc_shift = box_length / 2.0 if halfbox_shift else 0.0

    bins = 10

    # Transfer to GPU
    group_indices_d = cp.asarray(group_indices, dtype=cp.int64)

    # Extract group positions along centering direction
    group_pos_h = all_positions[group_indices, direction].copy()
    group_pos_d = cp.asarray(group_pos_h, dtype=cp.float64)

    # Compute initial histogram on GPU
    def compute_histogram_gpu(data_d: "cp.ndarray") -> "cp.ndarray":
        """Compute histogram using custom CUDA kernel."""
        counts_d = cp.zeros(bins, dtype=cp.int32)
        kernel = _get_histogram_kernel()

        threads = 256
        blocks = (len(data_d) + threads - 1) // threads

        kernel((blocks,), (threads,), (
            data_d,
            np.int32(len(data_d)),
            np.float64(range_min),
            np.float64(range_max),
            np.int32(bins),
            counts_d
        ))

        # Convert to density
        total_count = cp.sum(counts_d)
        bin_width = (range_max - range_min) / bins
        density_d = counts_d.astype(cp.float64) / (total_count * bin_width + 1e-10)
        return density_d

    density_d = compute_histogram_gpu(group_pos_d)
    density_h = cp.asnumpy(density_d)

    max_val = np.max(density_h)
    min_val = np.min(density_h)
    delta = min_val + (max_val - min_val) / 3.0

    total_shift = 0.0
    shift_pbc_kernel = _get_shift_pbc_kernel()

    threads = 256
    blocks_group = (n_group + threads - 1) // threads

    # Centering loop (GPU)
    for iteration in range(max_iter):
        # Check convergence
        if not (density_h[0] > delta or density_h[9] > delta):
            break

        # Check for failure
        total_shift += shift_increment
        if total_shift >= box_length:
            raise RuntimeError("Centering failure: maximum shift exceeded")

        # Apply shift to group positions with PBC (GPU kernel)
        shift_pbc_kernel((blocks_group,), (threads,), (
            group_pos_d,
            np.int32(n_group),
            np.float64(shift_increment),
            np.float64(box_length),
            np.float64(pbc_shift)
        ))

        # Recompute histogram
        density_d = compute_histogram_gpu(group_pos_d)
        density_h = cp.asnumpy(density_d)

    # Compute group mean (GPU)
    group_mean = float(cp.mean(group_pos_d).item())

    box_half = 0.0 if halfbox_shift else box_length / 2.0

    # Apply final shift to all atoms (GPU kernel)
    # Transfer all positions to GPU
    all_pos_dir_d = cp.asarray(all_positions[:, direction], dtype=cp.float64)

    final_center_kernel = _get_final_center_kernel()
    blocks_all = (n_atoms + threads - 1) // threads

    final_center_kernel((blocks_all,), (threads,), (
        all_pos_dir_d,
        np.int32(n_atoms),
        np.float64(total_shift),
        np.float64(group_mean),
        np.float64(box_half),
        np.float64(box_length),
        np.float64(pbc_shift)
    ))

    # Transfer back to CPU
    all_positions[:, direction] = cp.asnumpy(all_pos_dir_d)

    return total_shift


__all__ = ["center_gpu", "GPUBackendUnavailable"]
