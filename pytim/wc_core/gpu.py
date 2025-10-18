"""CuPy-backed helpers mirroring the C++ Willard–Chandler KDE routine."""

from __future__ import annotations

import numpy as np

try:
    import cupy as cp
except Exception:  # pragma: no cover - optional dependency
    cp = None

import time


class GPUBackendUnavailable(RuntimeError):
    """Raised when the CuPy backend cannot be used."""


def _ensure_cupy() -> None:
    if cp is None:
        raise GPUBackendUnavailable(
            "CuPy is not available. Install cupy-cudaXX or use surface_backend='cpu'.")


def _build_grid_cell_list(grid_d: "cp.ndarray",
                          box_d: "cp.ndarray",
                          cutoff: float) -> tuple["cp.ndarray", "cp.ndarray", "cp.ndarray", int, int, int]:
    cutoff_d = cp.asarray(cutoff, dtype=cp.float32)
    ncell_d = cp.maximum(1, cp.ceil(box_d / cutoff_d)).astype(cp.int32)
    cell_size_d = box_d / ncell_d.astype(cp.float32)
    nx = int(ncell_d[0].item())
    ny = int(ncell_d[1].item())
    nz = int(ncell_d[2].item())
    total_cells = nx * ny * nz

    def cell_index_coords(coords: "cp.ndarray", L: float, cs: float, dim: int) -> "cp.ndarray":
        wrapped = cp.mod(coords, L)
        idx = cp.floor(wrapped / cs).astype(cp.int32)
        return cp.minimum(idx, dim - 1)

    ix = cell_index_coords(grid_d[:, 0], box_d[0], cell_size_d[0], nx)
    iy = cell_index_coords(grid_d[:, 1], box_d[1], cell_size_d[1], ny)
    iz = cell_index_coords(grid_d[:, 2], box_d[2], cell_size_d[2], nz)

    cell_id = ix + nx * (iy + ny * iz)
    order = cp.argsort(cell_id)
    grid_idx_sorted = order.astype(cp.int32)
    sorted_cell_id = cell_id[order]

    counts = cp.bincount(sorted_cell_id, minlength=total_cells)
    cell_starts = cp.empty(total_cells + 1, dtype=cp.int64)
    cell_starts[0] = 0
    cell_starts[1:] = cp.cumsum(counts, dtype=cp.int64)
    return grid_idx_sorted, cell_starts, cell_size_d.astype(cp.float32), nx, ny, nz


_KERNEL_SOURCE = r"""
extern "C" __global__
void accumulate_gaussians_f32(
    const float* __restrict__ grid,
    const float* __restrict__ pos,
    const float* __restrict__ box_len,
    const float* __restrict__ half_box,
    const float* __restrict__ cell_size,
    const int nx, const int ny, const int nz,
    const long long* __restrict__ cell_starts,
    const int* __restrict__ grid_idx_sorted,
    const float cutoff_sq,
    const float scale,
    float* __restrict__ out,
    const int K)
{
    int p = blockDim.x * blockIdx.x + threadIdx.x;
    if (p >= K) return;

    const float px = pos[3 * p + 0];
    const float py = pos[3 * p + 1];
    const float pz = pos[3 * p + 2];

    auto cell_index = [&](float v, float L, float cs, int dim) -> int {
        float w = fmodf(v, L);
        if (w < 0.0f) w += L;
        int c = (int)floorf(w / cs);
        if (c >= dim) c = dim - 1;
        return c;
    };

    const int cx = cell_index(px, box_len[0], cell_size[0], nx);
    const int cy = cell_index(py, box_len[1], cell_size[1], ny);
    const int cz = cell_index(pz, box_len[2], cell_size[2], nz);

    auto wrap = [](int v, int dim) -> int {
        int m = v % dim;
        if (m < 0) m += dim;
        return m;
    };

    auto flat_index = [&](int ix, int iy, int iz) -> long long {
        return (long long)ix + (long long)nx * ((long long)iy + (long long)ny * (long long)iz);
    };

    for (int dx = -1; dx <= 1; ++dx) {
        const int nx_i = wrap(cx + dx, nx);
        for (int dy = -1; dy <= 1; ++dy) {
            const int ny_i = wrap(cy + dy, ny);
            for (int dz = -1; dz <= 1; ++dz) {
                const int nz_i = wrap(cz + dz, nz);
                const long long cid = flat_index(nx_i, ny_i, nz_i);
                const long long start = cell_starts[cid];
                const long long end = cell_starts[cid + 1];
                for (long long k = start; k < end; ++k) {
                    const int gi = grid_idx_sorted[k];
                    const float gx = grid[3 * gi + 0];
                    const float gy = grid[3 * gi + 1];
                    const float gz = grid[3 * gi + 2];

                    float dxv = gx - px;
                    float dyv = gy - py;
                    float dzv = gz - pz;

                    if (dxv > half_box[0]) dxv -= box_len[0];
                    if (dxv < -half_box[0]) dxv += box_len[0];
                    if (dyv > half_box[1]) dyv -= box_len[1];
                    if (dyv < -half_box[1]) dyv += box_len[1];
                    if (dzv > half_box[2]) dzv -= box_len[2];
                    if (dzv < -half_box[2]) dzv += box_len[2];

                    const float dist2 = dxv * dxv + dyv * dyv + dzv * dzv;
                    if (dist2 > cutoff_sq) continue;

                    const float dens = __expf(-dist2 / scale);
                    atomicAdd(&out[gi], dens);
                }
            }
        }
    }
}
"""

_RAW_KERNEL = None


def _get_kernel() -> "cp.RawKernel":
    global _RAW_KERNEL
    if _RAW_KERNEL is None:
        _RAW_KERNEL = cp.RawKernel(_KERNEL_SOURCE, "accumulate_gaussians_f32")
    return _RAW_KERNEL


def evaluate_density_gpu(pos: np.ndarray,
                         grid: np.ndarray,
                         sigma: float,
                         box: np.ndarray,
                         *,
                         chunk_size: int = 4096,
                         max_chunk_mem: int = 256 * 1024 * 1024) -> np.ndarray:
    """CuPy implementation of ``evaluate_pbc_fast_auto`` returning float64."""
    _ensure_cupy()

# Method 1: Manual timing with synchronization
    begin = cp.cuda.Event()
    stop = cp.cuda.Event()

    begin.record()
    pos_h = np.asarray(pos, dtype=np.float32, order="C")
    grid_arr = np.asarray(grid, dtype=np.float32, order="C")
    if grid_arr.ndim != 2:
        raise ValueError("grid must be a 2D array")

    if grid_arr.shape[0] == 3:
        grid_h = np.ascontiguousarray(grid_arr[::-1].T)
        n_points = grid_arr.shape[1]
    elif grid_arr.shape[1] == 3:
        grid_h = np.ascontiguousarray(grid_arr[:, ::-1])
        n_points = grid_arr.shape[0]
    else:
        raise ValueError("grid must have shape (3, N) or (N, 3)")
    box_h = np.asarray(box, dtype=np.float32).reshape(3)

    if pos_h.ndim != 2 or pos_h.shape[1] != 3:
        raise ValueError("pos must have shape (N, 3)")
    if grid_h.ndim != 2 or grid_h.shape[1] != 3:
        raise ValueError("grid must have shape (M, 3)")
    if np.any(box_h <= 0):
        raise ValueError("Box lengths must be positive")

    sigma32 = np.float32(sigma)
    cutoff = np.float32(2.5) * sigma32
    cutoff_sq = cutoff * cutoff
    scale = np.float32(2.0) * sigma32 * sigma32

    grid_d = cp.asarray(grid_h)
    box_d = cp.asarray(box_h)
    half_box_d = box_d * np.float32(0.5)

    grid_idx_sorted, cell_starts, cell_size_d, nx, ny, nz = _build_grid_cell_list(
        grid_d, box_d, float(cutoff))

    out_d = cp.zeros(grid_d.shape[0], dtype=cp.float32)
    kernel = _get_kernel()

    bytes_per_particle = 3 * 4
    if bytes_per_particle <= 0:
        effective_chunk = chunk_size
    else:
        max_by_mem = max(1, max_chunk_mem // bytes_per_particle)
        effective_chunk = int(max(1, min(chunk_size, max_by_mem)))

    threads = 128
    M = pos_h.shape[0]
    for start in range(0, M, effective_chunk):
        end = min(start + effective_chunk, M)
        ksize = end - start
        if ksize <= 0:
            continue
        pos_chunk_d = cp.asarray(pos_h[start:end], dtype=cp.float32, order="C")
        blocks = (ksize + threads - 1) // threads
        kernel((blocks,), (threads,), (
            grid_d.ravel(),
            pos_chunk_d.ravel(),
            box_d,
            half_box_d,
            cell_size_d,
            np.int32(nx),
            np.int32(ny),
            np.int32(nz),
            cell_starts,
            grid_idx_sorted,
            np.float32(cutoff_sq),
            np.float32(scale),
            out_d,
            np.int32(ksize),
        ))
    stop.record()
    stop.synchronize()
    elapsed_time = cp.cuda.get_elapsed_time(begin, stop)  # milliseconds
#    print(f"GPU time: {elapsed_time:.2f} ms")

    return cp.asnumpy(out_d[:n_points]).astype(np.float64, copy=False)


__all__ = ["evaluate_density_gpu", "GPUBackendUnavailable"]
