"""Core building blocks used by the Willard–Chandler interface implementation."""

from .grid import build_grid
from .density import density_map
from .surface import compute_surface

__all__ = ["build_grid", "density_map", "compute_surface"]

try:  # pragma: no cover - optional dependency
    from .gpu import evaluate_density_gpu, GPUBackendUnavailable

    __all__.extend(["evaluate_density_gpu", "GPUBackendUnavailable"])
except Exception:  # pylint: disable=broad-except
    # CuPy (or its dependencies) is not available; keep CPU-only API.
    pass
