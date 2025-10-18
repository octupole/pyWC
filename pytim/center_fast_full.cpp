#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>

/*
 * FULLY OPTIMIZED centering that bypasses all Python/MDAnalysis overhead.
 *
 * This version operates directly on the position arrays without calling
 * back into Python, eliminating the bottleneck from repeated MDAnalysis
 * operations in the centering loop.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

namespace {

using Array2d = py::array_t<double, py::array::c_style | py::array::forcecast>;
using Array1d = py::array_t<double, py::array::c_style | py::array::forcecast>;

/**
 * Rebox coordinates to [shift, shift + box_length)
 */
inline void rebox_coords(double* coords, std::size_t n, double box_length, double shift) {
    const double inv_box = 1.0 / box_length;

    #pragma omp parallel for
    for (std::size_t i = 0; i < n; ++i) {
        // Shift to [0, box_length)
        double val = coords[i] - shift;
        // Apply periodic boundary
        val -= box_length * std::floor(val * inv_box);
        // Shift back
        coords[i] = val + shift;
    }
}

/**
 * Fast 1D histogram for centering algorithm.
 */
void compute_histogram_density(const double* data, std::size_t n_data,
                                int bins, double range_min, double range_max,
                                double* density) {
    const double bin_width = (range_max - range_min) / bins;
    const double inv_bin_width = 1.0 / bin_width;

    std::fill_n(density, bins, 0.0);
    std::vector<std::size_t> counts(bins, 0);

    #pragma omp parallel
    {
        std::vector<std::size_t> local_counts(bins, 0);

        #pragma omp for nowait
        for (std::size_t i = 0; i < n_data; ++i) {
            const double val = data[i];
            if (val >= range_min && val < range_max) {
                const int bin_idx = static_cast<int>((val - range_min) * inv_bin_width);
                if (bin_idx >= 0 && bin_idx < bins) {
                    local_counts[bin_idx]++;
                }
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < bins; ++i) {
                counts[i] += local_counts[i];
            }
        }
    }

    std::size_t total_count = 0;
    for (int i = 0; i < bins; ++i) {
        total_count += counts[i];
    }

    if (total_count > 0) {
        const double norm = 1.0 / (static_cast<double>(total_count) * bin_width);
        for (int i = 0; i < bins; ++i) {
            density[i] = static_cast<double>(counts[i]) * norm;
        }
    }
}

/**
 * FULLY OPTIMIZED centering that operates directly on position arrays.
 *
 * This avoids all Python/MDAnalysis overhead by:
 * 1. Working directly on NumPy arrays
 * 2. Not calling back into Python during the loop
 * 3. Handling reboxing in C++
 *
 * @param all_positions Full position array (N_atoms, 3) - modified in-place
 * @param group_indices Indices of atoms in the group to center
 * @param direction Dimension to center (0=x, 1=y, 2=z)
 * @param box Box dimensions [Lx, Ly, Lz]
 * @param halfbox_shift If true, center to origin; if false, to box middle
 * @return Total shift applied
 */
double center_full_optimized(Array2d all_positions,
                             py::array_t<int64_t> group_indices,
                             int direction,
                             Array1d box,
                             bool halfbox_shift) {

    // Validate inputs
    if (all_positions.ndim() != 2 || all_positions.shape(1) != 3) {
        throw std::invalid_argument("all_positions must have shape (N, 3)");
    }
    if (group_indices.ndim() != 1) {
        throw std::invalid_argument("group_indices must be 1D");
    }
    if (box.ndim() != 1 || box.shape(0) != 3) {
        throw std::invalid_argument("box must have length 3");
    }
    if (direction < 0 || direction > 2) {
        throw std::invalid_argument("direction must be 0, 1, or 2");
    }

    const std::size_t n_atoms = static_cast<std::size_t>(all_positions.shape(0));
    const std::size_t n_group = static_cast<std::size_t>(group_indices.shape(0));

    auto pos_ptr = all_positions.mutable_unchecked<2>();
    auto idx_ptr = group_indices.unchecked<1>();
    auto box_ptr = box.unchecked<1>();

    const double box_length = box_ptr(direction);
    const double shift_increment = box_length / 100.0;

    const double range_min = halfbox_shift ? -box_length / 2.0 : 0.0;
    const double range_max = halfbox_shift ? box_length / 2.0 : box_length;
    const double pbc_shift = halfbox_shift ? box_length / 2.0 : 0.0;

    // Extract group positions along centering direction
    std::vector<double> group_pos(n_group);
    for (std::size_t i = 0; i < n_group; ++i) {
        const auto atom_idx = static_cast<std::size_t>(idx_ptr(i));
        if (atom_idx >= n_atoms) {
            throw std::out_of_range("group_indices contains invalid index");
        }
        group_pos[i] = pos_ptr(atom_idx, direction);
    }

    // Compute initial histogram
    std::vector<double> density(10);
    compute_histogram_density(group_pos.data(), n_group, 10, range_min, range_max, density.data());

    double max_val = *std::max_element(density.begin(), density.end());
    double min_val = *std::min_element(density.begin(), density.end());
    const double delta = min_val + (max_val - min_val) / 3.0;

    double total_shift = 0.0;
    const int max_iter = 100;

    // Release GIL for the main loop
    {
        py::gil_scoped_release release;

        for (int iteration = 0; iteration < max_iter; ++iteration) {
            // Check convergence
            if (!(density[0] > delta || density[9] > delta)) {
                break;
            }

            // Check for failure
            total_shift += shift_increment;
            if (total_shift >= box_length) {
                throw std::runtime_error("Centering failure: maximum shift exceeded");
            }

            // Apply shift to group positions ONLY (don't update universe yet)
            #pragma omp parallel for
            for (std::size_t i = 0; i < n_group; ++i) {
                group_pos[i] += shift_increment;
                // Apply PBC to group positions
                double val = group_pos[i] - pbc_shift;
                val -= box_length * std::floor(val / box_length);
                group_pos[i] = val + pbc_shift;
            }

            // Recompute histogram
            compute_histogram_density(group_pos.data(), n_group, 10, range_min, range_max, density.data());
        }

        // Final centering: compute mean and shift to target
        double group_mean = 0.0;
        for (std::size_t i = 0; i < n_group; ++i) {
            group_mean += group_pos[i];
        }
        group_mean /= static_cast<double>(n_group);

        const double box_half = halfbox_shift ? 0.0 : box_length / 2.0;
        const double final_shift = total_shift - group_mean + box_half;

        // Apply total shift to all atoms (total_shift from iterations + final centering)
        #pragma omp parallel for
        for (std::size_t i = 0; i < n_atoms; ++i) {
            double val = pos_ptr(i, direction) + total_shift;
            // Apply PBC
            val -= pbc_shift;
            val -= box_length * std::floor(val / box_length);
            val += pbc_shift;
            // Now apply the final centering shift
            pos_ptr(i, direction) = val - group_mean + box_half;
        }
    }

    return total_shift;
}

}  // namespace

PYBIND11_MODULE(center_fast_full, m) {
    m.doc() = "Fully optimized centering that bypasses Python/MDAnalysis overhead";

    m.def("center_full_optimized", &center_full_optimized,
          py::arg("all_positions"),
          py::arg("group_indices"),
          py::arg("direction"),
          py::arg("box"),
          py::arg("halfbox_shift"),
          R"pbdoc(
          Fully optimized centering with no Python callbacks.

          This version operates directly on position arrays without calling
          back into Python/MDAnalysis, eliminating the bottleneck from
          repeated universe.coord.positions updates in the loop.

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

          Returns
          -------
          total_shift : float
              Total shift applied

          Notes
          -----
          - all_positions is modified in-place
          - Much faster than the version that calls back into Python
          - Uses OpenMP for parallelization
          )pbdoc");
}
